import numpy as np
from scipy.optimize import root
from typing import Dict, List

from reaction import Reaction, ReactionNetwork
from kinetics import KineticsCalculator


class GolgiModel:
    """
    Solves the steady-state mass balance equations for 4 cascaded CSTRs.
    Steady-state equation: c_{i,j} = c_{i,j-1} + tau_j * r_{i,j}
    """
    def __init__(
        self,
        network: ReactionNetwork,
        kinetics: KineticsCalculator,
        tau: np.ndarray,
        enzyme_dist: List[Dict[str, float]],
        donor_conc: Dict[str, float],
        initial_feed: np.ndarray,
    ):
        self.network = network
        self.kinetics = kinetics
        self.tau = tau
        self.enzyme_dist = enzyme_dist          # [{enzyme_name: concentration}, ...] one dict per compartment
        self.donor_conc = donor_conc            # {cosubstrate name: concentration}
        self.initial_feed = initial_feed
        self.n_structs = len(network.structures)
        self.n_compartments = len(tau)

    def _calc_single_rate(
        self,
        rxn: Reaction,
        c_j: np.ndarray,
        competitor_sums: Dict[str, float],
        enzyme_conc_j: Dict[str, float],
    ) -> float:
        """
        Computes the rate of a single reaction:
        r = kf * [Et] * [donor] * [Pi] / (Kmi * (Kmd + [donor]) * (1 + sum(Pj)/Kmj))
        """
        enzyme = rxn.enzyme
        sub_id = rxn.substrate_id
        kmi = rxn.km_value
        kf = enzyme.kf
        kmd = enzyme.Kmd

        et = enzyme_conc_j.get(enzyme.enzyme_name, 0.0)
        donor = self.donor_conc.get(enzyme.cosubstrate, 1.0)
        pi = c_j[sub_id] if 0 <= sub_id < len(c_j) else 0.0
        cs = competitor_sums.get(enzyme.enzyme_name, 1.0)

        denom = kmi * (kmd + donor) * cs
        if denom <= 0 or et <= 0:
            return 0.0
        return kf * et * donor * pi / denom

    def _calculate_net_rates(
        self,
        c_j: np.ndarray,
        competitor_sums: Dict[str, float],
        enzyme_conc_j: Dict[str, float],
        compartment_idx: int,
    ) -> np.ndarray:
        """
        Computes the net production rate for each glycoform in the specified compartment.
        Caches all reaction rates first, then aggregates by structure using
        Glycoform._consuming_reactions / _producing_reactions.
        """
        rxn_rates = {}
        for rxn_id, rxn in self.network.reactions.items():
            rxn_rates[rxn_id] = self._calc_single_rate(
                rxn, c_j, competitor_sums, enzyme_conc_j
            )

        r = np.zeros(self.n_structs)
        for i, glycan in self.network.structures.items():
            for rxn_id in glycan._consuming_reactions:
                r[i] -= rxn_rates.get(rxn_id, 0.0)
            for rxn_id in glycan._producing_reactions:
                r[i] += rxn_rates.get(rxn_id, 0.0)

        return r

    def residual(self, x: np.ndarray) -> np.ndarray:
        """
        Computes the residual F(x).
        x : shape (n_structs * n_compartments,)
        """
        c = x.reshape(self.n_structs, self.n_compartments)
        F = np.zeros_like(c)

        for j in range(self.n_compartments):
            c_j = c[:, j]
            c_prev = c[:, j - 1] if j > 0 else self.initial_feed

            competitor_sums = self.kinetics.compute_all_competitor_sums(c_j)
            enzyme_conc_j = self.enzyme_dist[j]
            r = self._calculate_net_rates(c_j, competitor_sums, enzyme_conc_j, j)

            F[:, j] = c_j - (c_prev + self.tau[j] * r)

        return F.flatten()

    def solve(
        self,
        method: str = 'krylov',
        tol: float = 1e-6,
        max_iter: int = 100,
        verbose: bool = True,
        **kwargs,
    ) -> np.ndarray:
        """
        Solves the steady state using scipy.optimize.root.
        """
        x0 = np.zeros(self.n_structs * self.n_compartments)
        for j in range(self.n_compartments):
            x0[j * self.n_structs : (j + 1) * self.n_structs] = self.initial_feed

        options = {'maxiter': max_iter, **kwargs}
        if verbose:
            options['disp'] = True

        result = root(self.residual, x0, method=method, tol=tol, options=options)

        if result.success:
            c = result.x.reshape(self.n_structs, self.n_compartments)
            return np.maximum(c, 0)
        raise RuntimeError(f"Steady-state solver did not converge: {result.message}")

    def solve_damped(
        self,
        tol: float = 1e-6,
        max_iter: int = 200,
        damping: float = 0.3,
        verbose: bool = True,
    ) -> np.ndarray:
        """
        Damped Picard iteration (no Jacobian computation, memory-friendly).
        """
        x = np.zeros(self.n_structs * self.n_compartments)
        for j in range(self.n_compartments):
            x[j * self.n_structs : (j + 1) * self.n_structs] = self.initial_feed

        for it in range(max_iter):
            F = self.residual(x)
            norm = np.linalg.norm(F)
            if verbose and it % 10 == 0:
                print(f"Iter {it}: ||F|| = {norm:.6e}")
            if norm < tol:
                if verbose:
                    print(f"Converged at iteration {it}, ||F|| = {norm:.6e}")
                break
            x = x - damping * F
            x = np.maximum(x, 0)

        return x.reshape(self.n_structs, self.n_compartments)
