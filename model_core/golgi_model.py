import warnings

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import root
from typing import Dict, List, Tuple
from collections import defaultdict

from reaction import Reaction, ReactionNetwork
from kinetics import KineticsCalculator


class GolgiModel:
    """
    Solves the steady-state mass balance equations for 4 cascaded CSTRs.
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

        # ---- helpers for the sparse extended-system solver ----
        enzyme_names_set = sorted(
            {rxn.enzyme.enzyme_name for rxn in self.network.reactions.values()}
        )
        self.enzyme_names = enzyme_names_set
        self.n_enzymes = len(enzyme_names_set)
        self.enzyme_idx = {name: i for i, name in enumerate(enzyme_names_set)}

        # Per-enzyme list of reaction ids (deterministic ordering)
        self.enzyme_rxn_ids: Dict[str, List[int]] = defaultdict(list)
        for rxn_id, rxn in self.network.reactions.items():
            self.enzyme_rxn_ids[rxn.enzyme.enzyme_name].append(rxn_id)
        for name in self.enzyme_rxn_ids:
            self.enzyme_rxn_ids[name].sort()

        self.block_size = self.n_structs + self.n_enzymes
        self.n_total = self.block_size * self.n_compartments

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
        Legacy residual F(x) on the concentration-only variable layout.
        x : shape (n_structs * n_compartments,)
        """
        c = x.reshape(self.n_compartments, self.n_structs)
        F = np.zeros_like(c)

        for j in range(self.n_compartments):
            c_j = c[j]
            c_prev = c[j - 1] if j > 0 else self.initial_feed

            competitor_sums = self.kinetics.compute_all_competitor_sums(c_j)
            enzyme_conc_j = self.enzyme_dist[j]
            r = self._calculate_net_rates(c_j, competitor_sums, enzyme_conc_j, j)

            F[j] = c_j - (c_prev + self.tau[j] * r)

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
        Legacy wrapper around scipy.optimize.root.
        """
        x0 = np.zeros(self.n_structs * self.n_compartments)
        for j in range(self.n_compartments):
            x0[j * self.n_structs : (j + 1) * self.n_structs] = self.initial_feed

        options = {'maxiter': max_iter, **kwargs}
        if verbose:
            options['disp'] = True

        result = root(self.residual, x0, method=method, tol=tol, options=options)

        if result.success:
            c = result.x.reshape(self.n_compartments, self.n_structs).T
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
        Damped Picard iteration (no Jacobian, memory-friendly).
        """
        x = np.zeros(self.n_structs * self.n_compartments)
        for j in range(self.n_compartments):
            x[j * self.n_structs : (j + 1) * self.n_structs] = self.initial_feed

        for it in range(max_iter):
            F = self.residual(x)
            norm = np.linalg.norm(F)
            if verbose and it % 8 == 0:
                print(f"Iter {it}: ||F|| = {norm:.6e}")
            if norm < tol:
                if verbose:
                    print(f"Converged at iteration {it}, ||F|| = {norm:.6e}")
                break
            x = x - damping * F
            x = np.maximum(x, 0)

        return x.reshape(self.n_compartments, self.n_structs).T

    # =====================================================================
    # Constrained Newton-Raphson with auxiliary Michaelis-Menten denominator
    # variables and a sparse direct linear solver.
    #
    # Reformulation (Krambeck & Betenbaugh, 2005):
    #
    # Original reaction rate for reaction r with enzyme e and substrate s:
    #     r = alpha_r * c[s] / (1 + sum_{r' in enzyme_e} c[s_r'] / Km_r')
    #
    # where alpha_r = kf * Et_e[j] * Donor / (Km_r * (Kmd_e + Donor)) is
    # parameter-only. We promote the competitor sum for each enzyme e to an
    # auxiliary variable d_e and add one equation per enzyme per compartment
    # to define it. The reaction rate becomes
    #     r = alpha_r * c[s] / d_e
    # so each r depends on only one c and one d, keeping the Jacobian sparse.
    #
    # Variables (one block of size S+E per compartment, S=n_structs, E=n_enzymes):
    #     [c_1, ..., c_S, d_1, ..., d_E]      per compartment
    #
    # Equations (same block layout):
    #     F_c(i, j) = c[i, j] - c[i, j-1] - tau_j * r_net_i(c_j, d_j) = 0
    #     F_d(e, j) = d[e, j] - compute_all_competitor_sums(c_j)[e]   = 0
    # =====================================================================




    def _setup_sparse_solver(self):
        """Precompute reaction-indexed constants and the Jacobian sparsity
        pattern for the extended (c, d) system. Runs once per model.

        The pattern is fixed (topology only), so we build a canonical
        duplicate-summed CSC structure up front; each Newton iteration then
        only refills the values array.
        """
        S, E, nc = self.n_structs, self.n_enzymes, self.n_compartments
        n = S + E
        R = len(self.network.reactions)
        n_total = n * nc

        # Per-reaction constants indexed by reaction id.
        rxn_sub = np.empty(R, dtype=np.int64)
        rxn_enz = np.empty(R, dtype=np.int64)
        rxn_inv_km = np.empty(R)
        for r in range(R):
            rxn = self.network.reactions[r]
            rxn_sub[r] = rxn.substrate_id
            rxn_enz[r] = self.enzyme_idx[rxn.enzyme.enzyme_name]
            rxn_inv_km[r] = 1.0 / rxn.km_value if rxn.km_value > 0 else 0.0
        self._rxn_sub = rxn_sub
        self._rxn_enz = rxn_enz
        self._rxn_inv_km = rxn_inv_km

        # alpha[r, j] = kf * Et[e, j] * Donor / (Km_r * (Kmd_e + Donor)); constant w.r.t. x.
        alpha = np.zeros((R, nc))
        for r in range(R):
            rxn = self.network.reactions[r]
            enz = rxn.enzyme
            donor = self.donor_conc.get(enz.cosubstrate, 1.0)
            denom = rxn.km_value * (enz.Kmd + donor)
            if denom <= 0:
                continue
            factor = enz.kf * donor / denom
            for j in range(nc):
                alpha[r, j] = factor * self.enzyme_dist[j].get(enz.enzyme_name, 0.0)
        self._alpha = alpha

        # Consume-incidence matrix for the op0 diagonal term of the Jacobian.
        cols_r = np.arange(R)
        self._C_mat = sp.csr_matrix((np.ones(R), (rxn_sub, cols_r)), shape=(S, R))

        # ---- local pattern within one compartment block ----
        pat_row, pat_col, pat_op, pat_r, pat_e = [], [], [], [], []

        # op0: dF_c(i)/dc(i) diagonal, one per structure
        for i in range(S):
            pat_row.append(i); pat_col.append(i); pat_op.append(0)
            pat_r.append(-1); pat_e.append(-1)
        # op1: dF_c(i)/dc(s_r), one per producing reaction
        for i in range(S):
            for rxn_id in sorted(self.network.structures[i]._producing_reactions):
                pat_row.append(i); pat_col.append(self.network.reactions[rxn_id].substrate_id)
                pat_op.append(1); pat_r.append(rxn_id); pat_e.append(-1)
        # op2: dF_c(i)/dd(e), one per (structure, enzyme) pair that shares a reaction
        op2_pairs = []
        for i in range(S):
            enz_set = set()
            for rxn_id in self.network.structures[i]._producing_reactions:
                enz_set.add(rxn_enz[rxn_id])
            for rxn_id in self.network.structures[i]._consuming_reactions:
                enz_set.add(rxn_enz[rxn_id])
            for e_idx in sorted(enz_set):
                pat_row.append(i); pat_col.append(S + e_idx)
                pat_op.append(2); pat_r.append(-1); pat_e.append(e_idx)
                op2_pairs.append((i, e_idx))
        # op3: dF_d(e)/dc(s_r), one per reaction
        for e_idx, e_name in enumerate(self.enzyme_names):
            for rxn_id in self.enzyme_rxn_ids[e_name]:
                pat_row.append(S + e_idx); pat_col.append(self.network.reactions[rxn_id].substrate_id)
                pat_op.append(3); pat_r.append(rxn_id); pat_e.append(-1)
        # op4: d-equation diagonal, one per enzyme
        for e_idx in range(E):
            pat_row.append(S + e_idx); pat_col.append(S + e_idx)
            pat_op.append(4); pat_r.append(-1); pat_e.append(-1)

        L_local = len(pat_row)
        self._L_local = L_local
        b0 = (0, S)
        b1 = (S, S + R)
        b2 = (S + R, S + R + len(op2_pairs))
        b3 = (S + R + len(op2_pairs), S + R + len(op2_pairs) + R)
        b4 = (S + R + len(op2_pairs) + R, L_local)
        self._seg = (b0, b1, b2, b3, b4)
        pat_r = np.array(pat_r, dtype=np.int64)
        pat_e = np.array(pat_e, dtype=np.int64)
        self._op1_r = pat_r[b1[0]:b1[1]]
        self._op3_r = pat_r[b3[0]:b3[1]]
        self._op2_e = pat_e[b2[0]:b2[1]]

        # PE[k, r] = +1 if r produces op2_pairs[k][0] via enzyme op2_pairs[k][1],
        #            -1 if r consumes it via the same enzyme. Lets us form the
        #            op2 numerator (sum of signed rates per (i, e) pair) as one matvec.
        op2_map = {pair: k for k, pair in enumerate(op2_pairs)}
        pe_rows, pe_cols, pe_vals = [], [], []
        for r in range(R):
            e_idx = rxn_enz[r]
            k = op2_map.get((self.network.reactions[r].product_id, e_idx))
            if k is not None:
                pe_rows.append(k); pe_cols.append(r); pe_vals.append(1.0)
            k = op2_map.get((rxn_sub[r], e_idx))
            if k is not None:
                pe_rows.append(k); pe_cols.append(r); pe_vals.append(-1.0)
        self._PE_mat = sp.csr_matrix(
            (pe_vals, (pe_rows, pe_cols)), shape=(len(op2_pairs), R)
        )

        # ---- global pattern: tile per-compartment blocks + feed-forward coupling ----
        pat_row = np.array(pat_row, dtype=np.int64)
        pat_col = np.array(pat_col, dtype=np.int64)
        tile = np.tile(np.arange(L_local), nc)
        comp_base = np.repeat(np.arange(nc) * n, L_local)
        rows_block = comp_base + pat_row[tile]
        cols_block = comp_base + pat_col[tile]
        if nc > 1:
            coup_j = np.repeat(np.arange(1, nc), S)
            coup_i = np.tile(np.arange(S), nc - 1)
            rows = np.concatenate([rows_block, coup_j * n + coup_i])
            cols = np.concatenate([cols_block, (coup_j - 1) * n + coup_i])
        else:
            rows, cols = rows_block, cols_block

        # Canonical (duplicates-summed) CSC structure + pattern->canonical map.
        key = cols.astype(np.int64) * n_total + rows
        _, first_idx, inv = np.unique(key, return_index=True, return_inverse=True)
        self._canon_map = inv
        self._n_canon = len(first_idx)
        coo = sp.coo_matrix(
            (np.ones(self._n_canon), (rows[first_idx], cols[first_idx])),
            shape=(n_total, n_total),
        )
        csc = coo.tocsc()
        self._csc_indptr = csc.indptr
        self._csc_indices = csc.indices

    def _extended_residual(self, x: np.ndarray) -> np.ndarray:
        """Residual of the extended system.

        F_c reuses :meth:`_calculate_net_rates` with the d variables injected
        as the competitor sums (its rates divide by whatever dict it is given),
        so the reaction-rate code has a single implementation. F_d reuses
        :meth:`KineticsCalculator.compute_all_competitor_sums` directly.
        """
        S, nc = self.n_structs, self.n_compartments
        n = S + self.n_enzymes
        F = np.empty(n * nc)
        for j in range(nc):
            c_j = x[j * n : j * n + S]
            d_j = x[j * n + S : (j + 1) * n]
            c_prev = x[(j - 1) * n : (j - 1) * n + S] if j > 0 else self.initial_feed
            d_as_cs = {name: d_j[idx] for name, idx in self.enzyme_idx.items()}
            r_net = self._calculate_net_rates(c_j, d_as_cs, self.enzyme_dist[j], j)
            F[j * n : j * n + S] = c_j - c_prev - self.tau[j] * r_net
            cs = self.kinetics.compute_all_competitor_sums(c_j)
            F[j * n + S : (j + 1) * n] = [
                d_j[idx] - cs.get(name, 1.0) for name, idx in self.enzyme_idx.items()
            ]
        return F

    def _assemble_jacobian(self, x: np.ndarray) -> sp.csc_matrix:
        """Assemble the analytical Jacobian of the extended residual at x.

        Rate r = alpha_r * c[s_r] / d[e_r], so
            dr/dc[s_r] = alpha_r / d[e_r]
            dr/dd[e_r] = -r / d[e_r]
        and each F_c(i) touches only c(i), the substrates of reactions producing
        i, and the d(e) of enzymes acting on i; F_d(e) touches c(s_r) of its own
        reactions and d(e). That keeps the Jacobian sparse.
        """
        S, nc = self.n_structs, self.n_compartments
        n = S + self.n_enzymes
        vals = np.empty(nc * self._L_local + max(0, nc - 1) * S)
        for j in range(nc):
            c_j = x[j * n : j * n + S]
            d_j = x[j * n + S : (j + 1) * n]
            tau_j = self.tau[j]
            a_over_d = self._alpha[:, j] / d_j[self._rxn_enz]
            r_vec = a_over_d * c_j[self._rxn_sub]
            op2_numer = self._PE_mat @ r_vec
            consume = self._C_mat @ a_over_d
            v = np.empty(self._L_local)
            v[self._seg[0][0]:self._seg[0][1]] = 1.0 + tau_j * consume
            v[self._seg[1][0]:self._seg[1][1]] = -tau_j * a_over_d[self._op1_r]
            v[self._seg[2][0]:self._seg[2][1]] = tau_j * op2_numer / d_j[self._op2_e]
            v[self._seg[3][0]:self._seg[3][1]] = -self._rxn_inv_km[self._op3_r]
            v[self._seg[4][0]:self._seg[4][1]] = 1.0
            vals[j * self._L_local : (j + 1) * self._L_local] = v
        if nc > 1:
            vals[nc * self._L_local:] = -1.0
        vals_c = np.bincount(self._canon_map, weights=vals, minlength=self._n_canon)
        return sp.csc_matrix(
            (vals_c, self._csc_indices, self._csc_indptr),
            shape=(n * nc, n * nc),
        )

    def solve_sparse(
        self,
        tol: float = 1e-6,
        max_iter: int = 50,
        damping: float = 1.0,
        verbose: bool = True,
    ) -> np.ndarray:
        """Constrained Newton-Raphson on the extended (c, d) system.

        Modern replacement for the MA28-based solver of Krambeck & Betenbaugh
        (2005): SuperLU via :func:`scipy.sparse.linalg.spsolve` handles the
        sparse direct factorization, the analytical Jacobian's sparsity pattern
        is frozen and pre-assembled in CSC form, and non-negativity of c and
        the d >= 1 lower bound are enforced by a safeguarded step plus
        backtracking on ||F||_inf. Solves all compartments as one system, as
        the reference found sequential compartment solves unnecessary.
        """
        if not hasattr(self, "_csc_indptr"):
            self._setup_sparse_solver()
        S, nc = self.n_structs, self.n_compartments
        n = S + self.n_enzymes

        cs0 = self.kinetics.compute_all_competitor_sums(self.initial_feed)
        d0 = np.array([cs0.get(name, 1.0) for name in self.enzyme_names])
        x = np.empty(n * nc)
        for j in range(nc):
            x[j * n : j * n + S] = self.initial_feed
            x[j * n + S : (j + 1) * n] = d0

        for it in range(max_iter):
            F = self._extended_residual(x)
            norm = np.linalg.norm(F, ord=np.inf)
            if verbose:
                print(f"[solve_sparse] iter {it}: ||F||_inf = {norm:.4e}")
            if norm < tol:
                break

            J = self._assemble_jacobian(x)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                dx = spla.spsolve(J, -F)
            if not np.all(np.isfinite(dx)):
                if verbose:
                    print("[solve_sparse] singular Jacobian; stopping")
                break

            # Box-constrained Newton: take the full step, clip each variable
            # to its bound, then backtrack on the CLIPPED step. Clipping
            # per-variable (instead of shrinking a global scalar step) keeps
            # one nearly-active variable from freezing the entire update.
            x_trial = x + dx
            for j in range(nc):
                np.maximum(x_trial[j * n : j * n + S], 0.0, out=x_trial[j * n : j * n + S])
                np.maximum(x_trial[j * n + S : (j + 1) * n], 1.0, out=x_trial[j * n + S : (j + 1) * n])
            dx_clip = x_trial - x  # zeroed-out components sit at a bound

            step = damping
            x_new = x_trial
            for _ in range(25):
                if np.linalg.norm(self._extended_residual(x_new), ord=np.inf) < norm:
                    break
                step *= 0.5
                x_new = x + step * dx_clip
                for j in range(nc):
                    np.maximum(x_new[j * n : j * n + S], 0.0, out=x_new[j * n : j * n + S])
                    np.maximum(x_new[j * n + S : (j + 1) * n], 1.0, out=x_new[j * n + S : (j + 1) * n])
            else:
                # No acceptable step: stop rather than drift uphill.
                if verbose:
                    print("[solve_sparse] line search failed; stopping")
                break
            x = x_new

        c = np.empty((S, nc))
        for j in range(nc):
            c[:, j] = x[j * n : j * n + S]
        return np.maximum(c, 0)
