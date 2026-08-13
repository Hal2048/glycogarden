import numpy as np
from typing import Dict, List
from collections import defaultdict

from reaction import ReactionNetwork


class KineticsCalculator:
    """
    Computes the competitive sum in enzyme reaction kinetics: 1 + sum(c[j]/Kmj).
    """
    def __init__(self, reaction_network: ReactionNetwork):
        self.network = reaction_network
        self._enzyme_reactions: Dict[str, List[int]] = defaultdict(list)
        for rxn_id, rxn in self.network.reactions.items():
            self._enzyme_reactions[rxn.enzyme.enzyme_name].append(rxn_id)

    def compute_all_competitor_sums(self, concentrations: np.ndarray) -> Dict[str, float]:
        """
        Computes the competitive sum 1 + sum(c[j]/Kmj) for all enzymes in the current compartment.

        Parameters
        ----------
        concentrations : np.ndarray
            Concentrations of glycans in the current compartment, indexed by glycan ID.

        Returns
        -------
        Dict[str, float]
            Enzyme name -> competitive sum.
        """
        result = {}
        for enzyme_name, rxn_ids in self._enzyme_reactions.items():
            competitor_sum = 1.0
            for rxn_id in rxn_ids:
                rxn = self.network.reactions[rxn_id]
                sub_id = rxn.substrate_id
                km = rxn.km_value
                if 0 <= sub_id < len(concentrations) and km > 0:
                    competitor_sum += concentrations[sub_id] / km
            result[enzyme_name] = competitor_sum
        return result
