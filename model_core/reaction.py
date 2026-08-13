from typing import Dict, List

from glycoform import Glycoform
from enzyme import Enzyme


class Reaction:
    """
    Stores information for a single reaction.
    """
    __slots__ = ('id', 'enzyme', 'substrate_id', 'product_id', 'km_value')

    def __init__(self, rxn_id: int, enzyme: Enzyme, sub_id: int, prod_id: int, km: float):
        self.id = rxn_id
        self.enzyme = enzyme
        self.substrate_id = sub_id
        self.product_id = prod_id
        self.km_value = km


class ReactionNetwork:
    """
    Stores all glycan structures and reactions, and auto-generates the network based on rules.
    """
    def __init__(self):
        self.structures: Dict[int, Glycoform] = {}
        self.reactions: Dict[int, Reaction] = {}
        self._structure_to_id: Dict[Glycoform, int] = {}
        self._next_id = 0
        self._next_reaction_id = 0

    def add_structure(self, glycan: Glycoform) -> int:
        """Adds a new structure and returns its unique ID; returns existing ID if already present."""
        if glycan in self._structure_to_id:
            return self._structure_to_id[glycan]
        struct_id = self._next_id
        self._next_id += 1
        self.structures[struct_id] = glycan
        self._structure_to_id[glycan] = struct_id
        return struct_id

    def generate_network(self, initial_structures: List[Glycoform], enzymes: List[Enzyme]):
        """Iteratively generates the entire network from initial structures by applying all rules."""
        queue = [self.add_structure(s) for s in initial_structures]

        while queue:
            current_id = queue.pop(0)
            current_glycan = self.structures[current_id]
            if current_glycan.available == 0:
                continue
            for enzyme in enzymes:
                if enzyme.condition(current_glycan):
                    product_glycan = enzyme.product_struc(current_glycan)
                    product_id = self.add_structure(product_glycan)

                    km_val = enzyme.Km * enzyme.Km_adjust(current_glycan)

                    rxn_id = self._next_reaction_id
                    self._next_reaction_id += 1

                    reaction = Reaction(
                        rxn_id=rxn_id,
                        enzyme=enzyme,
                        sub_id=current_id,
                        prod_id=product_id,
                        km=km_val
                    )
                    self.reactions[rxn_id] = reaction

                    current_glycan._consuming_reactions.add(rxn_id)
                    product_glycan._producing_reactions.add(rxn_id)
                    enzyme.substrate_km_map[current_id] = km_val

                    if product_id not in queue:
                        queue.append(product_id)
