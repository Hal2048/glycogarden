import re
from typing import Callable, Dict, List

from glycoform import Glycoform


class Enzyme:
    """
    Represents an enzymatic reaction rule and its kinetic parameters.
    """
    def __init__(
        self,
        enzyme_name: str,
        condition: Callable[[Glycoform], bool],
        product_struc: Callable[[Glycoform], Glycoform],
        cosubstrate: str,
        Km: float,
        Kmd: float,
        kf: float,
        Km_adjust: Callable[[Glycoform], float] = None,
    ):
        self.enzyme_name = enzyme_name
        self.condition = condition
        self.product_struc = product_struc
        self.cosubstrate = cosubstrate
        self.Km_adjust = Km_adjust
        self.Km = Km
        self.Kmd = Kmd
        self.kf = kf
        self.substrate_km_map: Dict[int, float] = {}


def build_enzymes(enzyme_rules) -> List[Enzyme]:
    """Builds a list of Enzyme objects from ENZYME_RULES.

    Each rule's ``name`` may encode the specific reaction (e.g. ``ManI_M9``);
    this is mapped to the matching ``ENZYME_CONC`` key (e.g. ``ManI_9``) so
    that compartment concentrations can be looked up correctly.
    """
    enzymes = []
    for rule in enzyme_rules:
        conc_key = re.sub(r"_M(\d+)$", r"_\1", rule["name"])
        enzyme = Enzyme(
            enzyme_name=conc_key,
            condition=eval("lambda g: " + rule["condition"]),
            product_struc=eval("lambda g: " + rule["product_struc"]),
            cosubstrate=rule["cosubstrate"],
            Km=rule["Km"],
            Kmd=rule["Kmd"],
            kf=rule["kf"],
            Km_adjust=lambda g, r=rule: next(
                (adj[1] for adj in r.get("adjustments", [])
                 if eval("lambda g: " + adj[0])(g)), 1.0
            )
        )
        enzymes.append(enzyme)
    return enzymes
