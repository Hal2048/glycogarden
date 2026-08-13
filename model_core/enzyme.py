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


Enzyme_RULES = [
    {
        "name": "ManI",
        "condition": "g.man > 5",
        "product_struc": "Glycoform(g.man-1, g.fuc, g.gnb, g.br1, g.br2, g.br3, g.br4, g.gal, g.sia)",
        "cosubstrate": None,
        "Km": 100, "Kmd": 0, "kf": 1924,
        "adjustments": []
    },
    {
        "name": "ManII",
        "condition": "g.man > 3 and g.br4 == 1 and g.gnb == 0",
        "product_struc": "Glycoform(g.man-1, g.fuc, g.gnb, g.br1, g.br2, g.br3, g.br4, g.gal, g.sia)",
        "cosubstrate": None,
        "Km": 200, "Kmd": 0, "kf": 1924,
        "adjustments": [("g.man == 4", 0.50)]
    },
    {
        "name": "FucT",
        "condition": "g.fuc == 0 and g.br4 > 0 and g.gnb == 0 and g.gal == 0",
        "product_struc": "Glycoform(g.man, g.fuc+1, g.gnb, g.br1, g.br2, g.br3, g.br4, g.gal, g.sia)",
        "cosubstrate": "GDP-Fuc",
        "Km": 25, "Kmd": 46, "kf": 253,
        "adjustments": [("g.br2 > 0", 0.048)]
    },
    {
        "name": "GnTI",
        "condition": "g.br4 == 0 and g.man == 5",
        "product_struc": "Glycoform(g.man, g.fuc, g.gnb, g.br1, g.br2, g.br3, g.br4+1, g.gal, g.sia)",
        "cosubstrate": "UDP-GlcNAc",
        "Km": 260, "Kmd": 170, "kf": 990,
        "adjustments": [("g.br2 == 0", 5.0)]
    },
    {
        "name": "GnTII",
        "condition": "g.br2 == 0 and g.man == 3 and g.br4 == 1 and g.gnb == 0",
        "product_struc": "Glycoform(g.man, g.fuc, g.gnb, g.br1, g.br2+1, g.br3, g.br4, g.gal, g.sia)",
        "cosubstrate": "UDP-GlcNAc",
        "Km": 190, "Kmd": 3100, "kf": 607,
        "adjustments": []
    },
    {
        "name": "GnTIII",
        "condition": "g.gnb == 0 and g.br4 > 0 and g.gal == 0",
        "product_struc": "Glycoform(g.man, g.fuc, g.gnb+1, g.br1, g.br2, g.br3, g.br4, g.gal, g.sia)",
        "cosubstrate": "UDP-GlcNAc",
        "Km": 4000, "Kmd": 3100, "kf": 607,
        "adjustments": [("g.br2 > 0", 0.048)]
    },
    {
        "name": "GnTIV",
        "condition": "g.br3 == 0 and g.br4 == 1 and g.gnb == 0",
        "product_struc": "Glycoform(g.man, g.fuc, g.gnb, g.br1, g.br2, g.br3+1, g.br4, g.gal, g.sia)",
        "cosubstrate": "UDP-GlcNAc",
        "Km": 3400, "Kmd": 8300, "kf": 187,
        "adjustments": [("g.br2 == 0", 5), ("g.br2 > 1 or g.br1 > 1", 1.5), ("g.br1 > 0", 0.178)]
    },
    {
        "name": "GnTV",
        "condition": "g.br1 == 0 and g.br2 == 1 and g.gnb == 0",
        "product_struc": "Glycoform(g.man, g.fuc, g.gnb, g.br1+1, g.br2, g.br3, g.br4, g.gal, g.sia)",
        "cosubstrate": "UDP-GlcNAc",
        "Km": 130, "Kmd": 3500, "kf": 1410,
        "adjustments": [("g.br3 > 0", 0.692)]
    },
    {
        "name": "GnTE_Br1",
        "condition": "g.br1 == 2",
        "product_struc": "Glycoform(g.man, g.fuc, g.gnb, g.br1+2, g.br2, g.br3, g.br4, g.gal, g.sia)",
        "cosubstrate": "UDP-GlcNAc",
        "Km": 700, "Kmd": 55, "kf": 25,
        "adjustments": [("g.br2 == 0", 4)]
    },
    {
        "name": "GnTE_Br2",
        "condition": "g.br2 == 2",
        "product_struc": "Glycoform(g.man, g.fuc, g.gnb, g.br1, g.br2+2, g.br3, g.br4, g.gal, g.sia)",
        "cosubstrate": "UDP-GlcNAc",
        "Km": 700, "Kmd": 55, "kf": 25,
        "adjustments": [("g.br1 == 0", 4)]
    },
    {
        "name": "GnTE_Br3",
        "condition": "g.br3 == 2",
        "product_struc": "Glycoform(g.man, g.fuc, g.gnb, g.br1, g.br2, g.br3+2, g.br4, g.gal, g.sia)",
        "cosubstrate": "UDP-GlcNAc",
        "Km": 7000, "Kmd": 55, "kf": 25,
        "adjustments": [("g.br4 == 0", 4)]
    },
    {
        "name": "GnTE_Br4",
        "condition": "g.br4 == 2",
        "product_struc": "Glycoform(g.man, g.fuc, g.gnb, g.br1, g.br2, g.br3, g.br4+2, g.gal, g.sia)",
        "cosubstrate": "UDP-GlcNAc",
        "Km": 7000, "Kmd": 55, "kf": 25,
        "adjustments": [("g.br3 == 0", 4)]
    },
    {
        "name": "GalT_Br1",
        "condition": "g.br1 == 1 or g.br1 == 4",
        "product_struc": "Glycoform(g.man, g.fuc, g.gnb, g.br1+1, g.br2, g.br3, g.br4, g.gal+1, g.sia)",
        "cosubstrate": "UDP-Gal",
        "Km": 150, "Kmd": 0, "kf": 8712,
        "adjustments": [("g.gnb > 0 and g.br2 > 0", 3.62)]
    },
    {
        "name": "GalT_Br2",
        "condition": "g.br2 == 1 or g.br2 == 4",
        "product_struc": "Glycoform(g.man, g.fuc, g.gnb, g.br1, g.br2+1, g.br3, g.br4, g.gal+1, g.sia)",
        "cosubstrate": "UDP-Gal",
        "Km": 135, "Kmd": 0, "kf": 8712,
        "adjustments": [("g.gnb > 0 and g.br2 > 0", 3.62)]
    },
    {
        "name": "GalT_Br3",
        "condition": "g.br3 == 1 or g.br3 == 4",
        "product_struc": "Glycoform(g.man, g.fuc, g.gnb, g.br1, g.br2, g.br3+1, g.br4, g.gal+1, g.sia)",
        "cosubstrate": "UDP-Gal",
        "Km": 80, "Kmd": 0, "kf": 8712,
        "adjustments": [("g.gnb > 0 and g.br2 > 0", 3.62)]
    },
    {
        "name": "GalT_Br4",
        "condition": "g.br4 == 1 or g.br4 == 4",
        "product_struc": "Glycoform(g.man, g.fuc, g.gnb, g.br1, g.br2, g.br3, g.br4+1, g.gal+1, g.sia)",
        "cosubstrate": "UDP-Gal",
        "Km": 4000, "Kmd": 0, "kf": 8712,
        "adjustments": [("g.gnb > 0 and g.br2 > 0", 3.62)]
    },
    {
        "name": "SiaT_Br1",
        "condition": "g.br1 == 2 or g.br1 == 5",
        "product_struc": "Glycoform(g.man, g.fuc, g.gnb, g.br1+1, g.br2, g.br3, g.br4, g.gal, g.sia+1)",
        "cosubstrate": "CMP-NeuAc",
        "Km": 260, "Kmd": 57, "kf": 484,
        "adjustments": [("g.sia > 1", 5.0)]
    },
    {
        "name": "SiaT_Br2",
        "condition": "g.br2 == 2 or g.br2 == 5",
        "product_struc": "Glycoform(g.man, g.fuc, g.gnb, g.br1, g.br2+1, g.br3, g.br4, g.gal, g.sia+1)",
        "cosubstrate": "CMP-NeuAc",
        "Km": 260, "Kmd": 57, "kf": 484,
        "adjustments": [("g.sia > 1", 5.0)]
    },
    {
        "name": "SiaT_Br3",
        "condition": "g.br3 == 2 or g.br3 == 5",
        "product_struc": "Glycoform(g.man, g.fuc, g.gnb, g.br1, g.br2, g.br3+1, g.br4, g.gal, g.sia+1)",
        "cosubstrate": "CMP-NeuAc",
        "Km": 260, "Kmd": 57, "kf": 484,
        "adjustments": [("g.sia > 1", 5.0)]
    },
    {
        "name": "SiaT_Br4",
        "condition": "g.br4 == 2 or g.br4 == 5",
        "product_struc": "Glycoform(g.man, g.fuc, g.gnb, g.br1, g.br2, g.br3, g.br4+1, g.gal, g.sia+1)",
        "cosubstrate": "CMP-NeuAc",
        "Km": 260, "Kmd": 57, "kf": 484,
        "adjustments": [("g.sia > 1", 5.0)]
    },
]


def build_enzymes() -> List[Enzyme]:
    """Builds a list of Enzyme objects from Enzyme_RULES."""
    enzymes = []
    for rule in Enzyme_RULES:
        enzyme = Enzyme(
            enzyme_name=rule["name"],
            condition=eval("lambda g: " + rule["condition"]),
            product_struc=eval("lambda g: " + rule["product_struc"]),
            cosubstrate=rule["cosubstrate"],
            Km=rule["Km"],
            Kmd=rule["Kmd"],
            kf=rule["kf"],
            Km_adjust=lambda g, rule=rule: next(
                (adj[1] for adj in rule.get("adjustments", [])
                 if eval("lambda g: " + adj[0])(g)), 1.0
            )
        )
        enzymes.append(enzyme)
    return enzymes
