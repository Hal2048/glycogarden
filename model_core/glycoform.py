from collections import namedtuple


class Glycoform(namedtuple("Glycoform", ["man1", "man2", "man3", "fuc", "gnb", "br1", "br2", "br3", "br4", "gal", "sia", "br_o"])):
    """
    Represents an N-glycan structure with a 12-digit encoding.

    """

    def __init__(self, *args, **kwargs):
        self._consuming_reactions = set()
        self._producing_reactions = set()

    def __new__(cls, man1, man2, man3, fuc, gnb, br1, br2, br3, br4, gal, sia, br_o=0):
        if not all(isinstance(x, int) and x >= 0 for x in [man1, man2, man3, fuc, gnb, br1, br2, br3, br4, gal, sia, br_o]):
            raise ValueError("All parameters must be non-negative integers.")
        obj = super().__new__(cls, man1, man2, man3, fuc, gnb, br1, br2, br3, br4, gal, sia, br_o)
        obj.available = 0 if any(b > 6 for b in (br1, br2, br3, br4)) else 1
        return obj
