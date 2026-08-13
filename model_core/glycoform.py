from collections import namedtuple


class Glycoform(namedtuple("Glycoform", ["man", "fuc", "gnb", "br1", "br2", "br3", "br4", "gal", "sia"])):
    """
    Represents an N-glycan structure with a 9-digit encoding.
    Example: (3, 1, 1, 6, 5, 3, 1, 5, 2)
    """

    def __init__(self, *args, **kwargs):
        self._consuming_reactions = set()
        self._producing_reactions = set()

    def __new__(cls, man, fuc, gnb, br1, br2, br3, br4, gal, sia):
        if not all(isinstance(x, int) and x >= 0 for x in [man, fuc, gnb, br1, br2, br3, br4, gal, sia]):
            raise ValueError("All parameters must be non-negative integers.")
        obj = super().__new__(cls, man, fuc, gnb, br1, br2, br3, br4, gal, sia)
        obj.available = 0 if any(b > 6 for b in (br1, br2, br3, br4)) else 1
        return obj
