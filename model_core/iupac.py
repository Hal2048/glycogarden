"""Convert the 12-field Krambeck encoding to an IUPAC-style glycan name."""

from glycoform import Glycoform


GLYCOFORM_FIELDS = (
    "man1", "man2", "man3", "fuc", "gnb",
    "br1", "br2", "br3", "br4", "gal", "sia", "br_o",
)


def _coerce_code(glycoform_or_tuple):
    """
    Accept either a Glycoform (namedtuple with 12 fields) or a 12-tuple and
    return the canonical (man1, man2, man3, fuc, gnb, br1, br2, br3, br4,
    gal, sia, br_o) tuple.
    """
    if hasattr(glycoform_or_tuple, "_fields") and len(glycoform_or_tuple._fields) == 12:
        return tuple(getattr(glycoform_or_tuple, name) for name in GLYCOFORM_FIELDS)
    tup = tuple(glycoform_or_tuple)
    if len(tup) != 12:
        raise ValueError(
            f"Expected 12-tuple Glycoform, got {len(tup)} values: {tup}"
        )
    return tup


def krambeck_to_iupac(glycan_form: Glycoform):
    """
    将12位Krambeck编码转换为IUPAC表示法
    """

    glycan_code = _coerce_code(glycan_form)
    man1, man2, man3, fuc, gnb, br1, br2, br3, br4, gal, sia, br_o = _coerce_code(glycan_code)
    
    frag = {
        0:                                          "",
        1:                                   "(GlcNAc",
        2:                            "(Galb1-4GlcNAc",
        3:                  "(Neu5Aca2-3Galb1-4GlcNAc",
        4:                  "(GlcNAcb1-3Galb1-4GlcNAc",
        5:           "(Galb1-4GlcNAcb1-3Galb1-4GlcNAc",
        6: "(Neu5Aca2-3Galb1-4GlcNAcb1-3Galb1-4GlcNAc",
    }

    def expand(br, link):
        return f"{frag[br]}{link}" + ")" if br > 0 else ""

    # Upper arm: Br1 (b1-6) + Br2 (b1-2) -> Man(a1-6)
    if br1 > 0 or br2 > 0:
        upper = expand(br1, "b1-6") + expand(br2, "b1-2") + "Mana1-6"
    else:
        # man1/man2 表示 α1-6 主臂 / α1-3 分支上的 Man 数量
        main = "Mana1-2" * (man1 - 1) + "Mana1-6" if man1 > 0 else ""
        branch = ""
        if man2 == 2:
            branch = "(Mana1-2Mana1-3)" 
        elif man2 == 1:
            branch = "(Mana1-3)"
        upper = main + branch + "Mana1-6"

    lower = "Mana1-3"
    # Lower arm: Br3 (b1-4) + Br4 (b1-2) -> Man(a1-3)
    if br_o > 0:
        lower = "(" + "Mana1-6" * br_o + ")" + lower
    
    if br3 > 0 or br4 > 0:
        lower = expand(br3, "b1-4") + expand(br4, "b1-2")+ lower
    else:
        lower = "Mana1-2" * man3 + lower

    # Core
    core = "Manb1-4GlcNAcb1-4"
    if fuc:
        core += "(Fuca1-6)"
    core += "GlcNAcb1"

    # Bisecting GlcNAc hangs off the central Man
    if gnb:
        core = "(GlcNAcb1-4)" + core


    return upper + (f"({lower})" if lower else "") + core