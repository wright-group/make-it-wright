import makeitwright as mw
from makeitwright import datasets

andor = mw.andor
parse = mw.parsers.parse


def test_import_andor():
    """smokescreen to see if importing fails"""
    p = datasets.PL  
    filepath = p.parent
    filename = p.stem

    data1 = parse(filepath, objective="10", keywords=filename + ".asc")
    data2 = mw.parsers.fromAndorNeo(p)
    assert data1.variable_names == data2.variable_names == ("wl", "y")


if __name__ == "__main__":
    test_import_andor()
