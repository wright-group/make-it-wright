import makeitwright as mw
from makeitwright import datasets

andor = mw.andor
roi = mw.helpers.roi
parse = mw.parsers.parse
plot = mw.spectra.plot_spectra


def test_import_andor():
    p = datasets.PL  
    filepath = p.parent
    filename = p.stem

    data1 = parse(filepath, objective="10", keywords=filename + ".asc")
    data2 = mw.parsers.fromAndorNeo(p)
