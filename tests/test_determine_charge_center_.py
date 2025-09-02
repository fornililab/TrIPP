from tripp._determine_charge_center_ import determine_charge_center
import MDAnalysis as mda
from io import StringIO
import pytest
import numpy as np 
class TestDetermineChargeCenter:
    def pdb1():
        pdb = """
ATOM      1  N   ASP     1      44.260  49.990  17.910  1.00  0.00           N
ATOM      2  C   ASP     2      44.260  49.990  17.910  1.00  0.00           N
ATOM      3  O   LYS     3      43.080  50.770  17.580  1.00  0.00           O
"""
        selection = 'resid 2'
        return pdb, selection
    def pdb2():
        pdb = """
ATOM      1  OD1 ASPG    1      44.260  49.990  17.910  1.00  0.00           O
ATOM      2  OD2 ASPG    1      44.260  49.990  17.910  1.00  0.00           O
"""
        selection = 'resid 1'
        return pdb, selection
    @pytest.mark.parametrize("pdb, selection", [pdb1(), pdb2()])
    def test_determine_charge_center_fail(self, pdb, selection):
        pstream = mda.lib.util.NamedStream(StringIO(pdb),'tmp.pdb')
        universe = mda.Universe(pstream, format='PDB')
        pstream.reset()
        with pytest.raises(Exception):
            determine_charge_center(universe, selection)
        
        