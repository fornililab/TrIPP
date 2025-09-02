from tripp._create_mda_universe_ import create_mda_universe, create_propka_compatible_universe
from io import StringIO
import MDAnalysis as mda
import os
import pytest
import glob

def get_data_path(filename):
    return os.path.join(os.path.dirname(__file__), "data", filename)
def check_dir(path):
    if not os.path.isdir(path):
        # Make directory if not present
        os.makedirs(path) 
    else:
        # Remove files named .temp* in the directory before proceeding.
        [os.remove(file) for file in glob.glob(f'{path}/.temp*')]
    return path

topology_file = get_data_path('lyso.pdb')
trajectory_file = get_data_path('lyso_test.xtc')

class Test_create_mda_universe:
    def test_create_mda_universe(self):
        u = create_mda_universe(topology_file=topology_file,
                                trajectory_file=trajectory_file)
        assert u.atoms.n_atoms == 1960
        assert len(u.trajectory) == 30
        
    def test_resname_correction(self):
        
        pdb = """
ATOM      1  N   ASPH    1      44.260  49.990  17.910  1.00  0.00           N
ATOM      2  O   LYS     2      43.080  50.770  17.580  1.00  0.00           O
ATOM      3  OXT LYS     2      43.080  50.770  17.580  1.00  0.00           O
"""
        universe = mda.Universe(StringIO(pdb), format='PDB')
        corrected_universe = create_propka_compatible_universe(universe,
                                                               hetatm_resname=None,
                                                               custom_terminal_oxygens=None,
                                                               custom_resname_correction=None)
        assert corrected_universe.atoms.residues.resnames[0] == 'ASP'
        
        pdb = """
ATOM      1  N   ABCD    1      44.260  49.990  17.910  1.00  0.00           N
ATOM      2  O   LYS     2      43.080  50.770  17.580  1.00  0.00           O
ATOM      3  OXT LYS     2      43.080  50.770  17.580  1.00  0.00           O
"""
        pstream = mda.lib.util.NamedStream(StringIO(pdb),'tmp.pdb')
        universe = mda.Universe(pstream, format='PDB')
        pstream.reset()
        with pytest.raises(NameError):
            incorrect_universe = create_propka_compatible_universe(universe,
                                                               hetatm_resname=None,
                                                               custom_terminal_oxygens=None,
                                                               custom_resname_correction=None)

        corrected_universe = create_propka_compatible_universe(universe,
                                                        hetatm_resname=None,
                                                        custom_terminal_oxygens=None,
                                                        custom_resname_correction={'ABCD':'ASP'})
        assert corrected_universe.atoms.residues.resnames[0] == 'ASP'
    
    def test_hetatm_correction(self):
        pdb = """
ATOM      1  N   ABCD    1      44.260  49.990  17.910  1.00  0.00           N
ATOM      2  O   LYS     2      43.080  50.770  17.580  1.00  0.00           O
ATOM      3  OXT LYS     2      43.080  50.770  17.580  1.00  0.00           O
"""
        pstream = mda.lib.util.NamedStream(StringIO(pdb),'tmp.pdb')
        universe = mda.Universe(pstream, format='PDB')
        pstream.reset()
        corrected_universe = create_propka_compatible_universe(universe,
                                                        hetatm_resname='ABCD',
                                                        custom_terminal_oxygens=None,
                                                        custom_resname_correction=None)
        assert corrected_universe.atoms.residues.record_types[0] == 'HETATM'
        
        pdb = """
ATOM      1  N   ABCD    1      44.260  49.990  17.910  1.00  0.00           N
ATOM      2  N   DEF     2      44.260  49.990  17.910  1.00  0.00           N
ATOM      3  O   LYS     3      43.080  50.770  17.580  1.00  0.00           O
ATOM      4  OXT LYS     3      43.080  50.770  17.580  1.00  0.00           O
"""
        pstream = mda.lib.util.NamedStream(StringIO(pdb),'tmp.pdb')
        universe = mda.Universe(pstream, format='PDB')
        pstream.reset()
        corrected_universe = create_propka_compatible_universe(universe,
                                                        hetatm_resname=['ABCD','DEF'],
                                                        custom_terminal_oxygens=None,
                                                        custom_resname_correction=None)
        assert corrected_universe.atoms.residues.record_types[0] == 'HETATM'
        assert corrected_universe.atoms.residues.record_types[1] == 'HETATM'
        
        topology_file = get_data_path('myosin_elc_test.pdb')
        universe = mda.Universe(topology_file)
        with pytest.raises(NameError):
            incorrect_universe = create_propka_compatible_universe(universe,
                                                                   hetatm_resname=None,
                                                                   custom_terminal_oxygens=None,
                                                                   custom_resname_correction=None)
        corrected_universe = create_propka_compatible_universe(universe,
                                                                hetatm_resname=['ADP','PI2','MG'],
                                                                custom_terminal_oxygens=None,
                                                                custom_resname_correction=None)
        assert (corrected_universe.select_atoms('resname ADP PI2 MG').atoms.record_types == 'HETATM').all()
    
    def test_custom_terminal_oxygens(self):
        pdb = """
ATOM      1  N   ASP     1      44.260  49.990  17.910  1.00  0.00           N
ATOM      2  N   ASP     2      44.260  49.990  17.910  1.00  0.00           N
ATOM      3  OB1 LYS     3      43.080  50.770  17.580  1.00  0.00           O
ATOM      4  OB2 LYS     3      43.080  50.770  17.580  1.00  0.00           O
"""
        pstream = mda.lib.util.NamedStream(StringIO(pdb),'tmp.pdb')
        universe = mda.Universe(pstream, format='PDB')
        pstream.reset()
        with pytest.raises(NameError):
            incorrect_universe = create_propka_compatible_universe(universe,
                                                               hetatm_resname=None,
                                                               custom_terminal_oxygens=None,
                                                               custom_resname_correction=None)
        corrected_universe = create_propka_compatible_universe(universe,
                                                        hetatm_resname=None,
                                                        custom_terminal_oxygens=['OB1','OB2'],
                                                        custom_resname_correction=None)
    def test_myosin_elc_fail(self):
        topology_file = get_data_path('myosin_elc_test_nochB.pdb')
        trajectory_file = get_data_path('myosin_elc_test.xtc')
        with pytest.raises(Exception):
            universe = create_mda_universe(topology_file=topology_file,
                                trajectory_file=trajectory_file)
