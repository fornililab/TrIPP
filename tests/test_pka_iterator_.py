from tripp._create_mda_universe_ import create_mda_universe, create_propka_compatible_universe
from tripp._pka_iterator_ import pka_iterator
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
u = create_mda_universe(topology_file=topology_file,
                        trajectory_file=trajectory_file)
class Test_pka_iterator:
    def parm1():
        hetatm_resname = None
        custom_terminal_oxygens = None
        custom_resname_correction = None
        output_directory = check_dir('test_output/pka_iterator')
        trajectory_slice=[0,10]
        corrected_universe = create_propka_compatible_universe(u,
                                                               hetatm_resname=hetatm_resname,
                                                               custom_terminal_oxygens=custom_terminal_oxygens,
                                                               custom_resname_correction=custom_resname_correction)
        mutation_selections = None
        optargs = []
        return (trajectory_slice, corrected_universe, output_directory, mutation_selections, optargs)
    def parm2():
        hetatm_resname = None
        custom_terminal_oxygens = None
        custom_resname_correction = None
        output_directory = check_dir('test_output/pka_iterator')
        trajectory_slice=[0,10]
        corrected_universe = create_propka_compatible_universe(u,
                                                               hetatm_resname=hetatm_resname,
                                                               custom_terminal_oxygens=custom_terminal_oxygens,
                                                               custom_resname_correction=custom_resname_correction)
        mutation_selections = 'chainID A and resid 35 52'
        optargs = []
        return (trajectory_slice, corrected_universe, output_directory, mutation_selections, optargs)
    def parm3():
        hetatm_resname = None
        custom_terminal_oxygens = None
        custom_resname_correction = None
        output_directory = check_dir('test_output/pka_iterator')
        trajectory_slice=[0,10]
        corrected_universe = create_propka_compatible_universe(u,
                                                               hetatm_resname=hetatm_resname,
                                                               custom_terminal_oxygens=custom_terminal_oxygens,
                                                               custom_resname_correction=custom_resname_correction)
        mutation_selections = None
        optargs = ['-d']
        return (trajectory_slice, corrected_universe, output_directory, mutation_selections, optargs)
    @pytest.mark.parametrize("trajectory_slice, corrected_universe, output_directory, mutation_selections, optargs", [parm1(),parm2(),parm3()])
    def test_pka_iterator(self, trajectory_slice, corrected_universe, output_directory, mutation_selections, optargs):
        data, log_contents = pka_iterator(trajectory_slice=trajectory_slice,
                                        universe=corrected_universe,
                                        output_directory=output_directory,
                                        mutation_selections=mutation_selections,
                                        optargs=optargs)
        assert len(data) == trajectory_slice[1] - trajectory_slice[0]
        for key in ['residue_identifier_list', 'pka_list', 'buriedness_list']:
            assert key in data[0][0.0].keys()
        if mutation_selections is not None:
            mut_ag = corrected_universe.select_atoms(mutation_selections)
            assert (mut_ag.residues.resnames).all() == 'ALA'

    def test_pka_iterator_propka_Warning(self):
        """Testing that PROPKA warnings are captured in the log_contents."""
        pdb = """
ATOM      1  N   ASP     1      44.260  49.990  17.910  1.00  0.00           N
ATOM      2  N   ASP     2      44.260  49.990  17.910  1.00  0.00           N
ATOM      3  O   LYS     3      43.080  50.770  17.580  1.00  0.00           O
ATOM      4  OXT LYS     3      43.080  50.770  17.580  1.00  0.00           O
"""
        pstream = mda.lib.util.NamedStream(StringIO(pdb),'tmp.pdb')
        universe = mda.Universe(pstream, format='PDB')
        pstream.reset()
        output_directory = check_dir('test_output/pka_iterator_propka_warning')
        corrected_universe = create_propka_compatible_universe(universe,
                                                        hetatm_resname=None,
                                                        custom_terminal_oxygens=None,
                                                        custom_resname_correction=None)
        data, log_contents = pka_iterator(trajectory_slice=[0,1],
                                        universe=corrected_universe,
                                        output_directory=output_directory,
                                        mutation_selections=None,
                                        optargs=[])
        assert len(log_contents) > 0