import unittest
from lattice_mc.simulation import Simulation
from lattice_mc.lattice import Lattice
from lattice_mc.species import Species
from lattice_mc.atom import Atom
from unittest.mock import Mock, PropertyMock, patch
import numpy as np

class SimulationTestCase( unittest.TestCase ):
    """Test for Species class"""

    def test_simulation_is_initialised( self ):
        simulation = Simulation()
        self.assertEqual( simulation.lattice, None )
        self.assertEqual( simulation.number_of_atoms, None )
        self.assertEqual( simulation.number_of_equilibration_jumps, 0 )
        self.assertEqual( simulation.number_of_atoms, None )
        self.assertEqual( simulation.has_run, False )

    def test_reset( self ):
        simulation = Simulation()
        simulation.lattice = Mock( spec=Lattice )
        simulation.atoms = Mock( spec=Species )
        simulation.atoms.atoms = [ Mock( spec=Atom ), Mock( spec=Atom ) ]
        simulation.reset()
        simulation.lattice.reset.assert_called()
        simulation.atoms.atoms[0].reset.assert_called()
        simulation.atoms.atoms[1].reset.assert_called()

    def test_set_number_of_atoms( self ):
        simulation = Simulation()
        simulation.lattice = Mock( spec=Lattice )
        mock_atoms = [ Mock( spec=Atom ), Mock( spec=Atom ) ]
        simulation.lattice.populate_sites = Mock( return_value = mock_atoms )
        with patch( 'lattice_mc.species.Species' ) as mock_Species:
            mock_Species.return_value = 'some atoms'
            simulation.set_number_of_atoms( 3 )
            self.assertEqual( simulation.number_of_atoms, 3 )
            self.assertEqual( simulation.atoms, 'some atoms' )
            mock_Species.assert_called_with( mock_atoms )
    
    def test_set_number_of_jumps( self ):
        simulation = Simulation()
        simulation.set_number_of_jumps( 32 )
        self.assertEqual( simulation.number_of_jumps, 32 )
         
    def test_set_number_of_equilibration_jumps( self ):
        simulation = Simulation()
        simulation.set_number_of_equilibration_jumps( 32 )
        self.assertEqual( simulation.number_of_equilibration_jumps, 32 )
    
    def test_define_lattice_from_file( self ):
        with patch( 'lattice_mc.init_lattice.lattice_from_sites_file' ) as mock_lattice_from_file:
            mock_lattice_from_file.return_value = 'foo'
            simulation = Simulation()
            cell_lengths = np.array( [ 1.0, 2.0, 3.0 ] )
            simulation.define_lattice_from_file( 'filename', cell_lengths )
            self.assertEqual( simulation.lattice, 'foo' ) 
            mock_lattice_from_file.assert_called_with( 'filename', cell_lengths=cell_lengths )

    def test_set_nn_energy( self ):
        simulation = Simulation()
        simulation.lattice = Mock( spec=Lattice )
        simulation.lattice.set_nn_energy = Mock()
        simulation.set_nn_energy( 'foo' )
        simulation.lattice.set_nn_energy.assert_called_with( 'foo' )

    def test_set_cn_energies( self ):
        simulation = Simulation()
        simulation.lattice = Mock( spec=Lattice )
        simulation.lattice.set_cn_energies = Mock()
        simulation.set_cn_energies( 'foo' )
        simulation.lattice.set_cn_energies.assert_called_with( 'foo' )

    def test_set_site_energies( self ):
        simulation = Simulation()
        simulation.lattice = Mock( spec=Lattice )
        simulation.lattice.set_site_energies = Mock()
        simulation.set_site_energies( 'foo' )
        simulation.lattice.set_site_energies.assert_called_with( 'foo' )

    def test_run_fails_with_no_lattice( self ):
        simulation = Simulation()
        simulation.atoms = 'foo'
        simulation.number_of_jumps = 'bar'
        simulation.lattice = None
        with self.assertRaises( AttributeError ):
            simulation.run()
    
    def test_run_fails_with_no_atoms( self ):
        simulation = Simulation()
        simulation.atoms = None
        simulation.number_of_jumps = 'bar'
        simulation.lattice = 'foo'
        with self.assertRaises( AttributeError ):
            simulation.run()
        
    def test_run_fails_with_no_number_of_jumps( self ):
        simulation = Simulation()
        simulation.atoms = 'bar'
        simulation.number_of_jumps = None
        simulation.lattice = 'foo'
        with self.assertRaises( AttributeError ):
            simulation.run()

    def test_run_without_equilibration_steps( self ):
        simulation = Simulation()
        simulation.atoms = 'a'
        simulation.lattice = Mock( spec=Lattice )
        simulation.lattice.jump = Mock()
        simulation.number_of_jumps = 10
        simulation.run()
        self.assertEqual( simulation.lattice.jump.call_count, 10 )

    def test_run_for_time( self ):
        simulation = Simulation()
        simulation.atoms = 'a'
        simulation.lattice = Mock( spec=Lattice )
        simulation.lattice.time = 0.0
        def fake_jump_method():
            simulation.lattice.time += 1.0
        simulation.lattice.jump = fake_jump_method
        simulation.run( for_time=10.0 )
        self.assertEqual( simulation.lattice.time, 10.0 )
        self.assertEqual( simulation.number_of_jumps, 10 )

    def test_run_with_equilibration_steps( self ):
        simulation = Simulation()
        simulation.atoms = 'a'
        simulation.lattice = Mock( spec=Lattice )
        simulation.lattice.jump = Mock()
        simulation.reset = Mock()
        simulation.number_of_equilibration_jumps = 20
        simulation.number_of_jumps = 30
        simulation.run()
        self.assertEqual( simulation.lattice.jump.call_count, 20+30 )
        simulation.reset.assert_called()

class SimulationResultsTestCase( unittest.TestCase ):

    def setUp( self ):
        self.simulation = Simulation()
        self.simulation.has_run = True

    def test_tracer_correlation( self ):
        s = self.simulation
        s.atoms = Mock()
        s.atoms.sum_dr_squared = PropertyMock( return_value=10.0 )
        s.number_of_jumps = 5
        self.assertEqual( s.tracer_correlation, 2.0 )

    def test_new_tracer_correlation( self ):
        s = self.simulation
        s.atoms = Mock( spec=Species )
        s.atoms.tracer_correlation = Mock( return_value=3.3 )
        self.assertEqual( s.new_tracer_correlation, 3.3 )

    def test_tracer_diffusion_coefficient( self ):
        s = self.simulation
        s.atoms = Mock( spec=Species )
        s.atoms.sum_dr_squared = Mock( return_value=15.0 )
        s.number_of_atoms = 5
        s.lattice = Mock( spec=Lattice )
        s.lattice.time = 12.0
        self.assertEqual( s.tracer_diffusion_coefficient, 15.0 / ( 6.0 * 5.0 * 12.0 ) )

    def test_collective_correlation( self ):
        s = self.simulation
        s.atoms = Mock( spec=Species )
        s.atoms.collective_dr_squared = Mock( return_value=12.0 )
        s.number_of_jumps = 4
        self.assertEqual( s.collective_correlation, 3.0 )

    def test_new_collective_correlation( self ):
        s = self.simulation
        s.atoms = Mock( spec=Species )
        s.atoms.collective_correlation = Mock( return_value=3.0 )
        self.assertEqual( s.new_collective_correlation, 3.0 )

    def test_collective_diffusion_coefficient( self ):
        s = self.simulation
        s.atoms = Mock( spec=Species )
        s.atoms.collective_dr_squared = Mock( return_value = 36.0 )
        s.lattice = Mock( spec=Lattice )
        s.lattice.time = 2.0
        self.assertEqual( s.collective_diffusion_coefficient, 3.0 )

    def test_collective_diffusion_coefficient_per_atom( self ):
        s = self.simulation
        with patch( 'lattice_mc.simulation.Simulation.collective_diffusion_coefficient', new_callable=PropertyMock ) as mock_cdc:
            mock_cdc.return_value = 8.0
            s.number_of_atoms = 2
            self.assertEqual( s.collective_diffusion_coefficient_per_atom, 4.0 )

    def test_average_site_occupations( self ):
        s = self.simulation
        s.lattice = Mock( spec=Lattice )
        s.lattice.site_occupation_statistics = Mock( return_value='foo' )
        self.assertEqual( s.average_site_occupations, 'foo' )

    def test_setup_lookup_table( self ):
        s = self.simulation
        s.lattice = Mock( spec=Lattice )
        with patch( 'lattice_mc.lookup_table.LookupTable' ) as mock_lookup_table:
            mock_lookup_table.return_value = 'foo'
            s.setup_lookup_table() 
            self.assertEqual( s.lattice.jump_lookup_table, 'foo' )
            mock_lookup_table.assert_called_with( s.lattice, 'nearest-neighbour' )

    def test_setup_lookup_table_with_hamiltonian( self ):
        s = self.simulation
        s.lattice = Mock( spec=Lattice )
        with patch( 'lattice_mc.lookup_table.LookupTable' ) as mock_lookup_table:
            mock_lookup_table.return_value = 'foo'
            s.setup_lookup_table( hamiltonian='coordination_number' )
            self.assertEqual( s.lattice.jump_lookup_table, 'foo' )
            mock_lookup_table.assert_called_with( s.lattice, 'coordination_number' )

    def test_setup_lookup_table_with_non_implemented_hamiltonian( self ):
        s = self.simulation
        with self.assertRaises( ValueError ):
            s.setup_lookup_table( hamiltonian='bar' )

class SimulationNoResultsTestCase( unittest.TestCase ):

    def test_tracer_correlation( self ):
        s = Simulation()
        self.assertEqual( s.tracer_correlation, None )
    
    def test_new_tracer_correlation( self ):
        s = Simulation()
        self.assertEqual( s.new_tracer_correlation, None )

    def test_tracer_diffusion_coefficient( self ):
        s = Simulation()
        self.assertEqual( s.tracer_diffusion_coefficient, None )

    def test_collective_correlation( self ):
        s = Simulation()
        self.assertEqual( s.collective_correlation, None )

    def test_new_collective_correlation( self ):
        s = Simulation()
        self.assertEqual( s.new_collective_correlation, None )

    def test_collective_diffusion_coefficient( self ):
        s = Simulation()
        self.assertEqual( s.collective_diffusion_coefficient, None )

    def test_collective_diffusion_coefficient_per_atom( self ):
        s = Simulation()
        self.assertEqual( s.collective_diffusion_coefficient_per_atom, None )

if __name__ == '__main__':
    unittest.main() 
