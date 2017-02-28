import unittest
import lattice_mc
import random

class IntegrationTestCase( unittest.TestCase ):

    def test_simulation_runs_with_selected_sites( self ):
        a, b, c, = 4, 4, 4
        spacing = 1.0
        n_atoms = 32
        selected_sites = 'L'
        number_of_jumps = 10
        s = lattice_mc.Simulation()
        s.lattice = lattice_mc.init_lattice.cubic_lattice( a, b, c, spacing )
        s.lattice.transmute_sites( 'L', 'X', n_atoms )
        s.set_number_of_atoms( n_atoms, selected_sites=selected_sites )
        [ self.assertEqual( site.label, selected_sites ) for site in s.lattice.occupied_sites() ]
        s.setup_lookup_table()
        s.set_number_of_jumps( number_of_jumps )
        self.assertEqual( s.tracer_correlation, None )
        self.assertEqual( s.collective_correlation, None )
        self.assertEqual( s.tracer_diffusion_coefficient, None )
        self.assertEqual( s.collective_diffusion_coefficient_per_atom, None )
        self.assertEqual( s.average_site_occupations, None )
        s.run()
        self.assertNotEqual( s.tracer_correlation, None )
        self.assertNotEqual( s.collective_correlation, None )
        self.assertNotEqual( s.tracer_diffusion_coefficient, None )
        self.assertNotEqual( s.collective_diffusion_coefficient_per_atom, None )
        self.assertNotEqual( s.average_site_occupations, None )

if __name__ == '__main__':
    unittest.main()
