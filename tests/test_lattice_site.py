import unittest
from lattice_mc.lattice_site import Site
from unittest.mock import Mock

class SiteTestCase( unittest.TestCase ):
    """Test for Site class"""

    def setUp( self ):
        self.number = 'foo1'
        self.coordinates = 'foo2'
        self.neighbours = 'foo3'
        self.energy = 'foo4'
        self.label = 'foo5'
        self.cn_energies = { 'A' : [ 1.0, 2.0, 3.0 ], 'B' : [ 2.0, 4.0, 6.0 ] }
        self.site = Site( self.number, 
                          self.coordinates, 
                          self.neighbours, 
                          self.energy, 
                          self.label, 
                          cn_energies=self.cn_energies )
        self.neighbouring_sites = [ Mock( spec=Site, is_occupied=True, label='A' ),
                               Mock( spec=Site, is_occupied=False, label='A' ),
                               Mock( spec=Site, is_occupied=True, label='B' ),
                               Mock( spec=Site, is_occupied=False, label='B' ) ]
  
    def test_site_is_initialised( self ):
        self.assertEqual( self.site.number, self.number )
        self.assertEqual( self.site.index, Site.index - 1 )
        self.assertEqual( self.site.r, self.coordinates )
        self.assertEqual( self.site.neighbours, self.neighbours )
        self.assertEqual( self.site.p_neighbours, None )
        self.assertEqual( self.site.energy, self.energy )
        self.assertEqual( self.site.occupation, 0 )
        self.assertEqual( self.site.atom, None )
        self.assertEqual( self.site.is_occupied, False )
        self.assertEqual( self.site.label, self.label )
        self.assertEqual( self.site.time_occupied, 0.0 )
        self.assertEqual( self.site.cn_occupation_energies, self.cn_energies )

    def test_nn_occupation( self ):
        self.site.p_neighbours = self.neighbouring_sites
        self.assertEqual( self.site.nn_occupation(), 2 )

    def test_site_specific_nn_occupation( self ):
        self.site.p_neighbours = self.neighbouring_sites
        self.assertEqual( self.site.site_specific_nn_occupation(), { 'A' : 1, 'B' : 1 } )

    def test_site_specific_neighbours( self ):
        self.site.p_neighbours = self.neighbouring_sites
        self.assertEqual( self.site.site_specific_neighbours(), { 'A' : 2, 'B' : 2 } )

    def test_set_cn_occupation_energies( self ):
        self.site.set_cn_occupation_energies( 'foo7' )
        self.assertEqual( self.site.cn_occupation_energies, 'foo7' )

    def test_unchanged_cn_occupation_energy( self ):
        self.site.p_neighbours = self.neighbouring_sites
        self.site_specific_nn_occupation = Mock( return_value={ 'A' : 2, 'B' : 2 } )
        self.assertEqual( self.site.cn_occupation_energy(), 6.0 )

    def test_changed_cn_occupation_energy( self ):
        delta_occupation = { 'A' : 1, 'B' : -1 }
        self.site.p_neighbours = self.neighbouring_sites
        self.site_specific_nn_occupation = Mock( return_value={ 'A' : 2, 'B' : 2 } )
        self.assertEqual( self.site.cn_occupation_energy( delta_occupation=delta_occupation), 5.0 )

if __name__ == '__main__':
    unittest.main()
