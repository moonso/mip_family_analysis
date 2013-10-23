#!/usr/bin/env python
# encoding: utf-8
"""
test_family_parser.py

.ped and .fam always have 6 columns, these are

Family_ID 
Individual_ID 
Paternal_ID 
Maternal_ID 
Sex(1=male; 2=female; other=unknown)
Phenotype(-9 missing, 0 missing, 1 unaffected, 2 affected)


        self.family_type
        self.families
        
        def ped_parser(self, individual_line, line_count):
        def broad_parser(self, individual_line, line_count):

Created by Måns Magnusson on 2013-08-19.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
from Family_Analysis.utils import family_parser

class TestFamilyParser(object):
    """Test class for testing how the get_variant functions behave."""
    
    def setup_class(self):
        """Setup some variants"""
        ped_header = "#Fam_ID\tIDN(FDN-GDN-SDN(DS)\tFather_ID\tMother_ID\tSex\tPhenotype\n"
        my_sick_daughter = "30\t30-1-2A\t30-2-1U\t30-2-2U\t2\t2\n"
        my_healthy_father = "30\t30-2-1U\t0\t0\t1\t1\n"
        my_healthy_mother = "30\t30-2-2U\t0\t0\t2\t1\n"
        
        my_sick_son = "F59\t6402\t6404\t6405\t1\t2\n"
        my_sick_father = "F59\t6404\t0\t0\t1\t2\n"
        my_other_healthy_mother = "F59\t6405\t5\t6\t2\t1\n"
        
        path_to_test_get_family = os.path.dirname(os.path.abspath(__file__))
        self.ped_file = path_to_test_get_family + '/simple_family.ped'
        f = open(self.ped_file, 'w')
        f.write(ped_header)
        f.write(my_sick_daughter)
        f.write(my_healthy_mother)
        f.write(my_healthy_father)
        f.write(my_sick_son)
        f.write(my_other_healthy_mother)
        f.write(my_sick_father)
        f.close()
        self.my_family_parser = family_parser.FamilyParser(self.ped_file, 'ped')
        
    def test_family(self):
        """Check if the variant is stored with the correct id"""
        assert len(self.my_family_parser.families) == 2
        assert len(self.my_family_parser.families['30'].individuals) == 3
        assert len(self.my_family_parser.families['F59'].individuals) == 3

    def teardown_class(self):
        """Teardown the objects"""
        os.remove(self.ped_file)

def main():
    pass


if __name__ == '__main__':
    main()

