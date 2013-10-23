#!/usr/bin/env python
# encoding: utf-8
"""
test_indivdual.py

Test the individual class.

Created by Måns Magnusson on 2013-03-07.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
from Family_Analysis import individual
from Family_Analysis.variants import genotype


class TestIndividual(object):
    """Test class for testing how the individual class behave"""
    
    def setup_class(self):
        """Setup a simple family with family id 1, sick daughter id 1, healthy father id 2, healthy mother id 3"""
        self.daughter = individual.Individual(ind='1', family='1', mother='3', father='2', sex=2, phenotype=2)
        self.father = individual.Individual(ind='2', family='1', mother='0', father='0', sex=1, phenotype=1)
        self.mother = individual.Individual(ind='3', family='1', mother='0', father='0', sex=2, phenotype=1)
        self.daughter_genotypes = {}
        self.father_genotypes = {}
        self.mother_genotypes = {}
        self.daughter_genotypes['1_1_T_A'] = genotype.Genotype(GT='0/1')
        self.daughter_genotypes['1_3_A_C'] = genotype.Genotype(GT='1/1')
        self.father_genotypes['1_1_T_A'] = genotype.Genotype(GT='0/0')
        self.father_genotypes['1_3_A_C'] = genotype.Genotype(GT='0/1')
        self.mother_genotypes['1_1_T_A'] = genotype.Genotype(GT='./.')
        self.mother_genotypes['1_3_A_C'] = genotype.Genotype(GT='0/1')
    
    def test_daughter(self):
        """Test if the information about the daughter comes out correctly."""
        assert self.daughter.get_info() == {'ind_id':'1', 'family':'1', 'mother':'3', 'father':'2', 'sex':2, 'phenotype':2, 'phasing':False}
    
    def test_disease_status(self):
        """Test if the disease status is correct for an individual"""
        assert self.daughter.affected()
        assert not self.father.affected()
        assert not self.mother.affected()
    
    def test_add_genotype(self):
        """Test if the adding of genotypes works well."""
        for variant_id, genotype in self.daughter_genotypes.items():
            self.daughter.add_genotype(variant_id, genotype)
        assert len(self.daughter.genotypes) == 2
    
    def test_get_genotype_daughter(self):
        """Test if the correct genotype is returned when asked."""
        assert self.daughter.get_genotype('1_1_T_A').genotype == '0/1' 
    


def main():
    pass


if __name__ == '__main__':
    main()

