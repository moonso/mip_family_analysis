#!/usr/bin/env python
# encoding: utf-8
"""
test_family.py

Tests for the family class

Created by Måns Magnusson on 2013-03-13.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
from Family_Analysis import family,individual
from Family_Analysis.variants import genetic_variant, genotype



class TestFamily(object):
    """Test class for testing how the individual class behave"""
    
    def setup_class(self):
        """Setup a simple family with family id 1, sick daughter id 1, healthy father id 2, healthy mother id 3"""
        self.family = family.Family(family_id = '1') # Create a family 
        self.daughter = individual.Individual(ind = '1', family = '1', mother = '3', father = '2', sex = 2, phenotype = 2) # Create a sick daughter
        self.father = individual.Individual(ind = '2', family = '1', mother = '0', father = '0', sex = 1, phenotype = 1) # Create a healthy father
        self.mother = individual.Individual(ind = '3', family = '1', mother = '0', father = '0', sex = 2, phenotype = 1) # Create a healthy mother
        self.family.add_individual(self.daughter)
        self.family.add_individual(self.father)
        self.family.add_individual(self.mother)
        # genotypes_variant_1 = {'1': genotype.Genotype(GT='1/1'), '2': genotype.Genotype(GT='0/1'), '3': genotype.Genotype(GT='0/1')}
    
    def test_individuals(self):
        """Test if all individuals are at place"""
        assert self.daughter.ind_id in self.family.individuals
        assert self.mother.ind_id in self.family.individuals
        assert self.father.ind_id in self.family.individuals
        assert not '4' in self.family.individuals




def main():
    pass


if __name__ == '__main__':
    main()

