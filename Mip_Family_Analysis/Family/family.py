#!/usr/bin/env python
# encoding: utf-8
"""
family.py

Holds the meta information of a family and its individuals.

    - has a Individual
    - has a Variant

Attributes:

individuals DICT dictionary with family members on the form {<ind_id>:<Individual_obj>}
variants DICT dictionary with all the variants that exists in the family on the form {<var_id>:<Variant_obj>}

    
Methods:

    - print_individuals
    - print_all_variants
    - add_variant(attach a variant objekt and also a variant id to individuals)


Created by Måns Magnusson on 2012-11-01.
Copyright (c) 2012 __MyCompanyName__. All rights reserved.
"""

import sys
import os
from Mip_Family_Analysis.Utils.is_number import is_number
from collections import OrderedDict

class Family(object):
    """Base class for the family parsers."""
    def __init__(self, family_id, individuals = {}):
        super(Family, self).__init__()
        self.individuals = individuals # This is a dictionary on the form {ind_id:<Individual object>}
        self.family_id = family_id
        self.models_of_inheritance = [] # List of models of inheritance that should be prioritized.
        self.variants = OrderedDict() # Dictionary on the form {var_id:<Variant object>}
    
    def family_check(self):
        """Check if the family members break the structure of the family in any way, eg. nonexistent parent, wrong sex on parent..."""
        #TODO Make some tests for these
        pass
    
    def add_individual(self, individual_object):
        """Add an individual to the family."""
        self.individuals[individual_object.ind_id] = individual_object
    
    def add_variant(self, variant_object):
        """Adds the variants to the family."""
        self.variants[variant_object.variant_id] = variant_object
    
    def print_all_variants(self):
        """Print all original variants."""
        for variant in self.variants:
            print variant
            for individual in self.individuals:
                genotype = self.individuals[individual].genotypes[variant]
                print individual, genotype
    
    def __str__(self):
        """Print the family members of this family"""
        family = [individual for individual in self.individuals]
        return "\t".join(family)


def main():
    pass


if __name__ == '__main__':
    main()
