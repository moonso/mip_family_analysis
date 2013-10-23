#!/usr/bin/env python
# encoding: utf-8
"""
test_models_sick_father.py

Created by Måns Magnusson on 2013-03-13.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
from Mip_Family_Analysis.Family import family, individual
from Mip_Family_Analysis.Models import genetic_models
from Mip_Family_Analysis.Variants import genetic_variant, genotype


class TestModelsSickFather(object):
    """Test class for testing how the genetic models behave"""
    def setup_class(self):
        """Setup a simple family with family id 1, sick daughter id 1,
         healthy father id 2, healthy mother id 3"""
        # Setup family with sick kid, sick father and healthy mother:
        self.sick_father_family = family.Family(family_id = '2')
        sick_daugther = individual.Individual(ind='1', family='2',mother='3', father='2', sex=2, phenotype=2)
        sick_father = individual.Individual(ind='2', family='2',mother='0', father='0', sex=1, phenotype=2)
        healthy_mother = individual.Individual(ind='3', family='2',mother='0', father='0', sex=2, phenotype=1)

        #Setup variant with only autosomal dominant pattern
        self.dominant_variant = genetic_variant.Variant(chrom = '1', start = 5, stop=5, alternative = 'A',
                                            reference = 'C', identity = 'rs2230749')

        sick_daugther.add_genotype(self.dominant_variant.variant_id, genotype.Genotype(GT='0/1'))
        sick_father.add_genotype(self.dominant_variant.variant_id, genotype.Genotype(GT='0/1'))
        healthy_mother.add_genotype(self.dominant_variant.variant_id, genotype.Genotype(GT='0/0'))
        
        #Setup variant with only autosomal recessive pattern(Should not work here)
        self.recessive_variant = genetic_variant.Variant(chrom = '1', start = 10, stop=10, alternative = 'C',
                                            reference = 'T')

        sick_daugther.add_genotype(self.recessive_variant.variant_id, genotype.Genotype(GT='1/1'))
        sick_father.add_genotype(self.recessive_variant.variant_id, genotype.Genotype(GT='0/1'))
        healthy_mother.add_genotype(self.recessive_variant.variant_id, genotype.Genotype(GT='0/1'))

        self.sick_father_family.add_individual(sick_father)
        self.sick_father_family.add_individual(sick_daugther)
        self.sick_father_family.add_individual(healthy_mother)
                
        self.sick_father_family.add_variant(self.dominant_variant)
        self.sick_father_family.add_variant(self.recessive_variant)
        
        
        self.my_sick_father_model = genetic_models.genetic_models(self.sick_father_family)
    
    def test_dominant(self):
        """Check if the genetic models are followed for the heterozygote variant"""
        assert not self.sick_father_family.variants[self.dominant_variant.variant_id].ar
        assert not self.sick_father_family.variants[self.dominant_variant.variant_id].ar_dn
        assert self.sick_father_family.variants[self.dominant_variant.variant_id].ad
        assert not self.sick_father_family.variants[self.dominant_variant.variant_id].ad_dn
        assert not self.sick_father_family.variants[self.dominant_variant.variant_id].ar_comp
        assert not self.sick_father_family.variants[self.dominant_variant.variant_id].x_linked
        assert not self.sick_father_family.variants[self.dominant_variant.variant_id].x_linked_dn
        self.dominant_variant.check_models()
        assert 'AD' in self.dominant_variant.models
        assert len(self.dominant_variant.models) == 1
    
    def test_recessive(self):
        """Check if the genetic models are followed for the homozygote variant"""
        assert not self.sick_father_family.variants[self.recessive_variant.variant_id].ar
        assert not self.sick_father_family.variants[self.recessive_variant.variant_id].ar_dn
        assert not self.sick_father_family.variants[self.recessive_variant.variant_id].ad
        assert not self.sick_father_family.variants[self.recessive_variant.variant_id].ad_dn
        assert not self.sick_father_family.variants[self.recessive_variant.variant_id].ar_comp
        assert not self.sick_father_family.variants[self.recessive_variant.variant_id].x_linked
        assert not self.sick_father_family.variants[self.recessive_variant.variant_id].x_linked_dn
        self.recessive_variant.check_models()
        assert 'Na' in self.recessive_variant.models
        assert len(self.recessive_variant.models) == 1
    

def main():
    pass


if __name__ == '__main__':
    main()

