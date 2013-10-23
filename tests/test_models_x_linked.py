#!/usr/bin/env python
# encoding: utf-8
"""
test_models_x_linked.py

Test the so that the genetic models behave as suspected.


Created by Måns Magnusson on 2013-03-07.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
from Mip_Family_Analysis.Family import family, individual
from Mip_Family_Analysis.Models import genetic_models
from Mip_Family_Analysis.Variants import genetic_variant, genotype


class TestModelsXlinked(object):
    """Test class for testing how the genetic models behave with a recessive variant"""

    def setup_class(self):
        """Setup a simple family with family id 1, sick son id 1,
         healthy father id 2, healthy mother id 3"""
        # Setup family with sick kid, sick father and healthy mother:
        self.x_family = family.Family(family_id = '1')
        sick_son = individual.Individual(ind='1', family='1',mother='3', father='2', sex=1, phenotype=2)
        healthy_father = individual.Individual(ind='2', family='1',mother='0', father='0', sex=1, phenotype=1)
        healthy_mother = individual.Individual(ind='3', family='1',mother='0', father='0', sex=2, phenotype=1)

        #Setup variant with only autosomal dominant de novo pattern
        self.x_variant = genetic_variant.Variant(chrom = 'chrX', start = 5, stop=5, alternative = 'A',
                                            reference = 'C', identity = 'rs2230749')

        sick_son.add_genotype(self.x_variant.variant_id, genotype.Genotype(GT='0/1'))
        healthy_father.add_genotype(self.x_variant.variant_id, genotype.Genotype(GT='0/0'))
        healthy_mother.add_genotype(self.x_variant.variant_id, genotype.Genotype(GT='0/1'))
        
        #Setup variant with only autosomal recessive pattern
        self.x_dn_variant = genetic_variant.Variant(chrom = 'chrX', start = 10, stop=10, alternative = 'C',
                                            reference = 'T')

        sick_son.add_genotype(self.x_dn_variant.variant_id, genotype.Genotype(GT='0/1'))
        healthy_father.add_genotype(self.x_dn_variant.variant_id, genotype.Genotype(GT='0/0'))
        healthy_mother.add_genotype(self.x_dn_variant.variant_id, genotype.Genotype(GT='0/0'))

        self.x_family.add_individual(healthy_father)
        self.x_family.add_individual(sick_son)
        self.x_family.add_individual(healthy_mother)
                
        self.x_family.add_variant(self.x_dn_variant)
        self.x_family.add_variant(self.x_variant)
        
        
        self.my_healthy_father_model = genetic_models.genetic_models(self.x_family)
    
    def test_x_dn(self):
        """Check if the genetic models are followed for the heterozygote variant"""
        assert not self.x_family.variants[self.x_dn_variant.variant_id].ar
        assert not self.x_family.variants[self.x_dn_variant.variant_id].ar_dn
        assert not self.x_family.variants[self.x_dn_variant.variant_id].ad
        assert not self.x_family.variants[self.x_dn_variant.variant_id].ad_dn
        assert not self.x_family.variants[self.x_dn_variant.variant_id].ar_comp
        assert self.x_family.variants[self.x_dn_variant.variant_id].x_linked
        assert self.x_family.variants[self.x_dn_variant.variant_id].x_linked_dn
        self.x_dn_variant.check_models()
        assert 'X_denovo' in self.x_dn_variant.models
        assert len(self.x_dn_variant.models) == 2
    
    def test_x(self):
        """Check if the genetic models are followed for the homozygote variant"""
        assert not self.x_family.variants[self.x_variant.variant_id].ar
        assert not self.x_family.variants[self.x_variant.variant_id].ar_dn
        assert not self.x_family.variants[self.x_variant.variant_id].ad
        assert not self.x_family.variants[self.x_variant.variant_id].ad_dn
        assert not self.x_family.variants[self.x_variant.variant_id].ar_comp
        assert self.x_family.variants[self.x_variant.variant_id].x_linked
        assert not self.x_family.variants[self.x_variant.variant_id].x_linked_dn
        self.x_variant.check_models()
        assert 'X' in self.x_variant.models
        assert len(self.x_variant.models) == 1
    




def main():
    pass


if __name__ == '__main__':
    main()

