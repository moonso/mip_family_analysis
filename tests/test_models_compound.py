#!/usr/bin/env python
# encoding: utf-8
"""
test_models_compound.py

Test the so that the genetic models behave as suspected.


Created by Måns Magnusson on 2013-03-07.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
from Mip_Family_Analysis.Family import family, individual
from Mip_Family_Analysis.Models import genetic_models
from Mip_Family_Analysis.Variants import genetic_variant, genotype


class TestModelsCompound(object):
    """Test class for testing how the genetic models behave with a recessive variant"""

    def setup_class(self):
        """Setup a simple family with family id 1, sick son id 1,
         healthy father id 2, healthy mother id 3"""
        # Setup family with sick kid, sick father and healthy mother:
        self.recessive_family = family.Family(family_id = '1')
        sick_son = individual.Individual(ind='1', family='1',mother='3', father='2', sex=1, phenotype=2)
        healthy_father = individual.Individual(ind='2', family='1',mother='0', father='0', sex=1, phenotype=1)
        healthy_mother = individual.Individual(ind='3', family='1',mother='0', father='0', sex=2, phenotype=1)

        #Setup two variants with only autosomal recessive pattern
        self.recessive_comp_variant_1 = genetic_variant.Variant(chrom = '1', start = 5, stop=5, alternative = 'A',
                                            reference = 'C', identity = 'rs2230749', 
                                             all_info={'Ensemble_GeneID':'ENSG00000187634;'})

        self.recessive_comp_variant_2 = genetic_variant.Variant(chrom = '1', start = 10, stop=10, alternative = 'C',
                                            reference = 'T', identity = '.', 
                                             all_info={'Ensemble_GeneID':'ENSG00000187634;'})

        sick_son.add_genotype(self.recessive_comp_variant_1.variant_id, genotype.Genotype(GT='0/1'))
        healthy_father.add_genotype(self.recessive_comp_variant_1.variant_id, genotype.Genotype(GT='0/1'))
        healthy_mother.add_genotype(self.recessive_comp_variant_1.variant_id, genotype.Genotype(GT='0/0'))
        
        sick_son.add_genotype(self.recessive_comp_variant_2.variant_id, genotype.Genotype(GT='0/1'))
        healthy_father.add_genotype(self.recessive_comp_variant_2.variant_id, genotype.Genotype(GT='0/0'))
        healthy_mother.add_genotype(self.recessive_comp_variant_2.variant_id, genotype.Genotype(GT='0/1'))

        self.recessive_family.add_individual(healthy_father)
        self.recessive_family.add_individual(sick_son)
        self.recessive_family.add_individual(healthy_mother)
                
        self.recessive_family.add_variant(self.recessive_comp_variant_1)
        self.recessive_family.add_variant(self.recessive_comp_variant_2)
        
        
        self.my_healthy_father_model = genetic_models.genetic_models(self.recessive_family)
        for variant in self.recessive_family.variants:
            self.recessive_family.variants[variant].check_models()
        # assert True == False
    
    def test_recessive(self):
        """Check if the genetic models are followed for the heterozygote variant"""
        assert not self.recessive_family.variants[self.recessive_comp_variant_1.variant_id].ar
        assert not self.recessive_family.variants[self.recessive_comp_variant_1.variant_id].ar_dn
        assert not self.recessive_family.variants[self.recessive_comp_variant_1.variant_id].ad
        assert not self.recessive_family.variants[self.recessive_comp_variant_1.variant_id].ad_dn
        assert not self.recessive_family.variants[self.recessive_comp_variant_1.variant_id].x_linked
        assert not self.recessive_family.variants[self.recessive_comp_variant_1.variant_id].x_linked_dn
        assert not self.recessive_family.variants[self.recessive_comp_variant_2.variant_id].ar
        assert not self.recessive_family.variants[self.recessive_comp_variant_2.variant_id].ar_dn
        assert not self.recessive_family.variants[self.recessive_comp_variant_2.variant_id].ad
        assert not self.recessive_family.variants[self.recessive_comp_variant_2.variant_id].ad_dn
        assert not self.recessive_family.variants[self.recessive_comp_variant_2.variant_id].x_linked
        assert not self.recessive_family.variants[self.recessive_comp_variant_2.variant_id].x_linked_dn
        assert len(self.recessive_comp_variant_1.models) == 1
        assert len(self.recessive_comp_variant_1.models) == 1
        # assert True == False
    
    def test_compound(self):
        """Check if the variants are labeled as a compound pair"""
        assert self.recessive_family.variants[self.recessive_comp_variant_1.variant_id].ar_comp
        assert self.recessive_family.variants[self.recessive_comp_variant_2.variant_id].ar_comp
    
    # def test_recessive(self):
    #     """Check if the genetic models are followed for the homozygote variant"""
    #     assert self.recessive_family.variants[self.recessive_variant.variant_id].ar
    #     assert not self.recessive_family.variants[self.recessive_variant.variant_id].ar_dn
    #     assert not self.recessive_family.variants[self.recessive_variant.variant_id].ad
    #     assert not self.recessive_family.variants[self.recessive_variant.variant_id].ad_dn
    #     assert not self.recessive_family.variants[self.recessive_variant.variant_id].ar_comp
    #     assert not self.recessive_family.variants[self.recessive_variant.variant_id].x_linked
    #     assert not self.recessive_family.variants[self.recessive_variant.variant_id].x_linked_dn
    #     self.recessive_variant.check_models()
    #     assert 'AR_hom' in self.recessive_variant.models
    #     assert len(self.recessive_variant.models) == 1
    




def main():
    pass


if __name__ == '__main__':
    main()

