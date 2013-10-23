#!/usr/bin/env python
# encoding: utf-8
"""
test_models_no_model.py

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

        #Setup a variant where all are homozygote alternative
        self.homozygote_alternative = genetic_variant.Variant(chrom = '1', start = 5, stop=5, alternative = 'A', reference = 'C', 
                                                            identity = 'rs2230749', all_info={'Ensemble_GeneID':'ENSG00000187634;'})

        sick_son.add_genotype(self.homozygote_alternative.variant_id, genotype.Genotype(GT='1/1'))
        healthy_father.add_genotype(self.homozygote_alternative.variant_id, genotype.Genotype(GT='1/1'))
        healthy_mother.add_genotype(self.homozygote_alternative.variant_id, genotype.Genotype(GT='1/0'))

        #Setup a variant that is recessive de novo
        self.recessive_dn = genetic_variant.Variant(chrom = '1', start = 6, stop=6, alternative = 'C',reference = 'T', 
                                                    identity = '.', all_info={'Ensemble_GeneID':'ENSG00000187634;'})

        sick_son.add_genotype(self.recessive_dn.variant_id, genotype.Genotype(GT='1/1'))
        healthy_father.add_genotype(self.recessive_dn.variant_id, genotype.Genotype(GT='./.'))
        healthy_mother.add_genotype(self.recessive_dn.variant_id, genotype.Genotype(GT='0/1'))

        self.recessive_family.add_individual(healthy_father)
        self.recessive_family.add_individual(sick_son)
        self.recessive_family.add_individual(healthy_mother)
        
        self.recessive_family.add_variant(self.homozygote_alternative)        
        self.recessive_family.add_variant(self.recessive_dn)
        
        self.my_healthy_father_model = genetic_models.genetic_models(self.recessive_family)
        
        for variant in self.recessive_family.variants:
            self.recessive_family.variants[variant].check_models()
        # assert True == False
    
    def test_homozygote_alternative(self):
        """Check if the genetic models are followed for the heterozygote variant"""
        assert not self.recessive_family.variants[self.homozygote_alternative.variant_id].ar
        assert not self.recessive_family.variants[self.homozygote_alternative.variant_id].ar_dn
        assert not self.recessive_family.variants[self.homozygote_alternative.variant_id].ad
        assert not self.recessive_family.variants[self.homozygote_alternative.variant_id].ad_dn
        assert not self.recessive_family.variants[self.homozygote_alternative.variant_id].x_linked
        assert not self.recessive_family.variants[self.homozygote_alternative.variant_id].x_linked_dn
        assert len(self.homozygote_alternative.models) == 1
        assert 'Na' in self.homozygote_alternative.models
        # assert True == False
    
    def test_recessive_dn(self):
        """Check if the variants are labeled as a compound pair"""
        assert not self.recessive_family.variants[self.recessive_dn.variant_id].ar
        assert self.recessive_family.variants[self.recessive_dn.variant_id].ar_dn
        assert not self.recessive_family.variants[self.recessive_dn.variant_id].ad
        assert not self.recessive_family.variants[self.recessive_dn.variant_id].ad_dn
        assert not self.recessive_family.variants[self.recessive_dn.variant_id].x_linked
        assert not self.recessive_family.variants[self.recessive_dn.variant_id].x_linked_dn
        assert len(self.recessive_dn.models) == 1
        assert 'AR_hom_denovo' in self.recessive_dn.models
    


def main():
    pass


if __name__ == '__main__':
    main()

