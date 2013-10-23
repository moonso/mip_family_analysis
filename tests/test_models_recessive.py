#!/usr/bin/env python
# encoding: utf-8
"""
test_models_recessive.py

Test the so that the genetic models behave as suspected.


Created by Måns Magnusson on 2013-03-07.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
from Mip_Family_Analysis.Family import family, individual
from Mip_Family_Analysis.Models import genetic_models
from Mip_Family_Analysis.Variants import genetic_variant, genotype


class TestModelsRecessive(object):
    """Test class for testing how the genetic models behave with a recessive variant"""

    def setup_class(self):
        """Setup a simple family with family id 1, sick son id 1,
         healthy father id 2, healthy mother id 3"""
        # Setup family with sick kid, sick father and healthy mother:
        self.recessive_family = family.Family(family_id = '1')
        sick_son = individual.Individual(ind='1', family='1',mother='3', father='2', sex=1, phenotype=2)
        healthy_father = individual.Individual(ind='2', family='1',mother='0', father='0', sex=1, phenotype=1)
        healthy_mother = individual.Individual(ind='3', family='1',mother='0', father='0', sex=2, phenotype=1)

        #Setup variant with only autosomal dominant de novo pattern
        self.recessive_dn_variant = genetic_variant.Variant(chrom = '1', start = 5, stop=5, alternative = 'A',
                                            reference = 'C', identity = 'rs2230749')

        sick_son.add_genotype(self.recessive_dn_variant.variant_id, genotype.Genotype(GT='1/1'))
        healthy_father.add_genotype(self.recessive_dn_variant.variant_id, genotype.Genotype(GT='0/1'))
        healthy_mother.add_genotype(self.recessive_dn_variant.variant_id, genotype.Genotype(GT='0/0'))
        
        #Setup variant with only autosomal recessive pattern
        self.recessive_variant = genetic_variant.Variant(chrom = '1', start = 10, stop=10, alternative = 'C',
                                            reference = 'T')

        sick_son.add_genotype(self.recessive_variant.variant_id, genotype.Genotype(GT='1/1'))
        healthy_father.add_genotype(self.recessive_variant.variant_id, genotype.Genotype(GT='0/1'))
        healthy_mother.add_genotype(self.recessive_variant.variant_id, genotype.Genotype(GT='0/1'))

        #Setup potential recessive but does not follow any patterns
        self.almost_recessive_variant = genetic_variant.Variant(chrom = '1', start = 20, stop=20, alternative = 'C',
                                            reference = 'T')

        sick_son.add_genotype(self.almost_recessive_variant.variant_id, genotype.Genotype(GT='./.'))
        healthy_father.add_genotype(self.almost_recessive_variant.variant_id, genotype.Genotype(GT='0/1'))
        healthy_mother.add_genotype(self.almost_recessive_variant.variant_id, genotype.Genotype(GT='0/1'))

        self.recessive_family.add_individual(healthy_father)
        self.recessive_family.add_individual(sick_son)
        self.recessive_family.add_individual(healthy_mother)
                
        self.recessive_family.add_variant(self.recessive_dn_variant)
        self.recessive_family.add_variant(self.recessive_variant)
        self.recessive_family.add_variant(self.almost_recessive_variant)
        
        
        self.my_healthy_father_model = genetic_models.genetic_models(self.recessive_family)

    def test_almost_recessive(self):
        """Check if the genetic models are followed for the heterozygote variant"""
        self.almost_recessive_variant.check_models()
        assert self.recessive_family.variants[self.almost_recessive_variant.variant_id].ar
        assert not self.recessive_family.variants[self.almost_recessive_variant.variant_id].ar_dn
        assert not self.recessive_family.variants[self.almost_recessive_variant.variant_id].ad
        assert not self.recessive_family.variants[self.almost_recessive_variant.variant_id].ad_dn
        assert not self.recessive_family.variants[self.almost_recessive_variant.variant_id].ar_comp
        assert not self.recessive_family.variants[self.almost_recessive_variant.variant_id].x_linked
        assert not self.recessive_family.variants[self.almost_recessive_variant.variant_id].x_linked_dn
        assert 'AR_hom' in self.almost_recessive_variant.models
        assert len(self.recessive_dn_variant.models) == 1

    
    def test_recessive_dn(self):
        """Check if the genetic models are followed for the heterozygote variant"""
        assert not self.recessive_family.variants[self.recessive_dn_variant.variant_id].ar
        assert self.recessive_family.variants[self.recessive_dn_variant.variant_id].ar_dn
        assert not self.recessive_family.variants[self.recessive_dn_variant.variant_id].ad
        assert not self.recessive_family.variants[self.recessive_dn_variant.variant_id].ad_dn
        assert not self.recessive_family.variants[self.recessive_dn_variant.variant_id].ar_comp
        assert not self.recessive_family.variants[self.recessive_dn_variant.variant_id].x_linked
        assert not self.recessive_family.variants[self.recessive_dn_variant.variant_id].x_linked_dn
        self.recessive_dn_variant.check_models()
        assert 'AR_hom_denovo' in self.recessive_dn_variant.models
        assert len(self.recessive_dn_variant.models) == 1
    
    def test_recessive(self):
        """Check if the genetic models are followed for the homozygote variant"""
        assert self.recessive_family.variants[self.recessive_variant.variant_id].ar
        assert not self.recessive_family.variants[self.recessive_variant.variant_id].ar_dn
        assert not self.recessive_family.variants[self.recessive_variant.variant_id].ad
        assert not self.recessive_family.variants[self.recessive_variant.variant_id].ad_dn
        assert not self.recessive_family.variants[self.recessive_variant.variant_id].ar_comp
        assert not self.recessive_family.variants[self.recessive_variant.variant_id].x_linked
        assert not self.recessive_family.variants[self.recessive_variant.variant_id].x_linked_dn
        self.recessive_variant.check_models()
        assert 'AR_hom' in self.recessive_variant.models
        assert len(self.recessive_variant.models) == 1
    




def main():
    pass


if __name__ == '__main__':
    main()

