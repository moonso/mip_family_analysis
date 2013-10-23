#!/usr/bin/env python
# encoding: utf-8
"""
genetic_models.py

Genetic models take a family object with individuals and variants and annotates for each variant which models they follow in this family.

The following models are checked:

- Autosomal Dominant(AD)
- Autosomal Dominant De Novo(AD_DN)
- Autosomal Recessive(AR)
- Autosomal Recessive De Novo(AR_DN)
- Autosomal Recesive Compound.

In this model a variant must imply affected status, otherwise it can not be dominant. All sick has to be ay least heterozygote for the variant and all healthy can not have it.

We will assume that each individual that we have information about is present among the individual in self.family.individuals.


is the individual sick?

    - If the individual is homozygote alternative then AD/AD-denovo and AR/AR-denovo are ok

    - If the individual is is heterozygote then AD/AD-denovo are ok but AR/AR-denovo are not ok

    - If the individual is homozygote reference no model is ok

    - If the individual has no call all models are ok



is the individual healthy?

    - If the individual is homozygote alternative then no model is ok

    - If the individual is heterozygote then AR/AR-denove are ok but AD/AD-denovo are not ok

    - If the individual is homozygote referense all models are ok

    - If there is no call all models are ok



is there no known phenotype?

    - All models are ok for all variants



Created by Måns Magnusson on 2013-02-12.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
from Mip_Family_Analysis.Variants import genotype
from Mip_Family_Analysis.Utils import pair_generator


class genetic_models(object):
    """Holds the genetic models."""
    def __init__(self, family_object, gene_annotation = 'Ensembl', max_variants = 200):
        super(genetic_models, self).__init__()
        #Dictionary with <Gene ID> : [variant_id_1, variant_id_2, ...] for controlling the compound heterozygotes
        self.gene_variants = {}
        #Object with individuals and variants: 
        self.family = family_object
        print 'Variants:', self.family.variants
        #These are the genes with to many variants:
        genes_with_many_variants = []
        # For each variant in the family:
        for variant_id in self.family.variants:
                        
            # Only check X-linked for the variants in the X-chromosome:
            if self.family.variants[variant_id].chr == 'X':
                self.check_x_linked(variant_id)
            
            else:
                # Check the dominant model:
                self.check_dominant(variant_id)
                # Check the recessive model:
                self.check_recessive(variant_id)
            
            # Prepare for compounds by adding all variants for each gene:
            for gene_id in self.family.variants[variant_id].get_genes(gene_annotation):
                if gene_id in self.gene_variants:
                    self.gene_variants[gene_id].append(variant_id)
                elif gene_id != '':
                    self.gene_variants[gene_id] = [variant_id]
        counter = 0
        
        # Check all compound candidates:
        
        for gene_id in self.gene_variants:
        # If there are more than one variant in a gene we start to look for compounds.
            if len(self.gene_variants[gene_id]) > 1:
                if len(self.gene_variants[gene_id]) < max_variants:
                    self.check_compound(gene_id, self.gene_variants[gene_id])
                else:
                    genes_with_many_variants.append(gene_id)
        print 'Genes with more than', max_variants, 'variants: ', len(genes_with_many_variants) 
        for variant_id in self.family.variants:
            if len(self.family.variants[variant_id].ar_comp_genes) == 0:
                self.family.variants[variant_id].ar_comp = False
            else:
                for gene, variant_list in self.family.variants[variant_id].ar_comp_genes.items():
                    if len(variant_list) > 0:
                        self.family.variants[variant_id].ar_comp = True
            self.family.variants[variant_id].check_models()
    
    def check_x_linked(self, variant_id):
        """Check if the variant follows the x linked patter of inheritance in this family."""
        for individual_id in self.family.individuals:
            individual = self.family.individuals[individual_id]
            # Get the genotype for this variant for this individual
            genotype = individual.get_genotype(variant_id)
            
            # The case where the individual is healthy
            if individual.phenotype == 1:
                #The case where the individual is a male
                if individual.sex == 1:
                    if genotype.has_variant:
                    # If the individual is healthy, male and have a variation it can not be x-linked.
                        self.family.variants[variant_id].x_linked = False
                        self.family.variants[variant_id].x_linked_dn = False
                
                #The case where the individual is a female
                elif individual.sex == 2:
                    # If the individual is HEALTHY, female and is homozygote alternative it can not be x - linked.
                    if genotype.homo_alt:
                        self.family.variants[variant_id].x_linked = False
                        self.family.variants[variant_id].x_linked_dn = False
            
            # The case when the individual is sick
            elif individual.phenotype == 2:
                #If the individual is sick and homozygote ref it can not be                     x-linked
                if genotype.homo_ref:
                    self.family.variants[variant_id].x_linked = False
                    self.family.variants[variant_id].x_linked_dn = False
                elif genotype.has_variant:
                    self.check_parents('x_linked', individual_id, variant_id)
            # Else if phenotype is unknown we can not say anything about the model
    
    def check_dominant(self, variant_id):
        """Check if the variant follows the dominant pattern in this family."""
        for individual_id in self.family.individuals: 
            # Check in all individuals what genotypes that are in the trio based of the individual picked.
            individual = self.family.individuals[individual_id]    
            genotype = individual.get_genotype(variant_id) # Get the genotype for this variant for this individual
            if individual.phenotype == 1:# The case where the individual is healthy
                if genotype.has_variant:
                    # If the individual is healthy and have a variation on one or both alleles it can not be dominant.
                    self.family.variants[variant_id].ad = False
                    self.family.variants[variant_id].ad_dn = False
            elif individual.phenotype == 2:
                # The case when the individual is sick
                if genotype.homo_ref:
                    self.family.variants[variant_id].ad = False
                    self.family.variants[variant_id].ad_dn = False
                else: 
                # Now the ind is sick and have a variant ≠ ref, check parents for de novo
                    self.check_parents('dominant', individual_id, variant_id)
                # Else if phenotype is unknown we can not say anything about the model
    
    def check_recessive(self, variant_id):
        """Check if the variant follows the autosomal recessive pattern in this family."""
        for individual_id in self.family.individuals: 
            individual = self.family.individuals[individual_id]    
            genotype = individual.get_genotype(variant_id)
            # The case where the individual is healthy:
            if individual.phenotype == 1:
            # If the individual is healthy and homozygote alt the model is broken.
                if genotype.homo_alt:
                    self.family.variants[variant_id].ar = False
                    self.family.variants[variant_id].ar_dn = False
            elif individual.phenotype == 2:# The case when the individual is sick
            # In the case of a sick individual it must be homozygote alternative for compound heterozygote to be true.
            # Also, we can not exclude the model if no call.
                if genotype.homo_ref or genotype.heterozygote:
                    self.family.variants[variant_id].ar = False
                    self.family.variants[variant_id].ar_dn = False
                else:
                #Models are followed but we need to check the parents to see if de novo is followed or not.
                    self.check_parents('recessive', individual_id, variant_id)
    
    def check_compound(self, gene_id, list_of_variants):
        """Check which variants in the list that follow the compound heterozygous model. 
        We need to go through all variants and sort them into their corresponding genes 
        to see which that are candidates for compound heterozygotes first. 
        The cheapest way to store them are in a hash table. After this we need to go
         through all pairs, if both variants of a pair is found in a healthy individual
          the pair is not a deleterious compound heterozygote."""
        true_variants = []
        false_variants = []
        
        def add_false_variants(false_variants, variant_list):
            """Add the pairs that where found to be not true"""
            for variant_id in variant_list:
                if (variant_list[0],variant_list[1]) not in false_variants:
                    false_variants.append((variant_list[0],variant_list[1]))
                if (variant_list[1],variant_list[0]) not in false_variants:
                    false_variants.append((variant_list[1],variant_list[0]))
            return false_variants
        
        def add_true_variants(true_variants, variant_list):
            """Add variants to the genes that was found as pairs."""
            for variant_id in variant_list:
                if (variant_list[0],variant_list[1]) not in true_variants:
                    true_variants.append((variant_list[0],variant_list[1]))
                if (variant_list[1],variant_list[0]) not in true_variants:
                    true_variants.append((variant_list[1],variant_list[0]))
            return true_variants
        
        def add_comp_variants(gene_id, variant_pair):
            """Annotate the gene with the true compounds to each of the variants."""
            if gene_id in self.family.variants[variant_pair[0]].ar_comp_genes:
                if variant_pair[1] not in self.family.variants[variant_pair[0]].ar_comp_genes[gene_id]:
                    self.family.variants[variant_pair[0]].ar_comp_genes[gene_id].append(variant_pair[1])
            else:
                self.family.variants[variant_pair[0]].ar_comp_genes[gene_id] = [variant_pair[1]]
            if gene_id in self.family.variants[variant_pair[1]].ar_comp_genes:
                if variant_pair[0] not in self.family.variants[variant_pair[1]].ar_comp_genes[gene_id]:
                    self.family.variants[variant_pair[1]].ar_comp_genes[gene_id].append(variant_pair[0])
            else:
                self.family.variants[variant_pair[1]].ar_comp_genes[gene_id] = [variant_pair[0]]
        
        # Returns a generator with all possible pairs for this individual:
        my_pairs = pair_generator.Pair_Generator(list_of_variants)
        for pair in my_pairs.generate_pairs():
            variant_1 = pair[0]
            variant_2 = pair[1]
            for individual_id in self.family.individuals: # Check in all individuals what genotypes that are in the trio based of the individual picked.
                individual = self.family.individuals[individual_id]
                genotype_1 = individual.get_genotype(variant_1)
                genotype_2 = individual.get_genotype(variant_2)
                # If an individual has any of the variants homo alt, healthy or sick, it will not count as a compound:
                if genotype_1.homo_alt or genotype_2.homo_alt:
                    false_variants = add_false_variants(false_variants, [variant_1, variant_2])
                # We only need to check is the individual is sick and have both variants
                # otherwise we know compound is NOT true:
                elif genotype_1.heterozygote and genotype_2.heterozygote:# A sick individual must have the heterozygous pair.
                    if individual.phenotype == 2:# The case where the individual is affected
                        mother_id = individual.mother
                        father_id = individual.father
                        if mother_id == '0': # Non existing mother
                            mother_genotype_1 = Genotype()
                            mother_genotype_2 = Genotype()
                            mother_phenotype = 0
                        else:
                            mother_genotype_1 = self.family.individuals[mother_id].get_genotype(variant_1)
                            mother_genotype_2 = self.family.individuals[mother_id].get_genotype(variant_2)
                            mother_phenotype = self.family.individuals[mother_id].phenotype
                        if father_id == '0': # Non existing father
                            father_genotype_1 = genotype.Genotype()
                            father_genotype_2 = genotype.Genotype()
                            father_phenotype = 0
                        else:
                            father_genotype_1 = self.family.individuals[father_id].get_genotype(variant_1)
                            father_genotype_2 = self.family.individuals[father_id].get_genotype(variant_2)
                            father_phenotype = self.family.individuals[father_id].phenotype
                        # If a parent has both variants and is unaffected it can not be a compound.
                        # This will change when we get the phasing information.
                        if ((mother_genotype_1.heterozygote and mother_genotype_2.heterozygote and mother_phenotype == 1) or (father_genotype_1.heterozygote and father_genotype_2.heterozygote and father_phenotype == 1)):
                            false_variants = add_false_variants(false_variants, [variant_1, variant_2])
                        else:
                            true_variants = add_true_variants(true_variants, [variant_1, variant_2])
        for variant_pair in true_variants:
            if variant_pair not in false_variants:
                add_comp_variants(gene_id, variant_pair) 
    
    def check_parents(self, model, individual_id, variant_id):
        """Check if information in the parents can tell us if model is de novo or not. Model in ['recessive', 'compound', 'dominant']."""
        individual = self.family.individuals[individual_id]
        sex = individual.sex
        mother_id = individual.mother
        father_id = individual.father
        if mother_id == '0': # Non existing mother
            mother_genotype = genotype.Genotype()# This will initialise a genptype like './.' and nocall is True
            mother_phenotype = 0# This is the code for unknown phenotype
        else:
            mother_genotype = self.family.individuals[mother_id].get_genotype(variant_id)
            mother_phenotype = self.family.individuals[mother_id].phenotype
        if father_id == '0': # Non existing father
            father_genotype = genotype.Genotype()
            father_phenotype = 0
        else:
            father_genotype = self.family.individuals[father_id].get_genotype(variant_id)
            father_phenotype = self.family.individuals[father_id].phenotype
        if model == 'recessive':
            # If any of the parents doesent exist de novo will be true as the model is specifyed now
            if mother_id != '0' and father_id != '0':
            # If both parents have the variant or if one of the parents are homozygote the de novo model is NOT followed, otherwise de novo is true.
                if (mother_genotype.homo_alt or father_genotype.homo_alt) or (mother_genotype.has_variant and father_genotype.has_variant):
                    self.family.variants[variant_id].ar_dn = False
            if self.family.variants[variant_id].ar_dn:# If de novo is true then the it is only de novo
                    self.family.variants[variant_id].ar = False
                    
        elif model == 'dominant':
        # If one of the parents have the variant on any form the de novo model is NOT followed.
            if mother_genotype.has_variant or father_genotype.has_variant:
                self.family.variants[variant_id].ad_dn = False
            if self.family.variants[variant_id].ad_dn:# If variant is ad de novo then it is not ad
                self.family.variants[variant_id].ad = False
                
        elif model == 'x_linked':
            if sex == 1: #If male
            
                if mother_genotype.has_variant or father_genotype.has_variant:
                    self.family.variants[variant_id].x_linked_dn = False
            elif sex == 2: #If female
                if individual.genotypes[variant_id].heterozygote:
                    if (mother_genotype.homo_alt or father_genotype.homo_alt) or (mother_genotype.has_variant and father_genotype.has_variant):
                        self.family.variants[variant_id].x_linked_dn = False
                elif individual.genotypes[variant_id].homo_alt:
                    if mother_genotype.has_variant or father_genotype.has_variant:
                        self.family.variants[variant_id].x_linked_dn = False            
    




def main():
    pass


if __name__ == '__main__':
    main()

