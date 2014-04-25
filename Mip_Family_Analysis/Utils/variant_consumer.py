#!/usr/bin/env python
# encoding: utf-8
"""
variant_consumer.py

Class that takes a list of objects and return all unordered pairs as a generator.

If only one object? Raise Exception
 
Created by Måns Magnusson on 2013-03-01.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import multiprocessing
from pprint import pprint as pp

from Mip_Family_Analysis.Models import genetic_models, score_variants

class VariantConsumer(multiprocessing.Process):
    """Yeilds all unordered pairs from a list of objects as tuples, like (obj_1, obj_2)"""
    
    def __init__(self, task_queue, results_queue, family, verbosity = False):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.family = family
        self.results_queue = results_queue
        self.verbosity = verbosity
    
    def fix_variants(self, variant_batch):
        """Merge the variants into one dictionary, make shure that the compounds are treated right."""
        fixed_variants = {}
        for feature in variant_batch:
            for variant_id in variant_batch[feature]:
                if variant_id in fixed_variants:
                    # We need to add compound information from different features
                    if len(variant_batch[feature][variant_id]['Compounds']) > 0:
                        fixed_variants[variant_id]['Compounds'] = (
                         dict(list(variant_batch[feature][variant_id]['Compounds'].items()) +
                                    list(fixed_variants[variant_id]['Compounds'].items())))
                else:
                    fixed_variants[variant_id] = variant_batch[feature][variant_id]
        
        return fixed_variants
    
    def make_print_version(self, variant_dict):
        """Get the variants ready for printing"""
        for variant_id in variant_dict:
            model_list = []
            compounds_list = []
            #Remove the 'Genotypes' post since we will not need them for now
            variant_dict[variant_id].pop('Genotypes', 0)
            variant_dict[variant_id]['Rank_score'] = variant_dict[variant_id]['Individual_rank_score']
            
            #If there are compounds we add the compound scores to each pair
            if len(variant_dict[variant_id]['Compounds']) > 0:
                # Put the compound scores
                for compound_id in variant_dict[variant_id]['Compounds']:
                    compound_score = (int(variant_dict[variant_id]['Individual_rank_score']) + 
                                         int(variant_dict[compound_id]['Individual_rank_score']))
                    variant_dict[variant_id]['Compounds'][compound_id] = compound_score
                    compounds_list.append(compound_id+'='+str(compound_score))  
            
            for model in variant_dict[variant_id]['Inheritance_model']:
                if variant_dict[variant_id]['Inheritance_model'][model]:
                    model_list.append(model)
            
            if len(model_list) == 0:
                variant_dict[variant_id]['Inheritance_model'] = 'NA'
            else:
                variant_dict[variant_id]['Inheritance_model'] = ':'.join(model_list)
            
            #If AR_compound is the only valid model we need to check the rank score:
            if 'AR_comp' in model_list:
                if len(model_list) == 1:
                    variant_dict[variant_id]['Rank_score'] = min(variant_dict[variant_id]['Individual_rank_score'], 
                            max([value for value in variant_dict[variant_id]['Compounds'].values()]))
                
            if len(compounds_list) > 0:
                variant_dict[variant_id]['Compounds'] = ':'.join(compounds_list)
            else:
                variant_dict[variant_id]['Compounds'] = '-'
            
            variant_dict[variant_id]['Individual_rank_score'] = str( 
                                    variant_dict[variant_id]['Individual_rank_score'])                        
            variant_dict[variant_id]['Rank_score'] = str( 
                                    variant_dict[variant_id]['Rank_score'])                        
            
        return
    
    def run(self):
        """Run the consuming"""
        proc_name = self.name
        if self.verbosity:
            print('%s Starting!' % proc_name)
        while True:
            # A batch is a dictionary on the form {gene_1:{variant_id:variant_dict}, gene_2:{variant_id:variant_dict}}
            next_batch = self.task_queue.get()
            # if self.verbosity:
            #     if self.results_queue.full():
            #         print('Batch results queue Full! %s' % proc_name)
            #     if self.task_queue.full():
            #         print('Variant queue full! %s' % proc_name)
            if next_batch is None:
                self.task_queue.task_done()
                if self.verbosity:
                    print('%s: Exiting' % proc_name)
                break
            genetic_models.check_genetic_models(next_batch, self.family, self.verbosity, proc_name = proc_name)
            fixed_variants = self.fix_variants(next_batch)
            score_variants.score_variant(fixed_variants, self.family.models_of_inheritance)
            self.make_print_version(fixed_variants)
            
            self.results_queue.put(fixed_variants)
            self.task_queue.task_done()
        return
        
    

def main():
    pass

if __name__ == '__main__':
    main()