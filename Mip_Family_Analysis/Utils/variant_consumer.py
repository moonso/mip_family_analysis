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
            genetic_models.check_genetic_models(next_batch, self.family, self.verbosity, proc_name)
            fixed_variants = {}
            # # Make shore we only have one copy of each variant:
            for feature in next_batch:
                #Make one dictionary for each gene:
                for variant_id in next_batch[feature]:                        
                    if variant_id in fixed_variants:
                        # Add the compounds together
                        if len(next_batch[feature][variant_id]['Compounds']) > 0:
                            fixed_variants[variant_id]['Compounds'] = dict( 
                                next_batch[feature][variant_id]['Compounds'].items() +
                                fixed_variants[variant_id]['Compounds'].items())
                            #We do not want reference to itself as a compound:
                            fixed_variants[variant_id]['Compounds'].pop(variant_id, 0)
                            
                    else:
                        fixed_variants[variant_id] = next_batch[feature][variant_id]
            next_batch = {} # Free memory?
            score_variants.score_variant(fixed_variants, self.family.models_of_inheritance)
            
            # Fix the compound scores and make ready for printing:
            
            for variant_id in fixed_variants:
                # Remove 'Genotypes' since we don't need them anymore:
                fixed_variants[variant_id].pop('Genotypes', 0)
                # Set the rank score to individual rank score for now:
                fixed_variants[variant_id]['Rank_score'] = fixed_variants[variant_id]['Individual_rank_score']
                # if fixed_variants[variant_id]['Inheritance_model']['AR_compound'] and len(fixed_variants[variant_id]['Compounds']):
                #     print('disaster')
                #If there are copounds we add the compound scores to each pair
                if len(fixed_variants[variant_id]['Compounds']) > 0:
                    # Put the compound scores
                    for compound_id in fixed_variants[variant_id]['Compounds']:
                        compound_score = (fixed_variants[variant_id]['Individual_rank_score'] + 
                                             fixed_variants[compound_id]['Individual_rank_score'])
                        fixed_variants[variant_id]['Compounds'][compound_id] = compound_score
            
            for variant_id in fixed_variants:
                
                model_list = []
                for model in fixed_variants[variant_id]['Inheritance_model']:
                    if fixed_variants[variant_id]['Inheritance_model'][model]:
                        model_list.append(model)
                
                if len(model_list) == 0:
                    fixed_variants[variant_id]['Inheritance_model'] = 'NA'
                else:
                    fixed_variants[variant_id]['Inheritance_model'] = ':'.join(model_list)
                
                compounds_list = []                
                if 'AR_compound' in model_list:
                    # print('AR_compound %s' % str(fixed_variants[variant_id]['Compounds']))
                    if len(model_list) == 1:
                        fixed_variants[variant_id]['Rank_score'] = min(fixed_variants[variant_id]['Individual_rank_score'], 
                                max([value for value in fixed_variants[variant_id]['Compounds'].values()]))
                    
                    compounds_list = [comp_id + '=' + str(comp_value) for comp_id, comp_value in
                                         fixed_variants[variant_id]['Compounds'].items()]
                if len(compounds_list) > 0:
                    fixed_variants[variant_id]['Compounds'] = ':'.join(compounds_list)
                else:
                    fixed_variants[variant_id]['Compounds'] = '-'
                
                fixed_variants[variant_id]['Individual_rank_score'] = str( 
                                        fixed_variants[variant_id]['Individual_rank_score'])                        
                fixed_variants[variant_id]['Rank_score'] = str( 
                                        fixed_variants[variant_id]['Rank_score'])                        
                          
            self.results_queue.put(fixed_variants)
            self.task_queue.task_done()
        return
        
    

def main():
    pass

if __name__ == '__main__':
    main()