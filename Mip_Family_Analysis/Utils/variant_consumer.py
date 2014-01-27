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

from Mip_Family_Analysis.Models import genetic_models, score_variants

class VariantConsumer(multiprocessing.Process):
    """Yeilds all unordered pairs from a list of objects as tuples, like (obj_1, obj_2)"""
    
    def __init__(self, lock, task_queue, results_queue, family, verbosity = False):
        multiprocessing.Process.__init__(self)
        self.batch_queue = task_queue
        self.family = family
        self.lock = lock
        self.results = results_queue
        self.verbosity = verbosity
        
    
    def run(self):
        """Run the consuming"""
        proc_name = self.name
        while True:
            # A batch is a dictionary on the form {gene:{variant_id:variant}}
            next_batch = self.batch_queue.get()
            if next_batch is None:
                if self.verbosity:
                    print '%s: Exiting' % proc_name
                self.results.put('Done')
                self.batch_queue.task_done()
                break
            # print '%s: %s' % (proc_name, next_batch)
            genetic_models.check_genetic_models(next_batch, self.family)
            fixed_variants = {}
            for gene, variant_dict in next_batch.items():
                for variant_id, variant in variant_dict.items():
                    fixed_variants[variant_id] = variant
            
            for variant_id in fixed_variants:
                score_variants.score_variant(fixed_variants[variant_id], self.family.models_of_inheritance)
            
            for variant_id in fixed_variants:
                if len(fixed_variants[variant_id].ar_comp_variants) > 0:
                    for compound_id in fixed_variants[variant_id].ar_comp_variants:
                        compound_score = fixed_variants[variant_id].rank_score + fixed_variants[compound_id].rank_score
                        fixed_variants[variant_id].ar_comp_variants[compound_id] = compound_score            

            self.results.put(fixed_variants)
            self.batch_queue.task_done()
        return
        
    

def main():
    pass

if __name__ == '__main__':
    main()

