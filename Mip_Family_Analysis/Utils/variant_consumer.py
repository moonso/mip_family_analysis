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
    
    def __init__(self, task_queue, results_queue, lock, family, verbosity = False):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.family = family
        self.results_queue = results_queue
        self.verbosity = verbosity
        self.lock = lock
    
    def run(self):
        """Run the consuming"""
        proc_name = self.name
        while True:
            # A batch is a dictionary on the form {gene:{variant_id:(variant_dict, genotypes)}}
            next_batch = self.task_queue.get()
            if self.task_queue.empty():
                print 'Halllooooouuuu!!', proc_name
            if self.verbosity:
                if self.results_queue.full():
                    print 'Batch results queue Full!'
                if self.task_queue.full():
                    print 'Batch tasks queue Full!'
            if next_batch is None:
                self.task_queue.task_done()
                if self.verbosity:
                    print '%s: Exiting' % proc_name
                break
            # print '%s: %s' % (proc_name, next_batch)
            variant_batch = genetic_models.check_genetic_models(next_batch, self.family, self.verbosity, proc_name)
            fixed_variants = {}
            for gene in variant_batch:
                variant_dict = dict((variant_id, variant_info[0]) for variant_id, variant_info in variant_batch[gene].items())
                for variant_id, variant in variant_dict.items():
                    fixed_variants[variant_id] = variant
            fixed_variants = score_variants.score_variant(fixed_variants, self.family.models_of_inheritance)
            # with self.lock:
            #     for variant_id, variant in fixed_variants.items():
            #         print variant_id
            #         print variant
            # for variant_id in fixed_variants:
            #     if len(fixed_variants[variant_id].ar_comp_variants) > 0:
            #         for compound_id in fixed_variants[variant_id].ar_comp_variants:
            #             compound_score = fixed_variants[variant_id].rank_score + fixed_variants[compound_id].rank_score
            #             fixed_variants[variant_id].ar_comp_variants[compound_id] = compound_score
            
            self.results_queue.put(fixed_variants)
            self.task_queue.task_done()
        return
        
    

def main():
    pass

if __name__ == '__main__':
    main()