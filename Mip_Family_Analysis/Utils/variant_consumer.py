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
from mip_family_analysis.models import genetic_models

class VariantConsumer(multiprocessing.Process):
    """Yeilds all unordered pairs from a list of objects as tuples, like (obj_1, obj_2)"""
    
    def __init__(self, lock, task_queue, results_queue, family):
        multiprocessing.Process.__init__(self)
        self.batch_queue = task_queue
        self.results_queue = results_queue
        self.family = family
        self.lock = lock
        print 'HALLÅÅÅ!!!'
    
    def run(self):
        """Run the consuming"""
        proc_name = self.name
        while True:
            # A batch is a dictionary on the form {gene:{variant_id:variant}}
            next_batch = self.batch_queue.get()
            if next_batch is None:
                print '%s: Exiting' % proc_name
                self.batch_queue.task_done()
                break
            # print '%s: %s' % (proc_name, next_batch)
            genetic_models.check_genetic_models(self.family, next_batch[0], compound_check=next_batch[1])
            fixed_variants = {}
            for gene, variant_dict in next_batch.items():
                for variant_id, variant in variant_dict.items():
                    fixed_variant[variant_id] = variant
            with self.lock:
                for variant_id in fixed_variants:
                    print fixed_variant[variant_id]
            self.batch_queue.task_done()
            # self.results_queue.put(answer)
        return
        
    

def main():
    pass

if __name__ == '__main__':
    main()

