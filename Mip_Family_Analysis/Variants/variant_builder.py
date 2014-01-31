#!/usr/bin/env python
# encoding: utf-8
"""
variant_builder.py

Class that takes a variant line and put a variant object in a joined queue for future processing.

 
Created by Måns Magnusson on 2013-03-01.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
from datetime import datetime
import multiprocessing
if sys.version_info < (2, 7):
    from ordereddict import OrderedDict
else:
    from collections import OrderedDict

from Mip_Family_Analysis.Variants import genetic_variant, genotype
from Mip_Family_Analysis.Utils import get_genes

class VariantBuilder(multiprocessing.Process):
    """Yeilds all unordered pairs from a list of objects as tuples, like (obj_1, obj_2)"""
    
    def __init__(self, variant_file, variant_queue, head, verbosity = False):
        multiprocessing.Process.__init__(self)
        self.task_queue = variant_queue
        self.variant_file = variant_file
        # self.results_queue = task_queue , task_queue        
        self.individuals = head.individuals
        self.header_line = head.header
        self.verbosity = verbosity
        
    
    def run(self):
        """Run the consuming"""
        proc_name = self.name
        current_chrom = '1'
        beginning = True
        batch = {}
        start_parsing = datetime.now()
        chrom_start = start_parsing
        while True:
            # A batch list with the info of a variant
            next_variant_line = self.task_queue.get()
            # print self.results_queue.qsize()
            # if self.results_queue.full():
            #     print 'Variant queue Full!'
            if self.task_queue.full():
                print 'Variant tasks queue Full!'
            if next_variant_line is None:
                if self.verbosity:
                    print proc_name, 'Exiting!'
                    # print batch
                # self.results_queue.put(batch)
                # self.results_queue.put('Done')
                break
                
            variant = self.cmms_variant(next_variant_line)
            
            new_genes = variant.genes
            # If we look at the first variant, setup boundary conditions:
            if beginning:
                current_genes = new_genes
                beginning = False
                variant_batch = self.add_variant(batch, variant) # Add variant batch
            else:
                send = True
            
            #Check if we are in a space between genes:
                with self.lock:
                    print current_genes, new_genes
                if len(new_genes) == 0:
                    if len(current_genes) == 0:
                        send = False
            #If not check if we are in a consecutive region
                elif len(set.intersection(set(new_genes),set(current_genes))) > 0:
                    send = False
                if send:
                    # If there is an intergenetic region we do not look at the compounds.
                    # The tasks are tuples like (variant_list, bool(if compounds))
                    # print batch
                    # self.results_queue.put(batch)
                    current_genes = new_genes
                    batch = self.add_variant({}, variant)
                else:
                    current_genes = list(set(current_genes) | set(new_genes))
                    batch = self.add_variant(batch, variant) # Add variant batch
        return
    
    def add_variant(self, batch, variant):
        """Adds the variant to the proper gene(s) in the batch."""
        if len(variant.genes) == 0:
            if len(batch) == 0:
                batch['-'] = {variant.variant_id:variant}
            else:
                batch['-'][variant.variant_id] = variant
        for gene in variant.genes:
            if gene in batch:
                batch[gene][variant.variant_id] = variant
            else:
                batch[gene] = {variant.variant_id:variant}
        return batch
    
    
    def cmms_variant(self, splitted_variant_line):
        """Returns a variant object in the cmms format."""
        
        variant_info = OrderedDict()
        individual_genotypes = {} #DICT with {ind_id:{<genotype>}}
        counter = 0
        ensemble_entry = splitted_variant_line[5]
        hgnc_entry = splitted_variant_line[6]
        
        # These must be parsed separately
        hgnc_genes = get_genes.get_genes(hgnc_entry, 'HGNC')
        ensemble_genes = get_genes.get_genes(ensemble_entry, 'Ensemble')
    
        for entry in range(len(splitted_variant_line)):
            
            if 'IDN' in self.header_line[entry]:
                # Looks like IDN:11-1-2A
                individual = self.header_line[entry].split(':')[-1]
                if individual not in self.individuals:
                    raise SyntaxError('One of the individuals in the variant file \
                                is not in the ped file: %s' % individual)
            variant_info[self.header_line[entry]] = splitted_variant_line[entry]
                
        chrom = variant_info['Chromosome']
        start = variant_info['Variant_start']
        stop = variant_info['Variant_stop']
        alternative = variant_info['Alternative_allele']
        reference = variant_info['Reference_allele']
        identity = variant_info['Dbsnp_rs_nr']
        my_variant = genetic_variant.Variant(chrom , start, stop, reference, 
                        alternative, identity, genes=hgnc_genes, all_info=variant_info)
    
        # Add the genotypes to variant:
        
        for individual in self.individuals:
            genotype_arguments = {} # args for genotype class
            key = 'IDN:' + individual
            # gt_info looks like 11-1-2A:GT=0/1:PL=32,3,2:...
            for gt_info in variant_info[key].split(':')[1:]:
                value_pair = gt_info.split('=')
                genotype_arguments[value_pair[0]] = value_pair[-1]
            my_genotype = genotype.Genotype(GT=genotype_arguments.get('GT','./.'), 
                                            AD=genotype_arguments.get('AD','.,.'), 
                                            DP=genotype_arguments.get('DP','0'), 
                                            GQ=genotype_arguments.get('GQ','0'))
            my_variant.genotypes[individual] = my_genotype
   
        return my_variant
    
        
    

def main():
    pass

if __name__ == '__main__':
    main()

