#!/usr/bin/env python
# encoding: utf-8
"""
variant_parser.py


Parse a file with variant info, this can be a .vcf file, an annotated annovar file, 
a annotated .txt cmms file, a annotated .txt cmms_ranked .

Create a variant objects and a dictionary with individuals that have a dictionary with genotypes for each variant.

Created by Måns Magnusson on 2013-01-17.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import argparse
import shelve
from Mip_Family_Analysis.Variants import genetic_variant, genotype
from collections import OrderedDict

class VariantParser(object):
    """Parses a file with family info and creates a family object with individuals."""
    def __init__(self, infile, file_type):
       super(VariantParser, self).__init__()
       self.variant_type = file_type
       self.chrom_shelves = OrderedDict() # ordered dict with {<chr>: <path_to_shelve_with_variants>}
       self.header_lines = []
       self.individuals = []
       self.metadata = []
       chrom_change = False
       beginning = True
       
       with open(infile, 'r') as f:
           line_count = 0
           for line in f:
               line = line.rstrip()
               line_count += 1
               if chrom_change:
                   # Close the shelve since we are at a new chromosome:
                   current_shelve.close()
                   # Make a new shelve:
                   shelve_name = 'chrom_' + new_chrom + '.shelve'
                   current_shelve = shelve.open(shelve_name)
                   # Add the filename to our dictionary:
                   self.chrom_shelves[new_chrom] = shelve_name
                   # Add the first variant from the new chromosome:
                   current_shelve[my_variant.variant_id] = my_variant
                   # Restart:
                   chrom_change = False
                   current_chrom = new_chrom
                   
               if line[0] != '#' and len(line) > 1:
                   my_variant = self.cmms_variant(line)
                   
                   if beginning:
                       # Init
                       current_chrom = my_variant.chr
                       new_chrom = my_variant.chr
                       shelve_name = 'chrom_' + new_chrom + '.shelve'
                       current_shelve = shelve.open(shelve_name)
                       self.chrom_shelves[new_chrom] = shelve_name
                       current_shelve[my_variant.variant_id] = my_variant
                       beginning = False
                   else:    
                       new_chrom = my_variant.chr
                       if new_chrom == current_chrom:
                           current_shelve[my_variant.variant_id] = my_variant
                       else:
                           chrom_change = True
               elif line[:1] != '##':
                   # If necesary we can write something to parse the headers.
                   self.header_lines = line[1:].split()
                   self.metadata.append(line)
               elif line[:2] == '##':
                   self.metadata.append(line)
       current_shelve[my_variant.variant_id] = my_variant
       current_shelve.close()
    
    def cmms_variant(self, variant_line):
        """Returns a Variant objekt from the cmms variant format. Creates a list with Genotype objects that becomes a member of the variant."""        
        variant_line.rstrip()
        variant_info = OrderedDict()
        individual_genotypes = {} #DICT with {ind_id:{<genotype>}}
        counter = 0
        variant_line = variant_line.split('\t')
        
        for entry in range(len(variant_line)):
            
            if 'IDN' in self.header_lines[entry]:
                # Looks like IDN:11-1-2A
                individual = self.header_lines[entry].split(':')[-1]
                if individual not in self.individuals:
                    self.individuals.append(individual)
            variant_info[self.header_lines[entry]] = variant_line[entry]
                
        counter += 1
        chrom = variant_info['Chromosome']
        start = variant_info['Variant_start']
        stop = variant_info['Variant_stop']
        alternative = variant_info['Alternative_allele']
        reference = variant_info['Reference_allele']
        identity = variant_info['Dbsnp_rs_nr']
        my_variant = genetic_variant.Variant(chrom, start, stop, reference, alternative, identity, variant_info)
        
        # Add the genotypes to variant:
        
        for individual in self.individuals:
            genotype_arguments = {} # args for genotype class
            key = 'IDN:' + individual
            # gt_info looks like 11-1-2A:GT=0/1:PL=32,3,2:...
            for gt_info in variant_info[key].split(':')[1:]:
                value_pair = gt_info.split('=')
                genotype_arguments[value_pair[0]] = value_pair[-1]
            my_genotype = genotype.Genotype(GT=genotype_arguments.get('GT','./.'), AD=genotype_arguments.get('AD','.,.'), DP=genotype_arguments.get('DP','0'), GQ=genotype_arguments.get('GQ','0'))
            my_variant.genotypes[individual] = my_genotype
        return my_variant
        


def main():
    parser = argparse.ArgumentParser(description="Parse different kind of pedigree files.")
    parser.add_argument('variant_file', type=str, nargs=1 , help='A file with variant information.')
    args = parser.parse_args()
    infile = args.variant_file[0]
    file_type = 'cmms'
    my_parser = VariantParser(infile, file_type)
    number_of_variants = 0
    for chrom in my_parser.chrom_shelves:
        current_db = shelve.open(my_parser.chrom_shelves[chrom])
        print chrom, len(current_db)
        for variant in current_db:
            print current_db[variant]
        number_of_variants += len(current_db)
        current_db.close()
        os.remove(my_parser.chrom_shelves[chrom])
    print 'Number of shelved:', number_of_variants


if __name__ == '__main__':
    main()
