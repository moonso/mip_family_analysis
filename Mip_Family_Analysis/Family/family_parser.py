#!/usr/bin/env python
# encoding: utf-8
"""
get_family.py


Parse a file with family info, this can be a .ped file, a .fam, a .txt(CMMS style) 
file or a .txt(Broad style) file.
.ped and .fam always have 6 columns, these are

Family_ID Individual_ID Paternal_ID Maternal_ID 
Sex(1=male; 2=female; other=unknown) 
Phenotype(-9 missing, 0 missing, 1 unaffected, 2 affected)

The .txt allways have two columns on the form 

Individual_ID key=value

Where keys can be fid(=Family Id), mom, dad, sex, phenotype


If a pedigree file includes information about several families this must be taken care
 of by the parser by creating several family objects and then add information about the
  familymembers to theese familys. 

Create a family object and its family members from different types of input file
Created by Måns Magnusson on 2013-01-17.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import argparse
from Mip_Family_Analysis.Family import family, individual

class FamilyParser(object):
    """Parses a file with family info and creates a family object with individuals."""
    def __init__(self, infile, family_type):
        super(FamilyParser, self).__init__()
        self.family_type = family_type
        self.families = {}
        with open(infile, 'r') as f:
            line_count = 0
            for line in f:
                line_count += 1
                if line[0] != '#' and len(line) > 1:
                    individual_line = line.rstrip()
                    if family_type == 'cmms':
                        self.cmms_parser(individual_line, line_count)
                    elif family_type == 'fam':
                        self.ped_parser(individual_line, line_count)
                    elif family_type == 'ped':
                        self.ped_parser(individual_line, line_count)
                    elif family_type == 'broad':
                        self.broad_parser(individual_line, line_count)
                    my_counter = 1
    
    def ped_parser(self, individual_line, line_count):
        """Parse a .ped ped file."""
        if len(individual_line) < 6:
            print 'Malformed ped file, to few entrys on line ' + str(line_count)
            print individual_line
            sys.exit()
        line = individual_line.split()
        fam_id = line[0]
        if fam_id not in self.families:
            new_family = family.Family(family_id = fam_id, individuals = {})
            self.families[fam_id] = new_family
        ind = line[1]
        father = line[2]
        mother = line[3]
        sex = line[4]
        phenotype = line[5]
        my_individual = individual.Individual(ind, fam_id, mother, father, sex, phenotype)
        self.families[my_individual.family].add_individual(my_individual)
    
    def broad_parser(self, individual_line, line_count):
        """Parse a broad style .txt family file."""
        
        def parse_entry(entry):
            """Parses a entry from the .txt type of family file."""
            key        = ''
            value    = ''
            i        = 0
            while entry[i] != '=':
                key    += entry[i]
                i    += 1
                value    = entry[i+1:]
            return key, value
        
        if len(line) > 1:
            line =     individual_line.split()
            info = {} #Dictionary with info of individual
            info[ind] = line[0]
            other_info = line[1].split(';')
            for entry in other_info:
                key, value = parse_entry(entry)
                if key in ['fid','family', 'family_id']:
                    info[family] = value
                elif key in ['mom', 'mother']:
                    info[mother] = value
                elif key in ['dad', 'father']:
                    info[father] = value
                elif key in ['sex', 'gender']:
                    info[sex] = value
                elif key in ['phenotype', 'pheno']:
                    info[phenotype] = value
                else:
                    print 'Unknown parameter', key, value
                    sys.exit()
            if fam_id not in self.families:
                self.families['fam_id'] = family.Family(fam_id)
            self.families['fam_id'].add_individual(ind=ind_id, father=father, mother=mother, sex=sex, phenotype=pheno)
    

def main():
    parser = argparse.ArgumentParser(description="Parse different kind of pedigree files.")
    parser.add_argument('pedigree_file', type=str, nargs=1 , help='A file with pedigree information.')
    parser.add_argument('-ped', '--ped', action="store_true", help='Pedigree file is in ped format.')
    parser.add_argument('-fam', '--fam', action="store_true", help='Pedigree file is in fam format.')
    parser.add_argument('-broad', '--broad', action="store_true", help='Pedigree file is in broad format.')
    args = parser.parse_args()
    infile = args.pedigree_file[0]
    file_type = 'ped'
    if args.ped:
        file_type = 'ped'
    if args.fam:
        file_type = 'ped'
    if args.broad:
        file_type = 'broad'
    my_parser = FamilyParser(infile, file_type)
    print my_parser.families
    for family in my_parser.families:
        print my_parser.families[family]


if __name__ == '__main__':
    main()
