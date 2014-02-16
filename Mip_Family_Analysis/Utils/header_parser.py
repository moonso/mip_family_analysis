#!/usr/bin/env python
# encoding: utf-8
"""
header_parser.py


Parse the header of a variant file.

Create a variant objects and a dictionary with individuals that have a dictionary with genotypes for each variant.

Created by Måns Magnusson on 2013-01-17.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import re
import argparse
if sys.version_info < (2, 7):
    from ordereddict import OrderedDict
else:
    from collections import OrderedDict

from pprint import pprint as pp

### TODO make a proper vcf parser ###


class HeaderParser(object):
    """Parses a file with family info and creates a family object with individuals."""
    def __init__(self, infile):
        super(HeaderParser, self).__init__()
        self.metadata=OrderedDict()
        self.header=[]
        self.line_counter = 0
        self.individuals = []
        self.metadata_pattern = re.compile(r'''\#\#COLUMNNAME="(?P<colname>[^"]*)"
            (?P<info>.*)''', re.VERBOSE)
        with open(infile, 'rb') as f:
            for line in f:
                self.line_counter += 1
                line = line.rstrip()
                if line.startswith('##'):
                    match = self.metadata_pattern.match(line)
                    if not match:
                        raise SyntaxError("One of the metadata lines is malformed: %s" % line)
                    matches = [match.group('colname'), match.group('info')]
                    self.metadata[match.group('colname')] = line
                elif line.startswith('#'):
                    self.header = line[1:].split('\t')
                    for entry in self.header:
                        if entry[:3] == 'IDN':
                            self.individuals.append(entry.split(':')[1])
                        else:
                            self.check_header(entry)
                else:
                    break
        
    def add_metadata(self, column_name, data_type=None, version=None, description=None, dbname=None, delimiter='\t'):
        """Add metadata info to the header."""
        data_line = '##COLUMNAME='+'"'+ column_name +'"'
        if column_name not in self.metadata:
            if data_type:
                if data_type not in ['Float', 'String', 'Integer']:
                    raise SyntaxError("Type must be 'Float', 'String' or 'Integer'. You tried: %s" % data_type)
                data_line += delimiter + 'TYPE="' + data_type + '"'
            if version:
                data_line += delimiter + 'VERSION="' + version + '"'
            if description:
                data_line += delimiter + 'DESCRIPTION="' + description + '"'
            if dbname:
                data_line += delimiter + 'SCOUTHEADER="' + dbname + '"'
            self.metadata.pop(column_name, 0)
            self.metadata[column_name] = data_line
        return
    
    def add_header(self, header_name):
        """Add an entry to the header line. The entry must be specified among the metadata lines first."""
        self.check_header(header_name)
        if header_name not in self.header:
            self.header.append(header_name)
        return
    
    def get_headers_for_print(self):
        """Returns a list with the metadata lines on correct format."""
        lines_for_print = []
        for header in self.metadata:
            lines_for_print.append(self.metadata[header])
        lines_for_print.append('\t'.join(self.header))
        lines_for_print[-1] = '#' + lines_for_print[-1]
        return lines_for_print
    
    def check_header(self, entry):
        """All entrys in the header must be specified in the metadata lines."""
        if entry not in self.metadata:
            raise SyntaxError("Header entry must be described in the metadata lines. Entry: %s is not in metadata." % entry)
    


def main():
    parser = argparse.ArgumentParser(description="Parse different kind of pedigree files.")
    parser.add_argument('variant_file', type=str, nargs=1 , help='A file with variant information.')
    args = parser.parse_args()
    infile = args.variant_file[0]
    my_parser = HeaderParser(infile)
    my_parser.add_metadata('Individual_rank_score', data_type='Integer', description='This is the correct rank score if the variant only follows the AR_comp model.', dbname='Individual Rank Score')
    my_parser.add_header('Individual_rank_score')
    try:
        my_parser.add_header('Rank_score')
    except SyntaxError:
        print 'Must specify the new header in metadata first'
    else:
        for line in my_parser.get_headers_for_print():
            print line
        print my_parser.header
    # print '\t'.join(my_parser.header)
    # print my_parser.line_counter
    # print my_parser.individuals
    # for individual in my_parser.individuals:
    #     for genotype in my_parser.individuals[individual]:
    #         print individual, genotype, my_parser.individuals[individual][genotype]
    


if __name__ == '__main__':
    main()
