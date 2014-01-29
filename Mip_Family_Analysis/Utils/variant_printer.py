#!/usr/bin/env python
# encoding: utf-8
"""
variant_printer.py


Print the variants of a results queue to a file.


Created by Måns Magnusson on 2013-01-17.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import multiprocessing



class VariantPrinter(multiprocessing.Process):
    """docstring for VariantPrinter"""
    def __init__(self, results, outfile, number_of_consumers, verbosity=False):
        multiprocessing.Process.__init__(self)
        self.results = results
        self.outfile = outfile
        self.verbosity = verbosity
        print self.verbosity
        self.number_of_consumers = number_of_consumers
    
    def run(self):
        """Starts the printing"""
        # Print the results to a temporary file:
        number_of_finished = 0
        proc_name = self.name
        if self.verbosity:
            print proc_name ,'starting!'
        with open(self.outfile, 'w') as file_handle:
            while True:
                next_result = self.results.get()
                if type(next_result) == type('a'):
                    if next_result == 'Done':
                        number_of_finished += 1
                    if number_of_finished == self.number_of_consumers:
                        if self.verbosity:
                            print 'All variants printed!'
                        break
                else:
                    for variant_id in next_result:
                        file_handle.write('\t'.join(next_result[variant_id].get_cmms_variant())+'\n')
        return

def main():
    pass

if __name__ == '__main__':
    main()
