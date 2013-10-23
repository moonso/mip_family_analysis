#!/usr/bin/env python
# encoding: utf-8
"""
test_get_variant_vcf.py

Created by Måns Magnusson on 2013-06-18.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import StringIO
from Mip_Family_Analysis.Variants import variant_parser

class TestVariant(object):
    """Test class for testing how the get_variant functions behave."""
    
    def setup_class(self):
        """Setup some variants"""
        vcf_header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t6405\t6404\n"
        my_simple_variant = (
            "1\t209949474\trs41303065\tG\tA\t4448.64\tPASS\tAC=1;AF=0.50;AN=2;BaseQRankSum=1.391;DB;DP=74"
            ";Dels=0.00;FS=13.174;HRun=1;HaplotypeScore=2.7803;InbreedingCoeff=-0.1333;MQ=58.10;MQ0=0"
            ";MQRankSum=2.354;QD=13.65;ReadPosRankSum=0.642;VQSLOD=5.7742;culprit=HaplotypeScore"
            ";set=variant2\tGT:AD:DP:GQ:PL\t0/1:33,41:74:99:1240,0,984\t./.\n")
        my_double_variant = (
            "1\t209969979\trs5780538\tA\tC,T\t1556.66\tPASS\tAC=2,0;AF=0.50,0.00;"
            "AN=4;DB;DP=25;FS=0.000;HRun=0;MQ0=0;set=Intersection\tGT:AD:DP:GQ:PL\t"
            "0/1:10,1:11:48.24\t1/2:20,4:24:43.94\n")
        path_to_test_variant = os.path.dirname(os.path.abspath(__file__))
        self.vcf_file = path_to_test_variant + 'simple_variant.vcf'
        f = open(self.vcf_file, 'w')
        f.write(vcf_header)
        f.write(my_simple_variant)
        f.write(my_double_variant)
        f.close()
        self.my_variant_parser = variant_parser.VariantParser(self.vcf_file, 'vcf')
        
    def test_variants(self):
        """Check if the variant is stored with the correct id"""
        assert len(self.my_variant_parser.variants) == 3
        assert self.my_variant_parser.individuals['6405']['1_209949474_G_A'].genotype == '0/1'
        assert self.my_variant_parser.individuals['6404']['1_209949474_G_A'].genotype == './.'
        assert self.my_variant_parser.individuals['6404']['1_209969979_A_C'].genotype == '1/2'

    def teardown_class(self):
        """Teardown the objects"""
        os.remove(self.vcf_file)
    

    # def test_location(self):
    #     """Test if the genomic positions are correct for the variant."""
    #     print self.my_variant.get_variant()
    #     assert self.my_variant.get_variant() == {'chrom': '1', 'start': '19202926', 'stop': '19202927', 'ref': 'T', 'alt': 'C', 'identity': 'rs2230706'}
    # 



def main():
    pass


if __name__ == '__main__':
    main()

