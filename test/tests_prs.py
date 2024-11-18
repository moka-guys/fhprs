# This file will run when the docker image is being built - as specified in Dockerfile

import unittest
from unittest.mock import patch, mock_open
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from fh import PRS

class TestReadVCF(unittest.TestCase):

    @patch("builtins.open",side_effect=FileNotFoundError) # generating a mocked error
    def test_file_not_found(self, mock_open):
        # test that FileNotFoundError is raised and handled properly
        with self.assertRaises(FileNotFoundError) as context:
            PRS("non_existent.vcf")
        self.assertEqual(str(context.exception),"The specified VCF file 'non_existent.vcf' could not be found.")

        # ran test, OK

    @patch("builtins.open",side_effect=PermissionError) # generating a mocked error
    def test_permission_error(self, mock_open):
        with self.assertRaises(PermissionError) as context:
            PRS("restricted.vcf")
        self.assertEqual(str(context.exception), "Permission denied when trying to open the VCF file 'restricted.vcf'.")

        # ran test, OK

    @patch("builtins.open", side_effect=Exception("Unknown error")) # generating a mocked error
    def test_generic_exception(self, mock_open):
        # Test that a generic Exception is raised and handled properly
        with self.assertRaises(Exception) as context:
            PRS("error.vcf")
        self.assertEqual(str(context.exception), "An unexpected error occurred while opening the VCF file: Unknown error")

        # ran test, OK

    def test_read_complete_vcf(self):
        # Assume complete.vcf exists in the current test directory
        vcf_file = os.path.join(os.path.dirname(__file__), 'complete.vcf')
        
        # Create an instance of PRS with the complete VCF file
        prs = PRS(vcf_file)
        
        # Call _readGenotypes to process the file
        prs._readGenotypes()

        # Check if the genotypes are correctly extracted 
        self.assertIn("1:55038977", prs.genotypes)  # the first SNP in complete.vcf
        self.assertEqual(prs.genotypes["1:55038977"], "GA")  # The genotype of the first SNP in complete.vcf

        # ran test, OK

    def test_read_incomplete_vcf(self):
        # Assume incomplete.vcf exists in the current test directory
        vcf_file = os.path.join(os.path.dirname(__file__), 'incomplete.vcf')
        
        # Create an instance of PRS with the incomplete VCF file
        prs = PRS(vcf_file)
        
        # Call _readGenotypes to process the file
        prs._readGenotypes()

        # Check if the genotypes are correctly extracted 
        self.assertIn("1:55038977", prs.genotypes) # This SNP SHOULD be extracted
        self.assertEqual(prs.genotypes["1:55038977"], "GA")  # The genotype of the first SNP in complete.vcf
        self.assertNotIn("11:126374057", prs.genotypes)  # The missing SNP in incomplete.vcf - should NOT be present
        
        # ran test OK

    def test_score_complete_vcf(self):
        # Assuming complete.vcf exists in current test directory
        vcf_file = os.path.join(os.path.dirname(__file__), 'complete.vcf')
        
        # Create an instance of PRS with the complete VCF file
        prs = PRS(vcf_file)
        score_range = prs.scoreGenotypes()
        self.assertEqual(score_range[0],1.007) # checks that the lower score value is as we expect for complete.vcf from manual calculation
        self.assertEqual(score_range[1],1.007) # checks that the upper score value is the same as the lower score.

    def test_score_incomplete_vcf(self):
        # Assuming incomplete.vcf exists in current test directory
        vcf_file = os.path.join(os.path.dirname(__file__), 'incomplete.vcf')
        
        # Create an instance of PRS with the incomplete VCF file
        prs = PRS(vcf_file)
        score_range = prs.scoreGenotypes()
        self.assertEqual(score_range[0],0.907) # checks that the lower score value is as we expect for incomplete.vcf from manual calculation
        self.assertEqual(score_range[1],1.007) # checks that the upper score differs and is as we expect from manual calc.

    def test_risk_complete_vcf(self):
        # Assuming complete.vcf exists in current test directory
        vcf_file = os.path.join(os.path.dirname(__file__), 'complete.vcf')
        
        # Create an instance of PRS with the complete VCF file
        prs = PRS(vcf_file)
        risk = prs.risk()
        self.assertEqual(risk,["high-7"]) # checks that the risk category and decile are as expected for complete.vcf

    def test_risk_incomplete_vcf(self):
        # Assuming incomplete.vcf exists in current test directory
        vcf_file = os.path.join(os.path.dirname(__file__), 'incomplete.vcf')
        
        # Create an instance of PRS with the incomplete VCF file
        prs = PRS(vcf_file)
        risk = prs.risk()
        self.assertEqual(risk,["intermediate-5","high-7"]) # checks that the risk category and decile are as expected for incomplete.vcf

    def test_empty_vcf(self):
        # Assuming empty.vcf exists in current test directory
        vcf_file = os.path.join(os.path.dirname(__file__), 'empty.vcf')

        # Test that a ValueError Exception is raised and handled properly
        with self.assertRaises(Exception) as context:
            prs = PRS(vcf_file)
        self.assertEqual(str(context.exception), f"The data in the VCF file '{vcf_file}' is either missing or corrupted. Please check file.")

    def test_headers_only_vcf(self):
        # Assuming headers_only.vcf exists in current test directory
        vcf_file = os.path.join(os.path.dirname(__file__), 'headers_only.vcf')

        # Test that a ValueError Exception is raised and handled properly
        with self.assertRaises(Exception) as context:
            prs = PRS(vcf_file)
            
        self.assertEqual(str(context.exception), f"The data in the VCF file '{vcf_file}' is either missing or corrupted. Please check file.")

    def invalid_vcf(self):
        # Assuming invalid.vcf exists in current test directory
        vcf_file = os.path.join(os.path.dirname(__file__), 'invalid.vcf')

        # Test that a ValueError Exception is raised and handled properly
        with self.assertRaises(Exception) as context:
            prs = PRS(vcf_file)
            
        self.assertEqual(str(context.exception), f"The data in the VCF file '{vcf_file}' is either missing or corrupted. Please check file.")

 
        


 
        





# test a vcf with totally wild values idk

        
