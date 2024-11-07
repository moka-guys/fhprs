import unittest
from unittest.mock import patch, mock_open
from fh import PRS

class TestReadVCF(unittest.TestCase):

    @patch("builtins.open",side_effect=FileNotFoundError)
    def test_file_not_found(self,mock_open):
        # test that FileNotFoundError is raised and handled properly
        with self.assertRaises(FileNotFoundError) as context:
            PRS("non_existent.vcf")
        self.assertEqual(str(context.exception),"The specified VCF file 'non_existent.vcf' could not be found.")

        # ran test, OK

    @patch("builtins.open",side_effect=PermissionError)
    def test_permissio_error(self,mock_open):
        with self.assertRaises(PermissionError) as context:
            PRS("restricted.vcf")
        self.assertEqual(str(context.exception), "Permission denied when trying to open the VCF file 'restricted.vcf'.")

        # ran test, OK

    @patch("builtins.open", side_effect=Exception("Unknown error"))
    def test_generic_exception(self, mock_open):
        # Test that a generic Exception is raised and handled properly
        with self.assertRaises(Exception) as context:
            PRS("error.vcf")
        self.assertEqual(str(context.exception), "An unexpected error occurred while opening the VCF file: Unknown error")
