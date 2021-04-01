import sys
import os
import re
import csv
import vcf
from collections import defaultdict


'''generic flatMap function'''
flat_map = lambda f, xs: [y for ys in xs for y in f(ys)]


class PRS(object):
  """
  Class to easily calculate Polygenic Risk Score for FH.
  see Talmund et al 2013

  Attributes
  ----------
  SCORES : nested hashtable [location][genotype]
      Static risk scores for 12 SNPs from Talmud et al. 2013 (DOI:10.1016/S0140-6736(12)62127-8)
  RISKRANGES : tuple of tuples
      decile ranges for risk score and the associated risk level (according to Bristol)
  vcf_file : str
      VCF file path (uncompressed)
  locations : [str]
      SNP locations as extracted from the SCORES structure
  sample_index : int
      Index of sample for which score is calculated

  Methods
  -------
  _readGenotypes()
      Prints the animals name and what sound it makes
  scoreGenotypes()
      Calculates PRS and returns min and max value (which are identical if all SNPs have been genotyped)
  risk()
      Returns decile number for risk score and the associated risk term
  
  """
  SCORES = {
    "1:55504650": { #rs2479409 (PCSK9)
      "GG":0.104,
      "GA":0.052,
      "AA":0
    }, 
    "1:109818306":{ #rs629301 (CELSR2)
      "GG":0,
      "GT":0.15,
      "TT":0.3
    }, 
    "2:21263900":{ #rs1367117 (APOB)
      "GG":0,
      "GA":0.1,
      "AA":0.2
    }, 
    "2:44072576":{ #rs4299376 (ABCG8)
      "GG":0.142,
      "GT":0.071,
      "TT":0
    }, 
    "6:16127407":{ #rs3757354 (MYLIP) CORRECTED
      "CC":0.074,
      "CT":0.037,
      "TT":0
    }, 
    "6:26093141":{ #rs1800562 (HFE)
      "GG":0.114,
      "GA":0.057,
      "AA":0
    }, 
    "6:160578860":{ #rs1564348 (SLC22A1) CORRECTED
      "TT":0.028,
      "TC":0.014,
      "CC":0
    },
    "11:126243952":{ #rs11220462 (ST3GAL4)
      "GG":0,
      "GA":0.05,
      "AA":0.1
    }, 
    "14:24883887":{ #rs8017377 (NYNRIN)
      "GG":0,
      "GA":0.029,
      "AA":0.058
    },
    "19:11202306":{ #rs6511720 (LDLR)
      "GG":0.36,
      "GT":0.18,
      "TT":0
    },
    "19:45411941,19:45412079":{ #rs429358,rs7412
      "TTTT":-0.9,
      "TCTT":-0.4,
      "TTCT":-0.4,
      "TCCT":-0.2,
      "TTCC": 0,
      "CCTT": 0,
      "TCCC": 0.1,
      "CCCT": 0.1,
      "CCCC": 0.2
    }
  }
  RISKRANGES = (
		(-0.5,0.58,'low'),
		(0.58,0.73,'low'),
		(0.73,0.81,'low'),
		(0.81,0.88,'intermediate'),
		(0.88,0.93,'intermediate'),
		(0.93,0.98,'high'),
		(0.98,1.02,'high'),
		(1.02,1.08,'high'),
		(1.08,1.16,'high'),
		(1.16,1.46,'high')
  )

  def __init__(self, vcf_file, sample_index=0):
    """
    Creates class instance by extracting the Genotypes required for risk score calculation
    ----------
    vcf_file : str
        VCF file path (uncompressed)
    sample_index : str
        Index of sample for which score is calculated (defaults to first sample in VCF)
    """
    self.vcf_file = vcf_file
    self.locations = flat_map(lambda x: x.split(','), self.SCORES.keys()) 
    self.sample_index = sample_index
    self._readGenotypes()

  def _readGenotypes(self):
    """
    Extracts genotypes from VCF at position specified in SCORES hash table
    """
    self.genotypes = defaultdict(None)
    vcf_reader = vcf.Reader(open(self.vcf_file, 'r'))
    for record in vcf_reader:
      location = ':'.join([re.sub("^chr","",record.CHROM),str(record.POS)])
      if location in self.locations:
        # own GT extract function to ensure correct ordering of ALLELES
        gt_bases = []
        for num in sorted(record.samples[self.sample_index].gt_alleles):
          try:
            gt_bases.append(str(record.alleles[int(num)]).upper())
          except:
            break
        if len(gt_bases)==len(record.samples[self.sample_index].gt_alleles):
          try:
            self.genotypes[location] = ''.join(gt_bases)
          except:
            pass
        ## builtin methods dont guarantee Genotype order REF,ALT
        # genotype = record.samples[self.sample_index].gt_bases.upper()
        # self.genotypes[location] = gt_bases

  def scoreGenotypes(self):
    """
    Calculates risk score for extracted genotypes.
    Will return a range if not all SNPs were genotyped, or two indentical values otherwise
    """
    score_range = [0,0]
    for l,s in self.SCORES.items():
      locations = l.split(',')
      genotypes = list(map(lambda x: self.genotypes[x] if x in self.genotypes.keys() else None, locations))
      # lookup
      allele_scores = self.SCORES[l]
      try:
        alleles = ''.join(genotypes)
        score_range[0] += allele_scores[alleles]
        score_range[1] += allele_scores[alleles]
      except:
        score_range[0] += min(list(allele_scores.values()))
        score_range[1] += max(list(allele_scores.values()))
    return score_range
    
  def risk(self):
    """
    Returns risk decile and associated term
    """
    score_range = self.scoreGenotypes()
    risk_strings = []
    for score in score_range:
      decile, risk = 'NA', 'Out of range'
      for i,r in enumerate(self.RISKRANGES):
        if r[0] <= score < r[1]:
            decile = str(i+1)
            risk = r[2]
            break
      risk_strings.append('{} ({})'.format(decile,risk))
      if score_range[0] == score_range[1]:
        break
    return " - ".join(risk_strings)


if __name__ == "__main__":
  vcf_file = sys.argv[1]
  vcf = PRS(vcf_file)
  print(vcf.genotypes)
  print(vcf.scoreGenotypes())
  print(vcf.risk())
