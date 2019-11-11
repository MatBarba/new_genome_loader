#!env python3

from BCBio import GFF

import eHive
import json
import magic
import os
import subprocess as sp
import sys

from os.path import join as pj

from Bio import SeqIO

import gzip
import io
from math import floor

class LoadSequenceData(eHive.BaseRunnable):

    def param_defaults(self):
        return {
            'cs_order' : 'chunk,contig,non_ref_scaffold,scaffold,superscaffold,linkage_group,chromosome',
            'IUPAC' : 'RYKMSWBDHV',
        }

    def run(self):
        # params
        en_root = self.param_required("ensembl_root_dir")
        wd = self.param_required("work_dir")
      
        cs_order= self.param("cs_order")
        IUPAC = self.param("IUPAC")
        
        # 
        errors = []

        # initialize whatever
        fasta_raw = self.from_param("manifest_data", "fasta_dna")
        agps = self.from_param("manifest_data", "agp")
        genome = self.from_param("manifest_data", "genome")
        sra = self.from_param("manifest_data", "seq_region")
        
        # rename IUPAC to N symbols using sed
        os.makedirs(pj(wd,"fasta"), exist_ok=True)
        iupac_cmd = r'''{_cat} {_file} | sed -r '/^[^>]/ {{ s/[{_IUPAC}]+/N/g; s/{_iupac}/n/g }}' > {_out}'''.format(
            _cat = self.is_gz(fasta_raw) and "zcat" or "cat",
            _file = fasta_raw,
            _IUPAC = IUPAC.upper(),
            _iupac = IUPAC.lower(),
            _out = pj(wd, "fasta", "seq_no_iupac.fasta")
            )
        try:
            print("running %s" % (iupac_cmd), file = sys.stderr)
            res = sp.run(iupac_cmd, shell=True, check=True)
        except sp.CalledProcessError as e:
            errors += [ "Execution failed: {} ".format(e) ]

        # end
        # TODO: use catch raise catch instead
        if errors:
            raise Exception("Integrity test failed: " + str(errors))

    def is_gz(self, filename):
      try:   
          return magic.from_file(filename, mime=True) == 'application/x-gzip'
      except:
          ...
      return False 
 
    # TODO: add some metafunc setter getter
    def from_param(self, param, key):
        data = self.param_required(param)
        if key not in data:
            raise Exception("Missing required %s data: %s" % (param , key))
        return data[key]

