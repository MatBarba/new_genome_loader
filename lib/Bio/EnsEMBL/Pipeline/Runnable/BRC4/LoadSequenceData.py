#!env python3

from BCBio import GFF

import eHive
import json
import magic
import os
import subprocess as sp
import sys

from os.path import dirname, join as pj

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
        errors = []
        # params
        en_root = self.param_required("ensembl_root_dir")
        wd = self.param_required("work_dir")

        # initialize whatever
        genome = self.from_param("manifest_data", "genome")
        sra = self.from_param("manifest_data", "seq_region")
        
        # split into contigs
        # TODO: add to agps,.... 

        # rename IUPAC to N symbols using sed
        fasta_raw = self.from_param("manifest_data", "fasta_dna")
        fasta_clean = pj(wd, "fasta", "seq_no_iupac.fasta")
        try:
            self.remove_IUPAC(fasta_raw, fasta_clean)
        except sp.CalledProcessError as e:
            errors += [ "Execution failed: {} ".format(e) ]

        # order agp based on cs_order 
        cs_order_lst = map(lambda x: x.strip(), self.param("cs_order").split(","))
        cs_order = { e:i for i,e in enumerate(filter(lambda x: len(x)>0, cs_order_lst)) }

        agps = self.from_param("manifest_data", "agp")
        cs_used_set = frozenset(sum(map(lambda x: x.split("-"), agps.keys()),[]))
        cs_unknown = cs_used_set.difference(cs_order.keys()) 
        if (len(cs_unknown) > 0):
            raise Exception("Unknown coordinate system(s) %s" % {str(cs_unknown)})

        # when loading agp sort by:
        #   highest component (cmp) cs level
        #   lowest difference between cs ranks (asm - cmp)
        #   i.e: chromosome-scaffold
        agp_cs_pairs = list(map(lambda x: [x]+x.split("-"), agps.keys()))
        agp_levels = [ (x[0], cs_order[x[1]], cs_order[x[2]]) for x in agp_cs_pairs ]
        bad_agps = list(filter(lambda x: x[1] < x[2], agp_levels))
        if (len(bad_agps) > 0):
            raise Exception("component cs has higher order than assembled cs %s" % {str(bad_agps)})
         
        agp_levels_sorted = [ e[0] for e in sorted(agp_levels, key=lambda x:(-x[2], x[1]-x[2])) ]
        print(" ".join(agp_levels_sorted), file=sys.stderr)
        #prune agps
        # make pruned agps
        # use them for loading

        """
i.e., when chunk(1),contig(2),non_ref_scaffold(3),scaffold(4),superscaffold(5),linkage_group(6),chromosome(7)
chromosome_scaffold 7:4
scaffold_contig 4:2
chromosome_contig 7:2
contig_chunk 2:1
scaffold_chunk 4:1
chromosome_chunk 7:1
do not use agp mapping line, if src (lowest cs rank component is in the used list).
fix gaps. remove leading/trailing gaps.
don't include chains made of gaps only.
        """
        # end
        # TODO: use catch raise catch instead
        if errors:
            raise Exception("Integrity test failed: " + str(errors))

    # STAGES
    def remove_IUPAC(self, from_file, to_file):
        IUPAC = self.param("IUPAC")
        os.makedirs(dirname(to_file), exist_ok=True)
        cmd = r'''{_cat} {_file} | sed -r '/^[^>]/ {{ s/[{_IUPAC}]+/N/g; s/{_iupac}/n/g }}' > {_out}'''.format(
            _cat = self.is_gz(from_file) and "zcat" or "cat",
            _file = from_file,
            _IUPAC = IUPAC.upper(),
            _iupac = IUPAC.lower(),
            _out = to_file
        )
        print("running %s" % (cmd), file = sys.stderr)
        return sp.run(cmd, shell=True, check=True)

    # UTILS
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

