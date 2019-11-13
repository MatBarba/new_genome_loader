#!env python3

from BCBio import GFF

import eHive
import gzip
import json
import magic
import os
import subprocess as sp
import sys

from os.path import dirname, join as pj

from Bio import SeqIO

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
        #   i.e: chromosome-scaffold scaffold-chunk chromosome-chunk
        agp_cs_pairs = list(map(lambda x: [x]+x.split("-"), agps.keys()))
        agp_levels = [ (x[0], cs_order[x[1]], cs_order[x[2]]) for x in agp_cs_pairs ]
        bad_agps = list(filter(lambda x: x[1] < x[2], agp_levels))
        if (len(bad_agps) > 0):
            raise Exception("component cs has higher order than assembled cs %s" % {str(bad_agps)})
         
        agp_levels_sorted = [ e[0] for e in sorted(agp_levels, key=lambda x:(-x[2], x[1]-x[2])) ]
        print(" ".join(agp_levels_sorted), file=sys.stderr)

        cs_rank = {e:i for i,e in enumerate(sorted(cs_used_set,key=lambda x:-cs_order[x]), start=1)} 
        #prune agps
        agps_pruned_dir = pj(wd, "agps_pruned") 
        agps_pruned = {}
        used_components = set()
        for asm_cmp in agp_levels_sorted:
           agp_file_src = agps[asm_cmp] 
           agp_file_dst = pj(agps_pruned_dir, asm_cmp + ".agp")
           if self.agp_prune(agp_file_src, agp_file_dst, used_components) > 0:
               agps_pruned[asm_cmp] = agp_file_dst

        print(agps_pruned, file = sys.stderr)
        # for asm_cmp in filter(lambda k: k in agps_pruned, agp_levels_sorted):

        # end
        # TODO: use catch raise catch instead
        if errors:
            raise Exception("Loading sequence data failed: " + str(errors))

    # STAGES
    def agp_prune(self, from_file, to_file, used):
        # reomve used component 
        #   and GAPS as they are not used by 'ensembl-pipeline/scripts/load_agp.pl'
        os.makedirs(dirname(to_file), exist_ok=True)
        open_ = self.is_gz(from_file) and gzip.open or open 
        writes = 0
        with open_(from_file, "r") as src:
            with open(to_file, "w") as dst:
                for line in src:
                    fields = line.strip().split("\t")
                    ( asm_id, asm_start, asm_end, asm_part,
                      type_,
                      cmp_id, cmp_start, cmp_end, cmp_strand
                    ) = fields
                    if type_ in "NU" or cmp_id in used:
                        continue 
                    used.add(cmp_id)
                    print(line.strip(), file = dst)
                    writes += 1
        return writes


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

