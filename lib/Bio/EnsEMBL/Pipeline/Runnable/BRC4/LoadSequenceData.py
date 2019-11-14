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
            raise Exception("component cs has higher order than assembled cs %s" % (str(bad_agps)))
         
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

        asm_v = self.from_param("genome_data","assembly")["name"]

        sequence_rank = max(cs_rank.values())
        for (cs, rank) in sorted(cs_rank.items(), key=lambda p: -p[1]):
           logs = pj(wd, "load", "%02d_%s" %(rank, cs) )
           if (rank == sequence_rank): 
               self.load_cs_data(cs, rank, "fasta", asm_v, fasta_clean, logs, seq_level = True)
           else:
               useful_agps = list(filter(lambda x: cs in x, agps_pruned.keys()))
               if len(useful_agps) == 0:
                   raise Exception("non-seq_level cs %s has no agps to assemble it from" % (cs))
               for pair, agp_file_pruned in map(lambda k: (k, agps_pruned[k]), useful_agps):
                   if (not pair.startswith(cs+"-")):
                       continue
                   self.load_cs_data(cs, rank, pair, asm_v, agp_file_pruned, logs)

        # TODO
        # ! add parameter to turn pruning off
        # set ENA attribs
        # add synonims from json and synonims with version (remove from sr name) for nonENA
        # nullify sequence? versions in cs 
        # set attributes MT, circular, location, etc
        # ???

        # end
        # TODO: use catch raise catch instead
        if errors:
            raise Exception("Loading sequence data failed: " + str(errors))

    # STAGES
    def load_cs_data(self, cs, rank, pair, asm_v, src_file, log_pfx, seq_level = False):
        # NB load_seq_region.pl and load_agp.pl are not failing on parameter errors (0 exit code)
        os.makedirs(dirname(log_pfx), exist_ok=True)
        self.load_seq_region(cs, rank, asm_v, src_file, log_pfx, seq_level)
        if not seq_level: 
            self.load_agp(pair, asm_v, src_file, log_pfx)

    def load_seq_region(self, cs, rank, asm_v, src_file, log_pfx, seq_level = False):
        en_root = self.param_required("ensembl_root_dir")
        cmd = (r'''{_loader} {_db_string} -coord_system_version {_asm_v} -default_version ''' +
               r'''    -rank {_rank} -coord_system_name {_cs} {_sl_flag} -{_tag}_file {_file}''' +
               r'''     > {_log}.stdout 2> {_log}.stderr''').format(
            _loader = "perl %s" % (pj(en_root, r"ensembl-pipeline/scripts/load_seq_region.pl")),
            _db_string = self.db_string(),
            _asm_v = asm_v,
            _rank = rank,
            _cs = cs,
            _sl_flag = seq_level and "-sequence_level" or "",
            _tag = seq_level and "fasta" or "agp",
            _file = src_file,
            _log = "%s_seq" % (log_pfx),
        )
        print("running %s" % (cmd), file = sys.stderr)
        return sp.run(cmd, shell=True, check=True)
  
    def load_agp(self, pair, asm_v, src_file, log_pfx):
        en_root = self.param_required("ensembl_root_dir")
        (asm_n, cmp_n) = pair.strip().split("-")
        cmd = (r'''{_loader} {_db_string} -assembled_version {_asm_v} ''' +
               r'''    -assembled_name {_asm} -component_name {_cmp} ''' +
               r'''    -agp_file {_file} ''' +
               r'''    > {_log}.stdout 2> {_log}.stderr''').format(
            _loader = "perl %s" % (pj(en_root, r"ensembl-pipeline/scripts/load_agp.pl")),
            _db_string = self.db_string(),
            _asm_v = asm_v,
            _asm = asm_n,
            _cmp = cmp_n,
            _file = src_file,
            _log = "%s_agp_%s" % (log_pfx, pair.replace("-","_")),
        )
        print("running %s" % (cmd), file = sys.stderr)
        return sp.run(cmd, shell=True, check=True)

    def db_string(self):
        return "-dbhost {host_} -dbport {port_} -dbuser {user_} -dbpass {pass_} -dbname {dbname_} ".format(
            host_ = self.param("dbsrv_host"),
            port_ = self.param("dbsrv_port"),
            user_ = self.param("dbsrv_user"),
            pass_ = self.param("dbsrv_pass"),
            dbname_ = self.param("db_name")
        )


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

