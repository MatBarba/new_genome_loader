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
            'unversion_scaffolds' : 0,
            'versioned_sr_syn_src' : 'INSDC', # 50710
            'sr_syn_src' : 'ensembl_internal_synonym', # 50803
        }

    def run(self):
        errors = []
        # params
        en_root = self.param_required("ensembl_root_dir")
        wd = self.param_required("work_dir")

        # initialize whatever
        genome = self.from_param("manifest_data", "genome")
        sra = self.from_param("manifest_data", "seq_region")

        # TODO
        # split into contigs, add AGP
        # load data with no agps ??? m.b. create empty cs-cs agps

        # rename IUPAC to N symbols using sed
        fasta_raw = self.from_param("manifest_data", "fasta_dna")
        fasta_clean = pj(wd, "fasta", "seq_no_iupac.fasta")
        try:
            self.remove_IUPAC(fasta_raw, fasta_clean)
        except sp.CalledProcessError as e:
            errors += [ "Execution failed: {} ".format(e) ]

        agps = self.from_param("manifest_data", "agp")
        cs_order = self.coord_sys_order(self.param("cs_order"))
        cs_rank = self.used_cs_ranks(agps, cs_order)

        agps_pruned_dir = pj(wd, "agps_pruned")
        agps_pruned = self.prune_agps(agps, cs_order, agps_pruned_dir, self.param_bool("prune_agp"))

        self.load_seq_data(fasta_clean, agps_pruned, cs_rank, pj(wd, "load"))

        self.add_contig_ena_attrib(pj(wd, "load", "set_ena"))
        if self.param_bool("unversion_scaffolds"):
            self.unversion_scaffolds(cs_rank, pj(wd, "unversion_scaffolds"))

        seq_reg_meta = self.from_param("manifest_data", "seq_region")
        self.add_synonyms(seq_reg_meta, self.param_bool("unversion_scaffolds"))

        # TODO
        # nullify sequence? versions in cs
        # set attributes MT, circular, location, etc
        # ???

        # end
        # TODO: use catch raise catch instead
        if errors:
            raise Exception("Loading sequence data failed: " + str(errors))

    # STAGES
    def add_synonyms(self, meta_file, unversioned = False):
        pass

    def unversion_scaffolds(self, cs_rank, logs):
        # non-versioned syns for contigs, versioned for the rest
        seq_cs, max_rank = max([ (c, r) for c, r in cs_rank.items()], key = lambda k: k[1])
        for cs in cs_rank:
            if cs == seq_cs:
                xdb = self.param("sr_syn_src")
                self.copy_sr_name_to_syn(cs, xdb, pj(logs, "cp2syn", cs))
                self.sr_name_unversion(cs, "seq_region_synonym", "synonym", pj(logs, "unv_srs", cs))
            else:
                xdb = self.param("versioned_sr_syn_src")
                self.copy_sr_name_to_syn(cs, xdb, pj(logs, "cp2syn", cs))
                self.sr_name_unversion(cs, "seq_region", "name", pj(logs, "unv_sr", cs))

    def run_sql_req(self, sql, log_pfx):
        os.makedirs(dirname(log_pfx), exist_ok=True)
        en_root = self.param_required("ensembl_root_dir")

        cmd = r'''{_dbcmd} -url "{_srv}{_dbname}" -sql '{_sql}' > {_out} 2> {_err}'''.format(
            _dbcmd = 'perl %s/ensembl-hive/scripts/db_cmd.pl' %(en_root),
            _srv = self.param("dbsrv_url"),
            _dbname = self.param("db_name"),
            _sql = sql,
            _out = log_pfx + ".stdout",
            _err = log_pfx + ".stderr"
        )
        print("running %s" % (cmd), file = sys.stderr)
        return sp.run(cmd, shell=True, check=True)


    def add_contig_ena_attrib(self, log_pfx):
        # Add ENA attrib for contigs (no sequence_level checks -- just cs name)
        #   (see ensembl-datacheck/lib/Bio/EnsEMBL/DataCheck/Checks/SeqRegionNamesINSDC.pm)
        sql = r'''insert into seq_region_attrib (seq_region_id, attrib_type_id, value)
                select
                  sr.seq_region_id, at.attrib_type_id, "ENA"
                from
                  seq_region sr, coord_system cs, attrib_type at
                where   sr.coord_system_id = cs.coord_system_id
                    and cs.name = "contig"
                    and at.code = "external_db"
              ;'''
        return self.run_sql_req(sql, log_pfx)

    def copy_sr_name_to_syn(self, cs, x_db, log_pfx):
        asm_v = self.from_param("genome_data","assembly")["name"]
        sql = r'''insert into seq_region_synonym (seq_region_id, synonym, external_db_id)
                  select
                      sr.seq_region_id, sr.name, xdb.external_db_id
                  from
                     seq_region sr, external_db xdb, coord_system cs
                  where   xdb.db_name = "%s"
                      and sr.coord_system_id = cs.coord_system_id
                      and cs.name = "%s"
                      and cs.version = "%s"
                      and sr.name like "%%._"
                ;''' % (x_db, cs, asm_v)
        return self.run_sql_req(sql, log_pfx)


    def sr_name_unversion(self, cs, tbl, fld, log_pfx):
        # select synonym, substr(synonym,  1, locate(".", synonym, length(synonym)-2)-1)
        #     from seq_region_synonym  where synonym like "%._"
        asm_v = self.from_param("genome_data","assembly")["name"]
        sql = r'''update {_tbl} t, seq_region sr, coord_system cs
                    set
                      t.{_fld} = substr(t.{_fld},  1, locate(".", t.{_fld}, length(t.{_fld})-2)-1)
                    where t.{_fld} like "%._"
                      and t.seq_region_id = sr.seq_region_id
                      and sr.coord_system_id = cs.coord_system_id
                      and cs.name = "{_cs}"
                      and cs.version = "{_asm_v}"
                ;'''.format(
                    _tbl = tbl,
                    _fld = fld,
                    _cs = cs,
                    _asm_v = asm_v
                )
        return self.run_sql_req(sql, log_pfx)

    def coord_sys_order(self, cs_order_str):
        cs_order_lst = map(lambda x: x.strip(), cs_order_str.split(","))
        return { e:i for i,e in enumerate(filter(lambda x: len(x)>0, cs_order_lst)) }

    def used_cs_ranks(self, agps, cs_order):
        cs_used_set = frozenset(sum(map(lambda x: x.split("-"), agps.keys()),[]))
        cs_unknown = cs_used_set.difference(cs_order.keys())
        if (len(cs_unknown) > 0):
            raise Exception("Unknown coordinate system(s) %s" % {str(cs_unknown)})
        return { e:i for i,e in enumerate(sorted(cs_used_set,key=lambda x:-cs_order[x]), start=1) }

    def prune_agps(self, agps, cs_order, agps_pruned_dir, pruning = True):
        # when loading agp sort by:
        #   highest component (cmp) cs level (lowest rank)
        #   lowest difference between cs ranks (asm - cmp)
        #   i.e: chromosome-scaffold scaffold-chunk chromosome-chunk
        agp_cs_pairs = list(map(lambda x: [x]+x.split("-"), agps.keys()))
        agp_levels = [ (x[0], cs_order[x[1]], cs_order[x[2]]) for x in agp_cs_pairs ]
        bad_agps = list(filter(lambda x: x[1] < x[2], agp_levels))
        if (len(bad_agps) > 0):
            raise Exception("component cs has higher order than assembled cs %s" % (str(bad_agps)))

        agp_levels_sorted = [ e[0] for e in sorted(agp_levels, key=lambda x:(-x[2], x[1]-x[2])) ]

        #prune agps
        agps_pruned = {}
        used_components = set()
        if not pruning:
            used_components = None
        for asm_cmp in agp_levels_sorted:
            agp_file_src = agps[asm_cmp]
            agp_file_dst = pj(agps_pruned_dir, asm_cmp + ".agp")
            if self.agp_prune(agp_file_src, agp_file_dst, used_components) > 0:
                agps_pruned[asm_cmp] = agp_file_dst
        return agps_pruned


    def load_seq_data(self, fasta, agps, cs_rank, log_pfx):
        asm_v = self.from_param("genome_data","assembly")["name"]

        sequence_rank = max(cs_rank.values())
        for (cs, rank) in sorted(cs_rank.items(), key=lambda p: -p[1]):
           logs = pj(log_pfx, "%02d_%s" %(rank, cs) )
           if (rank == sequence_rank):
               self.load_cs_data(cs, rank, "fasta", asm_v, fasta, logs, seq_level = True)
           else:
               useful_agps = list(filter(lambda x: cs in x, agps.keys()))
               if len(useful_agps) == 0:
                   raise Exception("non-seq_level cs %s has no agps to assemble it from" % (cs))
               load_region = True
               for pair, agp_file_pruned in map(lambda k: (k, agps[k]), useful_agps):
                   if (not pair.startswith(cs+"-")):
                       continue
                   self.load_cs_data(cs, rank, pair, asm_v, agp_file_pruned, logs, load_region)
                   load_region = False

    def load_cs_data(self,
                     cs, rank, pair, asm_v,
                     src_file, log_pfx,
                     load_region = True, seq_level = False):
        # NB load_seq_region.pl and load_agp.pl are not failing on parameter errors (0 exit code)
        os.makedirs(dirname(log_pfx), exist_ok=True)
        if load_region:
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


    def agp_prune(self, from_file, to_file, used = None):
        # reomve used component
        #   and GAPS as they are not used by 'ensembl-pipeline/scripts/load_agp.pl'
        os.makedirs(dirname(to_file), exist_ok=True)
        open_ = self.is_gz(from_file) and gzip.open or open
        if used is None:
            cmd = r'''{_cat} {_file} > {_out}'''.format(
                _cat = self.is_gz(from_file) and "zcat" or "cat",
                _file = from_file,
                _out = to_file
            )
            print("running %s" % (cmd), file = sys.stderr)
            sp.run(cmd, shell=True, check=True)
            return 1
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

    def param_bool(self, param):
        val = self.param(param)
        return bool(val) and "0" != val

