#!env python3

import eHive
import os
import subprocess as sp
import sys

from os.path import dirname, join as pj

class FillMetadata(eHive.BaseRunnable):

    def param_defaults(self):
        return {
            'division' : '',
            'ignore' : [],
            'copy' : {},
        }


    def run(self):
        # params
        wd = self.param_required("work_dir")
        division = self.param("division")
        ignore = self.param("ignore")
        copy = self.param("copy")
        genome_data = self.param_required("genome_data")

        # update division
        if "species" in genome_data:
            if "division" not in genome_data["species"]:
                genome_data["species"]["division"] = division

        # flattern and dump
        flat = self.flattern(genome_data, ignore)
        flat += list(map( lambda p: ( copy[p[0]], p[1] ), filter(lambda x: x[0] in copy, flat) ))
        flat = list(filter(lambda x: x[0] not in ignore, flat))
        if len(flat) <= 0:
            return

        # dump
        os.makedirs(wd, exist_ok=True)
        sql_file = pj(wd, "insert_meta.sql")
        with open(sql_file, "w") as sql:
            print("insert ignore into meta (species_id, meta_key, meta_value) values", file=sql)
            c = ""
            for _key, _val in flat:
                if isinstance(_val, bool):
                    _val = int(_val)
                if isinstance(_val, str):
                    _val = '"%s"' % (_val)
                print ('%s (1, "%s", %s)' % (c, str(_key), str(_val)), file = sql)
                c = ","
            print(";", file=sql)
        # run insert sql
        self.run_sql_req(sql_file, log_pfx = sql_file, from_file = True)


    def flattern(self, data, ignore_lst, pfx = None):
        # vectors with values only
        if isinstance(data, list):
            return [ ( pfx, v ) for v in list ]
        elif isinstance(data, dict):
            ignore_lst = frozenset(ignore_lst)
            res = []
            for (k, v) in data.items():
                new_pfx = ".".join(filter(lambda p: p != None, [pfx, k]))
                res.append(self.flattern(v, ignore_lst, new_pfx))
            return sum(res, [])
        return data != None and [ (pfx, data) ] or []


    def run_sql_req(self, sql, log_pfx, from_file = False):
        os.makedirs(dirname(log_pfx), exist_ok=True)
        en_root = self.param_required("ensembl_root_dir")

        sql_option = r''' -sql '{_sql}' '''.format(_sql = sql)
        if from_file:
            sql_option = r''' < '{_sql}' '''.format(_sql = sql)

        cmd = r'''{_dbcmd} -url "{_srv}{_dbname}" {_sql_option} > {_out} 2> {_err}'''.format(
            _dbcmd = 'perl %s/ensembl-hive/scripts/db_cmd.pl' %(en_root),
            _srv = self.param("dbsrv_url"),
            _dbname = self.param("db_name"),
            _sql_option = sql_option,
            _out = log_pfx + ".stdout",
            _err = log_pfx + ".stderr"
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
