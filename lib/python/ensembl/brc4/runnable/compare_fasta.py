#!env python3

import eHive
import gzip
import json
import io
import re
import sys
from functools import partial

from Bio import SeqIO
from os import path
from ensembl.brc4.runnable.parser import Parser

class SeqGroup():
    def __init__(self, sequence, identifier=None):
        self.sequence = sequence
        self.length = len(self.sequence)
        self.ids = []
        if identifier: self.add_id(identifier)
        self.count = len(self.ids)
    
    def __str__(self):
        return ", ".join(self.ids)
    
    def add_id(self, identifier):
        self.ids.append(identifier)
        self.count = len(self.ids)

class compare_fasta(eHive.BaseRunnable):

    def param_defaults(self):
        return {
        }

    def run(self):
        report = self.param_required("report")
        fasta1 = self.param_required("fasta1")
        fasta2 = self.param_required("fasta2")
        map_dna_path = self.param_required("seq_regions")
        output_dir = self.param_required("output_dir")
        species = self.param_required("species")
        name = self.param_required("comparison_name")

        map_dna = self.get_map(map_dna_path)
        seq1 = self.get_fasta(fasta1, map_dna)
        seq2 = self.get_fasta(fasta2, map_dna)

        (stats, diffs, seq_map) = self.compare_seqs(seq1, seq2)

        # Print mapping to a file (add report data)
        map_file = output_dir + "/" + species + "_" + name + ".map"
        self.print_map(seq_map, map_file, report)

        # Print full list of results in a file
        output_file = output_dir + "/" + species + "_" + name + ".log"
        print("Write results in %s" % output_file)
        with open(output_file, "w") as out_fh:
            for line in diffs:
                out_fh.write(line + "\n")

        # Print the stats separately
        out = {
            "species": species,
            "stats": stats
        }
        self.dataflow(out, 2)

    def print_map(self, seq_map, map_file, report_file):

        report_parser = Parser()
        report_seq = report_parser.get_report_regions(report_file)
        report = self.add_report_to_map(seq_map, report_seq)

        print("Write map in %s" % map_file)
        with open(map_file, "w") as out_fh:
            out_fh.write(json.dumps(report, sort_keys=True, indent=4))

    def add_report_to_map(self, seq_map, report_seq):

        accession_version = r'\.\d+$'
        report = []
        for insdc_name, old_name in seq_map.items():
            if insdc_name not in report_seq:
                raise Exception("No INSDC %s found in report" % insdc_name)
            else:
                seqr = report_seq[insdc_name]
                seqr["name"] = old_name
                seqr["EBI_seq_region_name"] = old_name
                brc4_name = insdc_name
                brc4_name = re.sub(accession_version, '', brc4_name)
                seqr["BRC4_seq_region_name"] = brc4_name
                syns = [{
                    "source": "INSDC",
                    "name": insdc_name
                }]
                seqr["synonyms"] = syns
                report.append(seqr)
        return report

    def get_map(self, map_path):

        print("Read file %s" % map_path)
        data = self.get_json(map_path)

        map_dna = {}

        for seqr in data:
            if "synonyms" in seqr:
                for syn in seqr["synonyms"]:
                    if syn["name"] == "INSDC":
                        map_dna[name] = syn["value"]

        return map_dna

    def get_json(self, json_path):
        with open(json_path) as json_file:
            return json.load(json_file)
        
    def build_seq_dict(self, seqs):
        """Build a seq dict taking duplicates into account"""
        
        seqs_dict = dict()
        
        for name, seq in seqs.items():
            if seq in seqs_dict:
                seqs_dict[seq].add_id(name)
            else:
                seqs_dict[seq] = SeqGroup(seq, name)
        
        return seqs_dict
                

    def get_fasta(self, fasta_path, map_dna):

        print("Read file %s" % fasta_path)
        sequences = {}
        _open = partial(gzip.open, mode='rt') if fasta_path.endswith(
            '.gz') else open
        with _open(fasta_path) as fasta_fh:
            for rec in SeqIO.parse(fasta_fh, "fasta"):
                name = rec.id
                if name in map_dna:
                    name = map_dna[name]
                sequences[name] = re.sub(r"[^CGTA]", "N", str(rec.seq.upper()))
        return sequences

    def compare_seqs(self, seq1, seq2):
        comp = []
        diff = abs(len(seq1) - len(seq2))
        stats = {
            "accession": None,
            "seq_count_1": len(seq1),
            "seq_count_2": len(seq2),
            "num_diff_seq": diff,
            "common": 0,
            "only1": 0,
            "only2": 0,
            "max_only1": 0,
            "max_only2": 0,
            "only1_200": 0,
            "only1_1000": 0,
            "only2_200": 0,
            "only2_1000": 0,
            "other_locations": 0,
            "summary": None,
            "organellar_summary": None
        }

        value = None  # variable used for summary
        org_value = "no_organelles_present"   # variable used for organellar_summary

        # Compare number of sequences
        if len(seq1) != len(seq2):
            comp.append("WARNING: Different number of sequences: %d vs %d" % (
                len(seq1), len(seq2)))
            value = "mismatch"
            org_value = "identical"
        else:
            comp.append("Same number of sequences: %d" % len(seq1))

        # Compare sequences
        seqs1 = self.build_seq_dict(seq1)
        seqs2 = self.build_seq_dict(seq2)
        
        common, group_comp = self.find_common_groups(seqs1, seqs2)
        comp += group_comp
        only1 = {seq: group for seq, group in seqs1.items() if not seq in seqs2}
        only2 = {seq: group for seq, group in seqs2.items() if not seq in seqs1}

        if len(seq1) == len(seq2) and len(common) == len(seq1):
            if len(only1) == 0 and len(only2) == 0:
                value = "identical"
        else:
            value = "mismatch"

        # comparing the organellar sequences
        report = self.param_required("report")
        report_parser = Parser()
        report_seq = report_parser.get_report_regions(report)
        report = self.add_report_to_map(common, report_seq)
        map_dna_path = self.param_required("seq_regions")
        data = self.get_json(map_dna_path)
        org = []
        org1 = []
        org2 = []
        org3 = []
        # get the accession
        fasta1 = self.param_required("fasta1")
        get_acc = fasta1
        acc = get_acc.split('/')
        accession = acc[4]
        # Gathering data from the INSDC report file and storing it into a list
        for k, v in report_seq.items():
            if "location" in v:
                if v['location'] not in ("chromosome", "nuclear_chromosome", "linkage_group"):
                    x = v['location']
                    org1.append(x)
                    org.append(k)
            else:
                pass

        # Gathering data from Seq_json file and storing it into a list
        for j in data:
            for k, v in j.items():
                if "location" in k:
                    if v not in ("chromosome", "nuclear_chromosome", "linkage_group"):
                        y = j['EBI_seq_region_name']
                        org3.append(v)
                        org2.append(y)
                else:
                    pass

        # checking for INSDC names for organellar sequences found in core
        for i in range(0, len(org2)):
            for p, l in report_seq.items():
                if l['name'] == org2[i]:
                    org2[i] = p
        for index, i in enumerate(org2):
            if i not in org:
                org.append(i)
                org1.append(org3[index])
            else:
                print("you have duplicates")

        # checking for multiple entries of organellar seq
        myList = [i.split('.')[0] for i in org]
        myList1 = [j[:-1] for j in myList]  # similar accession
        a = list(set(myList1))

        # comparing organellar sequences with common, only1 and only2
        count = 0
        for index, i in enumerate(org):
            if i == 'na':
                comp.append("MISSING accession in the report (na)")
            else:
                if i in common.values():
                    count = count+1
                    comp.append("%s (both) in location: %s" % (i, org1[index]))
                    if value == "identical" and count > 0:
                        org_value = "organellar_present"
                elif i in only1.values():
                    count = count+1
                    comp.append("%s (only1) in  location: %s" %
                                (i, org1[index]))
                    org_value = "unknown_with_organellar"
                else:
                    count = count+1
                    comp.append("%s (only2) in location: %s" %
                                (i, org1[index]))
                    org_value = "unknown_with_organellar"

        # if the mistmatch is due to added organell
        if len(seq1) > len(seq2):
            greater_len = len(seq1)
        else:
            greater_len = len(seq2)
        diff_common = greater_len - len(common)

        if diff != 0:
            if diff == count and diff_common == count:
                org_value = "organellar_present"

        # checking if multiple entries of organellar chromosomes are present
        if len(myList1) != len(a):
            org_value = "WARNING:Multiple_entry"

        # updating the stats
        stats["accession"] = accession
        stats["common"] = len(common)
        stats["only1"] = len(only1)
        stats["only2"] = len(only2)
        stats["other_locations"] = count
        stats["summary"] = value
        stats["organellar_summary"] = org_value

        if len(only1) > 0 or len(only2) > 0:
            comp.append("\nCommon sequences: %d" % len(common))

            if only1:
                stats["max_only1"] = max(len(seq) for seq in only1.keys())

                # Only list sequences where the length is > 200
                mini = {seq: name for seq, name in only1.items()
                        if len(seq) <= 200}
                maxi = {seq: name for seq, name in only1.items()
                        if len(seq) > 200}

                if mini and len(mini) > 3000:
                    comp.append(
                        "WARNING: Ignoring %d sequences from 1 with length <= 200" % len(mini))
                    only1 = maxi

            if only1:
                # Only list sequences where the length is > 1000
                mini = {seq: name for seq, name in only1.items()
                        if len(seq) <= 1000}
                maxi = {seq: name for seq, name in only1.items()
                        if len(seq) > 1000}

                if mini and len(mini) > 3000:
                    comp.append(
                        "WARNING: Ignoring %d sequences from 1 with length <= 1000" % len(mini))
                    only1 = maxi

            if only1:
                total = sum([len(seq) for seq in only1.keys()])
                comp.append("WARNING: Sequences only in 1: %d (%d)" %
                            (len(only1), total))
                only_seq1 = {name: len(seq) for seq, name in only1.items()}
                for name, length in sorted(only_seq1.items(), key=lambda x: x[1]):
                    comp.append("\tOnly in 1: %s (%d)" % (name, length))

            if only2:
                stats["max_only2"] = max(len(seq) for seq in only2.keys())

                # Only list sequences where the length is > 200
                mini = {seq: name for seq, name in only2.items()
                        if len(seq) <= 200}
                maxi = {seq: name for seq, name in only2.items()
                        if len(seq) > 200}

                if mini and len(mini) > 3000:
                    comp.append(
                        "WARNING: Ignoring %d sequences from 2 with length <= 200" % len(mini))
                    only2 = maxi

            if only2:
                # Only list sequences where the length is > 1000
                mini = {seq: name for seq, name in only2.items()
                        if len(seq) <= 1000}
                maxi = {seq: name for seq, name in only2.items()
                        if len(seq) > 1000}

                if mini and len(mini) > 3000:
                    comp.append(
                        "WARNING: Ignoring %d sequences from 2 with length <= 1000" % len(mini))
                    only2 = maxi

            if only2:
                total = sum([len(seq) for seq in only2.keys()])
                comp.append("WARNING: Sequences only in 2: %d (%d)" %
                            (len(only2), total))
                only_seq2 = {name: len(seq) for seq, name in only2.items()}
                for name, length in sorted(only_seq2.items(), key=lambda x: x[1]):
                    comp.append("\tOnly in 2: %s (%d)" % (name, length))

        return (stats, comp, common)

    def find_common_groups(self, seqs1, seqs2):
        
        print(len(seqs1))
        print(len(seqs2))
        
        comp = []
        common = {}
        for seq1, group1 in seqs1.items():
            if seq1 in seqs2:
                group2 = seqs2[seq1]
                
                # Check that the 2 groups have the same number of sequences
                if group1.count == group2.count:
                    if group1.count == 1:
                        common[group1.ids[0]] = group2.ids[0]
                    else:
                        comp.append(f"Matched 2 identical groups of sequences: {group1} and {group2}")
                        possible_id2 = " OR ".join(group2.ids)
                        for id1 in group1.ids:
                            common[id1] = possible_id2
                            
                else:
                    comp.append(f"Matched 2 different groups of sequences ({group1.count} vs {group2.count}): {group1} and {group2}")
        
        print(len(common))
        return common, comp
