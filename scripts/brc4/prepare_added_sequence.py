#!env python3

import re, sys, argparse
import os, json

from Bio import SeqIO
from Bio import GenBank
import tempfile

class FormattedFilesGenerator():
    """
    Contains a parser to load data from a file, and output a set of files that follow our schema for input into a core database
    """
    
    locations = {
            "mitochondrion" : "mitochondrial_chromosome",
            "apicoplast" : "apicoplast_chromosome",
            "chloroplast" : "chloroplast_chromosome",
            "chromoplast" : "chromoplast_chromosome",
            "cyanelle" : "cyanelle_chromosome",
            "leucoplast" : "leucoplast_chromosome",
            }
    
    def __init__(self, output_dir):
        self.output_dir = output_dir
        try:
            os.mkdir(self.output_dir)
        except FileExistsError:
            pass
            
        
        self.genome_json = os.path.join(self.output_dir, 'genome.json')
        self.seq_regions_json = os.path.join(self.output_dir, 'seq_regions.json')
        self.dna_fasta = os.path.join(self.output_dir, 'dna.fasta')
        self.pep_fasta = os.path.join(self.output_dir, 'pep.fasta')
        self.genes_gff = os.path.join(self.output_dir, 'genes.gff')

        self.seqs = []
        
    def parse_genbank(self, gb_file):
        """
        Load a sequence from a Genbank file
        """
        
        organella = self._get_organella(gb_file)
        
        with open(gb_file, "r") as gbh:
            for record in SeqIO.parse(gbh, "genbank"):
                # We don't want the record description (especially for the fasta file)
                record.description = ""
                record.organelle = None
                if record.id in organella:
                    record.organelle = organella[record.id]
                self.seqs.append(record)
            
            self._write_genome_json()
            self._write_seq_regions_json()
            self._write_dna_fasta()
    
    def _get_organella(self, gb_file):
        """
        Retrive the organelle from the genbank file, using the specific GenBank object,
        because SeqIO does not support this field
        """
        organella = {}
        with open(gb_file, "r") as gbh:
            for record in GenBank.parse(gbh):
                accession = record.version
                for q in record.features[0].qualifiers:
                    if q.key == "/organelle=":
                        organelle = q.value.replace('"', '')
                        organella[record.version] = organelle
        return organella
    
    def _write_dna_fasta(self):
        with open(self.dna_fasta, "w") as fasta_fh:
            SeqIO.write(self.seqs, fasta_fh, "fasta")
    
    def _write_seq_regions_json(self):
        json_array = []
        
        for seq in self.seqs:
            seq_obj = {
                    "name" : seq.id,
                    "coord_system_level" : "chromosome",
                    "circular" : (seq.annotations["topology"] == "circular"),
                    "codon_table" : 1,
                    "length" : len(seq.seq),
                    }
            if seq.organelle:
                seq_obj["location"] = self._prepare_location(seq.organelle)
                print(f"Warning: '{seq.organelle}' is an organelle: make sure to change the codon table number in {self.seq_regions_json} manually if it is not the standard codon table")
            
            # Additional attributes for Ensembl
            seq_obj["added_seq"] = {
                    "accession" : seq.id,
                    "asm" : {
                        "pr_nam" : "",
                        "pr_url" : "",
                        }
                    }
            if not seq_obj["added_seq"]["asm"]["pr_nam"]:
                print(f"Warning: please add the relevant provider name for the assembly in {self.seq_regions_json}")
            if not seq_obj["added_seq"]["asm"]["pr_url"]:
                print(f"Warning: please add the relevant provider url for the assembly in {self.seq_regions_json}")

            # Additional attributes for gene set, if any
            # TODO

            json_array.append(seq_obj)
        with open(self.seq_regions_json, "w") as seq_fh:
            seq_fh.write(json.dumps(json_array, indent=4))

    def _prepare_location(self, organelle):
        """
        Given an organelle name, returns the SO term corresponding to its location
        """
        if organelle in self.locations:
            return self.locations[organelle]
        else:
            raise Exception(f"Unkown organelle: {location}")

    def _write_genome_json(self):
        """
        Write a draft for the genome json file
        Only the production_name is needed, but the rest of the fields need to be given
        for the validation of the json file
        """
        genome_data = {
                "species" : {
                    "production_name" : "",
                    "taxonomy_id" : 0,
                    },
                "assembly" : {
                    "accession" : "GCA0000000",
                    "version" : 1
                    }
                }

        if not genome_data["species"]["production_name"]:
            print(f"Warning: please add the relevant production_name for this genome in {self.genome_json}")
            
                
        
        with open(self.genome_json, "w") as genome_fh:
            genome_fh.write(json.dumps(genome_data, indent=4))

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Prepare a genbank file for an additional sequence')
    
    parser.add_argument('--gb', type=str, required=True,
                help='Genbank file')
    parser.add_argument('--output_dir', type=str, required=True,
                help='Output dir')


    args = parser.parse_args()
    formatter = FormattedFilesGenerator(args.output_dir)
    formatter.parse_genbank(args.gb)

if __name__ == "__main__":
    main()
