#!env python3

from redminelib import Redmine
import argparse
import os, json, re, time
import requests
import xml.etree.ElementTree as ET
 
url = 'https://redmine.apidb.org'
default_fields = dict(
        status_name = 'Data Processing (EBI)',
        cf_17 = "Data Processing (EBI)",
        )
insdc_pattern = "^GC[AF]_\d{9}(\.\d+)?$"
accession_api_url = "https://www.ebi.ac.uk/ena/browser/api/xml/%s"
veupathdb_id = 1976

def retrieve_genomes(redmine, output_dir, build=None):
    """
    Get genomes metadata from Redmine, store them in json files.
    Each issue/new_genome is stored as one file in the output dir
    """
    issues = get_all_genomes(redmine, build)
    if not issues:
        print("No files to create")
        return
    
    # Create the output dir
    try:
        os.mkdir(output_dir)
    except:
        pass
    
    for issue in issues:
        genome_structure = parse_genome(issue)
        if not genome_structure:
            print("Skipped issue %d (%s). Not enough metadata." % (issue.id, issue.subject))
            continue

        try:
            organism = genome_structure["BRC4"]["organism_abbrev"]
            organism_file = output_dir + "/" + organism + ".json"
            f = open(organism_file, "w")
            json.dump(genome_structure, f, indent=True)
            f.close()
        except Exception as error:
            print("Skipped issue %d (%s). %s." % (issue.id, issue.subject, error))
            pass

def get_all_genomes(redmine, build=None):
    """
    Query Redmine to get all new genomes, with or without genes
    """
    genomes_with_genes = get_issues(redmine, "Genome sequence and Annotation", build)
    genomes_without_genes = get_issues(redmine, "Assembled genome sequence without annotation", build)
#    genomes_without_genes = []
    print("%d issues for new genomes with genes found" % len(genomes_with_genes))
    print("%d issues for new genomes without genes found" % len(genomes_without_genes))
    
    issues = genomes_with_genes + genomes_without_genes
    
    print("%d issues found" % len(issues))
    return issues
 
def parse_genome(issue, warn=True):
    """
    Extract genome metadata from a Redmine issue
    Return a nested dict
    """
    
    customs = get_custom_fields(issue)
    genome = {
            "BRC4": {},
            "species": {},
            "assembly": {},
            "genebuild": {},
            }
    
    # Get GCA accession
    if "GCA number" in customs:
        accession = customs["GCA number"]["value"]
        accession = check_accession(accession)
        if not accession:
            return
        genome["assembly"]["accession"] = accession
    else:
        print("No accession for issue %d (%s)" % (issue.id, issue.subject))

    # Get BRC4 component
    if "Component DB" in customs:
        components = customs["Component DB"]["value"]
        if len(components) == 1:
            genome["BRC4"]["component"] = components[0]
        elif len(components) > 1:
            raise Exception("More than 1 component for new genome " + str(issue.id))
    else:
        print("No component for issue %d (%s)" % (issue.id, issue.subject))

    # Get Organism abbrev
    try:
        abbrev = customs["Organism Abbreviation"]["value"]
        if abbrev:
            genome["BRC4"]["organism_abbrev"] = abbrev
        else:
            if warn:
                print("No organism abbrev could be found for %s" % issue.id)
            return
    except:
        if warn:
            print("Can't get organism abbrev for %s" % issue.id)

    # Warn to get GFF2Load
    try:
        gff_path = customs["GFF 2 Load"]["value"]
        if gff_path and warn:
            print("GFF2Load: separate gff file for %s: %s (issue %d)" % (genome["BRC4"]["organism_abbrev"], gff_path, issue.id))
    except:
        pass

    # Warn for replacement
    try:
        if customs["Replacement genome?"]["value"] == "Yes" and warn:
            print("Replacement: the organism %s is a replacement (issue %d)" % (genome["BRC4"]["organism_abbrev"], issue.id))
    except:
        pass

    return genome

def check_accession(accession):
    """
    Check the accession string format
    Returns a cleaned up accession
    """
    accession = accession.strip()
    
    # Remove the url if it's in one
    # There might even be a trailing url
    accession = re.sub(r"^.+/([^/]+)/?", r"\1", accession)
    
    if re.match(insdc_pattern, accession):
        return accession
    else:
        return

def get_custom_fields(issue):
    """
    Put all Redmine custom fields in a dict instead of an array
    Return a dict
    """
    
    cfs = {}
    for c in issue.custom_fields:
        cfs[c["name"]] = c
    return cfs

def get_issues(redmine, datatype, build=None):
    """
    Retrieve all issue for new genomes, be they with or without gene sets
    Return a Redmine ResourceSet
    """
    
    other_fields = { "cf_94" : datatype }
    if build:
        version_id = get_version_id(redmine, build)
        other_fields["fixed_version_id"] = version_id

    return list(get_ebi_issues(redmine, other_fields))

def get_version_id(redmine, build):
    """
    Given a build number, get the version id for it
    """
    versions = redmine.version.filter(project_id=veupathdb_id)
    version_name = "Build " + str(build)
    version_id = [version.id for version in versions if version.name == version_name]
    return version_id
    
def get_ebi_issues(redmine, other_fields=dict()):
    """
    Get EBI issues from Redmine, add other fields if provided
    Return a Redmine ResourceSet
    """

    # Other fields replace the keys that already exist in default_fields
    search_fields = { **default_fields, **other_fields }
    
    return redmine.issue.filter(**search_fields)

def add_genome_organism_abbrev(redmine, build, species_list_path, update=False):
    """
    Retrieve genome issues, get the corresponding INSDC accession, and generate an organism abbrev
    Put the abbrev back to the Redmine ticket
    By default do not update the ticket
    """
    
    # Get the list of existing organism_abbrev
    current_orgs = get_current_orgs(species_list_path)
    
    # Get the genomes issues
    issues = get_all_genomes(redmine, build)
    if not issues:
        print("No Redmine tickets to update")
        return
    
    # Get the taxonomy for each issue
    for issue in issues:
        time.sleep(0.1)
        genome = parse_genome(issue, warn=False)
        custom = get_custom_fields(issue)
        try:
            accession = custom["GCA number"]["value"]
            full_name = custom["Experimental Organisms"]["value"]
            organism_abbrev = make_organism_abbrev(full_name)
            print("New organism_abbrev for %s: %s (from %s)" % (issue.id, organism_abbrev, full_name))
            if organism_abbrev in current_orgs:
                print("New organism name already exists for %s: %s (from %s)" % (issue.id, organism_abbrev, full_name))
                continue
        except:
            print("Could not generate an organism_abbrev for issue %d (%s)" % (issue.id, issue.subject))
            continue

def get_current_orgs(list_path):
    orgs = dict()
    if not list_path: return orgs

    # Column that contains the organism_abbrev
    org_column = 3

    with open(list_path, "r") as orgs_fh:
        print("Try to read %s" % list_path)
        for line in orgs_fh:
            line.strip()
            cols = line.split("\t")
            orgs[cols[org_column]] = True
    print("%d known organism_abbrevs" % len(orgs))
    
    return orgs 
    
def make_organism_abbrev(name):
    
    items = name.split(" ")
    genus = items[0]
    species = items[1]
    strain_abbrev = "".join(items[2:])
    
    genus = re.sub("^[|]$", "", genus)
    strain_abbrev = re.sub(r"isolate|strain|breed|str\.|subspecies", "", strain_abbrev, re.IGNORECASE)
    strain_abbrev = re.sub(r"[\/\(\)]", "", strain_abbrev)
    
    organism_abbrev = genus[0].lower() + species[0:3] + strain_abbrev
    return organism_abbrev
    

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Retrieve metadata from Redmine')
    
    parser.add_argument('--key', type=str, required=True,
                help='Redmine authentification key')
    parser.add_argument('--output_dir', type=str,
                help='Output_dir')
    # Choice
    parser.add_argument('--get', choices=['genomes', 'rnaseq', 'dnaseq', 'organism_abbrev'], required=True,
                help='Get genomes, rnaseq, or dnaseq issues')
    # Optional
    parser.add_argument('--build', type=int,
                help='Restrict to a given build')
    parser.add_argument('--update_redmine', type=int,
                help='Actually update Redmine for the organism_abbrev (dry run by default)')
    parser.add_argument('--species_list', type=str,
                help='A list of existing organism_abbrevs to ensure uniqueness')
    args = parser.parse_args()
    
    # Start Redmine API
    redmine = Redmine(url, key=args.key)
    
    # Choose which data to retrieve
    if args.get == 'genomes':
        if args.output_dir:
            retrieve_genomes(redmine, args.output_dir, args.build)
        else:
            print("Need --output_dir")
            return
    elif args.get == 'rnaseq':
        print("RNA-Seq Redmine retrieval to be implemented")
    elif args.get == 'dnaseq':
        print("DNA-Seq Redmine retrieval to be implemented")
    elif args.get == 'organism_abbrev':
        add_genome_organism_abbrev(redmine, build=args.build, update=args.update_redmine, species_list_path=args.species_list)
    else:
        print("Need to say what data you want to --get: genomes? rnaseq? dnaseq? organism_abbrev?")

if __name__ == "__main__":
    main()
