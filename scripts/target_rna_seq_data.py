###########################
# target_rna_seq_data.py  #
###########################
"""
Adds data to a reference server that goes with the UCSC RNA recompute
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import os
import shutil
import json
import google.protobuf.text_format as tf
import pysam
import numpy
import utils
import sys
import generate_gff3_db
import tempfile
import zipfile
import csv
import datetime
import string
import sqlite3
import urllib
from collections import defaultdict
import collections

utils.ga4ghImportGlue()

# We need to turn off QA because of the import glue
import ga4gh.datarepo as datarepo  # NOQA
import ga4gh.datamodel.references as references  # NOQA
import ga4gh.datamodel.datasets as datasets  # NOQA
import ga4gh.datamodel.variants as variants  # NOQA
import ga4gh.datamodel.reads as reads  # NOQA
import ga4gh.datamodel.ontologies as ontologies  # NOQA
import ga4gh.datamodel.sequence_annotations as sequenceAnnotations  # NOQA
import ga4gh.datamodel.bio_metadata as biodata  # NOQA
import ga4gh.datamodel.rna_quantification as rnaQuantifications  # NOQA
import ga4gh.repo.rnaseq2ga as rnaseq2ga


# save_files_locally()
# Requires wget
def save_files_locally(data):
    print("Gonna download {} indexes!".format(len(data)))
    for row in data:
        download_url = make_address(row['name'], os.path.basename(row['indexUrl']))
        os.system("wget {}".format(download_url))

# EBI's ftp into amazon url's
def make_address(dirname, filename):
    base = "http://"
    url = "s3.amazonaws.com/1000genomes/"
    return base + url + "phase3/data/" + dirname + "/alignment/" + filename


# parses the tcga rna quantification data in a tsv file 
# with target data and returns individual and biosample dictionaries
def parse_file_target(filename):
    ontology_dict = populate_ontology_dict('target_ontology_terms.json')
    bio_samples = []
    individuals = []
    print("Loading biodata tsv from target")
    with open(filename, 'r') as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter=str("\t"), quoting=csv.QUOTE_NONE)
        id_map = {}
        for row in reader:
            description = "{}".format(row['study'])
            info = {}
            for key in row:
                info[key] = [row[key]]
            biosample = {
                "name": row['analysis_id'],
                "description": description,
                "disease": ontology_dict[row['disease_name']],
                "created": row['uploaded'],
                "updated": row['modified'],
                "info": info
            }
            if row['filename'] == 'male':
                sex = {
                    "id": "PATO:0020001",
                    "term": "male genotypic sex",
                    "sourceName": "PATO",
                    "sourceVersion": "2015-11-18"
                }
            elif row['filename'] == 'female':
                sex = {
                    "id": "PATO:0020002",
                    "term": "female genotypic sex",
                    "sourceName": "PATO",
                    "sourceVersion": "2015-11-18"
                }
            else:
                sex = None
            individual = {
                "name": row['participant_id'],
                "description": description,
                "created": row['uploaded'],
                "updated": row['modified'],
                "species": {
                    "term": row['disease_name'],
                    "id": "NCBITaxon:9606",
                    "sourceName": "http://purl.obolibrary.org/obo",
                    "sourceVersion": "2016-02-02"},
                "sex": sex,
                "info": info
            }
            bio_samples.append(biosample)
            individuals.append(individual)

    return individuals, bio_samples


# populates a dictionary of ontology terms from a json file
# made from ontobee.org and a files  ontology terms 
def populate_ontology_dict(filename):
    json_file = open(filename)
    ontology_dict = json_file.read()
    ontology_dict = json.loads(ontology_dict)
    return (ontology_dict)


# parse list of biosamples and filter them distinctly
def filter_data(data):
    print("filtering data")
    bio_list = []
    count = 0
    for k in range(len(data)):
        d = data[k]
        name = d['name']
        distinct_name = filter(lambda x: x['name'] == name, data)
        for j, item in enumerate(distinct_name):
            if (j > 0): item['name'] = item['name'] + '-' + (j > 0) * str(j)
    print("filtering on data complete")
    return data


# disease_ontology
def initialize_disease_ontology(filename):
    ontology_dict = {}
    print("Loading biodata tsv from tcga")
    empty_dict = {}
    with open(filename, 'r') as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter=str("\t"), quoting=csv.QUOTE_NONE)
        for row in reader:
            disease_name = row['disease_name']
            ontology_dict[disease_name] = empty_dict
    print(ontology_dict)
    return ontology_dict


def add_individual_messages(new_individuals, individuals, dataset, repo):
    for individual in individuals:
        new_individual = biodata.Individual(dataset, individual['name'])
        new_individual.populateFromJson(json.dumps(individual))
        repo.insertIndividual(new_individual)
        new_individuals.append(new_individual)


def add_bio_sample_messages(new_bio_samples, bio_samples, dataset, repo, sample):
    new_bio_samples_list = []
    for bio_sample in bio_samples:
        new_bio_sample = biodata.BioSample(dataset, bio_sample['name'])
        new_bio_sample.populateFromJson(json.dumps(bio_sample))
        individual_id = bio_sample['info']
        new_bio_sample.setIndividualId(str(individual_id['analysis_id']))
        repo.insertBioSample(new_bio_sample)
        new_bio_samples_list.append(new_bio_sample)
    return new_bio_samples_list

# populates database relations with data from each person
# in the RnaQuantificationSet, RnaQuantification, and ExpressionLevel  
@utils.Timed()
def main():

    tsv_location_target = 'target.tsv'

    individuals_target, bio_samples_target = parse_file_target(tsv_location_target)

    individuals_target = filter_data(individuals_target)
    bio_samples_target = filter_data(bio_samples_target)
    repoPath = os.path.join("target_repo.db")
    repo = datarepo.SqlDataRepository(repoPath)
    reference_set = references.HtslibReferenceSet("NCBI37")
    # initialize the repo database if it still exists 
    if (os.path.isfile("repo.db") == True): os.system("rm repo.db")
    repo.open("w")
    repo.initialise()
    dataset = datasets.Dataset("rna-quantifications")
    dataset.setDescription("Quantifications of RNA from different data sources.")
    repo.insertDataset(dataset)
    reference_set_path = 'hg38.fa'
    reference_set = references.HtslibReferenceSet("GRCh38")
    reference_set.populateFromFile(reference_set_path)
    reference_set.setDescription("GRCh38: the hg38 assembly of the human genome")
    reference_set.setNcbiTaxonId(9606)
    for reference in reference_set.getReferences():
        reference.setNcbiTaxonId(9606)
    repo.insertReferenceSet(reference_set)

    print ("Inserting target individuals")
    new_individuals_target = []
    add_individual_messages(new_individuals_target, individuals_target, dataset, repo)

    print ("Inserting target biosamples")
    new_bio_samples_target = []
    new_bio_samples_list = add_bio_sample_messages(new_bio_samples_target, bio_samples_target, dataset, repo, 'target' )


    for bio_sample in new_bio_samples_list:
        dataset.addBioSample(bio_sample)

    seq_ontology = ontologies.Ontology("so-xp")
    ontology_file_path = 'so-xp-simple.obo'
    seq_ontology.populateFromFile(ontology_file_path)
    seq_ontology._id = "so-xp"
    repo.insertOntology(seq_ontology)
    repo.addOntology(seq_ontology)

    gencode_file_path = 'gencode.v25lift37.annotation.gff3.db'
    gencode = sequenceAnnotations.Gff3DbFeatureSet(dataset, "gencode_v25lift37")
    gencode.setOntology(seq_ontology)
    gencode.populateFromFile(gencode_file_path)
    gencode.setReferenceSet(reference_set)
    dataset.addFeatureSet(gencode)
    repo.insertFeatureSet(gencode)
    rna_base = 'rna_data'
    kallisto_quant_location = os.path.join(rna_base, 'sqlite/target_kallisto_rnaseq.db')
    rsem_gene_quant_location = os.path.join(rna_base, 'sqlite/target_rsem_gene_rnaseq.db')
    rsem_transcript_location = os.path.join(rna_base, 'sqlite/target_rsem_transcript_rnaseq.db')
    store_kallisto = rnaseq2ga.RnaSqliteStore(kallisto_quant_location)
    store_kallisto.createTables()
    store_rsem_gene = rnaseq2ga.RnaSqliteStore(rsem_gene_quant_location)
    store_rsem_gene.createTables()
    store_rsem_transcript = rnaseq2ga.RnaSqliteStore(rsem_transcript_location)
    store_rsem_transcript.createTables()

    directory_contents = os.listdir('rna_data/target')
    count = 1 
    for directory in directory_contents:
        if (directory is not None):
            print("Adding RNA for " + directory + " individual: " + str(count) + "/" + str(len(directory_contents)))
            inner_directory_contents = os.listdir('rna_data/target/' + directory)
            for direct in inner_directory_contents:
                tool_location = 'rna_data/target/' + directory + '/' + direct
                if direct == 'Kallisto':
                    direct = 'kallisto'
                    filename_location = os.path.join(tool_location, u'{0}.abundance.tsv'.format(str(directory)))
                    rnaseq2ga.rnaseq2ga(
                        filename_location, kallisto_quant_location, directory,
                        direct, dataset=dataset, featureType='transcript',
                        description='RNA seq quantification data from the target project using ' + direct,
                        programs=None,
                        featureSetNames='gencode_v25lift37',
                        readGroupSetNames=None,
                        bioSampleId=dataset.getBioSampleByName(directory).getLocalId())
                    continue
                else:
                    direct = 'rsem'
                    filename_location = os.path.join(tool_location, u'{0}.rsem_genes.results'.format(str(directory)))
                    rnaseq2ga.rnaseq2ga(
                        filename_location, rsem_gene_quant_location, directory,
                        direct, dataset=dataset, featureType='gene',
                        description='RNA seq quantification data from the target project using ' + direct,
                        programs=None,
                        featureSetNames='gencode_v25lift37',
                        readGroupSetNames=None,
                        bioSampleId=dataset.getBioSampleByName(directory).getLocalId())
                    filename_location = os.path.join(tool_location, u'{0}.rsem_isoforms.results'.format(str(directory)))
                    rnaseq2ga.rnaseq2ga(
                        filename_location, rsem_transcript_location, directory,
                        direct, dataset=dataset, featureType='transcript',
                        description='RNA seq quantification data from the target project using ' + direct,
                        programs=None,
                        featureSetNames='gencode_v25lift37',
                        readGroupSetNames=None,
                        bioSampleId=dataset.getBioSampleByName(directory).getLocalId())
                    continue
        count += 1
    rnaQuantificationSet = rnaQuantifications.SqliteRnaQuantificationSet(dataset,"TARGET RNA quantifications with Kallisto and RSEM")
    rnaQuantificationSet.setReferenceSet(reference_set)
    rnaQuantificationSet.populateFromFile(kallisto_quant_location)
    rnaQuantificationSet.populateFromFile(rsem_gene_quant_location)
    rnaQuantificationSet.populateFromFile(rsem_transcript_location)
    repo.insertRnaQuantificationSet(rnaQuantificationSet)
    repo.commit()
    print("output written to {}!".format(repoPath))

if __name__ == "__main__":
    main()
