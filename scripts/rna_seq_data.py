####################
# rna_seq_data.py  #
####################
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
import ga4gh.datamodel.rna_quantification as rnaQuantifications #NOQA

# save_files_locally()
# Requires wget
def save_files_locally(data):
  print("Gonna download {} indexes!".format(len(data)))
  for row in data:
      download_url = make_address(row['name'], os.path.basename(row['indexUrl'])) 
      os.system("wget {}".format(download_url)) 
	

# converts the data in the merged json, which describes the locations of BAMs in
# EBI's ftp into amazon url's
def make_address( dirname, filename ):
    base = "http://"
    url = "s3.amazonaws.com/1000genomes/"
    return base + url + "phase3/data/" + dirname + "/alignment/" + filename  

# parses the gtex rna quantification data in a tsv file and returns
# individual and biosample dictionaries
def parse_file_gtex(filename):
  ontology_dict = populate_ontology_dict('gtex_ontology_terms.json')
  bio_samples = []
  individuals = []
  print("Loading biodata tsv from gtex")
  with open(filename, 'r') as tsvfile:
      reader = csv.DictReader(tsvfile,delimiter=str("\t"), quoting=csv.QUOTE_NONE)
      for row in reader:
        description = "{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}".format(
          # 1) e.g. SAMN01994019  
          row['Comment[BioSD_Sample]'],
          # 2) e.g. SRS408659
          row['Comment[ENA_SAMPLE]'],
          # 3) e.g. 931478
          row['Comment[dbGAP_SAMPLE]'],
          # 4) e.g. 590922
          row['Comment[gap_subject_id]'],
          # e.g. Homoe Sapiens  
          row['Characteristics[organism]'],
          # e.g. GTEX-R3RS
          row['Characteristics[individual]'],
          # 5) tissue location (e.g. Adrenal Gland) 
          row['Comment[original body site annotation]'],
          # 6) cell line 
          row['Characteristics[cell line]'],
          # 7) organism part (e.g. adrenal gland)
          row['Characteristics[organism part]'],
          # 8) clinical info 
          row['Characteristics[clinical information]'],
          # 9) Protocal REF (e.g. P-MTAB-41088)
          row['Protocol REF'],
          # 10) Extract Name (e.g. Solexa-135894)
          row['Extract Name'],
          # 11) Material type (e.g. total RNA)
          row['Material type'],
          # 12) e.g. sequence library (e.g. cDNA)
          row['Comment[LIBRARY_SELECTION]'],
          # 13) e.g. transcriptomic
          row['Comment[LIBRARY_SOURCE]'],
          # 14) e.g. RNA-Seq
          row['Comment[LIBRARY_STRATEGY]'],
          # 15) e.g. PAIRED 
          row['Comment[LIBRARY_LAYOUT]'],
          # 16) e.g. (note there are 2 protocol REF fields)
          row['Protocol REF'],
          # 17) e.g. BI
          row['Performer'],
          # 18) e.g. SRX261374
          row['Assay Name'],
          # 19) e.g. sequencing assay
          row['Technology type'],
          # 20) e.g. SRX261374
          row['Comment[ENA_EXPERIMENT]'],
          # 21) e.g. Illumina HiSeq 2000
          row['Comment[INSTRUMENT_MODEL]'],
          # 22) Scan Name (file name, e.g. SRR817649_1.fastq.gz)
          row['Scan Name'],
          # 23) e.g. SRR817649
          row['Comment[ENA_RUN]'],
          # 24) submitted file name
          row['Comment[SUBMITTED_FILE_NAME]'],
          # 25) e.g. HG19_Broad_variant
          row['Comment[assembly]'],
          # 26) e.g. adrenal gland
          row['Factor Value[organism part]'])
        info = {}
        for key in row:
          info[key] = [row[key]]
        # TODO update to use schemas
        disease = ontology_dict[row['Characteristics[disease]']]
        if ( bool(disease) == False): disease = None
        biosample = {
             "name": row['Source Name'],
             "description": description,
             "disease": disease,  # Ontology term
             "created": datetime.datetime.now().isoformat(),
             "updated": datetime.datetime.now().isoformat(),
             "info": info
        }
        if row['Characteristics[sex]'] == 'male':
           sex = {
               "id": "PATO:0020001",
               "term": "male genotypic sex",
               "sourceName": "PATO",
               "sourceVersion": "2015-11-18"
        }
        elif row['Characteristics[sex]'] == 'female':
          sex = {
            "id": "PATO:0020002",
            "term": "female genotypic sex",
            "sourceName": "PATO",
            "sourceVersion": "2015-11-18"
          }
        else:
          sex = None
        individual = {
               "name": row['Characteristics[individual]'],
               "description": description,
               "created": datetime.datetime.now().isoformat(),
               "updated": datetime.datetime.now().isoformat(),
               "species": {
                   "term": row['Characteristics[organism]'],
                   "id": "NCBITaxon:9606",
                   "sourceName": "http://purl.obolibrary.org/obo",
                   "sourceVersion": "2016-02-02"},
               "sex": sex,
               "info": info
        }
        bio_samples.append(biosample)
        individuals.append(individual)
  return individuals, bio_samples 

# parses the tcga rna quantification data in a tsv file 
# with target data and returns individual and biosample dictionaries
def parse_file_target(filename):
  ontology_dict = populate_ontology_dict('target_ontology_terms.json')
  bio_samples = []
  individuals = []
  print("Loading biodata tsv from target")
  with open(filename, 'r') as tsvfile:
    reader = csv.DictReader(tsvfile,delimiter=str("\t"), quoting=csv.QUOTE_NONE)
    id_map ={} 
    for row in reader:
      description = "{}".format(row['study'])
      info = {}
      for key in row:
        info[key] = [row[key]]
      biosample = {
           "name": row['sample_type_name'],
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
      else:sex = None
      individual = {
        "name" : row['participant_id'],
        "description": description,
        "created" : row['uploaded'],
        "updated" : row['modified'],
        "species": {
            "term": row['disease_name'],
            "id": "NCBITaxon:9606",
            "sourceName": "http://purl.obolibrary.org/obo",
            "sourceVersion": "2016-02-02"},
        "sex": sex,
        "info" : info 
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

# parses the tcga rna quantification data in a tsv file and returns
# individual and biosample dictionaries
def parse_file_tcga(filename):
  ontology_dict = populate_ontology_dict('tcga_ontology_terms.json')
  bio_samples = []
  individuals = []
  print("Loading biodata tsv from tcga")
  with open(filename, 'r') as tsvfile:
      reader = csv.DictReader(tsvfile,delimiter=str("\t"), quoting=csv.QUOTE_NONE)
      id_map ={}
      for row in reader:
        description = "{}, {}".format(
        row['barcode'],
		row['disease_name'])
        info = {}
        for key in row:
          info[key] = [row[key]]
        # TODO update to use schemas
        #if (row['sample_id'] == row['participant_id']):
        #  print ("{} {}".formatrow['sample_id'], row['participant_id'])
        #  assert(id_map.get(row['sample_id'] == row['participant_id']))
        #else:
        #  id_map[row['sample_id']] = row['participant_id']	
        biosample = { 
             "name": row['sample_id'],
             "description": description,
             "disease": ontology_dict[row['disease_name']],  # Ontology term
             "created": row['published'],
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
               "created": row['published'],
               "updated": row['modified'],
			   "species": {
                   "term": row['disease'],
                   "id": row['sample_id'],
                   "sourceName": "http://purl.obolibrary.org/obo/doid.owl",
                   "sourceVersion": "25:03:2016 16:27"},
               "sex": sex,
               "info": info
        }
        bio_samples.append(biosample)
        individuals.append(individual)
  return individuals, bio_samples

# parse list of biosamples and filter them distinctly
def filter_data(data):
    print ("filtering data")
    bio_list = []
    count = 0
    for k in range(len(data)):
       d = data[k]
       name = d['name']
       distinct_name = filter(lambda x: x['name'] == name, data)
       for j,item in enumerate(distinct_name):
            if ( j>0):item['name'] = item['name'] + '-' + (j>0) * str(j)
    print ("filtering on data complete")
    return data
        
# disease_ontology
def initialize_disease_ontology(filename):
  ontology_dict = {}
  print("Loading biodata tsv from tcga")
  empty_dict = {}
  with open(filename, 'r') as tsvfile:
      reader = csv.DictReader(tsvfile,delimiter=str("\t"), quoting=csv.QUOTE_NONE)
      for row in reader:
        disease_name = row['disease_name']	
        ontology_dict[disease_name] = empty_dict
  print (ontology_dict)
  return ontology_dict

def add_individual_messages(new_individuals, individuals, dataset, repo):
    for individual in individuals:
      new_individual = biodata.Individual(dataset, individual['name'])
      new_individual.populateFromJson(json.dumps(individual))
      repo.insertIndividual(new_individual)
      new_individuals.append(new_individual)

def add_bio_sample_messages(new_bio_samples, bio_samples, dataset, repo, sample ):
    for bio_sample in bio_samples:
        new_bio_sample = biodata.BioSample(dataset, bio_sample['name'])
        new_bio_sample.populateFromJson(json.dumps(bio_sample))
        individual_id = bio_sample['info']
        if ( sample == 'gtex' ):new_bio_sample.setIndividualId(str(individual_id['Characteristics[individual]']))
        else: new_bio_sample.setIndividualId(str(individual_id['participant_id']))
        repo.insertBioSample(new_bio_sample)
     
 
# populates database relations with data from each person
# in the RnaQuantificationSet, RnaQuantification, and ExpressionLevel  
@utils.Timed()
def main():
    tsv_location_gtex = 'gtex.tsv'
    tsv_location_tcga = 'tcga.tsv'
    tsv_location_target = 'target.tsv'
    individuals_gtex, bio_samples_gtex = parse_file_gtex(tsv_location_gtex)
    individuals_tcga, bio_samples_tcga = parse_file_tcga(tsv_location_tcga)
    individuals_target, bio_samples_target = parse_file_target(tsv_location_target)
    individuals_gtex = filter_data(individuals_gtex)
    bio_samples_gtex = filter_data(bio_samples_gtex)
    individuals_tcga = filter_data(individuals_tcga)
    bio_samples_tcga = filter_data(bio_samples_tcga)
    individuals_target = filter_data(individuals_target)
    bio_samples_target = filter_data(bio_samples_target)
    repoPath = os.path.join("repo.db")
    repo = datarepo.SqlDataRepository(repoPath)
    # initialize the repo database if it still exists 
    if ( os.path.isfile("repo.db") == True ):os.system( "rm repo.db" )
    repo.open("w")
    repo.initialise()
    dataset = datasets.Dataset("1kgenomes")
    dataset.setDescription("Variants from the 1000 Genomes project and GENCODE genes annotations")
    repo.insertDataset(dataset) 
    print("Inserting gtex individuals")
    new_individuals_gtex = []
    add_individual_messages(new_individuals_gtex, individuals_gtex, dataset, repo)
    print ("Inserting tcga individuals")
    new_individuals_tcga = []
    add_individual_messages(new_individuals_tcga, individuals_tcga, dataset, repo)
    print ("Inserting target individuals")
    new_individuals_target = []
    add_individual_messages(new_individuals_target, individuals_target, dataset, repo)
    print ("Inserting gtex biosamples")
    new_bio_samples_gtex = []
    add_bio_sample_messages(new_bio_samples_gtex, bio_samples_gtex, dataset, repo, 'gtex' )
    print ("Inserting tcga biosamples")
    new_bio_samples_tcga = []
    add_bio_sample_messages(new_bio_samples_tcga, bio_samples_tcga, dataset, repo, 'tcga' )
    print ("Inserting target biosamples")
    new_bio_samples_target = []
    add_bio_sample_messages(new_bio_samples_target, bio_samples_target, dataset, repo, 'target' )
    repo.commit()
    print ( "output written to {}!".format(repoPath))			

if __name__ == "__main__":
    main()
	#test_ontospy()
