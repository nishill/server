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
import ontospy
utils.ga4ghImportGlue()

# We need to turn off QA because of the import glue
import ga4gh.datarepo as datarepo  # NOQA
import ga4gh.datamodel.references as references  # NOQA
import ga4gh.datamodel.datasets as datasets  # NOQA
import ga4gh.datamodel.variants as variants  # NOQA
import ga4gh.datamodel.reads as reads  # NOQA
import ga4gh.datamodel.ontologies as ontologies  # NOQA
import ga4gh.datamodel.sequenceAnnotations as sequenceAnnotations  # NOQA
import ga4gh.datamodel.bio_metadata as biodata  # NOQA
import ga4gh.datamodel.rna_quantification as rnaQuantifications #NOQA


def test_ontospy():
	print ("testing")
	g = ontospy.Graph("http://xmlns.com/foaf/0.1/")
	#g = ontospy.Graph()


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
        biosample = {
             "name": row['Source Name'],
             "description": description,
             "disease": row['Characteristics[disease]'],  # Ontology term
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

# populates a dictionary of ontology terms from a json file
# made from ontobee.org and tcga ontology terms 
def populate_ontology_dict():
  #tsv_location_tcga = 'tcga.tsv'
  #ontology_dict = initialize_disease_ontology(tsv_location_tcga)   
  json_file = open('out.json')
  ontology_dict = json_file.read()
  ontology_dict = json.loads(ontology_dict)
  return (ontology_dict)

#the tcga rna quantification data in a tsv file and returns
# individual and biosample dictionaries
def parse_file_tcga(filename):
  ontology_dict = populate_ontology_dict()
  bio_samples = []
  individuals = []
  print("Loading biodata tsv from tcga")
  with open(filename, 'r') as tsvfile:
      reader = csv.DictReader(tsvfile,delimiter=str("\t"), quoting=csv.QUOTE_NONE)
      #id_map ={}
      for row in reader:
        description = "{}{}".format(
        row['barcode'],
		row['disease_name'])
        info = {}
        for key in row:
          info[key] = [row[key]]
        # TODO update to use schemas
        #if (row['sample_id'] == row['participant_id']):
	#		print ("{} {}".formatrow['sample_id'], row['participant_id'])
          
         #assert(id_map.get(row['sample_id'] == row['participant_id']))
        #else:
        #  id_map[row['sample_id']] = row['participant_id']	
        biosample = { 
             "name": row['sample_id'],
             "description": description,
             # diagnosis --> SampleCharacteristics, or SampleFeatures, or SampleOntologies 
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
def filter_biosamples_tcga(biosamples):
    print ("filtering")
    bio_list = []
    count = 0
    for k in range(len(biosamples)):
       bs = biosamples[k]
       name = bs['name']
       distinct_name = filter(lambda x: x['name'] == name, biosamples)
       for j,item in enumerate(distinct_name):
            item['name'] = item['name'] + '-' + (j>0) * str(j)
    return biosamples 

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

  return ontology_dict

# specification_list
def specification_list():
	spec_list = [
	# Uterine Corpus Endometrioid Carcinoma 
    {"id": "DOID:2871","term":"endometrial carcinoma" ,"source_name":"Disease Ontology","source_version":"25:03:2016 16:27"},
	# Adrenocortical carcinoma
    {"id": "DOID:3948" ,"term":"Adrenal cortical carcinoma" ,"source_name":"Disease Ontology" ,"source_version":"25:03:2016 16:27" },
    # Thymoma 
    {"id": "DOID:3275","term":"thymoma" ,"source_name":"Disease Ontology" ,"source_version":"25:03:2016 16:27" },
    # Ovarian serous cystadenocarcinoma 
	{"id":"DOID:5746" ,"term":"ovarian serous cystadenocarcinoma" ,"source_name":"Disease Ontology" ,"source_version": "25:03:2016 16:27" },
	# Kidney renal clear cell carcinoma
	{"id":"DOID:4467" ,"term":" clear cell kidney carcinoma" ,"source_name":"Disease Ontology" ,"source_version":"25:03:2016 16:27" },
	# Thyroid carcinoma 
	{"id":"DOID:3963" ,"term":"Thyroid carcinoma" ,"source_name": "Disease Ontology" ,"source_version":"25:03:2016 16:27" },
	# Cervical squamous cell carcinoma and endocervical adenocarcinoma
	{"id": "DOID:3744","term":"Cervical squamous cell carcinoma and endocervical adenocarcinoma" ,"source_name":"Disease Ontology" ,"source_version":"25:03:2016 16:27"},
	# Uveal Melanoma 
	{"id":"DOID:6039" ,"term":"melanoma of Uvea" ,"source_name":"Disease Ontology" ,"source_version":"25:03:2016 16:27" },
	# Stomach adenocarcinoma 
	{"id":"DOID:3717", "term":"Stomach adenocarcinoma", "source_name":"Disease Ontology" ,"source_version":"25:03:2016 16:27" },
	# Lung squamous cell carcinoma 
	{"id":"DOID:3907" , "term":"squamous cell carcinoma of lung", "source_name":"Disease Ontology", "source_version":"25:03:2016 16:27" },
	# Colon adenocarcinoma
	{"id":"DOID:234" ,"term":"adenocarcinoma of the colon" ,"source_name":"Disease Ontology" ,"source_version":"25:03:2016 16:27" },
	# Liver hepatocellular carcinoma
	{"id":"DOID:686" ,"term":"Liver and Intrahepatic bile duct carcinoma" ,"source_name":"Disease Ontology" ,"source_version":"25:03:2016 16:27" },
	# Esophageal carcinoma
	{"id":"DOID:1107" ,"term":"cancer of esophagus" ,"source_name":"Disease Ontology" ,"source_version":"25:03:2016 16:27" },
	# Mesothelioma
	{"id":"DOID:4486" ,"term":"mixed Mesothelioma" ,"source_name":"Disease Ontology" ,"source_version":"25:03:2016 16:27" },
	# Head and Neck squamous cell carcinoma
	{"id":"DOID:5520" ,"term":"Head and Neck squamous cell carcinoma" ,"source_name":"Disease Ontology" ,"source_version":"25:03:2016 16:27" },
	# Pheochromocytoma and Paraganglioma 
	{"id":"DOID:0050892" ,"term":"Pheochromocytoma and Paraganglioma" ,"source_name":"Disease Ontology" ,"source_version":"25:03:2016 16:27" },
	# Breast invasive carcinoma
	{"id":"DOID:3610" ,"term":"Invasive mucinous breast carcinoma" ,"source_name":"Disease Ontology" ,"source_version":"25:03:2016 16:27"},
	# Uterine Carcinosarcoma
	{"id":"DOID:6171" ,"term":"mixed mullerian sarcoma of uterus" ,"source_name":"Disease Ontology" ,"source_version":"25:03:2016 16:27" },
	# Rectum adenocarcinoma
	{"id": "DOID:1996","term":"adenocarcinoma of rectum" ,"source_name":"Disease Ontology" ,"source_version":"25:03:2016 16:27" },
	# Skin Cutaneous Melanoma
	{"id": "DOID:8923","term":"skin melanoma" ,"source_name":"Disease Ontology" ,"source_version":"25:03:2016 16:27" },
	# Prostate adenocarcinoma
	{"id":"DOID:2526" ,"term":"Prostate adenocarcinoma" ,"source_name":"Disease Ontology" ,"source_version": "25:03:2016 16:27"},
	# Kidney Chromophobe
	{"id":"DOID:4471" ,"term":"Chromophobe carcinoma of kidney" ,"source_name":"Disease Ontology" ,"source_version": "25:03:2016 16:27"},
	# Pancreatic adenocarcinoma
	{"id": "DOID:4074" ,"term":"adenocarcinoma of the pancreas" ,"source_name":"Disease Ontology" ,"source_version": "25:03:2016 16:27"},
	# Kidney renal papillary cell carcinoma
	{"id": "DOID:4465/DOID:8063" ,"term": "Papillary renal cell carcinoma" ,"source_name":"Disease Ontology" ,"source_version": "25:03:2016 16:27"},
	# Brain Lower Grade Glioma
	{"id":"DOID:0060108" ,"term":"lower grade glioma" ,"source_name":"Disease Ontology" ,"source_version": "25:03:2016 16:27"},
	# Lymphoid Neoplasm Diffuse Large B-cell Lymphoma
	{"id":"DOID:0050745" ,"term":"diffuse large B-cell lymphoma" ,"source_name":"Disease Ontology" ,"source_version": "25:03:2016 16:27"},
	# Controls 
	{"id":"" ,"term":"" ,"source_name":"Disease Ontology" ,"source_version": "25:03:2016 16:27"},
	# Lung adenocarcinoma
	{"id":"DOID:3910" ,"term":"adenocarcinoma of lung" ,"source_name":"Disease Ontology" ,"source_version": "25:03:2016 16:27"},
	# Bladder Urothelial Carcinoma 
	{"id":"DOID:4006" ,"term":"urothelial bladder carcinoma" ,"source_name":"Disease Ontology" ,"source_version": "25:03:2016 16:27"},
	# Cholangiocarcinoma 
	{"id": "DOID:4947","term":"cholangiosarcoma" ,"source_name":"Disease Ontology" ,"source_version": "25:03:2016 16:27"},
	# Testicular Germ Cell Tumors
	{"id":"DOID:5557" ,"term":"germ cell tumor of testis" ,"source_name":"Disease Ontology" ,"source_version": "25:03:2016 16:27"},
	# Glioblastoma multiforme
	{"id":"DOID:3068" ,"term":"Glioblastoma multiforme" ,"source_name":"Disease Ontology" ,"source_version": "25:03:2016 16:27"},
	# Sarcoma
	{"id":"DOID:1115" ,"term":"tumor of soft tissue and skeleton" ,"source_name":"Disease Ontology" ,"source_version": "25:03:2016 16:27"} ]	
	return spec_list


# populates database relations with data from each person
# in the RnaQuantificationSet, RnaQuantification, and ExpressionLevel  
@utils.Timed()
def main():
    #index_list_path = 'merged.json'
    #download_indexes = False
    tsv_location_tcga = 'tcga.tsv'
    #individuals_gtex, bio_samples_gtex = parse_file_gtex(tsv_location_gtex)
    individuals_tcga, bio_samples_tcga = parse_file_tcga(tsv_location_tcga)
    #ontology_dict = populate_ontology_dict()
    #print (ontology_dict)
	#filter_biosamples_tcga(individuals_tcga)
    #filter_biosamples_tcga(bio_samples_tcga)
    print (bio_samples_tcga)
    #repoPath = os.path.join("repo.db")
    #repo = datarepo.SqlDataRepository(repoPath)
    #if ( os.path.isfile("repo.db") == True ):
    #    os.system( "rm repo.db" )
    #repo.open("w")
    #repo.initialise()
    #dataset = datasets.Dataset("1kgenomes")
    #dataset.setDescription("Variants from the 1000 Genomes project and GENCODE genes annotations")
    #repo.insertDataset(dataset) 
    #print("Inserting tcga individuals")
    #new_individuals= []
    #for individual in individuals_tcga:
    #  new_individual = biodata.Individual(dataset, individual['name'])
    #  new_individual.populateFromJson(json.dumps(individual))
    #  repo.insertIndividual(new_individual)
    #  new_individuals.append(new_individual)



	#print("Load list of read group sets")
    #with open (index_list_path) as merged:
    #    data = json.load(merged)
    #    print("Found {} read group sets".format(len(data)))
    #    if download_indexes:
    #      save_files_locally(data)
    #    # TODO might have to do something smart about pointing to different index locations
    #    for row in data:
    #        print("Adding {}".format(row['name']))
    #        download_url = make_address(row['name'], os.path.basename(row['dataUrl']))
    #        name = row['name']
    #        read_group_set = reads.HtslibReadGroupSet(dataset, name)
    #        read_group_set.populateFromFile(download_url, row['indexUrl'])
    #        # could optimize by storing biodata in a name:value dictionary
    #        for read_group in read_group_set.getReadGroups():
    #          for bio_sample in new_bio_samples:
    #              if bio_sample.getLocalId() == read_group.getSampleName():
    #                  read_group.setBioSampleId(bio_sample.getId())
    #        read_group_set.setReferenceSet(reference_set)
    #        repo.insertReadGroupSet(read_group_set)
    
    
    #repo.commit()
    #print ( "database filled!")			

if __name__ == "__main__":
    main()
	#test_ontospy()
