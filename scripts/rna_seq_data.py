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

# parses the tcga rna quantification data in a tsv file and returns
# individual and biosample dictionaries
def parse_file_tcga(filename):
  bio_samples = []
  individuals = []
  print("Loading biodata tsv from tcga")
  with open(filename, 'r') as tsvfile:
      reader = csv.DictReader(tsvfile,delimiter=str("\t"), quoting=csv.QUOTE_NONE)
      for row in reader:
        description = "{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}".format(
	    # 1) e.g. TCGA 
        row['study'],
        # 2) e.g. TCGA-DQ-5630-01A-01R-1873-07  
        row['barcode'],
        # 3) e.g. HNSC
        row['disease'],
		# 4) e.g. Head and Neck squamous cell carcinoma
		row['disease_name'],
		# 5) e.g. TP
		row['sample_type'],
		# 6) e.g. 01
		row['sample_type_name'],
		# 7) e.g. RNA
		row['analyte_type'],
		# 8) e.g. RNA-Seq
		row['library_type'],
		# 9) e.g. UNC-LCCC
		row['center'],
		# 10) e.g. UNC-LCCC 
		row['center_name'],
		# 11) e.g. ILLUMINA 
		row['platform'],
		# 12) e.g. Illumina
		row['platform_name'],
		# 13) e.g. unaligned
		row['assembly'],
		# 14) e.g. UNCID_2200798.13c839d1-f77b-4bde-9bd5-8d5165d94346.111122_UNC15-SN850_0148_AC05D7ACXX_8_GGCTAC.tar.gz
		row['filename'],
		# 15) e.g. 5747943025
		row['files_size'],
		# 16) e.g. d81203da215be180128f260b788900b5
		row['checksum'],
		# 17) e.g. 5557a728-1827-4aff-b28b-f004d835f9d6
		row['analysis_id'],
		# 18) e.g. 13c839d1-f77b-4bde-9bd5-8d5165d94346 
		row['aliquot_id'],
		# 19) e.g. 91a712b5-a724-48d0-ba06-4f6857683463
		row['participant_id'],
		# 20) e.g. a0cdb208-e66f-4ad7-86fa-5077145a4a68
		row['sample_id'],
		# 21) e.g. A5 
		row['tss_id'],
		# 22) e.g.  SRS128873 
		row['sample_accession'],
		# 23) e.g. 2013-09-27
		row['published'],
		# 24) e.g. 2013-09-25
		row['uploaded'],
		# 25) e.g. 2013-09-27
		row['modified'],
		# 26) e.g. Live
		row['state'],
		# 27) ...
        row['reason'])
        info = {}
        for key in row:
          info[key] = [row[key]]
        # TODO update to use schemas
        biosample = { 
             "name": row['sample_id'],
             "description": description,
             "disease": row['disease'],  # Ontology term
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

# parse list of individuals and filter them distinctly
#def filter_individuals_tcga(individuals):

#	print ("filtering")

	
# main():
# populates database relations with data from each person
# in the RnaQuantificationSet, RnaQuantification, and ExpressionLevel  
@utils.Timed()
def main():
    #index_list_path = 'merged.json'
    #download_indexes = False
    tsv_location_gtex = 'gtex.tsv'
    individuals_gtex, bio_samples_gtex = parse_file_gtex(tsv_location_gtex)
    tsv_location_tcga = 'tcga.tsv'
    individuals_tcga, bio_samples_tcga = parse_file_tcga(tsv_location_tcga)
    #print (individuals_tcga)
	#repoPath = os.path.join("repo.db")
    #repo = datarepo.SqlDataRepository(repoPath)
    #if ( os.path.isfile("repo.db") == True ):
    #    os.system( "rm repo.db" )
    #repo.open("w")
    #repo.initialise()
    #dataset = datasets.Dataset("1kgenomes")
    #dataset.setDescription("Variants from the 1000 Genomes project and GENCODE genes annotations")
    #repo.insertDataset(dataset) 
    #print("Inserting biosamples")
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
