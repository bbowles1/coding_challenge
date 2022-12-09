#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 15:43:36 2022

Goal: Extract sample annotations from input VCF file

Annotate variants with:
    1. Depth of sequencing coverage at the site of variation ('depth')
    2. Number of reads supporting the variant ('var_reads')
    3. Percentage of the reads supporting the variant versus reference ('var_percent')
    4. Use VEP to obtain:
        A. Gene ('gene_id' and 'gene_symbol')
        B. Type of variant ('variant_class')
        C. Variant effect ('effect' and 'hgvsc')
    5. The minor allele frequency of each variant ('minor allele freq')

@author: bbowles
"""

import pandas as pd
import numpy as np
import sys
import os
import argparse

# parse inputs, if config path is provided
parser = argparse.ArgumentParser(description='Annotate a VCF file.')
parser.add_argument("-i", "--input", type=str, required=True,
                    help='(Required) Path to input VCF file.')
parser.add_argument("-o", "--output", type=str, default='annot.tsv', nargs='?',
                    help='Path to output annotation file (annot.tsv by default).')
parser.add_argument("-v", "--verbose", type=str, default='N', const='Y', nargs='?',
                   help='Display sumary stats for depth, variant read count, and variant read percentage? Y/N')
parser.add_argument("-t", "--test", type=str, default='N', const='Y', nargs='?',
                   help='Annotate a random test subset of the input VCF? Y/N')
parser.add_argument("-c", "--cosmic", type=str, default='N', const='Y', nargs='?',
                   help='Annotate variants with COSMIC IDs? Y/N')

# print help if no arguments provided
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()

if args.verbose.lower() not in ['y','n']:
    raise Exception("Invalid input for 'verbose' option. Please provide a Y or N.")
if args.test.lower() not in ['y','n']:
    raise Exception("Invalid input for 'test' option. Please provide a Y or N.")
if args.cosmic.lower() not in ['y','n']:
    raise Exception("Invalid input for 'cosmic' option. Please provide a Y or N.")

# set up module import
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
#sys.path.append(os.path.abspath(SCRIPT_DIR))
sys.path.append(os.path.abspath(os.path.sep.join([SCRIPT_DIR, "modules"])))
import tempustools as tt


# import VCF
df = tt.VCF_to_df(args.input)

if args.test.lower() == 'y':
    df = df.sample(n=250)
elif args.test.lower() not in ['y', 'n']:
    raise Exception("Invalid input for 'test' option. Please provide a Y or N.")


# explode 'sample' annotations to new columns
df[ df.FORMAT.unique()[0].split(':') ] = df['sample'].str.split(':', expand=True)

# explode the delimited INFO column to separate df columns
colnames = [i.split('=')[0] for i in df.INFO.iloc[0].split(';')]
df[colnames] = df.INFO.apply(lambda row: ';'.join([i.split('=')[1] for i in row.split(';')])).str.split(';',expand=True)

# Expand multiallelic ccalls
df = tt.expand_multiallelic(df)

if args.verbose.lower() == 'y':
    # get depth of sequencing coverage at each site
    # this is already contained in the TC column
    print('\nSummary of sequencing depth for input variants:')
    print(df.TC.describe().to_string())
    
    
    # get supporting reads for each variant
    print('\nSummary of supporting reads for input variants:')
    print(df.TR.describe().to_string())
    
    # get VAF for each input variant
    print('\nSummary of supporting read percetange:')
    print((df.TR.astype(int) / df.TC.astype(int)).describe().to_string())
    print('\n')

# save depth, read count, read percent as cols in VCF
df['var_perc'] = df.TR.astype(int) / df.TC.astype(int)
df.rename(columns={'TC':'depth', 'TR':'var_reads'}, inplace=True)





# format an API query column in our data
df['apiq'] = df.apply(lambda row: ' '.join([str(i) for i in [row.CHROM, ' ', row.POS, " .", row.REF, row.ALT, ".", ".", "."]]), axis=1)

# run VEP query
decoded = tt.vepquery(df, 200)   


# take json output and unpack annotations
decoded = decoded.explode('transcript_consequences')
    
# extract transcript consequences
decoded = tt.annot_extract(decoded, 'transcript_consequences')

# explode values to new column
# variant type is listed as variant_class
decoded['gene_symbol'] = decoded.annot.apply(pd.Series)['gene_symbol']
decoded['gene_id'] = decoded.annot.apply(pd.Series)['gene_id']
decoded['effect'] = decoded.annot.apply(pd.Series)['consequence_terms']
decoded['hgvsc'] = decoded.annot.apply(pd.Series)['hgvsc']


# get cols of interest
decoded = decoded[['apiq', 'variant_class', 'gene_symbol', 'gene_id', 'effect', 'hgvsc']]


# run some quick checks to make sure df sizes are not changing
if decoded[['apiq','gene_id','variant_class', 'gene_symbol']].drop_duplicates().shape[0] != decoded[['apiq','gene_id']].drop_duplicates().shape[0]:
    raise Exception('Unanticipated additional rows in gene-level information')

# groupby variant info (apiq) and gene_id, then
# compress hgvsc and effect annotations into a delimited string
hgvsc = decoded.loc[decoded.hgvsc.notna()].groupby(['apiq','gene_id','variant_class', 'gene_symbol'])['hgvsc'].unique().apply(lambda x: ','.join((x))).reset_index()[['apiq','gene_id','hgvsc']]
decoded = decoded.groupby(['apiq','gene_id','variant_class', 'gene_symbol'])['effect'].unique().apply(lambda x: ','.join((x))).reset_index()

# merge hgvsc back into our slimmed dataframe
decoded = decoded.merge(hgvsc, on=['apiq','gene_id'], how='left')

# handle duplicates in effect and hgvsc cells (I believe caused by output from VEP, not concatenation process in pandas)
decoded.loc[:, 'effect'] = decoded.effect.apply(lambda x: ','.join(np.unique(x.split(','))))
decoded.loc[decoded.hgvsc.notna(), 'hgvsc'] = decoded.loc[decoded.hgvsc.notna()].hgvsc.apply(lambda x: ','.join(np.unique(x.split(','))))


# fill nan
decoded.effect.fillna('intergenic', inplace=True)
decoded.fillna('.', inplace=True)

# merge back to original df
decoded.loc[:, 'apiq'] = decoded.apiq.str.replace(' ','')
df.loc[:, 'apiq'] = df.apiq.str.replace(' ','')
df = df.merge(decoded, on='apiq').drop(columns='apiq')

# select cols of interest
cols = ['CHROM','POS', 'REF', 'ALT', 'depth', 'var_reads', 'var_perc',
        'variant_class', 'gene_symbol', 'gene_id', 'effect', 'hgvsc']
df = df[cols]

print('VEP region-matching complete.')





print('Querying allele frquency.')

# get test data
df['apiq'] = df.hgvsc.str.split(',').str[0]

# query MAF information
maf = tt.mafquery(df[['CHROM','apiq']].drop_duplicates())

if args.cosmic.lower() == 'y':
    # find cosmic entries in 'colocated variants' col
    maf['cosmic'] = maf.colocated_variants
    maf.loc[maf.cosmic.notna(), 'cosmic'] = maf.loc[maf.cosmic.notna()].cosmic.apply(lambda x: [j for j in [
        i['id'] for i in x if 'id' in i] if 'COSV' in j])

    # format column entries, grab correct ID from cols with multiple COSMIC entries    
    maf.loc[maf.cosmic.str.len()>1, 'cosmic'] = maf.loc[maf.cosmic.str.len()>1].cosmic.str[0]
    mapper = maf.loc[maf.cosmic.astype(bool) & maf.cosmic.notna()][['apiq','cosmic']]
    mapper.loc[:, 'cosmic'] = mapper.cosmic.str[0]
    mapper = mapper[['apiq','cosmic']].dropna()

    # map COSMIC entries back to df 
    mapper = dict(zip(mapper.apiq, mapper.cosmic))
    df['cosmic'] = df.apiq.map(mapper)

# drop entries within co_vars that do not contain frequencies info
maf.loc[maf.colocated_variants.notna(), 'colocated_variants'] = maf.colocated_variants.dropna().apply(
    lambda x: [i for i in x if 'frequencies' in i.keys()])

# stack all frequency data
freq = maf.set_index('apiq',drop=True).colocated_variants.apply(pd.Series).stack().apply(pd.Series).minor_allele_freq.reset_index()[['apiq', 'minor_allele_freq']]
freq = freq.dropna().drop_duplicates()

# merge MAF into df, fill NaN entries
maf = maf.merge(freq, how='left', on='apiq')
maf = maf[['apiq', 'minor_allele_freq']].drop_duplicates()

# merge back to df, fill nan
df = df.merge(maf, how='left', on='apiq')
df.loc[:, 'minor_allele_freq'] = df.minor_allele_freq.fillna(0)
df.loc[:, 'cosmic'] = df.cosmic.fillna('.')

# subset to output cols
if args.cosmic.lower() == 'y':
    df = df[ cols+['minor_allele_freq','cosmic'] ]
else:
    df = df[ cols+['minor_allele_freq'] ]

# sort data by CHROM and POS
df = tt.sort(df)

# save output
print('All annotations completed!')
df.to_csv(args.output, sep='\t', index=False)



