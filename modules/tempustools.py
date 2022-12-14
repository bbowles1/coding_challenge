#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 07:25:49 2022

@author: bbowles
"""

import pandas as pd
import numpy as np
import requests
import json
import sys
from math import ceil

def VCF_to_df(input_path):
    
    if type(input_path) != str:
        print('Please provide a string file path as input.')
    
    # get header
    with open(input_path, "r") as file:
    
        for line in file.readlines():
            if '#CHROM' in line:
                header = line.replace('\n','').replace('#','').split('\t')
                break
    
    dtypes = {'CHROM':str,
              'POS':int,
              'ID':str,
              'REF':str,
              'ALT':str,
              'QUAL':float,
              'FILTER':str,
              'INFO':str,
              'FORMAT':str}
    
    # import overlap file as dataframe
    df = pd.read_csv(input_path, sep='\t', comment='#', names=header, dtype=dtypes)
    
    return(df)




def expand_multiallelic(df):
    
    # takes a pandas dataframe with VCF 4.0 columns
    # explodes multiallelic annotations to new dataframe rows
    
    if df.loc[df.REF.str.contains(',') | df.ALT.str.contains(',')].empty:
        print('No rows to explode!')
        return(df)

    else:
        
        # create empty dataframe to hold exploded VCF info
        expanded = pd.DataFrame(columns=df.columns)
        
        for column in ['REF','ALT']:
        
            # iterate through multiallelic annotations of fixed size (n=2,3,... annotations in a column)
            for length in df[column].str.split(',').str.len().unique():
                            
                if length > 1:
                    
                    # subset data
                    multi = df.loc[df[column].str.contains(',')]
                    multi = multi.loc[multi[column].str.split(',').str.len()==length]
                    
                    # duplicate index N times
                    prevdex = multi.index.repeat(length)
            
                    def boolcheck(array):
                        # check if a series contains only 1 unique value
                        # return the single unique value if so
                        if array.size==1:
                            return(array[0])
                        else:
                            return(False)
                            
                    # select all columns with incorrect number of delimiters - do not explode these
                    icols = [i for i in multi.columns if boolcheck(multi[i].astype(str).str.count(',').unique())!=length-1 ]
                    icols = icols + [i for i in ['FORMAT', 'INFO', 'GL'] if i not in icols]
                    dcols = [i for i in df.columns if i not in icols]
            
                    # set non-delimited columns as index
                    multi.set_index(icols, inplace=True)
            
                    # split all remaining columns apart
                    multi = pd.DataFrame( [multi[col].str.split(',') for col in multi.columns] ).transpose()
                    multi = multi.apply(pd.Series.explode).reset_index()
                    multi.set_index(prevdex, inplace=True)
            
                    # concat exploded columns to expanded dataframe
                    expanded = pd.concat([expanded, multi])
                    
    # merge all exploded variants back to original dataframe, and sort
    df = pd.concat( [ df.loc[~( df.REF.str.contains(',') | df.ALT.str.contains(','))], expanded ])
    df.sort_index(inplace=True)
    
    return(df)



def chunker(seq, size):
    # chunk data into n sizes
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))


def collapse_list(series):
    # collapse lists in a pandas series to a single str
    return ( ','.join(
        (series.apply(pd.Series).stack().reset_index(drop=True).dropna().unique())))
    
def checktype(series):
    # check types in a pandas series, return dtype if there is only one dtype
    types = series.apply(type).unique()
    if types.size == 1:
        return(types[0])
    else:
        return(object)
    
    
# now vectorize it
def annot_extract(df, target_cell, newname='annot'):
    
    # extract nested json annotations from a target cell
    # return copy of dataframee with an added 'newname' col
            
    # apply json_norm to cell
    col = df.loc[df[target_cell].notna()][target_cell].apply(pd.json_normalize)
    
    # use list comprehension + lambda fn to obtain json info in a single series
    annot = pd.Series (
        col.apply( lambda cell:
        {i : (';'.join(cell[i].astype(str).unique()) ) 
         if checktype( cell[i])!=list else collapse_list(cell[i]) 
         for i in cell.columns} ) )
    
    # add info to new column, propagate nan values
    df[newname] = np.nan
    df.loc[df[target_cell].notna(), newname] = annot
    
    return(df)
            
    
def vepquery(df, size=200):
    
    nchunks = ceil(df.shape[0] / size)
    n=1
    decoded=pd.DataFrame()
    
    for chunk in chunker(df, size):
        
        print('Querying VEP: batch', n, '/', nchunks )
    
        # format a query of our variants
        query = chunk.apiq.to_list()
        query = json.dumps(
            {"variants":query}
            )
        
        params = {'transcript_version':'true', 'variant_class':'true', 'hgvs':'true'}
        
        # query the VEP API
        server = "https://grch37.rest.ensembl.org"
        ext = "/vep/homo_sapiens/region"
        headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
        r = requests.post(server+ext, headers=headers, data=query, params=params)
        if not r.ok:
          r.raise_for_status()
          sys.exit()
         
        # decode output
        if decoded.empty:
            
            # decode to pandas df
            decoded = pd.DataFrame(r.json()).rename(columns={'input':'apiq'})
            
        else:
            decoded = pd.concat([decoded,
                                 pd.DataFrame(r.json()).rename(columns={'input':'apiq'})])
            
        n+=1
        
    return (decoded)



def mafquery(df, size=200):
    
    nchunks = ceil(df.shape[0] / size)
    n=1
    decoded=pd.DataFrame()
    
    for chunk in chunker(df, size):
        
        print('Querying MAF: batch', n, '/', nchunks )
    
        # format a query of our variants
        query = list(chunk.apiq.unique())
        query = json.dumps(
            {"hgvs_notations":query}
            )
        
        # query the VEP API
        server = "https://grch37.rest.ensembl.org"
        ext = "/vep/homo_sapiens/hgvs"
        headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
        r = requests.post(server+ext, headers=headers, data=query)
        if not r.ok:
          r.raise_for_status()
          sys.exit()
         
        # decode output
        if decoded.empty:
            
            # decode to pandas df
            decoded = pd.DataFrame(r.json()).rename(columns={'input':'apiq'})
            
        else:
            decoded = pd.concat([decoded,
                                 pd.DataFrame(r.json()).rename(columns={'input':'apiq'})])
            
        n+=1
        
    return (decoded)


def sort(input_df): # takes BED and VCF type files
    
    # sorts an input df with CHROM, POS columns. Does not work for BED files.
    if all([i in input_df.columns for i in ['CHROM','POS']]):
    
        # flag for CHROM formatting
        chrflag = False
        
        # check POS format
        input_df.loc[:, 'POS'] = input_df.POS.astype(int)
            
        if input_df.CHROM.str.contains('chr').any():
            chrflag = True
        
        # strip chr format from input_df
        input_df.loc[:, 'CHROM'] = input_df.CHROM.astype(str).str.replace('chr','')
        
        # create a list of all 23 chromosome, impose order on input df
        chrom_list = [str(i) for i in (list(range(1,23)))]
        chrom_list.append('X')
        chrom_list.append('Y')
        index_order = [i for i in chrom_list if i in input_df.CHROM.astype(str).str.strip('chr').to_list()]
        input_df.loc[:, 'CHROM'] = input_df.CHROM.astype(str)
        input_df.set_index('CHROM', inplace=True)
        input_df = input_df.loc[index_order, :]
        
        if chrflag == True:
            # pop out chrom index, reset to "chr#" format
            input_df.reset_index(inplace=True)
            input_df.loc[:, 'CHROM'] = 'chr'+ input_df.CHROM
        else:
            input_df.reset_index(inplace=True)
    
        # use groupby to sort by POS
        input_df = input_df.groupby('CHROM').apply(lambda x: x.sort_values(by='POS', ascending=True)).reset_index(drop=True)
        
        return(input_df)
    
    # sorts an input df with CHROM, POS columns. Does not work for VCF files.
    elif all([i in input_df.columns for i in ['chrom','chromStart', 'chromEnd']]):
            
        # flag for chrom formatting
        chrflag = False
        
        # check POS format
        input_df.loc[:, 'chromStart'] = input_df.chromStart.astype(int)
        input_df.loc[:, 'chromEnd'] = input_df.chromEnd.astype(int)
            
        if input_df.chrom.str.contains('chr').any():
            chrflag = True
            
        # strip chr format from input_df
        input_df.loc[:, 'chrom'] = input_df.chrom.astype(str).str.replace('chr','')
        
        # create a list of all 23 chromosomes
        chrom_list = [str(i) for i in (list(range(1,23)))]
        chrom_list.append('X')
        chrom_list.append('Y')
        index_order = [i for i in chrom_list if i in input_df.chrom.astype(str).str.strip('chr').to_list()]
        input_df.loc[:, 'chrom'] = input_df.chrom.astype(str)
        input_df.set_index('chrom', inplace=True)
        input_df = input_df.loc[index_order, :]
        
        if chrflag == True:
            # pop out chrom index, reset to "chr#" format
            input_df.reset_index(inplace=True)
            input_df.loc[:, 'chrom'] = 'chr'+ input_df.chrom
        else:
            input_df.reset_index(inplace=True)
        
        # use groupby to sort by chromStart, chromEnd
        input_df = input_df.groupby('chrom').apply(lambda x: x.sort_values(by=(['chromStart','chromEnd']), ascending=True)).reset_index(drop=True)
        
        return(input_df)
    
    else:
        print('Could not detect input format - missing required VCF (chrom, POS) or BED (chrom, chromStart, chromEnd) columns!')
        return(input_df)


