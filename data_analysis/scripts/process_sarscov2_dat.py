import argparse
import pandas as pd
import numpy as np
import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import string
import itertools
import os
import matplotlib.gridspec as gridspec
from collections import Counter

def process_vcf(vcf_path, filter_allow):
    sample_name = \
        'CoV_' + vcf_path.split('/')[-1].split('_')[1]
    with open(vcf_path, 'r') as fp:
        for line in fp:
            if line[0:2] != '##':
                header=line.split('\n')[0].split('\t')
                break
    vcf = \
        pd.read_csv(vcf_path, sep='\t', header=None, comment='#')
    vcf.columns = header
    vcf['SAMPLE'] = sample_name
    vcf['FILTER_PASS'] = \
        vcf['FILTER'].isin(filter_allow)
    # todo make this better
    vcf['AF'] = \
        vcf.loc[:,'INFO'].str.split(';').apply(lambda k: 
            [i.split('=')[1] for i in k if i.split('=')[0]=='AF'][0]) 
    vcf['DP4'] = \
        vcf.loc[:,'INFO'].str.split(';').apply(lambda k: 
            [sum([int(j) for j in i.split('=')[1].split(',')]) 
            for i in k if i.split('=')[0]=='DP4'][0])
    vcf['AF'] = vcf.loc[:,'AF'].astype(float)
    vcf['TYPE'] = \
        np.where((vcf['REF'].str.len()==1) & 
            (vcf['REF'].str.len()==1), 'SNP','INDEL')
    vcf = vcf[['SAMPLE', '#CHROM', 
        'POS','REF', 'ALT', 'AF', 'DP4', 
        'FILTER_PASS', 'TYPE']]
    return(sample_name, vcf)


def get_pair_dat(dat_1, dat_2):
    pair_dat = dat_1.merge(dat_2, 
        left_on=['POS', 'REF', 'ALT'], 
        right_on=['POS', 'REF', 'ALT'], 
        how='outer')
    pair_dat.loc[:,['AF_x', 'AF_y']] = \
        pair_dat.loc[:,['AF_x', 'AF_y']].fillna(0)
    pair_dat.loc[:,['AF_x', 'AF_y']] = \
        pair_dat.loc[:,['AF_x', 'AF_y']].astype(float)
    pair_dat['SAMPLE_x'] = pair_dat.SAMPLE_x.dropna().iloc[0]
    pair_dat['SAMPLE_y'] = pair_dat.SAMPLE_y.dropna().iloc[0]
    return(pair_dat)



def get_transmission_pairs(dat):
    import ast
    dat['donor_pairs'] = \
        dat['donor_pairs'].apply(lambda k: ast.literal_eval(k))
    dat['recip_pairs'] = \
        dat['recip_pairs'].apply(lambda k: ast.literal_eval(k))
    # key is transmission pair ID, value is donor 
    donor_dict = {}
    # key is transmission pair ID, value is recipient
    recip_dict = {}
    for idx, row in dat.iterrows():
        for item in row['donor_pairs']:
            donor_dict[item] = row['sample_name']
        for item in row['recip_pairs']:
            recip_dict[item] = row['sample_name']
    if set(donor_dict.keys()) != set(recip_dict.keys()):
        raise Exception('incomplete transmission pair data in metadata')
    pair_dict = {}
    for key,value in donor_dict.items():
        pair_dict[key] = (value, recip_dict[key])
    return(pair_dict)

     
def get_clonal_diffs(pair_dat, threshold=0.05):
    return(pair_dat[((pair_dat['AF_x'] < threshold) & 
            (pair_dat['AF_y'] > 1-threshold)) | 
        ((pair_dat['AF_x'] > 1-threshold) & 
            (pair_dat['AF_y'] < threshold))])




def run():
    parser = argparse.ArgumentParser()
    # input files
    # todo outdir currently doesn't get used by anything
    parser.add_argument('--metadata')
    parser.add_argument('--vcfDir', 
        help='directory with all vcf files')
    parser.add_argument('--filterAllow',
        default=['PASS', 'min_af_0.010000'],
        nargs='+',
        help='which vcf filter strings to allow')
    args = parser.parse_args()
    # reads in metadata
    metadata = pd.read_csv(args.metadata, sep=',')
    transmission_pairs = get_transmission_pairs(metadata)
    metadata['date'] = pd.to_datetime(metadata['date'], format='%m/%d/%y', errors='coerce')

    vcf_dat = {}
    for path in glob.glob(args.vcfDir):
        sample, sample_dat = \
            process_vcf(path, args.filterAllow)
        vcf_dat[sample] = sample_dat
    
    all_pair_dat = pd.concat([
        get_pair_dat(vcf_dat[i[1][0]], vcf_dat[i[1][1]])\
            [['SAMPLE_x', 'SAMPLE_y', 'POS', 'REF', 'ALT', 'AF_x', 'AF_y', 'FILTER_PASS_x', 'FILTER_PASS_y']]
        for i in transmission_pairs.items()])

    all_a = []
    for t  in [0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01, 0.005, 0.0025]:
        all_a.append(all_pair_dat.query('REF != ALT & FILTER_PASS_x != False & FILTER_PASS_y != False').reset_index(drop=True).assign(
                    clonal = lambda k: ((k['AF_x'] >= 1-t) & (k['AF_y'] <= t)) | 
                        ((k['AF_x'] <= t) & (k['AF_y'] >= 1-t))).\
                groupby(['SAMPLE_x', 'SAMPLE_y']).\
                agg({'clonal': sum}).reset_index().assign(t=t))

    pd.concat(all_a).pivot(index=['SAMPLE_x', 'SAMPLE_y'], columns='t', values='clonal').\
        reset_index().\
        rename({'SAMPLE_x': 'ID1', 'SAMPLE_y':'ID2'}).to_csv('data/sars2_data_processed.csv')


if __name__ == "__main__":
    run()


