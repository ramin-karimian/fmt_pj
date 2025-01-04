import pickle
import pandas as pd
import os
import numpy as np


def normalize_bysample(df,sample2count):
    for cl in df.index:
        if cl in ["S35","S36","S37","S38"]: continue
        df.loc[cl,:] = df.loc[cl,:] / sample2count[cl]
    return df


def remove_single_double_features(df):
    ## remove singletons and doubletons
    l = []
    for cl in df.columns:
        s = df[cl].apply(lambda x : 1 if x >0 else 0).sum()
        if s <=2:
            l.append(cl)
    df.drop(columns= l, inplace = True)
    return df, l

def remove_low_abundance_taxa(df,conf,remove_below_half=True,th=0.001,th_sample=0.1):
    # samples should be columns
    df = df.T
    for i in df.columns:
        df[i] = df[i]/np.sum(df[i])

    rmvdict = {}
    below_half = {}

    sample_count = 0
    for sample in df.columns:
        if "S" not in sample: continue
        sample_count +=1
        for i in df.index:
            if df[sample][i] <= th:
                feature = i
                if feature not in rmvdict: rmvdict[feature]=0
                rmvdict[feature] +=1
            if df[sample][i] <= 0:
                if feature not in below_half: below_half[feature]=0
                below_half[feature] +=1

    rmvlist = [k for k,v in rmvdict.items() if v>=int(th_sample*sample_count)]
    below_half_list = [k for k,v in below_half.items() if v>=0.5*sample_count]
    if remove_below_half : rmvlist.extend(below_half_list)
    return  df.T , set(rmvlist)


def extract_genes_list(card_data_folder):
    genesList = []
    geneIdsList = []
    genesDict = {}
    geneFamilyList = []
    for dr in os.listdir(card_data_folder):
        if "S" not in dr: continue
        if dr in ["S35","S36","S37","S38"]: continue
        for f in os.listdir(f"{card_data_folder}/{dr}"):
            if ".gene" not in f: continue
            filepath = f"{card_data_folder}/{dr}/{f}"
            df = pd.read_csv(filepath,sep='\t')
            glist, idlist = df['ARO Term'].values, df['ARO Accession'].values
            genesList.extend(glist)
            geneIdsList.extend(idlist)
            geneFamilyList.extend(df['AMR Gene Family'].values)
            genesDict, genesFamilyDict = genes_dict(genesDict,glist, geneFamilyList,df)

    genesList , geneIdsList, geneFamilyList = list(set(genesList)), list(set(geneIdsList)), list(set(geneFamilyList))
    return genesList, geneIdsList, geneFamilyList, genesDict, genesFamilyDict


def genes_dict(genesDict,glist, genesFamily, df):
    genesFamilyDict = {}
    for g in glist:
        if g not in genesDict.keys():
            genesDict[g]= {}
            genesDict[g]["AMR Gene Family"] = df[df['ARO Term']==g]['AMR Gene Family'].values[0]
            genesDict[g][ 'Reference Model Type'] = df[df['ARO Term']==g][ 'Reference Model Type'].values[0]
            genesDict[g]['Reference DB'] = df[df['ARO Term']==g]['Reference DB'].values[0]
            genesDict[g]['Drug Class'] = df[df['ARO Term']==g]['Drug Class'].values[0]
            genesDict[g]['Resistance Mechanism'] = df[df['ARO Term']==g]['Resistance Mechanism'].values[0]
            genesDict[g]['ARO Accession'] = df[df['ARO Term']==g]['ARO Accession'].values[0]
            genesDict[g]['ID'] = df[df['ARO Term']==g]['ARO Accession'].values[0]
            genesDict[g]['Resistomes & Variants: Observed Pathogen(s)'] = df[df['ARO Term']==g]['Resistomes & Variants: Observed Pathogen(s)'].values[0]
        else:
            assert genesDict[g]["AMR Gene Family"] == df[df['ARO Term']==g]['AMR Gene Family'].values[0] and genesDict[g]['Drug Class'] == df[df['ARO Term']==g]['Drug Class'].values[0] and  genesDict[g]['Resistance Mechanism'] == df[df['ARO Term']==g]['Resistance Mechanism'].values[0] and  genesDict[g]['ARO Accession'] == df[df['ARO Term']==g]['ARO Accession'].values[0] , f"g = {g} error attributes for "
    for i in range(len(glist)):
        if genesFamily[i] not in genesFamilyDict:
            genesFamilyDict[genesFamily[i]] = {}
            genesFamilyDict[genesFamily[i]]['ID'] = len(genesFamilyDict)+1
    return genesDict , genesFamilyDict

def genefamily_data_convertor(df,conf):
    family_df = pd.DataFrame(columns=[conf['readstype']])
    for i,mf in df.groupby('AMR Gene Family'):
        count = mf[conf['readstype']].sum()
        family_df.loc[i]= [count]
    return family_df


def convert_allele_to_gene_data(df,conf):

    # collist = ['ARO Term',"ARO Accession",'Reference Model Type','Reference DB','All Mapped Reads','AMR Gene Family']
    collist = [conf['readstype'],'AMR Gene Family']

    df.drop(index = [x for x in df.index if df['Reference DB'][x]!="CARD"],inplace = True)
    new_df = pd.DataFrame(columns = collist  )
    for i,mf in df.groupby("ARO Term"):
        new_df.loc[i] = [mf[conf['readstype']].sum(),mf['AMR Gene Family'].values[0]]

    return new_df


def create_summed_table(table,genesDict,by):

    fea_list = list(set([y.strip()  for x in pd.DataFrame(genesDict).T[by].values for y in x.split(";")]))

    sumtable = pd.DataFrame(index = table.index , columns = fea_list )
    sumtable = sumtable.fillna(0)

    for col in table.columns:
        ds = [x.strip() for x in genesDict[col][by].split(";")]
        for d in ds:
            sumtable[d] += table[col]

    return sumtable, fea_list



if __name__=="__main__":

    conf = {
        'qs':25,
        'minlen':77,
        'headcrop':10,
        'method':['bowtie2','kma'][0],
        "db" : ["CARD","Resistomes & Variants"][0],
        'normal':[False,True][1],
        'genefamily': ["gene","gnflor",'gndrg'][0], # do not use gnefml anymore,
        "remove_below_half":[True,False][1], # True removes almost everything
        'normalization type':["normsamp","normsampfea"][0],
        'reads':['all','complete'][0],
        'th_abundance': 0.001, # 0.001 below 0.1% relative abundance
        'th_num_samples': 1, # 0.1 at least 10% of the samples
        }

    conf['readstype'] = "All Mapped Reads" if conf['reads']=='all' else "Completely Mapped Reads"

    trimmed_folder = f'kneaddata_trimmed_v3_{conf["qs"]}_{conf["minlen"]}_{conf["headcrop"]}'

    card_data_folder = f'card_dir_{conf["method"]}_{trimmed_folder}'
    conf['method'] = f"{conf['method']}-{conf['normalization type']}-{conf['reads']}-{conf['th_abundance']*100}-{conf['th_num_samples']*100}{'_belowhalf' if conf['remove_below_half'] else ''}"
    output_folder = f"{trimmed_folder}_rgi-{conf['method']}-{conf['genefamily']}"
    print(output_folder)
    if f"{output_folder}" not in os.listdir("outputs"):
        os.mkdir(f"outputs/{output_folder}")
        os.mkdir(f"outputs/{output_folder}/figures")
        os.mkdir(f"outputs/{output_folder}/figures/diversities")
        os.mkdir(f"outputs/{output_folder}/figures/pca")
        os.mkdir(f"outputs/{output_folder}/figures/pcoa")
        os.mkdir(f"outputs/{output_folder}/figures/boxplots")
        os.mkdir(f"outputs/{output_folder}/network")


    with open(f"outputs/sample2count.pkl","rb") as f:
        sample2count = pickle.load(f)

    genesList, geneIdsList, geneFamilyList, genesDict, genesFamilyDict = extract_genes_list(card_data_folder)

    table = pd.DataFrame(columns= genesList)
    family_table = pd.DataFrame(columns=geneFamilyList)
    inds = []
    for dr in os.listdir(card_data_folder):
        if "S" not in dr: continue
        for f in os.listdir(f"{card_data_folder}/{dr}"):
            if ".allele" not in f or "json"  in f: continue
            filepath = f"{card_data_folder}/{dr}/{f}"
            id = f.split("_")[1].split(".")[0]

            df_allele = pd.read_csv(filepath,sep='\t')
            df = convert_allele_to_gene_data(df_allele.copy(),conf)   # this should ve changed
            family_df = genefamily_data_convertor(df,conf)
            table = table.append(df[conf['readstype']], ignore_index=True)
            family_table = family_table.append(family_df[conf['readstype']], ignore_index=True)
            inds.append(id)


    table.index = inds
    table = table.fillna(0)
    family_table.index = inds
    family_table = family_table.fillna(0)

    ## create drug_class table
    drugtable , _ = create_summed_table(table,genesDict,by='Drug Class')
    family_table_summed , _ = create_summed_table(table,genesDict,by='AMR Gene Family')


    if conf['genefamily']=='gene': finaltable = table.copy()
    elif conf['genefamily']=='gnflor': finaltable = family_table_summed.copy()
    elif conf['genefamily']=='gndrg': finaltable = drugtable.copy()

    ## remove singletons (already removed) and doubletons
    _ , l = remove_single_double_features(finaltable.copy())
    finaltable.drop(columns=l,inplace = True)
    ## removing low abundance markers or families

    tmp_table , rmvlist = remove_low_abundance_taxa(finaltable.copy(),conf.copy(),remove_below_half=conf['remove_below_half'],th=conf['th_abundance'],th_sample=conf['th_num_samples'])
    finaltable.drop(columns=rmvlist,inplace = True)


    if conf['normal']:
        finaltable = normalize_bysample(finaltable.copy(),sample2count)

    with open(f"outputs/{output_folder}/{conf['genefamily']}_table_rgi-{conf['method']}-{conf['genefamily']}.pkl",'wb') as f:
        pickle.dump([finaltable,None],f)

    with open(f"outputs/{output_folder}/genesDict.pkl",'wb') as f:
        pickle.dump([genesDict,genesFamilyDict],f)
