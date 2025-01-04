import pandas as pd
import pickle
import os



def taxa_mapper(df):
    taxid2rank = {}
    taxid2name = {}
    taxid2clade_name = {}
    name2taxid = {}
    name2clade_name = {}
    name2clade_taxid = {}
    feat_dit = {}
    clade2rank = {}
    cladename2taxid = {}

    for i in range(len(df)):
        taxid2rank[df['taxid'][i]] = df['taxlevel'][i]
        taxid2name[df['taxid'][i]] = df['taxname'][i]
        taxid2clade_name[df['taxid'][i]] = df['clade_name'][i]
        name2taxid[df['taxname'][i]] = df['taxid'][i]
        name2clade_name[df['taxname'][i]] = df['clade_name'][i]
        name2clade_taxid[df['taxname'][i]] = df['clade_taxid'][i]
        clade2rank[df['clade_name'][i]] = df['taxlevel'][i]
        cladename2taxid[df['clade_name'][i]] = df['taxid'][i]
        feat_dit[id] = {}
        feat_dit[id]['taxName'] = df['taxname'][i]
        feat_dit[id]['taxRank'] = df['taxlevel'][i]
    return taxid2rank , name2taxid, name2clade_name, name2clade_taxid, clade2rank, feat_dit, cladename2taxid, taxid2name, taxid2clade_name

def seperate_by_taxRank(df , clade2rank ):
    df['taxlevel'] = [clade2rank[x] for x in df.index ]
    df_total = {}
    for i,x in df.groupby('taxlevel'):
        x.drop(columns = ['taxlevel'],inplace=True)
        df_total[i] =  x
    return df_total , df

def separate_by_organism(df,conf):
    rmvlist = []
    for i in range(len(df)):
        clade = df.index[i]
        kingdom = clade.split("|")[0]
        if conf['mapper'][kingdom] != conf['organism']:
            rmvlist.append(df.index[i])
    df.drop(index = rmvlist,inplace = True)
    return df, rmvlist

def remove_single_double_features(df, double = True):
    ## remove singletons and doubletons
    if double: num = 2
    else: num =1
    l = []
    df = df.T
    ## now taxa are columns
    for cl in df.columns:
        s = df[cl].apply(lambda x : 1 if x >0 else 0).sum()
        if s <=num:
            l.append(cl)
    df.drop(columns= l, inplace = True)
    return df, l

def remove_low_abundance_taxa(df,total,conf,remove_below_half=True,th=0.1,th_sample=1):
    # samples as columns
    rmvdict = {}
    below_half = {}
    remove_below_half = conf['remove_below_half']
    df = df/total

    ### remove singletone and soubletone taxa
    _ , l = remove_single_double_features(df.copy(),double = conf['double'])
    df.drop(index= l, inplace = True)

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
    rmvlist.extend(l)
    below_half_list = [k for k,v in below_half.items() if v>=0.5*sample_count]
    if remove_below_half : rmvlist.extend(below_half_list)
    return  df , set(rmvlist)


if __name__ == "__main__":
    org2id = {"bac": "k__Bacteria",  'vir':"k__Viruses"  }

    conf = {
            'level' : ['S','G','F','O','C','P','D'][0],
            'normal':'normalized',
            'qs':25,
            'minlen':77,
            'headcrop':10,
            'db':"metaphlanDB",
            "organism": ["bac", 'vir'][0],
            "mapper":{"k__Bacteria":'bac',"k__Archaea":'archaea',"k__Viruses":'vir',"UNKNOWN":"UNKNOWN"},
            "remove_below_half":[True,False][0],
            'normalization type':["normsamp"][0],
            'th_abundance':0.0001, # 0.001 below 0.1% relative abundance
            'th_num_samples': 1, # 0.1 at least 10% of the samples
            'double':True
            }
    conf["organismID"] = org2id[conf['organism']]

    conf['db'] = f"metaphlanDB_{conf['organism']}_{conf['normalization type']}-{conf['th_abundance']*100}-{conf['th_num_samples']*100}{'_belowhalf' if conf['remove_below_half'] else ''}"

    conf['trimmed_folder'] = f'kneaddata_trimmed_v3_{conf["qs"]}_{conf["minlen"]}_{conf["headcrop"]}'

    metaphlan_folder = f'metaphlan_{conf["trimmed_folder"]}'
    datafolder = f"{conf['trimmed_folder']}_{conf['db']}"

    if f"{datafolder}" not in os.listdir("outputs"):
        os.mkdir(f"outputs/{datafolder}")
        os.mkdir(f"outputs/{datafolder}/figures")
        os.mkdir(f"outputs/{datafolder}/figures/diversities")
        os.mkdir(f"outputs/{datafolder}/figures/pca")
        os.mkdir(f"outputs/{datafolder}/figures/pcoa")
        os.mkdir(f"outputs/{datafolder}/figures/boxplots")
        os.mkdir(f"outputs/{datafolder}/network")


    outputpath = f"outputs/{datafolder}/all_tables.xlsx"
    writer = pd.ExcelWriter(outputpath)
    lns = open(f"outputs/metaphlan_merged_abundance_table_{metaphlan_folder}.txt").readlines()
    df = pd.DataFrame(columns=lns[1].split("\t"))
    for i in range(2,len(lns)):
        lns[i] = lns[i].replace("\n",'')
        df.loc[i-2] = lns[i].split("\t")

    mapper = {1:"K",2:"P",3:"C",4:"O",5:"F",6:"G",7:'S'}
    df['taxname'] = df['clade_name'].apply(lambda x: x.split("|")[-1])
    df['taxid'] = df['clade_taxid'].apply(lambda x: x.split("|")[-1])
    df['taxlevel'] = df['clade_name'].apply(lambda x: mapper[len(x.split("|"))])

    ## rearrange the order of the columns
    cols=[]
    for cl in df.columns:
        if "_S" in cl:
            cols.append(cl.split("_")[1])
        else:
            cols.append(cl)
    df.columns = cols
    df.to_excel(writer,sheet_name="data")


    taxid2rank , name2taxid, name2clade_name, name2clade_taxid, clade2rank, feat_dit, cladename2taxid, taxid2name, taxid2clade_name = taxa_mapper(df.copy())

    df.index = df['clade_name'].values
    df.drop(columns=['clade_name','clade_taxid','taxname','taxid','taxlevel'],inplace = True)
    df = df.astype("float")
    total =  -1*(df.loc['UNKNOWN']- 100)
    df.drop( index = ["UNKNOWN"], inplace = True)

    df.to_excel(writer,sheet_name="not_normalized")

    ## remove low abundance taxa and S/D (Sample normalization is done by metaphlan naturally since metaphlan gives percentage values it is only about rescaling the values and nothing more
    _ ,rmvlist1 = remove_low_abundance_taxa(df.copy(),total,conf,remove_below_half=conf['remove_below_half'],th=conf['th_abundance'],th_sample=conf['th_num_samples'])
    for i in df.columns:
        df[i] = df[i].apply(lambda x: float(x))

    df.drop(index = rmvlist1,inplace = True)

    df.to_excel(writer,sheet_name="Normalized")


    ## Separate organism of interest
    df2 , rmvlist2 = separate_by_organism(df.copy(),conf)

    df.to_excel(writer,sheet_name=f"{conf['organism']}_alllevels_Normalized")

    ## Separate different taxa levels (ranks)
    df_total, _ = seperate_by_taxRank(df2.copy(), clade2rank)

    for lvl in ['S','G','F','O','C','P','K']:
        df_total[lvl].to_excel(writer,sheet_name=f"{conf['organism']}_{lvl}_Normalized")

    with open(f"outputs/{datafolder}/table_{conf['db']}_normalized.pkl",'wb') as f:
        pickle.dump([df_total,df2],f)

    with open(f"outputs/{datafolder}/mapper_{metaphlan_folder}_normalized.pkl",'wb') as f:
        pickle.dump([taxid2rank , name2taxid, name2clade_name, name2clade_taxid, cladename2taxid,taxid2name, taxid2clade_name],f)

    with open(f"outputs/{datafolder}/feat_dict_{conf['trimmed_folder']}.pkl",'wb') as f:
        pickle.dump(feat_dit ,f)

    writer.save()
    writer.close()
