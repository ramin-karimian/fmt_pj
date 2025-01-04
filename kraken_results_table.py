import re
import pickle
from src.ncbi_funcs import *
from src.utils import *

def mapper(df):
    id2rank = {}
    name2taxid = {}
    taxid2name = {}
    feat_dit = {}
    for i in range(len(df)):
        id2rank[df['taxaID'][i]] = df['rank'][i]
        name2taxid[df['name'][i]] = df['taxaID'][i]
        taxid2name[df['taxaID'][i]] = df['name'][i]
        id = df['taxaID'][i]
        feat_dit[id] = {}
        feat_dit[id]['taxName'] = df['name'][i]
        feat_dit[id]['taxRank'] = df['rank'][i]
    return id2rank, name2taxid, taxid2name, feat_dit

def normalize_bysample(df):
    # df columns should be samples
    for s in df.columns:
        if "S" not in s: continue
        total_count = int(df[df['name']=='unclassified'][s]) + int(df[df['name']=='root'][s])
        df[s] = df[s].apply(lambda x: int(x))
        df[s] = df[s] / total_count
    return df

def remove_single_double_features(df):
    ## remove singletons and doubletons
    l = []
    df.drop(columns=['name','rank','taxaID'], inplace =True)
    df = df.T
    ## now taxa are columns
    for cl in df.columns:
        s = df[cl].apply(lambda x : 1 if x >0 else 0).sum()
        if s <=2:
            l.append(cl)
    df.drop(columns= l, inplace = True)
    return df, l

def remove_low_abundance_taxa(df,conf,th=0.001,th_sample=0.1):
    remove_below_half = conf['remove_below_half']
    for i in df.columns:
        if i in ['rank','name']: continue
        df[i] = df[i]/df[i][1]

    df.drop(index = [i for i in df.index if df['name'][i] in ["unclassified","root"]],columns=['rank','name'], inplace = True)
    rmvdict = {}
    below_half = {}

    sample_count = 0
    for sample in df.columns:
        if "_S" not in sample: continue
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
    return  df , set(rmvlist)

def separate_by_organism(df, conf, notfoundlist, name2id):
    rmvlist = []
    for i in df.index:
        id = df['taxaID'][i]
        name = df['name'][i]
        name2id[name] = id
        try:
            lineage = ncbi.get_lineage(df['taxaID'][i])
        except:
            lineage = []
            notfoundlist.append((df['name'][i], df['taxaID'][i]))

        if lineage == None: lineage = []
        if conf['organism']=="protozoa":
            state = 'y'
            for idd in conf['organismID'] :
                if idd in lineage:
                    state = 'n'
                    break
            if state=='y':  rmvlist.append(i)
        else:
            if conf['organismID'] not in lineage:
                rmvlist.append(i)
    df = df.drop(index=rmvlist)
    return df, rmvlist, notfoundlist, name2id

def seperate_by_taxRank(df, conf):
    notfoundlist = []
    name2id = {}

    df_total = {}
    for i, x in df.groupby('rank'):
        if i in ['U']: continue
        if re.search(f"\d",i): continue
        x, _, notfoundlist, name2id = separate_by_organism(x.copy(), conf, notfoundlist, name2id)
        x.drop(columns=['rank', 'name', 'taxaID'],inplace = True)
        x.columns = [cl.split("_")[1] for cl in x.columns]
        df_total[i] = x
        if True in x.index.duplicated():
            print("Duplicated name (WITHIN SAME RANK) due to the difference in taxaids but same taxa name")
    notfoundlist_df = pd.DataFrame(df,index=[j[1] for j in notfoundlist])
    df_total['ncbinotfound']=notfoundlist_df
    df.drop(columns=['taxaID'], inplace=True)
    return df_total, df, notfoundlist,notfoundlist_df, name2id


if __name__ == "__main__":
    org2id = {"bac": 2,   'vir': 10239 } ## in order to separate the organisms, IDS are forom NCBI
    conf = {
        'qs': 25,
        'minlen': 77,
        'headcrop': 10,
        'k-mer_length': 35,
        "organism": ["bac", 'vir'][-1],
        'normal':[False,True][1],
        "remove_below_half":[True,False][0], # True removes almost everything
        'normalization type':["normsamp"][0],
        'th_abundance':0.0001,
        'th_num_samples': 1,

    }
    conf["organismID"] = org2id[conf['organism']]

    with open("outputs/patients_dict.pkl", 'rb') as f:
        sampleid_info = pickle.load(f)



    db = f"krk{conf['k-mer_length']}_{conf['organism']}_{conf['normalization type']}-{conf['th_abundance']*100}-{conf['th_num_samples']*100}{'_belowhalf' if conf['remove_below_half'] else ''}"

    trimmed_folder = f'kneaddata_trimmed_v3_{conf["qs"]}_{conf["minlen"]}_{conf["headcrop"]}'

    kraken_folder = f'kraken_{trimmed_folder}_{db}'

    ## Constructing output folders
    if f"{trimmed_folder}_{db}" not in os.listdir("outputs"):
        os.mkdir(f"outputs/{trimmed_folder}_{db}")
        os.mkdir(f"outputs/{trimmed_folder}_{db}/figures")
        os.mkdir(f"outputs/{trimmed_folder}_{db}/figures/diversities")
        os.mkdir(f"outputs/{trimmed_folder}_{db}/figures/pca")
        os.mkdir(f"outputs/{trimmed_folder}_{db}/figures/pcoa")
        os.mkdir(f"outputs/{trimmed_folder}_{db}/figures/boxplots")
        os.mkdir(f"outputs/{trimmed_folder}_{db}/network")

    tablepath = f"outputs/{trimmed_folder}_{db}/table_{db}.xlsx"

    writer = pd.ExcelWriter(tablepath)

    table_taxon = ''
    table_clade = ''
    table_normalized = ''
    sample2count = {}

    sheetnames = ['#reads_clade', '#reads_taxon']
    for f in os.listdir(f'kraken_{trimmed_folder}_kraken_db_{conf["k-mer_length"]}'):
        filepath = f"kraken_{trimmed_folder}_kraken_db_{conf['k-mer_length']}/{f}"
        if ".report" not in f: continue
        id = f.split("_")[0] + "_" + f.split("_")[1]

        df = pd.read_csv(filepath, names=['percentage', '#reads_clade', '#reads_taxon', 'rank', 'taxaID', 'name'],
                         sep='\t', lineterminator='\n')
        df = df.sort_values(by='taxaID')
        df['name'] = df['name'].apply(lambda x: x.replace("\r", ''))
        df['name'] = df['name'].apply(lambda x: x.strip())
        df.index = range(len(df))
        if len(table_clade) == 0:
            table_taxon = pd.DataFrame(df, columns=['rank', 'taxaID', 'name'])
            table_clade = pd.DataFrame(df, columns=['rank', 'taxaID', 'name'])

        table_taxon[id] = df["#reads_taxon"].values
        table_clade[id] = df["#reads_clade"].values
        total_paired_reads = df[df['name'] == "root"]['#reads_clade'].values[0] + \
                             df[df['name'] == "unclassified"]['#reads_clade'].values[0]
        sample2count[id.split("_")[1]] = total_paired_reads

    table_taxon.to_excel(writer, sheet_name="taxon - not normal")
    table_clade.to_excel(writer, sheet_name="clade - not normal")

    ## Creating mapping dictionaries

    id2rank, name2taxid, taxid2name, feat_dit = mapper(table_clade.copy())
    table_clade.index = table_clade['taxaID'].values

    table_clade_notnormal = table_clade.copy()
    ##  normalize by sample
    if conf['normal']:
        table_clade = normalize_bysample(table_clade)
    else:
        for s in table_clade.columns:
            if "S" not in s: continue
            table_clade[s] = table_clade[s].apply(lambda x: int(x))
    ## remove singletons and doubletons
    _, l0 = remove_single_double_features(table_clade.copy())
    table_clade.drop(index=l0,inplace = True)


    ## Remove the low abundace taxa (lower than 0.1% for at least 10% of the samples)
    _ , rmvlist1 = remove_low_abundance_taxa(table_clade.copy(),conf,th=conf['th_abundance'] if conf['organism']!='fun' else conf['th_abundance']/100,th_sample=conf['th_num_samples']) # th is for feature abundance (percentage) in each sample, th_sample means it is for at least 10 % of amples
    table_clade.drop(index = rmvlist1,inplace = True)

    table_clade.to_excel(writer, sheet_name="Normalized")

    ## Separate taxa levels
    df_total_clade, df_clade, notfoundlist,notfoundlist_df, name2id = seperate_by_taxRank(table_clade.copy(), conf)


    for lvl in ['S','G','F','O','C','P','K']:
        df_total_clade[lvl].to_excel(writer,sheet_name=f"{conf['organism']}_{lvl}_Normalized")

    writer.save()
    writer.close()

    with open("outputs/name2id.pkl", 'wb') as f:
        pickle.dump(name2id, f)


    with open(f"{tablepath[:-5]}_{'normalized' if conf['normal'] else 'non-normal'}.pkl", "wb") as f:
        pickle.dump([df_total_clade, df_clade], f)
    with open(f"outputs/{trimmed_folder}_{db}/mapper_{kraken_folder}_{'normalized' if conf['normal'] else 'non-normal'}.pkl","wb") as f:
        pickle.dump([id2rank, name2taxid, taxid2name], f)
    with open(f"outputs/sample2count.pkl", "wb") as f:
        pickle.dump(sample2count, f)
    with open(f"outputs/{trimmed_folder}_{db}/feat_dict_{trimmed_folder}.pkl", 'wb') as f:
        pickle.dump(feat_dit, f)

