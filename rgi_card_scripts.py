import os
import time
import numpy as np

def exe_cmd(cmd):
    t1 = time.time()
    print("---------------------------------------------------")
    print(f"executing : \n{cmd}")
    os.system(cmd)
    t2 = time.time()
    print(f"\nTook {np.around(t2-t1,decimals=2)} seconds\n\n")
    print("---------------------------------------------------\n\n\n\n")

if __name__=="__main__":
    conf = {
        'threads':3,
        'aligner':"bowtie2"
    }
    trimmed_folder = 'kneaddata_trimmed_v2_test_25_75_10'
    card_folder = f'card_dir_test_{trimmed_folder}_{conf["aligner"]}'

    if card_folder not in os.listdir():
        os.mkdir(card_folder)

    files = {}
    for fo in os.listdir(trimmed_folder):
        for f in os.listdir(f"{trimmed_folder}/{fo}"):
            filepath = f"{trimmed_folder}/{fo}/{f}"
            if ".fastq.gz" not in f: continue
            if "L001" in f or "L002" in f: continue
            if "_kneaddata_paired_" not in f: continue

            id = f.split("_")[1]
            if id not in files:
                files[id] = []
                files[id].append(f)
            else:
                files[id].append(f)

    lg = open(f"logs/rgi_inscript_logger_{card_folder}.txt","w")
    for id in files.keys():
        lg.write(f"started for id: {id} \n")
        if id not in os.listdir(card_folder):
            os.mkdir(f"{card_folder}/{id}")
        os.chdir(f"{card_folder}/{id}")
        for f in files[id]:
            if "_paired_1" in f : f1 = f
            elif "_paired_2" in f: f2 = f
        lg.write(f"Before cmd id : {id}\n")
        cmd = f"rgi bwt --read_one ../../{trimmed_folder}/{id}/{f1} --read_two ../../{trimmed_folder}/{id}/{f2} --aligner {conf['aligner']} --output_file {'_'.join(f1.split('_')[:2])} --threads {conf['threads']} --clean --include_wildcard >> {'_'.join(f1.split('_')[:2])}_logfile.txt 2>&1"
        lg.write(f"After cmd id : {id} \n\n")
        # not sorted results
        exe_cmd(cmd)
        os.chdir(f"../../")
    lg.close()

