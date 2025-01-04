import os

def exe_cmd(cmd):
    print(f"executing : \n{cmd}")
    os.system(cmd)
    print("\n\n")

if __name__=="__main__":
    print('\n\n\nStep 1\n')

    conf = {
        "qs":25,
        'threads':40,
        'minlen':77,
        'headcrop':10
    }

    db = "kraken_db" ## path to db



    trimmed_folder = f'kneaddata_trimmed_v3_{conf["qs"]}_{conf["minlen"]}_{conf["headcrop"]}'
    kraken_folder = f'kraken_{trimmed_folder}_{db}'

    if kraken_folder not in os.listdir():
        os.mkdir(kraken_folder)

    files = {}
    for fo in os.listdir(trimmed_folder):
        for f in os.listdir(f"{trimmed_folder}/{fo}"):
            filepath = f"{trimmed_folder}/{fo}/{f}"
            if ".fastq.gz" not in f: continue
            if "L001" in f or "L002" in f: continue
            if "_kneaddata_paired_" not in f: continue
            print(f"{filepath}\n")
            id = f.split("_")[1]
            if id not in files:
                files[id] = []
                files[id].append(f)
            else:
                files[id].append(f)

    for id in files.keys():
        for f in files[id]:
            if "_paired_1" in f : f1 = f
            elif "_paired_2" in f: f2 = f
        cmd = f"kraken2 --use-names --threads {conf['threads']} --db {db} --report-zero-count --report " \
              f"{kraken_folder}/{f1.replace('.fastq.gz','').replace('_R1','').replace('_paired_1','')}.report --gzip-compressed --paired " \
              f"{trimmed_folder}/{id}/{f1} {trimmed_folder}/{id}/{f2} > {kraken_folder}/{f1.replace('.fastq.gz','').replace('_R1','').replace('_paired_1','')}.kraken"

        exe_cmd(cmd)
