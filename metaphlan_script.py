import os

def exe_cmd(cmd):
    print(f"executing : \n{cmd}")
    os.system(cmd)
    print("\n\n")

if __name__=="__main__":
    print('\n\n\nStep 1\n')
    conf = {
        'threads':5,
    }
    data_path = f""
    read1 = f""
    read2 = f""

    db = f"k2_pluspf_20210517_UPDATE"
    trimmed_folder = 'kneaddata_trimmed_v1'
    trimmed_fastqc_folder = 'kneaddata_trimmed_v1_fastqc'
    metaphlan_folder = f'metaphlan_{trimmed_folder}_{db}'

    if metaphlan_folder not in os.listdir():
        os.mkdir(metaphlan_folder)

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
        cmd = f"metaphlan -t rel_ab_w_read_stats --bowtie2db ../../data/metaphlan_databases/ {trimmed_folder}/{id}/{f1},{trimmed_folder}/{id}/{f2}  --bowtie2out {metaphlan_folder}/{f1.replace('.fastq.gz','').replace('_R1','').replace('_paired_1','')}.bowtie2.bz2 --nproc {conf['threads']}  --input_type fastq --add_viruses --sample_id_key {id} -o {metaphlan_folder}/{f1.replace('.fastq.gz','').replace('_R1','').replace('_paired_1','')}.txt"
        exe_cmd(cmd)

    cmd2 = f"merge_metaphlan_tables.py {metaphlan_folder}/*_kneaddata.txt > outputs/metaphlan_merged_abundance_table_{metaphlan_folder}.txt"
    exe_cmd(cmd2)



