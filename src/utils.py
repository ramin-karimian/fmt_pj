import os

import pandas as pd
import math


def split(df,sampleid_info,patient):
    x = pd.DataFrame(sampleid_info)
    x = x.T
    ind = x[x['patient']==patient].index
    xf = pd.DataFrame(df,columns = ind)
    return xf

def split_v2(df,sampleid_info,patient):
    x = pd.DataFrame(sampleid_info)
    x = x.T
    ind = x[x['patient']!=patient].index
    xf = pd.DataFrame(df,columns = ind)

    ind = x[x['patient']=='EDM'].index
    xf = xf.drop(columns = ind)
    return xf

def print_progess(i,ln,lim):
    if i%lim == 0:
        print(f"{i}/{ln}\n")



def millify(n):
    millnames = ['','K','M','B','T']
    n = float(n)
    millidx = max(0,min(len(millnames)-1,
                        int(math.floor(0 if n == 0 else math.log10(abs(n))/3))))

    return '{:.0f}{}'.format(n / 10**(3 * millidx), millnames[millidx])

def make_dir(path,dir):
    if dir not in os.listdir(path):
        os.mkdir(f"{path}/{dir}")

def make_directory(folder,dir):
    if dir not in os.listdir(folder):
        os.mkdir(f"{folder}/{dir}")
