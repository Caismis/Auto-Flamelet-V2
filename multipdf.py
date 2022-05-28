from funclib import *
from multiprocessing import Pool, cpu_count
import numpy as np
import platform
from tqdm import tqdm

def multipdf(input_data, switch):
    dataseries = np.load(file='dataseries.npy', allow_pickle=True)
    nzeta = input_data[13]
    if switch is False:
        nsp = len(dataseries[0]) - 3
        pbar = tqdm(total=len(dataseries), bar_format='{l_bar}>>>{bar}<<<', colour='RED')
    else:
        nsp = len(dataseries[0]) - 5
        pbar = tqdm(total=len(dataseries), bar_format='{l_bar}>>>{bar}<<<', colour='BLUE')
    pbar.set_description('Progress')
    update = lambda *args: pbar.update()
    rslt = []
    if platform.system() == 'Windows':
        pool = Pool(8)
    else:
        pool = Pool(cpu_count())
    bundle = [(dataseries[i], nsp, i, nzeta) for i in range(len(dataseries))]
    for package in bundle:
        r = pool.apply_async(presume, (package, switch, ), callback=update)
        # r = pool.apply_async(lim_presume, (package, switch, ), callback=update)
        rslt.append(r)
    for r in rslt:
        r.wait()
    pool.close()
    pool.join()
    uniset= [r.get() for r in rslt]
    if switch:
        return uniset[::-1]
    else:
        return uniset
    

def presume(package, switch):
    datas = package[0]
    nsp = package[1]
    nzeta = package[3]
    metadata = []
    zetas = np.linspace(0, 0.99, nzeta, endpoint=True)
    for z in zetas:
        intmassset = []
        inttemps = pdfinter(z, datas[nsp], datas[nsp + 2])      
        if switch:
            intrhos = pdfinter(z, datas[nsp + 4], datas[nsp + 2])
            intvar = pdfinter(z, datas[nsp + 3], datas[nsp + 2])
        else:
            intrhos = pdfinter(z, datas[nsp + 1], datas[nsp + 2])
        for isp in range(nsp):
            intmass = pdfinter(z, datas[isp], datas[nsp + 2])
            intmassset.append(intmass)
        if switch:
            metadata.append([intmassset, inttemps, intrhos, datas[nsp + 2], intvar])
        else:
            metadata.append([intmassset, inttemps, intrhos, datas[nsp + 2]])
    return metadata

