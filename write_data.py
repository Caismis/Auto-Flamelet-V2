from sympy import cse
from funclib import stdwrt
import numpy as np
from scipy.interpolate import griddata

def write_data(input_data, uniset, misc):
    print('---Write data to file---')
    dissrates = np.load(file='dissrates.npy', allow_pickle=True)
    outname = input_data[9]
    nzeta = input_data[13]
    npoint = input_data[12]
    zetas = np.linspace(0, 0.99, nzeta, endpoint=True)
    optf = open(outname, 'w')
    optf.write('Number of Flamelets' + '\n')
    optf.write(str(len(uniset)) + '\n')
    optf.write('Number of Zetas' + '\n')
    optf.write(str(nzeta) + '\n')
    optf.write('Number of Points\n')
    optf.write(str(npoint) + '\n')
    optf.write('Number of Species\n')
    optf.write(str(len(misc[1])) + '\n')
    for i in range(len(uniset)):
        optf.write('Dissipation rate\n')
        optf.write('{:<.6e}\n'.format(dissrates[i]))
        for j in range(nzeta):
            optf.write('Zeta\n')
            optf.write('{:<7.2f}\n'.format(zetas[j]))
            for k in range(len(uniset[i][j])):
                if k == 0:
                    for isp in range(misc[0]):
                        optf.write(misc[1][isp] + '\n')
                        stdwrt(uniset[i][j][k][isp], optf)
                elif k == 1:
                    optf.write('Temperature\n')
                    stdwrt(uniset[i][j][k], optf)
                elif k == 2:
                    optf.write('Density\n')
                    stdwrt(uniset[i][j][k], optf)
                elif k == 3:
                    optf.write('Z\n')
                    stdwrt(uniset[i][j][k], optf)
    print('---Completed---')
    optf.close()


def write_remap(input_data, uniset, misc):
    print('---Write data to file---')
    cseries = np.linspace(0, misc[2], len(uniset))
    optf = open('clist.d', 'w')
    optf.write('Number of C' + '\n')
    optf.write(str(len(uniset)) + '\n')
    optf.write('C value' + '\n')
    stdwrt(cseries, optf)
    optf.close()
    outname = input_data[9]
    nzeta = input_data[13]
    npoint = input_data[12]
    zetas = np.linspace(0, 1, nzeta, endpoint=False)
    normzs = np.linspace(0, 1, npoint)
    mapset = remapper(uniset, misc[0], npoint, nzeta)
    optf = open(outname, 'w')
    optf.write('Max Value of C' + '\n')
    optf.write('{:<.6e}\n'.format(misc[2]))
    optf.write('Number of Z' + '\n')
    optf.write(str(npoint) + '\n')
    optf.write('Number of Zetas' + '\n')
    optf.write(str(nzeta) + '\n')
    optf.write('Number of Points\n')
    optf.write(str(len(uniset)) + '\n')
    optf.write('Number of Species\n')
    optf.write(str(len(misc[1])) + '\n')
    for i in range(npoint):
        optf.write('Z\n')
        optf.write('{:<.6e}\n'.format(normzs[i]))
        for j in range(nzeta):
            optf.write('Zeta\n')
            optf.write('{:<7.2f}\n'.format(zetas[j]))
            currentcs = mapset[0][i][j]
            for k in range(len(mapset)):
                if k == 1:
                    for isp in range(misc[0]):
                        optf.write(misc[1][isp] + '\n')
                        # stdwrt(mapset[k][i][j][isp], optf)
                        stdwrt(cinterpol(currentcs, mapset[k][i][j][isp], cseries), optf)
                elif k == 2:
                    optf.write('Temperature\n')
                    # stdwrt(mapset[k][i][j], optf)
                    stdwrt(cinterpol(currentcs, mapset[k][i][j], cseries), optf)
                elif k == 0:
                    optf.write('C\n')
                    # stdwrt(mapset[k][i][j], optf)
                    stdwrt(cseries, optf)
                elif k == 3:
                    optf.write('Rate\n')
                    # stdwrt(mapset[k][i][j], optf)
                    stdwrt(cinterpol(currentcs, mapset[k][i][j], cseries), optf)
    print('---Completed---')
    optf.close()


def remapper(uniset, nsp, npoint, nzeta):
    print('---Remapping Data---')
    lambdas = len(uniset)
    normzs = np.linspace(0, 1, npoint)
    remapc = []
    #Progress Variable
    remapt = []
    #Temperature
    remapr = []
    #Progress Rate
    remaps = []
    #Species Fraction
    for i in range(len(normzs)):
        listt = []
        listc = []
        listr = []
        lists = []
        for j in range(nzeta):
            slicet =  np.empty(lambdas)
            slicec =  np.empty(lambdas)
            slicer =  np.empty(lambdas)
            slices =  []
            for isp in range(nsp):
                slices.append(np.empty(lambdas))
            for k in range(lambdas):
                slicec[k] = uniset[k][j][2][i]
            cindex = np.argsort(slicec)
            slicec = np.sort(slicec)
            for k in range(lambdas):
                slicet[k] = uniset[cindex[k]][j][1][i]
                slicer[k] = uniset[cindex[k]][j][4][i]
                for isp in range(nsp):
                    slices[isp][k] = uniset[cindex[k]][j][0][isp][i]
            listt.append(slicet)
            listc.append(slicec)
            listr.append(slicer)
            lists.append(slices)
        remapt.append(listt)
        remapc.append(listc)
        remapr.append(listr)
        remaps.append(lists)
    return [remapc, remaps, remapt, remapr]

def cinterpol(clist, nlist, stdc):
    outlist = []
    for c in stdc:
        if c <= clist[0]:
            outlist.append(nlist[0])
        elif c >= clist[-1]:
            outlist.append(nlist[-1])
        else:
            for i in range(len(clist)):
                if c >= clist[i] and c < clist[i + 1]:
                    interout = ((nlist[i + 1] - nlist[i])/(clist[i + 1] - clist[i]))*(c - clist[i + 1]) + nlist[i + 1]
                    outlist.append(interout)
                    break
    return outlist