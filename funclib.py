import math
import numpy as np
import cantera as ct
import scipy.integrate as integrate

def stdwrt(sourcelist, wrtfile):
    lines = len(sourcelist)//5
    remains = len(sourcelist)%5
    output = ''
    for line in range(lines):
        thisline = ''
        for item in sourcelist[line*5:(line + 1)*5]:
            if(abs(item) <= 10e-12):
                item = 0
            thisline += np.format_float_scientific(item, precision=9, exp_digits=2, unique=False) + ' '
        output += thisline + '\n'
    if remains != 0:
        for item in sourcelist[(line + 1)*5:]:
            if(abs(item) <= 10e-12):
                item = 0
            output += np.format_float_scientific(item, precision=9, exp_digits=2, unique=False) + ' '
        output += '\n'
    wrtfile.write(output)



def pdfinter(zeta, sourcelist, zs):
    varzs = np.zeros(len(sourcelist))
    pdfs = np.zeros(len(sourcelist))
    palphas = np.zeros(len(sourcelist))
    pbetas = np.zeros(len(sourcelist))
    points = 250
    hzs = np.zeros(points)
    helpzs = np.zeros(points)
    hzs = np.zeros(points)
    deltas = np.zeros(points)
    datalist = sourcelist.copy()
    varzs = zeta*(zs*(1.0 - zs))
    for i in range(len(zs)):
        # varzs[i] = (zeta)*(zs[i]*(1 - zs[i]))
        if varzs[i] >= 1e-7:
            palphas[i] = zs[i]*((zs[i]*(1.0 - zs[i]))/varzs[i] - 1)
            pbetas[i] = (1 - zs[i])*((zs[i]*(1.0 - zs[i]))/varzs[i] - 1)
            if palphas[i] > 500:
                pbetas[i] = pbetas[i]/palphas[i] * 500
                palphas[i] = 500
            if pbetas[i] > 500:
                palphas[i] = palphas[i]/pbetas[i] * 500
                pbetas[i] = 500
            
            if (palphas[i] > 1) and (pbetas[i] > 1):
                zmax = 0
                n1 = 0
                n2 = 0
                for j in range(len(zs)):
                    try:
                        pdfs[j] = (zs[j]**(palphas[i] - 1))*((1.0 - zs[j])**(pbetas[i] - 1))*min(math.gamma(min(palphas[i] + pbetas[i], 20)), 1e17)/(min(math.gamma(min(20, palphas[i])), 1e17)* \
                            min(math.gamma(min(20, pbetas[i])), 1e17))
                    except:
                        print(palphas[i], pbetas[i], zs[j])
                        raise OverflowError
                    if pdfs[j] > pdfs[max(j - 1, 0)]:
                        zmax = zs[j]
                if palphas[i]/pbetas[i] <= 0.5:
                    n1 = int(0.2*points)
                    n2 = int(0.8*points) + 1 
                elif palphas[i]/pbetas[i] >= 2:
                    n2 = int(0.2*points) + 1
                    n1 = int(0.8*points)
                else:
                    n1 = int(0.5*points)
                    n2 = int(0.5*points) + 1 
                ex1 = 0.9
                deltas[0] = (1.0 - ex1)/(1.0 - ex1**(n1 - 1))
                for j in range(n1):
                    deltas[j] = (ex1**j)*deltas[0]
                    hzs[j + 1] = hzs[j] + deltas[j]
                for j in range(1, n1):
                    hzs[j] *= zmax
                ex2 = 1.1
                deltas[0] = (1.0 - ex2)/(1.0 - ex2**(n2 - 1))
                for j in range(n2 - 1):
                    deltas[j] = (ex2**j)*deltas[0]
                    helpzs[j + 1] = helpzs[j] + deltas[j]
                for j in range(n2):
                    helpzs[j] *= (1.0 - zmax)
                for j in range(n2 - 1):
                    hzs[n1 + j] = hzs[n1 - 1] + helpzs[j]
                for j in range(points):
                    hzs[j] /= hzs[points - 1]
                pdfs = np.zeros(len(hzs))
                for j in range(len(hzs)):
                    pdfs[j] = (hzs[j]**(palphas[i] - 1))*((1.0 - hzs[j])**(pbetas[i] - 1))*min(math.gamma(min(20, palphas[i] + pbetas[i])), 1e17)/(min(math.gamma(min(20, palphas[i])), 1e17)* \
                        min(math.gamma(min(20, pbetas[i])), 1e17))
            elif (palphas[i] <= 1) and (pbetas[i] > 1):
                n1 = 0
                n2 = 0
                if palphas[i]/pbetas[i] > 0.5:
                    zmax = 0.5
                    ex1 = 1.1
                    n1 = int(0.7*points)
                    deltas[0] = (1.0 - ex1)/(1.0 - ex1**(n1 - 1))
                    for j in range(n1):
                        deltas[j] = (ex1**j)*deltas[0]
                        hzs[j + 1] = hzs[j] + deltas[j]
                    for j in range(1, n1):
                        hzs[j] *= zmax
                    ex2 = 1.1
                    n2 = int(0.3*points) + 1
                    deltas[0] = (1.0 - ex2)/(1.0 - ex2**(n2 - 1))
                    for j in range(n2 - 1):
                        deltas[j] = (ex2**j)*deltas[0]
                        helpzs[j + 1] = helpzs[j] + deltas[j]
                    for j in range(n2):
                        helpzs[j] *= (1.0 - zmax)
                    for j in range(n2 - 1):
                        hzs[n1 + j] = hzs[n1 - 1] + helpzs[j]
                else:
                    ex2 = 1.05
                    deltas[0] = (1.0 - ex2)/(1.0 - ex2**(points - 1))
                    for j in range(points - 1):
                        deltas[j] = (ex2**j)*deltas[0]
                        hzs[j + 1] = hzs[j] + deltas[j]

                for j in range(points):
                    hzs[j] /= hzs[points - 1]
                pdfs = np.zeros(len(hzs))
                for j in range(1, len(hzs)):
                    pdfs[j] = (hzs[j]**(palphas[i] - 1))*((1.0 - hzs[j])**(pbetas[i] - 1))*min(math.gamma(min(20, palphas[i] + pbetas[i])), 1e17)/(min(math.gamma(min(20, palphas[i])), 1e17)* \
                        min(math.gamma(min(20, pbetas[i])), 1e17))
                pdfs[0] = 1.5*pdfs[1]/palphas[i]
            elif (palphas[i] > 1) and (pbetas[i] <= 1):
                n1 = 0
                n2 = 0
                if palphas[i]/pbetas[i] < 2:
                    zmax = 0.5
                    ex1 = 1.1
                    n1 = int(0.3*points)
                    deltas[0] = (1.0 - ex1)/(1.0 - ex1**(n1 - 1))
                    for j in range(n1):
                        deltas[j] = (ex1**j)*deltas[0]
                        hzs[j + 1] = hzs[j] + deltas[j]
                    for j in range(1, n1):
                        hzs[j] *= zmax
                    ex2 = 0.9
                    n2 = int(0.7*points) + 1
                    deltas[0] = (1.0 - ex2)/(1.0 - ex2**(n2 - 1))
                    for j in range(n2 - 1):
                        deltas[j] = (ex2**j)*deltas[0]
                        helpzs[j + 1] = helpzs[j] + deltas[j]
                    for j in range(n2):
                        helpzs[j] *= (1.0 - zmax)
                    for j in range(n2 - 1):
                        hzs[n1 + j] = hzs[n1 - 1] + helpzs[j]
                else:
                    ex1 = 0.95
                    deltas[0] = (1.0 - ex1)/(1.0 - ex1**(points - 1))
                    for j in range(points - 1):
                        deltas[j] = (ex1**j)*deltas[0]
                        hzs[j + 1] = hzs[j] + deltas[j]
                
                for j in range(points):
                    hzs[j] /= hzs[points - 1]
                pdfs = np.zeros(len(hzs))
                for j in range(len(hzs) - 1):
                    pdfs[j] = (hzs[j]**(palphas[i] - 1))*((1.0 - hzs[j])**(pbetas[i] - 1))*min(math.gamma(min(20, palphas[i] + pbetas[i])), 1e17)/(min(math.gamma(min(20, palphas[i])), 1e17)* \
                        min(math.gamma(min(20, pbetas[i])), 1e17))
                pdfs[points - 1] = 1.5*pdfs[points - 2]/pbetas[i]
            elif (palphas[i] <= 1) and (pbetas[i] <= 1):
                n1 = 0
                n2 = 0
                zmax = 0.5
                ex1 = 1.1
                n1 = int(0.5*points)
                deltas[0] = (1.0 - ex1)/(1.0 - ex1**(n1 - 1))
                for j in range(n1):
                    deltas[j] = (ex1**j)*deltas[0]
                    hzs[j + 1] = hzs[j] + deltas[j]
                for j in range(1, n1):
                    hzs[j] *= zmax
                ex2 = 0.9
                n2 = int(0.5*points) + 1
                deltas[0] = (1.0 - ex2)/(1.0 - ex2**(n2 - 1))
                for j in range(n2 - 1):
                    deltas[j] = (ex2**j)*deltas[0]
                    helpzs[j + 1] = helpzs[j] + deltas[j]
                for j in range(n2):
                    helpzs[j] *= (1.0 - zmax)
                for j in range(n2 - 1):
                    hzs[n1 + j] = hzs[n1 - 1] + helpzs[j]
                for j in range(points):
                    hzs[j] /= hzs[points - 1]
                pdfs = np.zeros(len(hzs))
                for j in range(1, len(hzs) - 1):
                    pdfs[j] = (hzs[j]**(palphas[i] - 1))*((1.0 - hzs[j])**(pbetas[i] - 1))*min(math.gamma(min(20, palphas[i] + pbetas[i])), 1e17)/(min(math.gamma(min(20, palphas[i])), 1e17)* \
                        min(math.gamma(min(20, pbetas[i])), 1e17))
                pdfs[points - 1] = 1.5*pdfs[points - 2]/pbetas[i]
                pdfs[0] = 1.5*pdfs[1]/palphas[i]
            intpdf = 0
            for j in range(1, len(hzs)):
                intpdf += (hzs[j - 1] - hzs[j])*(pdfs[j - 1] + pdfs[j])/2
            hts = np.zeros(points)
            intt = 0.0
            hts[0] = sourcelist[0]
            hts[-1] = sourcelist[-1]
            for k in range(1, len(hzs) - 1):
                ubz = 0
                for l in range(len(zs)):
                    ubz = l
                    if hzs[k] < zs[l]:
                        break
                lbz = ubz - 1
                hts[k] = (sourcelist[ubz] - sourcelist[lbz])/max(zs[ubz] - zs[lbz], 1e-14) * (hzs[k] - zs[lbz]) + sourcelist[lbz]
                intt += (hzs[k - 1] - hzs[k])*(hts[k - 1]*pdfs[k - 1] + hts[k]*pdfs[k])/(2.0*intpdf)
            intt += (hzs[-2] - hzs[-1])*(hts[-2]*pdfs[-2] + hts[-1]*pdfs[-1])/(2.0 * intpdf)
            if (i != 0) and (i != len(zs) - 1):
                datalist[i] = intt
    return datalist

def non_increasing(L):
    return all(x>=y for x, y in zip(L, L[1:]))

def massprobe(press, tempf, tempo, mdot, compf, compo, mech, width, npoint):
    print('---Flow rate probe started---')
    left_enable = True
    initial_grid = width*np.linspace(0.0, 1.0, npoint)
    gas = ct.Solution(mech, width=width)
    flame = ct.CounterflowDiffusionFlame(gas, grid=initial_grid)
    ampl = 1.0
    oldampl = ampl
    flame.P = press
    flame.fuel_inlet.X = compf
    flame.fuel_inlet.T = tempf
    flame.oxidizer_inlet.X = compo
    flame.oxidizer_inlet.T = tempo
    logstart = 0
    factor = 1.1
    air_multi = 8
    while True:
        mflux = mdot*ampl
        flame.fuel_inlet.mdot = mflux
        flame.oxidizer_inlet.mdot = mflux*air_multi
        if logstart == 0:
            flame.set_initial_guess()
        flame.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2, prune=0.03)
        try:
            flame.solve(loglevel=0)
            mixfracset = flame.mixture_fraction('Bilger')
        except ct.CanteraError:
            if factor > 1.01:
                ampl = oldampl
                factor -= 0.01
                ampl *= factor
            continue
        o2rate = flame.net_production_rates[gas.species_index('O2')]
        for j in range(len(o2rate)):
            if abs(o2rate[j]) <= 1e-4:
                o2rate[j] = 0
        left = o2rate[:np.argmin(o2rate)]
        if (not (np.sort(left) == left[::-1]).all()) or o2rate[-1] != 0:
            print('Good initail value.')
            print(max(flame.T))
            guess = flame.to_solution_array()
            oldtmax = max(flame.T)
            break
        else:
            if left_enable:
                print('Good initail value.')
                print(max(flame.T))
                guess = flame.to_solution_array()
                oldtmax = max(flame.T)
                break
            else:
                print('No need for left, minimal is {}'.format(mflux))
                left_val = mflux
                guess = flame.to_solution_array()
                oldtmax = max(flame.T)
                break
        oldampl = ampl
        ampl *= factor
    print('---Initail value verified---')
    print('---To the minimum---')
    #Left
    if left_enable:
        factor = 0.6
        ampl = 1.0
        while True:
            try:
                flame.set_initial_guess()
                oldampl = ampl
                ampl *= factor
                mflux = mdot*ampl
                flame.fuel_inlet.mdot = mflux
                flame.oxidizer_inlet.mdot = mflux*air_multi
                flame.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2, prune=0.03)
                flame.radiation_enabled = False
                mixfracset = flame.mixture_fraction('Bilger')
                flame.solve(loglevel=0, auto=True)
                mixfracset = flame.mixture_fraction('Bilger')
                if max(flame.T) < oldtmax:
                    if factor >= 0.8:
                        print('Ended...')
                        break
                    else:
                        ampl = oldampl
                        factor *= 1.25
                else:
                    if factor >= 0.8:
                        print('Ended...')
                        break
                    print('Success, mass flow rate={}'.format(mflux))
                    print('Maximum Temperature={:7.2f}'.format(max(flame.T)))
                    oldtmax = max(flame.T)
                    left_val = mflux
                    last_flame = flame.to_solution_array()
            except ct.CanteraError:
                if factor >= 0.8:
                    print('Ended...')
                    break
                else:
                    ampl = oldampl
                    factor *= 1.25
    print('---To the maximum---')
    #right
    factor = 5
    ampl = 1.0
    first = True
    while True:
        if first:
            flame.set_initial_guess(guess)
            first = False
        else:
            flame.set_initial_guess(last_flame)
        oldampl = ampl
        ampl *= factor
        mflux = mdot*ampl
        flame.fuel_inlet.mdot = mflux
        flame.oxidizer_inlet.mdot = mflux*air_multi
        flame.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2, prune=0.03)
        try:
            flame.solve(loglevel=0)
            mixfracset = flame.mixture_fraction('Bilger')
        except ct.CanteraError:
            factor *= 0.95
            if factor <= 1:
                print('Ended at {}'.format(mflux))
                extinct_val = mflux
                break
            continue
        if max(flame.T) <= max(tempf, tempo) + 50:
            ampl = oldampl
            factor *= 0.95
            if factor <= 1:
                print('Ended at {}'.format(mflux))
                extinct_val = mflux
                break
        else:
            print('Success, mass flow rate={}'.format(mflux))
            print('Maximum Temperature={:7.2f}'.format(max(flame.T)))
            right_val = mflux
            last_flame = flame.to_solution_array()
    print('Min flow rate   = {}'.format(left_val))
    print('Max flow rate   = {}'.format(right_val))
    print('Extinction rate = {}'.format(extinct_val))
    print('---Flow rate probe finished---')
    return (left_val, right_val, extinct_val)

def arrayfilter(inlist):
    outlist = np.zeros(len(inlist))
    for i in range(len(inlist)):
        if abs(inlist[i]) < 1e-8:
            outlist[i] = 0.0
        elif abs(1 - inlist[i]) < 1e-8:
            outlist[i] = 1.0
        else:
            outlist[i] = inlist[i]
    return outlist

def posfilter(inlist):
    outlist = np.zeros(len(inlist))
    for i in range(len(inlist)):
        if inlist[i] < 1e-8:
            outlist[i] = 0.0
        elif abs(1 - inlist[i]) < 1e-8:
            outlist[i] = 1.0
        else:
            outlist[i] = inlist[i]
    return outlist

# New LIM method for beta-PDF
def lim_method(favre_var, fnorm, raw_var, raw_val):
    result = 0
    var_span = len(raw_var)
    alpha = favre_var*(1/fnorm - 1)
    beta = (1 - favre_var)*(1/fnorm - 1)
    bigb_upbase = np.zeros(var_span)
    bigb_upbase[1:-1] = (raw_var[1:-1]**(alpha))*((1 - raw_var[1:-1])**(beta - 1))
    bigb_up = integrate.trapezoid(bigb_upbase, raw_var)
    bigb_downbase = np.zeros(var_span)
    bigb_downbase[1:-1] = (raw_var[1:-1]**(alpha - 1))*((1 - raw_var[1:-1])**(beta - 1))
    bigb_down = integrate.trapezoid(bigb_downbase, raw_var)
    i_down = np.empty(var_span)
    i_up = np.empty(var_span)
    for i in range(var_span):
        i_down[i] = integrate.trapezoid(bigb_downbase[:i + 1], raw_var[:i + 1])/bigb_down
        i_up[i] = i_3 = integrate.trapezoid(bigb_upbase[:i + 1], raw_var[:i + 1])/bigb_up
    for i in range(1, var_span):
        a = (raw_val[i - 1]*raw_var[i] - raw_val[i]*raw_var[i - 1])/(raw_var[i] - raw_var[i - 1])
        b = (raw_val[i] - raw_val[i - 1])/(raw_var[i] - raw_var[i - 1])
        i_1 = i_down[i]
        i_2 = i_down[i - 1]
        i_3 = i_up[i]
        i_4 = i_up[i - 1]
        val = a*(i_1 - i_2) + b*(bigb_up/bigb_down)*(i_3 - i_4)
        result += val
    return result

def group_lim(raw_var, raw_val, fnorm):
    if fnorm == 0:
        return raw_val
    else:
        pdf_list = np.empty(len(raw_var))
        for i in range(len(raw_var)):
            pdf_list[i] = lim_method(raw_var[i], fnorm, raw_var, raw_val)
        return pdf_list