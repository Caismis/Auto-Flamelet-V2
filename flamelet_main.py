import cantera as ct
import numpy as np
import scipy.special as ss
from scipy import interpolate
from funclib import *
# Advanced Parameters
massread = False
masssave = True


def flamelet_main(input_data, switch, save):
    flamenum = input_data[0]
    zst = input_data[1]
    press = input_data[2]
    tempf = input_data[3]
    tempo = input_data[4]
    mdot = input_data[5]
    compf = input_data[6]
    compo = input_data[7]
    mech = input_data[8]
    moleout = input_data[10]
    width = input_data[11]
    npoint = input_data[12]
    air_multi = 1
    
    
    if switch:
        cdef = input_data[14]
    if massread:
        critical_rates = np.load(file='mass_range.npy', allow_pickle=True)
    else:
        critical_rates = massprobe(press, tempf, tempo, mdot, compf, compo, mech, width, npoint, air_multi)
    if masssave:
        np.save(file='mass_range.npy', arr=critical_rates)
    initial_grid = width*np.linspace(0.0, 1.0, npoint)
    gas = ct.Solution(mech, width=width)
    nsp = gas.n_species
    species_names = gas.species_names
    flame = ct.CounterflowDiffusionFlame(gas, grid=initial_grid)
    ampl = 1.0
    flame.P = press
    flame.fuel_inlet.X = compf
    flame.fuel_inlet.T = tempf
    flame.oxidizer_inlet.X = compo
    flame.oxidizer_inlet.T = tempo
    dataseries = []
    dissrates = []
    tmaxs = []
    flowrates = []
    flameout = False
    logstart = 0
    i = 0
    factor = 1.1
    minrate = critical_rates[0]
    maxrate = critical_rates[1]
    exrate = critical_rates[2]
    massflow = np.logspace(np.log10(minrate), np.log10(maxrate), flamenum, endpoint=False)
    # massflow = np.logspace(np.log10(minrate), np.log10(maxrate), flamenum - 10, endpoint=True)
    # massflow = np.append(massflow, np.linspace(maxrate, exrate, 3, endpoint=True))
    massflow = np.append(massflow, exrate)
    print('---Stable branch started---')
    massflowlen = len(massflow)
    i = 0
    while i < massflowlen:
        flux = massflow[i]
        dataset = []
        flame.fuel_inlet.mdot = flux
        flame.oxidizer_inlet.mdot = flux*air_multi
        if logstart == 0:
            flame.set_initial_guess()
        else:
            flame.set_initial_guess(flameout)
        flame.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2, prune=0.03)
        while True:
            try:
                flame.solve(loglevel=0, auto=True)
                if max(flame.T) < max(tempf, tempo) + 50 and flux != exrate:
                    print('Warning: Unexpected extinction at ' + str(flux))
                    print('Reducing Flow Rate...')
                    massflow = np.insert(massflow, i + 1, massflow[i])
                    massflowlen += 1
                    massflow[i] = (flux + massflow[i - 1])/2
                    flux = massflow[i]
                    flame.set_initial_guess(flameout)
                    flame.fuel_inlet.mdot = massflow[i]
                    flame.oxidizer_inlet.mdot = massflow[i]*air_multi
                    continue
                else:
                    break
            except ct.CanteraError:
                print('Warning: Unexpected error at ' + str(flux))
                print('Reducing Flow Rate...')
                massflow[i] = (flux + massflow[i - 1])/2
                flux = massflow[i]
                flame.set_initial_guess(flameout)
                flame.fuel_inlet.mdot = massflow[i]
                flame.oxidizer_inlet.mdot = massflow[i]*air_multi
                continue
        tempset = flame.T
        tmax = max(tempset)
        rhoset = flame.density
        a = flame.strain_rate('potential_flow_oxidizer')
        dissrate = (2*a/np.pi)*np.exp(-2.0*(ss.erfcinv(2.0*zst)**2))
        mixfracset = flame.mixture_fraction('Bilger')
        mixfracset = posfilter(mixfracset)
        if min(mixfracset) != 0:
            print(mixfracset)
        if switch:
            progressvar = np.zeros(len(mixfracset))
            varrate = np.zeros(len(mixfracset))
            # for num in range(len(mixfracset)):
            #     progressvar[num] = flame.Y[gas.species_index('CO2')][num] + flame.Y[gas.species_index('H2O')][num] + flame.Y[gas.species_index('CO')][num]
            #     varrate[num] = flame.net_production_rates[gas.species_index('CO2')][num] + flame.net_production_rates[gas.species_index('H2O')][num] + \
            #     flame.net_production_rates[gas.species_index('CO')][num]
            #     varrate[num] /= flame.density_mole[num]
            for item in cdef:
                progressvar += flame.Y[gas.species_index(item.strip())]
                varrate += flame.net_production_rates[gas.species_index(item.strip())]*gas.molecular_weights[gas.species_index(item.strip())]/flame.density
            # progressvar = flame.Y[gas.species_index('CO2')] + flame.Y[gas.species_index('H2O')] + flame.Y[gas.species_index('CO')]
            progressvar = posfilter(progressvar)
            # varrate /= flux*2
        for isp in range(nsp):
            dataset.append(posfilter(list(reversed(flame.Y[isp]))))
        dataset.append(list(reversed(tempset)))
        dataset.append(list(reversed(rhoset)))
        dataset.append(list(reversed(mixfracset)))
        if switch:
            dataset.append(list(reversed(varrate)))
            dataset.append(list(reversed(progressvar)))
        o2rate = flame.net_production_rates[gas.species_index('O2')]
        for j in range(len(o2rate)):
            if abs(o2rate[j]) <= 1e-4:
                o2rate[j] = 0
        left = o2rate[:np.argmin(o2rate)]
        if logstart == 0:
            # if (not (np.sort(left) == left[::-1]).all()) or o2rate[-1] != 0:
            #     print('Thick Flame at dissipation rate={:7.4f}, discard...'.format(dissrate))
            # else:
                if logstart == 0:
                    logstart = i
                    print('Start log...')
                print('Success, Max Temperature={:7.2f}, Dissipation Rate={:7.4f}'.format(tmax, dissrate))
                print('Flow rate={:7.4f}'.format(flux))
                dataseries.append(dataset)
                dissrates.append(dissrate)
                tmaxs.append(tmax)
                flowrates.append(flux)
                flameout = flame.to_solution_array()
        else:
            if tmax < max(tempf, tempo) + 50:
                    print('Extinct, Max Temperature={:7.2f}, Dissipation Rate={:7.4f}'.format(tmax, dissrate))
                    print('Flow rate={:7.4f}'.format(flux))
                    dissrates.append(dissrate)
                    if not switch:
                        dataseries.append(dataset)
                    else:
                        extinctdata = dataset
                    tmaxs.append(tmax)
                    flowrates.append(flux)
                    if switch:
                        extinctfunc = interpolate.interp1d(flame.grid, flame.T)
                    break
            else:
                print('Success, Max Temperature={:7.2f}, Dissipation Rate={:7.4f}'.format(tmax, dissrate))
                print('Flow rate={:7.4f}'.format(flux))
                dataseries.append(dataset)
                dissrates.append(dissrate)
                tmaxs.append(tmax)
                flowrates.append(flux)
                flameout = flame.to_solution_array()
                outgrid = flame.grid
                flameoutt = flame.T
        i += 1
        ampl *= factor
    normzs = np.linspace(0, 1, npoint)
    tmax_marker = tmaxs.index(max(tmaxs))
    dissrate_stable = dissrates[tmax_marker]
    print('---Stable branch finished---')

    dissrates = dissrates[tmax_marker:]
    dataseries = dataseries[tmax_marker:]
    tmaxs = tmaxs[tmax_marker:]

    if switch:
        print('---Unstable branch started---')
        max_attempt = 30
        unstable_reached = False
        dataset = []
        for i in range(0, 6):
            # left_bound = -0.8
            # right_bound = 0.8
            left_bound = 1
            right_bound = 0
            if unstable_reached == False:
                for attemp in range(max_attempt):
                    while True:
                        flux = flowrates[-2]
                        middle = (left_bound + right_bound)/2
                        deviate = i
                        flame.from_solution_array(flameout)
                        flame.fuel_inlet.mdot = flux
                        flame.oxidizer_inlet.mdot = flux*air_multi
                        newt = list(flame.T.copy())
                        newg = flame.grid.copy()
                        # desrate = (newt[-1] - newt[0])/(newg[-1] - newg[0])
                        # desb = newt[-1] - desrate*newg[-1]
                        # newt = flameoutt.copy()
                        oldt = flame.T.copy()
                        for j in range(len(newt)):
                            # extinctt = extinctfunc(newg[j])
                            # delta = newt[j] - extinctt
                            # newt[j] -= delta*middle
                            # if newt[j] < extinctt:
                            #     newt[j] = extinctt
                            extinctt = extinctfunc(newg[j])
                            newt[j] = middle*newt[j] + (1 - middle)*extinctt
                        for j in range(deviate):
                            newt.pop(-1)
                            # newt.pop(0)
                        for j in range(deviate):
                            newt.insert(0, newt[0])
                            # newt.append(newt[0])
                        newt = np.array(newt)
                        flame.set_profile('T', newg/width, newt)
                        flame.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2, prune=0.03)
                        try:
                            flame.solve(loglevel=0, auto=True)
                            break
                        except ct.CanteraError:
                            flowrates[-3] = (flowrates[-2] + flowrates[-3])/2
                            print('Error Occurred...')
                            continue
                    tempset = flame.T
                    if max(oldt) > max(tempset) + 5 and max(tempset) > 1000:
                        print('UNSTABLE, reference temperature = {:7.2f}, current temperature = {:7.2f}'.format(max(oldt), max(tempset)))
                        unstable_start = flame.to_solution_array()
                        unstable_reached = True
                        break
                    elif max(tempset) < 1000:
                        right_bound = middle
                        print('Extinct, reference temperature = {:7.2f}, current temperature = {:7.2f}'.format(max(oldt), max(tempset)))
                        # print('Shape factor: {}'.format(1 - middle))
                        print('Shape factor: {}'.format(middle))
                    else:
                        left_bound = middle
                        print('Stable, reference temperature = {:7.2f}, current temperature = {:7.2f}'.format(max(oldt), max(tempset)))
                        # print('Shape factor: {}'.format(1 - middle))
                        print('Shape factor: {}'.format(middle))
                if not unstable_reached:
                    print('Max attempts reached, increasing deviation...')
            else:
                break
        if unstable_reached:
            print('Converged...')
            a = flame.strain_rate('potential_flow_oxidizer')
            dissrate = (2*a/np.pi)*np.exp(-2.0*(ss.erfcinv(2.0*zst)**2))
            tmax = max(tempset)
            dissrates.append(dissrate)
            tmaxs.append(tmax)
            mixfracset = flame.mixture_fraction('Bilger')
            rhoset = flame.density
            progressvar = np.zeros(len(mixfracset))
            varrate = np.zeros(len(mixfracset))
            mixfracset = posfilter(mixfracset)
            # for num in range(len(mixfracset)):
            #     progressvar[num] = flame.Y[gas.species_index('CO2')][num] + flame.Y[gas.species_index('H2O')][num] + flame.Y[gas.species_index('CO')][num]
            #     varrate[num] = flame.net_production_rates[gas.species_index('CO2')][num] + flame.net_production_rates[gas.species_index('H2O')][num] + \
            #     flame.net_production_rates[gas.species_index('CO')][num]
            #     varrate[num] /= flame.density_mole[num]
            for item in cdef:
                progressvar += flame.Y[gas.species_index(item.strip())]
                varrate += flame.net_production_rates[gas.species_index(item.strip())]*gas.molecular_weights[gas.species_index(item.strip())]/flame.density
            progressvar = posfilter(progressvar)
            # varrate /= flux*2
            for isp in range(nsp):
                dataset.append(posfilter(list(reversed(flame.Y[isp]))))
            dataset.append(list(reversed(tempset)))
            dataset.append(list(reversed(rhoset)))
            dataset.append(list(reversed(mixfracset)))
            dataset.append(list(reversed(varrate)))
            dataset.append(list(reversed(progressvar)))
            dataseries.append(dataset)
        else:
            print('Unable to catch unstable branch.')
            raise Exception('The branch cannot be catched, check the code.')

        newflowrates = flowrates[:-2]
        newflowrates = newflowrates[::-1]
        start = True
        for ii in range(len(newflowrates)):
            deviate = 0
            max_attempt = 10
            unstable_reached = False
            dataset = []
            for i in range(0, 3):
                left_bound = -0.4
                right_bound = 0.4
                if unstable_reached == False:
                    for attemp in range(max_attempt):
                        while True:
                            flux = newflowrates[ii]
                            middle = (left_bound + right_bound)/2
                            deviate = i
                            if start:
                                flame.from_solution_array(unstable_start)
                            else:
                                flame.from_solution_array(last_flame)
                            flame.fuel_inlet.mdot = flux
                            flame.oxidizer_inlet.mdot = flux*air_multi
                            newt = list(flame.T.copy())
                            newg = flame.grid.copy()
                            desrate = (newt[-1] - newt[0])/(newg[-1] - newg[0])
                            desb = newt[-1] - desrate*newg[-1]
                            oldt = flame.T.copy()
                            for j in range(len(newt)):
                                extinctt = extinctfunc(newg[j])
                                delta = newt[j] - extinctt
                                newt[j] -= delta*middle
                                if newt[j] < extinctt:
                                    newt[j] = extinctt
                            for j in range(deviate):
                                newt.pop(-1)
                            for j in range(deviate):
                                newt.insert(0, newt[0])
                            newt = np.array(newt)
                            flame.set_profile('T', newg/width, newt)
                            # for j in range(len(newt)):
                            #     extinctt = newg[j]*desrate + desb
                            #     delta = newt[j] - extinctt
                            #     newt[j] -= delta*middle
                            #     if newt[j] < extinctt:
                            #         newt[j] = extinctt
                            # for j in range(deviate):
                            #     newt.pop(-1)
                            # for j in range(deviate):
                            #     newt.insert(0, newg[0]*desrate + desb)
                            # newt = np.array(newt)
                            # flame.set_profile('T', flame.grid/width, newt)
                            flame.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2, prune=0.03)
                            try:
                                flame.solve(loglevel=0, auto=True)
                                break
                            except ct.CanteraError:
                                print('Error Occurred...')
                                break
                        tempset = flame.T
                        if abs(max(oldt) - max(tempset)) <= 400 and max(tempset) > max(tempo, tempf) + 50:
                            print('UNSTABLE, reference temperature = {:7.2f}, current temperature = {:7.2f}'.format(max(oldt), max(tempset)))
                            unstable_reached = True
                            break
                        elif max(tempset) < max(tempo, tempf) + 50:
                            right_bound = middle
                            print('Extinct, reference temperature = {:7.2f}, current temperature = {:7.2f}'.format(max(oldt), max(tempset)))
                            print('Shape factor: {}'.format(1 - middle))
                        else:
                            left_bound = middle
                            print('Stable, reference temperature = {:7.2f}, current temperature = {:7.2f}'.format(max(oldt), max(tempset)))
                            print('Shape factor: {}'.format(1 - middle))
                    if not unstable_reached:
                        print('Max attempts reached, increasing deviation...')
                else:
                    break
            if unstable_reached:
                if dissrate < dissrate_stable:
                    break
                dataset = []
                start = False
                print('Converged...')
                last_flame = flame.to_solution_array()
                a = flame.strain_rate('potential_flow_oxidizer')
                dissrate = (2*a/np.pi)*np.exp(-2.0*(ss.erfcinv(2.0*zst)**2))
                tmax = max(tempset)
                dissrates.append(dissrate)
                tmaxs.append(tmax)
                mixfracset = flame.mixture_fraction('Bilger')
                rhoset = flame.density
                progressvar = np.zeros(len(mixfracset))
                varrate = np.zeros(len(mixfracset))
                mixfracset = posfilter(mixfracset)
                # for num in range(len(mixfracset)):
                #     progressvar[num] = flame.Y[gas.species_index('CO2')][num] + flame.Y[gas.species_index('H2O')][num] + flame.Y[gas.species_index('CO')][num]
                #     varrate[num] = flame.net_production_rates[gas.species_index('CO2')][num] + flame.net_production_rates[gas.species_index('H2O')][num] + \
                #     flame.net_production_rates[gas.species_index('CO')][num]
                #     varrate[num] /= flame.density_mole[num]
                for item in cdef:
                    progressvar += flame.Y[gas.species_index(item.strip())]
                    varrate += flame.net_production_rates[gas.species_index(item.strip())]*gas.molecular_weights[gas.species_index(item.strip())]/flame.density
                # progressvar = flame.Y[gas.species_index('CO2')] + flame.Y[gas.species_index('H2O')] + flame.Y[gas.species_index('CO')]
                # varrate /= flux*2
                for isp in range(nsp):
                    dataset.append(posfilter(list(reversed(flame.Y[isp]))))
                dataset.append(list(reversed(tempset)))
                dataset.append(list(reversed(rhoset)))
                dataset.append(list(reversed(mixfracset)))
                dataset.append(list(reversed(varrate)))
                dataset.append(list(reversed(progressvar)))
                dataseries.append(dataset)
            else:
                print('Unable to catch unstable branch, some data might lost.')
                break
        dataseries.append(extinctdata)
        print('---Unstable branch finished---')
    dataseries = np.array(dataseries, dtype=object)
    if switch:
        csts = np.empty(len(dissrates))
    for i in range(len(dissrates)):
        tempfunc = interpolate.interp1d((dataseries[i][nsp + 2]), (dataseries[i][nsp]))
        dataseries[i][nsp] = tempfunc(normzs)
        rhofunc = interpolate.interp1d((dataseries[i][nsp + 2]), (dataseries[i][nsp + 1]))
        dataseries[i][nsp + 1] = rhofunc(normzs)
        if switch:
            profunc = interpolate.interp1d((dataseries[i][nsp + 2]), (dataseries[i][nsp + 4]))
            varfunc = interpolate.interp1d((dataseries[i][nsp + 2]), (dataseries[i][nsp + 3]))
            dataseries[i][nsp + 4] = profunc(normzs)
            dataseries[i][nsp + 3] = varfunc(normzs)
            csts[i] = (max(dataseries[i][nsp + 4]))
        for isp in range(nsp):
            massfunc = interpolate.interp1d((dataseries[i][nsp + 2]), (dataseries[i][isp]))
            dataseries[i][isp] = massfunc(normzs)
        dataseries[i][nsp + 2] = normzs
    print('---Creating mole weights file---')
    fi = open(moleout, 'w')
    for i in range(len(gas.species_names)):
        fi.write(np.format_float_scientific(gas.molecular_weights[i], precision=9, exp_digits=2, unique=False) + '\n')
    fi.close()
    print('---', moleout, ' generation complete!---')
    if switch:
        misc_data = np.array([nsp, species_names, max(csts)], dtype=object)
    else:
        misc_data = np.array([nsp, species_names], dtype=object)
    if save:
        np.save(file='dataseries.npy', arr=dataseries)
        if switch:
            np.save(file='csts.npy', arr=csts)
            np.save(file='dissrates.npy', arr=dissrates)
        else:
            np.save(file='dissrates.npy', arr=dissrates)
        np.save(file='misc.npy', arr=misc_data)
        print('---Data has been saved as *.npy---')
    