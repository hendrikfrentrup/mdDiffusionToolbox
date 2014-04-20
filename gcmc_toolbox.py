#!/usr/bin/env python
#==============================================
#
#          'GCMC analysis toolbox'
#               H. Frentrup
#       Mon 31 Mar 2014 15:14:25 BST 
#
#==============================================

import numpy as np
import matplotlib.pyplot as pl


def logFileList(thisDir, searchString="log.run"):
    from os import listdir
    all_files=listdir(thisDir)

    log_files=[]
    for i in range(len(all_files)):
        if not all_files[i].find(searchString)==-1:
            log_files.append(all_files[i])
    log_files.sort()
    return log_files


def read_gcmc_log(filename, dtype):
    inFile = open(filename, 'r')

    data = [[] for dummy in xrange(len(dtype))]

    for line in inFile:
        if line.split() and line.split()[0]=="Step":
            break
            
    for line in inFile:
        if line.split() and line.split()[0]=="Loop":
            break
            
        fields = line.split()
        for i, number in enumerate(fields):
            data[i].append(number)
    
    for i in xrange(len(dtype)):
        data[i] = np.cast[dtype[i]](data[i])

    eq_data = np.rec.array(data, dtype=dtype) 

    data = [[] for dummy in xrange(len(dtype))]

    for line in inFile:
        if line.split() and line.split()[0]=="Step":
            break
            
    for line in inFile:
        if line.split() and line.split()[0]=="Loop":
            break
        fields = line.split()
        for i, number in enumerate(fields):
            data[i].append(number)

    for i in xrange(len(dtype)):
        data[i] = np.cast[dtype[i]](data[i])

    mc_data = np.rec.array(data, dtype=dtype) 

    data = [[] for dummy in xrange(len(dtype))]

    for line in inFile:
        if line.split() and line.split()[0]=="Step":
            break
            
    for line in inFile:
        if line.split() and line.split()[0]=="Loop":
            break
        fields = line.split()
        for i, number in enumerate(fields):
            data[i].append(number)

    for i in xrange(len(dtype)):
        data[i] = np.cast[dtype[i]](data[i])

    run_data = np.rec.array(data, dtype=dtype) 

    inFile.close()
    return eq_data, mc_data, run_data
    
    
def plotEqTimeSeries(fig, r, s, xlabel=r'Time $t$ / $\sigma (m/\epsilon)^{1/2}$', ylabel=r'Density / $\sigma^{-3}$'):
    pl.figure(fig)
    pl.plot(r,s, 'r:')
    pl.xlabel(xlabel)
    pl.ylabel(ylabel)


def plotRunTimeSeries(fig, u, v, color='b', xlabel=r'Time $t$ / $\sigma (m/\epsilon)^{1/2}$', ylabel=r'Density / $\sigma^{-3}$', s_eff=10):
    # statitical inefficiency needs to be determined from block_avg !!!
    from math import sqrt
    pl.figure(fig)
    avg=v.copy()
    avg[:]=v.mean()
    err=v.std()/sqrt(len(v)/s_eff)
    pl.plot(u,v, color+':')
    pl.plot(u,avg, 'k-')
    pl.plot(u,avg+err, 'k--')
    pl.plot(u,avg-err, 'k--')
    pl.xlabel(xlabel)
    pl.ylabel(ylabel)    
    


def avgNtoP(N, N_p, mu_p, mu_eos, p_eos):
    from scipy import interpolate
    import numpy as np
    offset_mu_b=-21.0 #mu_b.min()
    offset_mu_eos=-8.276427954090115 #np.interp(rho_b.min(), rho_eos, mu_eos)
    tck = interpolate.splrep(np.log(N_p), mu_p, s=0.01, k=5)
    spl_lnN = np.linspace(np.log(N).min(),np.log(N).max(),200)
    spl_mu = interpolate.splev(spl_lnN,tck,der=0)
    muAtN = np.interp(np.log(N), spl_lnN, spl_mu)
    pAtN = np.interp(muAtN-offset_mu_b,mu_eos-offset_mu_eos,p_eos)
    return pAtN    
    
    
    
    
    
    
    
def block_avg(run_stats):
    from scipy import stats
    t_b = []
    var_i=[]
    blocked_stats=[]
    # calculate averaged blocks
    for i in range(1,len(run_stats)/2+1):
        blocked_stats=[]
        for j in range(0,len(run_stats),i):
            if j+i <= len(run_stats):
#                if i>330:
#                    print 'b=',i,' j=', j,'-', j+i, 'n=', len(run_stats)
#                    print run_stats[j:j+i],'=>',sum(run_stats[j:j+i])/len(run_stats[j:j+i])
                blocked_stats.append(sum(run_stats[j:j+i])/len(run_stats[j:j+i]))
#            else:
#                break
#        if blocked_stats==[]:
#            import pdb; pdb.set_trace();
            
        t_b.append(i)
        var_i.append(np.array(blocked_stats).var())
#        if i>329:
#            import pdb; pdb.set_trace();
#            print i, np.array(blocked_stats).var()

    return range(len(run_stats)/2), var_i
