import numpy as np
import matplotlib.pyplot as pl

def getPoreIsotherm(thisDir):
    from math import sqrt
    from gcmc_toolbox import read_gcmc_log
    from gcmc_toolbox import logFileList
    from gcmc_toolbox import plotEqTimeSeries
    from gcmc_toolbox import plotRunTimeSeries

    log_files=logFileList(thisDir)
    log_files.reverse()

    mu=[]
    Np=[]
    dNp=[]
    rho_nom=[]
    drho_nom=[]

    s_eff=10 # statitical inefficiency    
    N_mu_id=234
       
    col_descr = np.dtype([('step', 'int32'), ('T', 'float64'), ('P', 'float64'), ('rho', 'float64'), ('N_f', 'float64'), \
                          ('N_t', 'float64'), ('T_f', 'float64'), ('T_s', 'float64'), ('potEng', 'float64'), ('kinEng', 'float64')])

    for i in range(len(log_files)):
        eq, mc, md = read_gcmc_log(thisDir+log_files[i], col_descr)

        plotEqTimeSeries(N_mu_id, eq.step, eq.N_f, ylabel=r'Particles in pore / $<N>$')
        plotRunTimeSeries(N_mu_id, mc.step, mc.N_f, ylabel=r'Particles in pore / $<N>$')
        plotRunTimeSeries(N_mu_id, md.step, md.N_f, color='c', ylabel=r'Particles in pore / $<N>$')

        #######################################################
        # !!!! CAREFUL With the mu values here str->float !!!!
        #######################################################        

        mu.append(log_files[i].split('mu')[1])
        Np.append(mc.N_f.mean())
        dNp.append(mc.N_f.std()/sqrt(len(mc.N_f)/s_eff))
        rho_nom.append(mc.rho.mean())
        drho_nom.append(mc.rho.std()/sqrt(len(mc.rho)/s_eff))
    
    return np.array(mu, dtype='float64'), np.array(Np, dtype='float64'), np.array(dNp, dtype='float64'),\
           np.array(rho_nom, dtype='float64'), np.array(drho_nom, dtype='float64')



## Analyze GCMC simulations
thatDir='/home/hf210/molsim/lammps_run/slit_pore/gcmc/eam_pore/'
systems= ['WCA_Hp2.5', 'LJ_Hp2.5', 'E2_Hp2.5', 'WCA_Hp5.0', 'LJ_Hp5.0', 'E2_Hp5.0']

for s in systems:
    mu_p, N_p, dN_p, rho_p, drho_p = getPoreIsotherm(thatDir+s+'/')

    ## Save data
    data=np.transpose((mu_p, N_p, dN_p))
    data.sort(axis=0)
    np.savetxt(thatDir+s+'/results.dat', data, fmt='   %6.8f'*len(data[0]))



