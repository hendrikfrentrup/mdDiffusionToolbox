import numpy as np
import matplotlib.pyplot as pl
#import plotToolboxThesis as box



def getBulkIsotherm(thisDir):
    from math import sqrt
    from gcmc_toolbox import read_gcmc_log
    from gcmc_toolbox import logFileList
    from gcmc_toolbox import plotEqTimeSeries
    from gcmc_toolbox import plotRunTimeSeries

    log_files=logFileList(thisDir)
    log_files.reverse()

    mu=[]
    rho=[]
    drho=[]
    p=[]
    dp=[]

    s_eff=10
    mu_rho_id=123
    plot_id=124
    
    col_descr = np.dtype([('step', 'int32'), ('T', 'float64'), ('P', 'float64'), ('rho', 'float64'), \
                          ('N', 'float64'), ('potEng', 'float64'), ('kinEng', 'float64')])
 
    for i in range(len(log_files)):
        eq, mc, md = read_gcmc_log(thisDir+log_files[i], col_descr)

#        plotEqTimeSeries(mu_rho_id, eq.step, eq.rho)
#        plotRunTimeSeries(mu_rho_id, mc.step, mc.rho)
#        plotRunTimeSeries(mu_rho_id, md.step, md.rho, color='c')

#        plotEqTimeSeries(p_rho_id, eq.step, eq.P, ylabel=r'Pressure $P$ / $\epsilon\sigma^{-3}$')
#        plotRunTimeSeries(p_rho_id, mc.step, mc.P, ylabel=r'Pressure $P$ / $\epsilon\sigma^{-3}$')
#        plotRunTimeSeries(p_rho_id, md.step, md.P, color='c', ylabel=r'Pressure $P$ / $\epsilon\sigma^{-3}$')
        
        #######################################################
        # !!!! CAREFUL With the mu values here str->float !!!!
        #######################################################        
        mu.append(log_files[i].split('mu')[1])
        rho.append(mc.rho.mean())
        drho.append(mc.rho.std()/sqrt(len(mc.rho)/s_eff))
        p.append(mc.P.mean())
        dp.append(mc.P.std()/sqrt(len(mc.P)/s_eff))

    return ( np.array(mu, dtype='float64'), np.array(rho, dtype='float64'), np.array(drho, dtype='float64') ,\
                                              np.array(p, dtype='float64'),   np.array(dp, dtype='float64'))


def getPressIsotherm(thisDir):
    from math import sqrt
    from gcmc_toolbox import read_gcmc_log
    from gcmc_toolbox import logFileList
    from gcmc_toolbox import plotEqTimeSeries
    from gcmc_toolbox import plotRunTimeSeries

    log_files=logFileList(thisDir, searchString='log.nvt')

    rho=[]
    p=[]
    dp=[]

    s_eff=10
    plot_id=765
    
    col_descr = np.dtype([('step', 'int32'), ('T', 'float64'), ('P', 'float64'), \
                          ('potEng', 'float64'), ('kinEng', 'float64'), ('rho', 'float64'), ('N', 'float64')])
 
    for i in range(len(log_files)):
        eq, md, empty = read_gcmc_log(thisDir+log_files[i], col_descr)

#        plotEqTimeSeries(plot_id, eq.step, eq.P, ylabel=r'Pressure $P$ / $\epsilon\sigma^{-3}$')
#        plotRunTimeSeries(plot_id, md.step, md.P, ylabel=r'Pressure $P$ / $\epsilon\sigma^{-3}$')

        rho.append(md.rho.mean())
        p.append(md.P.mean())
        dp.append(md.P.std()/sqrt(len(md.P)/s_eff))

    return ( np.array(rho, dtype='float64'), np.array(p, dtype='float64'),   np.array(dp, dtype='float64'))




## Analyze GCMC simulations
thatDir='/home/hf210/molsim/lammps_run/slit_pore/gcmc/bulk/'
mu_b, rho_b, drho_b, p_b, dp_b = getBulkIsotherm(thatDir+'mu-18/')

## Save data
data=np.transpose((mu_b,rho_b,drho_b,p_b,dp_b))
data.sort(axis=0)
np.savetxt(thatDir+'mu-18/'+'results.dat', data, fmt='   %6.8f'*len(data[0]))
