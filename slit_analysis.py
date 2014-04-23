#!/usr/bin/env python
#==============================================
#
#'Analyse output from slit pore EMD simulations'
#                H. Frentrup
#                2012-08-03
#
#==============================================
import pdb;

def readEMD_slitPore(filename, sample=10, corr=500000, run=5000000):
    import numpy as np    
    
    # sample frequency needs to be set?
    
    inFile = open(filename, "r")

    runLength = int(run) # length of the entire simulation
    corrLength = int(corr)  # length of the correlation
    nFreq = int(sample) # sample frequency

    nSamples = runLength/corrLength  # number of samples taken
    dtSample = corrLength/nFreq  # length of sample

    t=[]

    data=[[], [], [], []]
    #read first two lines and # of particles
    line = next(inFile).split()
    next(inFile)  # skip first line

    while line:
        try:
            line = next(inFile).split()
            t.append(line[0])
            data[0].append(line[1])
            data[1].append(line[2])
            data[2].append(line[3])
            data[3].append(line[4])
        except StopIteration:
            line=[]
    inFile.close()

    t.insert(0,'0.0')
    t.pop()
    data[0].insert(0, data[0].pop())
    data[1].insert(0, data[1].pop())
    data[2].insert(0, data[2].pop())
    data[3].insert(0, data[3].pop())

    time = np.array(t, dtype = 'float64')
    numData = np.array(data, dtype = 'float64')
            
    accData = np.zeros((4,dtSample-1))

    for i in np.array(range(nSamples))*dtSample:
        accData[0] += numData[0][i:i+dtSample-1]
        accData[1] += numData[1][i:i+dtSample-1]
        accData[2] += numData[2][i:i+dtSample-1]
        accData[3] += numData[3][i:i+dtSample-1]
    accData /= len(np.array(range(nSamples))*dtSample)
    
    return (time[0:dtSample-1], accData[0], accData[1], accData[2], accData[3])


def intVACF(time, vacf):
    import scipy.integrate
    return scipy.integrate.trapz(vacf, time)

    
def plotVACF(time, vacf, n=1):
    import pylab as pl
    pl.figure(n)
    pl.plot(time, vacf)
    pl.xlabel(r'Time $t$ / $\sigma (m/\epsilon)^{1/2}$')
    pl.ylabel(r'Velocity Auto Correlation Function / $\epsilon$/m')

    
def plotMSD(time, msd, n=1):    
    import pylab as pl
    pl.figure(n)
    pl.plot(time, msd)
    pl.xlabel(r'Time $t$ / $\sigma (m/\epsilon)^{1/2}$')
    pl.ylabel(r'Mean square displacement / $\sigma^{2}$')
    
    
def fitMSD(time, msd):     
    from scipy import polyfit
    lfb = round(len(msd)/4)-1 # Lower Fit Bound
    ufb = len(msd)-1 # Upper Fit Bound
    (a, b) = polyfit(time[lfb:ufb],msd[lfb:ufb],1)
    return (a, b)


def plotFitMSD(time, msd, fit_coeff, n=1):     
    import pylab as pl
    from scipy import polyval

    (a, b) = fit_coeff

    pl.figure(n)    
    fitData = polyval([a, b],time)
    pl.plot(time, fitData,':')
    pl.plot(time, msd, '-')
    pl.xlabel(r'Time $t$ / $\sigma (m/\epsilon)^{1/2}$')
    pl.ylabel(r'Mean square displacement / $\sigma^{2}$')






def readEMD_COM(filename):
    import numpy as np    
    inFile = open(filename, "r")
    t=[]
    data=[[], [], []]
    line = next(inFile).split() # skip first line

    while line:
        try:
            line = next(inFile).split()
            t.append(line[0])
            data[0].append(line[1])
            data[1].append(line[2])
            data[2].append(line[3])
        except StopIteration:
            line=[]
    inFile.close()

    time = np.array(t, dtype = 'float64')
    numData = np.array(data, dtype = 'float64')

    return (time, numData[0], numData[1], numData[2])


def calcDcfromCOMarray(com, corrLength=50000):
    from sys import stdout
    from sys import stdout
    import numpy as np
    corrL=corrLength
    
    counter=0    # initialize  counters for # of samples
    sampleCounter = np.zeros(int(corrL))
    time = np.zeros(int(corrL))    # initialize time and collective MSD vector
    collMSD = np.zeros((3, int(corrL)))

    for i,t in enumerate(com[0]):
        r = (com[1][i], com[2][i], com[3][i])
        if (counter==0 or counter>=corrL):
            stdout.write('\rReached t = '+ '{0:5.4f}'.format(t) + ' (trun = ('+'{0:5.4f}'.format(com[0][-1])+')')
            stdout.flush()
            counter = 0   # set counter to zero     
            t0 = t        # take first sample
            r0 = r
            sampleCounter[counter] += 1
        else: # take a normal sample
            time[counter] = t-t0
            collMSD[0][counter] += (r[0]-r0[0])**2
            collMSD[1][counter] += (r[1]-r0[1])**2
            collMSD[2][counter] += (r[2]-r0[2])**2
            sampleCounter[counter] += 1
        counter += 1
        
    collMSD[0] /= sampleCounter
    collMSD[1] /= sampleCounter
    collMSD[2] /= sampleCounter
    
    stdout.write('\nDone at t = '+ '{0:5.4f}'.format(t) + ' (trun = ('+'{0:5.4f}'.format(com[0][-1])+')\n')
    stdout.flush()    
    return (time, collMSD[0], collMSD[1], collMSD[2])


def otherDcfromCOMarray(com, corrLength=50000):
    from sys import stdout
    from sys import stdout
    import numpy as np
    corrL=corrLength
    
    counter=0    # initialize  counters for # of samples
    sampleCounter = []
 
    # initialize time and collective MSD vector
    tOrigins = []
    rOrigins = []
    time= []
    collMSD = [[], [], []]
        
    tOrigins.append(com[0][0]) # take initial time origin
    rOrigins.append( (com[1][0], com[2][0], com[3][0]) ) # take initial vector origin
 
    for i,t in enumerate(com[0]):
        r = np.array([com[1][i], com[2][i], com[3][i]])
        
        if (counter>=corrL):
            stdout.write('\rReached t = '+ '{0:5.4f}'.format(t) + ' (trun = ('+'{0:5.4f}'.format(com[0][-1])+')')
            stdout.flush()
            counter=0   # set counter to zero     
            tOrigins.append(t)  # take first sample
            rOrigins.append(r)  # ------"-------

        # take a normal sample
        time.append(t-tOrigins[0])
        disp=r-rOrigins[0]
        
        collMSD[0].append(disp[0]*disp[0])
        collMSD[1].append(disp[1]*disp[1])
        collMSD[2].append(disp[2]*disp[2])
        sampleCounter.append(1)
        for i in xrange(1,len(tOrigins)):
            t0=tOrigins[i]
            r0=rOrigins[i]
            index=len(time)-i*corrL
            disp=r-r0
            collMSD[0][index]+=disp[0]*disp[0]
            collMSD[1][index]+=disp[1]*disp[1]
            collMSD[2][index]+=disp[2]*disp[2]
            sampleCounter[index]+=1
            #pdb.set_trace()
        counter+=1

    #pdb.set_trace()        
    collMSD[0] /= np.array(sampleCounter)
    collMSD[1] /= np.array(sampleCounter)
    collMSD[2] /= np.array(sampleCounter)
    
    stdout.write('Reached t = '+ '{0:5.4f}'.format(t) + ' (trun = ('+'{0:5.4f}'.format(com[0][-1])+')\n')
    stdout.flush()      
    return (np.array(time), collMSD[0], collMSD[1], collMSD[2]) 


def cMSDarray(com):
    import numpy as np
    # initialize time and collective MSD vector
    time = np.zeros(len(com[0]))
    collMSD = np.zeros( (3, len(com[0])) )

    # set origin
    t0 = com[0][0]        
    r0 = np.array([com[1][0], com[2][0], com[3][0]])
    for i,t in enumerate(com[0]):
        time[i] = t-t0
        collMSD[0][i] += (com[1][i]-r0[0])**2
        collMSD[1][i] += (com[2][i]-r0[1])**2
        collMSD[2][i] += (com[3][i]-r0[2])**2
    return (time, collMSD[0], collMSD[1], collMSD[2])










def yaDcfromCOMarray(com, corrLength=20000):
    import numpy as np
    corrL=corrLength
    
    counter=0    # initialize  counters for # of samples
    sampleCounter = []
 
    # initialize time and collective MSD vector
    tOrigins = []
    rOrigins = []
    time= []
    collMSD = [[], [], []]
        
    tOrigins.append(com[0][0]) # take initial time origin
    rOrigins.append( (com[1][0], com[2][0], com[3][0]) ) # take initial vector origin
 
    for i,t in enumerate(com[0]):
        r = np.array([com[1][i], com[2][i], com[3][i]])
        
        if (counter>=corrL):   
            print( '\rReached t = '+ '{0:5.4f}'.format(t) + ' (trun = ('+'{0:5.4f}'.format(com[0][-1])+')' )
            counter=0   # set counter to zero     
            tOrigins.append(t)  # take first sample
            rOrigins.append(r)  # ------"-------

        # take a normal sample
        time.append(t-tOrigins[0])
        disp=r-rOrigins[0]
        pdb.set_trace()
        collMSD[0].append(disp[0]*disp[0])
        collMSD[1].append(disp[1]*disp[1])
        collMSD[2].append(disp[2]*disp[2])
        sampleCounter.append(1)
        for i in xrange(1,len(tOrigins)):
            t0=tOrigins[i]
            r0=rOrigins[i]
            index=len(time)-i*corrL
            disp=r-r0
            collMSD[0][index]+=disp[0]*disp[0]
            collMSD[1][index]+=disp[1]*disp[1]
            collMSD[2][index]+=disp[2]*disp[2]
            sampleCounter[index]+=1
        counter+=1
        
    collMSD[0] /= np.array(sampleCounter)
    collMSD[1] /= np.array(sampleCounter)
    collMSD[2] /= np.array(sampleCounter)
        
    return (np.array(time), collMSD[0], collMSD[1], collMSD[2]) 





def cVCFarray(cvel):
    import numpy as np
    # initialize time and collective MSD vector
    time = np.zeros(len(cvel[0]))
    collVCF = np.zeros( (3, len(cvel[0])) )

    # set origin
    t0 = cvel[0][0]        
    v0 = np.array([cvel[1][0], cvel[2][0], cvel[3][0]])
    for i,t in enumerate(cvel[0]):
        time[i] = t-t0
        collVCF[0][i] += (cvel[1][i]*v0[0])
        collVCF[1][i] += (cvel[2][i]*v0[1])
        collVCF[2][i] += (cvel[3][i]*v0[2])
    return (time, collVCF[0], collVCF[1], collVCF[2])


   
   
   
    


def calcDcfromCOM(filename, corrLength=10000):
    import numpy as np

    corrL=corrLength
    inFile = open(filename, "r")
    line = next(inFile) # skip first line

    # initialize time and counter for # of samples
    time = np.zeros(int(corrL))
    sampleCounter = np.zeros(int(corrL))
    counter=0
    # initialize collective MSD vector
    collMSD = np.zeros((3, int(corrL)))

    for line in inFile:
        # read current profile
        data = np.array(line.split())
        t=data[0].astype(float)
        r=np.array([data[1], data[2], data[3]], dtype = 'float64')

        if (counter==0 or counter>=corrL):
            print 'reached t=', t
            counter=0   # set counter to zero     
            t0=t        # take first sample
            r0=r
            sampleCounter[counter]+=1
        else: # take a normal sample
            time[counter]=t-t0
            disp=r-r0
            collMSD[0][counter]+=disp[0]*disp[0]
            collMSD[1][counter]+=disp[1]*disp[1]
            collMSD[2][counter]+=disp[2]*disp[2]
            sampleCounter[counter]+=1
        counter+=1
    
    inFile.close()

    collMSD[0]/=sampleCounter
    collMSD[1]/=sampleCounter
    collMSD[2]/=sampleCounter
        
    return (time, collMSD[0], collMSD[1], collMSD[2])


def otherDcfromCOM(filename, corrLength=5000):
    import numpy as np

    corrL=corrLength
    inFile = open(filename, "r")
    line = next(inFile) # skip first line

    # initialize time and counter for # of samples
    counter=0
    sampleCounter = []

    tOrigins = []
    rOrigins = []
    time= []
    collMSD = [[], [], []]

    # take first sample
    line = next(inFile)
    data = np.array(line.split())
    t=t=data[0].astype(float)
    r=np.array([data[1], data[2], data[3]], dtype = 'float64')
    tOrigins.append(t) # take initial time origin
    rOrigins.append(r) # take initial vector origin

    for line in inFile:
        # read current current line
        data = np.array(line.split())
        t=t=data[0].astype(float)
        r=np.array([data[1], data[2], data[3]], dtype = 'float64')

        if (counter>=corrL):   
            print 'reached t=', t
            counter=0   # set counter to zero     
            tOrigins.append(t)  # take first sample
            rOrigins.append(r)  # ------"-------

        # take a normal sample
        time.append(t-tOrigins[0])
        disp=r-rOrigins[0]
        collMSD[0].append(disp[0]*disp[0])
        collMSD[1].append(disp[1]*disp[1])
        collMSD[2].append(disp[2]*disp[2])
        sampleCounter.append(1)
        for i in range(1,len(tOrigins)):
            t0=tOrigins[i]
            r0=rOrigins[i]
            index=len(time)-i*corrL
            disp=r-r0
            collMSD[0][index]+=disp[0]*disp[0]
            collMSD[1][index]+=disp[1]*disp[1]
            collMSD[2][index]+=disp[2]*disp[2]
            sampleCounter[index]+=1
        counter+=1
    
    inFile.close()
    collMSD[0]/=np.array(sampleCounter)
    collMSD[1]/=np.array(sampleCounter)
    collMSD[2]/=np.array(sampleCounter)
        
    return (np.array(time), collMSD[0], collMSD[1], collMSD[2]) 
 
 
 
 
 
 
 
    
    
def calcDcfromCVF(filename, corrLength=10000):
    import numpy as np

    corrL=corrLength
    inFile = open(filename, "r")
    line = next(inFile) # skip first line

    # initialize time and counter for # of samples
    time = np.zeros(int(corrL))
    sampleCounter = np.zeros(int(corrL))
    counter=0
    # initialize collective MSD vector
    collVCF = np.zeros((3, int(corrL)))

    for line in inFile:
        # read current profile
        data = np.array(line.split())
        t=data[0].astype(float)
        v=np.array([data[1], data[2], data[3]], dtype = 'float64')

        if (counter==0 or counter>=corrL):
            print 'reached t=', t
            counter=0   # set counter to zero     
            t0=t        # take first sample
            v0=v
            sampleCounter[counter]+=1
        else: # take a normal sample
            time[counter]=t-t0
            vDotv=v*v0
            collVCF[0][counter]+=vDotv[0]
            collVCF[1][counter]+=vDotv[1]
            collVCF[2][counter]+=vDotv[2]
            sampleCounter[counter]+=1
        counter+=1
    
    inFile.close()

    collVCF[0]/=sampleCounter
    collVCF[1]/=sampleCounter
    collVCF[2]/=sampleCounter
        
    return (time, collVCF[0], collVCF[1], collVCF[2])    
    
    



# takes ages and not properly tried, although seems to work ok
def contDcfromCOM(filename, corrLength=5000):
    import numpy as np

    corrL=corrLength
    inFile = open(filename, "r")
    line = next(inFile) # skip first line

    # initialize time and counter for # of samples
    counter=0
    sampleCounter = []

    tOrigins = []
    rOrigins = []
    time= []
    collMSD = [[], [], []]

    for line in inFile:
        # read current current line
        data = np.array(line.split())
        t=t=data[0].astype(float)
        r=np.array([data[1], data[2], data[3]], dtype = 'float64')
        tOrigins.append(t) # take initial time origin
        rOrigins.append(r) # take initial vector origin

        # take a normal sample
        time.append(t-tOrigins[0])
        disp=r-rOrigins[0]
        collMSD[0].append(disp[0]*disp[0])
        collMSD[1].append(disp[1]*disp[1])
        collMSD[2].append(disp[2]*disp[2])
        sampleCounter.append(1)
        for i in range(1,len(tOrigins)):
            t0=tOrigins[i]
            r0=rOrigins[i]
            index=len(time)-i
            disp=r-r0
            collMSD[0][index]+=disp[0]*disp[0]
            collMSD[1][index]+=disp[1]*disp[1]
            collMSD[2][index]+=disp[2]*disp[2]
            sampleCounter[index]+=1
        counter+=1
        if (counter%corrL==0):
            print 'reached t=', t
    
    inFile.close()
    collMSD[0]/=np.array(sampleCounter)
    collMSD[1]/=np.array(sampleCounter)
    collMSD[2]/=np.array(sampleCounter)
        
    return (np.array(time), collMSD[0], collMSD[1], collMSD[2])

# also takes ages and not properly tried
def com2Dc(time,com):
    import numpy as np

    counter=0
    sampleCounter = [] 
    collMSD = [[], [], []]
    for j in range(len(time)):
        # read current current line
        t=time[j]
        r=np.array([com[0][j], com[1][j], com[2][j]])
        
        # take a normal sample
        disp=r-np.array([com[0][0], com[1][0], com[2][0]])
        collMSD[0].append(disp[0]*disp[0])
        collMSD[1].append(disp[1]*disp[1])
        collMSD[2].append(disp[2]*disp[2])
        sampleCounter.append(1)
        for i in range(1,j):
            t0=time[i]
            r0=np.array([com[0][i], com[1][i], com[2][i]])
            index=j-i
            disp=r-r0
            collMSD[0][index]+=disp[0]*disp[0]
            collMSD[1][index]+=disp[1]*disp[1]
            collMSD[2][index]+=disp[2]*disp[2]
            sampleCounter[index]+=1
        counter+=1
        if (counter%5000==0):
            print 'reached t=', t

    pdb.set_trace()
    collMSD[0]/=np.array(sampleCounter)
    collMSD[1]/=np.array(sampleCounter)
    collMSD[2]/=np.array(sampleCounter)
        
    return (np.array(time), collMSD[0], collMSD[1], collMSD[2])



# this has the wrong 'data output format', so rewritten above
def oldDcfromCOM(filename, corrLength=5000):
    import numpy as np

    corrL=corrLength
    inFile = open(filename, "r")
    line = next(inFile) # skip first line

    # initialize time and counter for # of samples
    time = np.zeros((int(corrL), 1))
    sampleCounter = np.zeros((int(corrL), 1))
    counter=0
    # initialize collective MSD vector
    collMSD = np.zeros((int(corrL), 3))

    for line in inFile:
        # read current profile
        data = np.array(line.strip('()\n ').split())
        t=np.array(data[0].astype(float)) #*dt
        r=np.array(data[1:4].astype(float))

        if (counter==0 or counter>=corrL):
            print counter, t
            counter=0   # set counter to zero     
            t0=t        # take first sample
            r0=r
            # set displacement to zero
            disp=np.array([0.0,0.0,0.0])
            sampleCounter[counter]+=1
        else: # take a normal sample
            disp+=r-r0
            time[counter]=t-t0
            collMSD[counter,:]+=disp*disp
            r0=r
            sampleCounter[counter]+=1
        counter+=1
    inFile.close()
    
    for i in range(sampleCounter.shape[0]):
        collMSD[i,:]=collMSD[i,:]/sampleCounter[i]
        
    return (time, collMSD[:,0], collMSD[:,1], collMSD[:,2])

