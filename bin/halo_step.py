import sys
import numpy as np
import ConfigParser
import matplotlib.pyplot as plt
from halo_makeDerivs import clusterNum
from time import time as timer

def main(argv):
    # Read Config
    starttime = timer()
    iniFile = "input/halo_makeDerivs.ini"
    Config = ConfigParser.SafeConfigParser()
    Config.optionxform = str
    Config.read(iniFile)
    paramList = []
    fparams = {}
    cosmo = {}
    stepSizes = {}
    fparams['hmf_model'] = Config.get('general','hmf_model')
    fparams['exp_name'] = Config.get('general','exp_name')
    for (key, val) in Config.items('hmf'):
        if ',' in val:
            param, step = val.split(',')
            paramList.append(key)
            fparams[key] = float(param)
            stepSizes[key] = float(step)
        else:
            fparams[key] = float(val)
    # Make a separate list for cosmology to add to massfunction
    for (key, val) in Config.items('cosmo'):
        if ',' in val:
            param, step = val.split(',')
            paramList.append(key)
            cosmo[key] = float(param)
            stepSizes[key] = float(step)
        else:
            cosmo[key] = float(val)
    fparams['cosmo'] = cosmo

    for paramName in paramList:
        #Make range for each parameter
        #First test: x2 the range, x0.01 the stepsize
        if paramName in cosmo:
            start = fparams['cosmo'][paramName] - stepSizes[paramName]
            end   = fparams['cosmo'][paramName] + stepSizes[paramName]
        else:
            start = fparams[paramName] - stepSizes[paramName]
            end   = fparams[paramName] + stepSizes[paramName]
        width = stepSizes[paramName]*0.01
        paramRange = np.arange(start,end+width,width)
        for paramVal in paramRange: 
            if paramName in cosmo:
                params = fparams.copy()
                params['cosmo'][paramName] = paramVal
            else:
                params = fparams.copy()
                params[paramName] = paramVal
            print paramName,paramVal
            N = clusterNum(params)
            np.savetxt("output/step/"+fparams['exp_name']+'_'+fparams['hmf_model']+"_"+paramName+"_"+str(paramVal)+".csv",N,delimiter=",")

        #----------------------------
        endtime = timer()
        print "Time elapsed: ",endtime-starttime
if (__name__ == '__main__'):
    main(sys.argv[1:])
