import sys
import numpy as np
import matplotlib.pyplot as plt
import ConfigParser
from hmf.hmf import MassFunction as mf
def getFig5(params):
    # Get massfunction for Tinker 08
    h = mf(Mmin=params['Mmin'],Mmax=params['Mmax'],dlog10m=params['dlog10m'],delta_c=params['delta_c'],delta_h=params['delta_h'],z=params['z'],hmf_model=params['hmf_model'])
    data = np.empty((len(h.M),2))
    data[:,0] = np.log10(h.M)
    data[:,1] = np.log10(h.M**2/h.mean_density*h.dndm)
    return data
def main(argv):

    # Read general Config                                                                                                                           
    iniFile = 'input/halo_massfunc.ini'
    Config = ConfigParser.SafeConfigParser()
    Config.optionxform = str
    Config.read(iniFile)
    genparams={}
    paramList=[]
    for (key, val) in Config.items('value'):
        if ',' in val:
            genparams[key] = [float(x) for x in val.split(',')]
        else:
            genparams[key] = float(val)
        paramList.append(key)
    for (key, val) in Config.items('name'):
        genparams[key] = val

    print genparams
    # Copy parameters array
    params = genparams.copy()

    # Get (x,y) data to plot Fig 5
    for key in paramList:
        if isinstance(genparams[key],list):
            print "Varying parameter ",key
            fig = plt.figure()
            ax = fig.add_subplot(111)
            for val in genparams[key]:
                params[key] = val
                data = getFig5(params)
                ax.plot(data[:,0],data[:,1],label=key+'='+str(val))
            ax.set_xlabel('log[M]')
            ax.set_ylabel('log[$(M^2/\\bar{\\rho_m})dn/dM$]')
            #ax.set_xlim([-0.65,0.5])
            ax.set_ylim(-3.6,-0.8)
            plt.legend(loc='lower left')
            plt.grid()
            figname = 'output/HaloMassFunc_'+key+"_"+params['hmf_model']+'.png'
            fig.savefig(figname, format='png')
            plt.close('all')
    
if (__name__ == "__main__"):
    main(sys.argv[1:])
