import sys
import numpy as np
import ConfigParser
import scipy
from hmf.hmf import MassFunction as mf
from scipy import constants as const
from astropy import units as u
import matplotlib.pyplot as plt
def plotter(x,y):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x,y)
    plt.grid()
    plt.show()
    return 0

def clusterNum(params,verbose=True):
    # Arrange pixels
    zrange = np.arange(params['zmin']+params['dz'],params['zmax']+params['dz'],params['dz'])

    N = np.zeros([len(zrange),len(zrange)])
    i = 0
    # Looping over all pixels
    for z in zrange:
        # Get Halo Mass Function 
        halo=mf(Mmin=params['Mmin'],Mmax=params['Mmax'],dlog10m=params['dlog10m'],delta_c=params['delta_c'],delta_h=params['delta_h'],hmf_model=params['hmf_model'],z=z)
        # .cosmo_params: to change parameters, .cosmo: modified object, .cosmo_model: default object (H0 = 67.74....)
        halo.cosmo_params = params['cosmo']
      
        # Get Mass Selection Function - unitless
        P = (1./2.)*scipy.special.erfc( (np.log(10**params['Mth'])-np.log(halo.M)) / (np.sqrt(2)*params['sigmalnM']) )
        P =  np.nan_to_num(P)
        
        # Plot Mass selection function
        #plotter(np.log10(halo.M),P)
        
        # Calculate number density of clusters in ith pixel - h^3/Mpc^3
        ni = 0.
        for k in range(len(P)):
            ni += halo.dlog10m*P[k]*halo.dndlog10m[k]
        ni = ni*halo.cosmo.h**3 # Convert to 1/Mpc^3 with little h

        # Calculate volume of ith pixel - Mpc^3
        # Formula: V_i=\frac{c}{H(z)}\chi^2\Delta\Omega\Delta z
        c = const.c/1000.0*u.km/u.s
        H = halo.cosmo.H(z)
        X = halo.cosmo.comoving_transverse_distance(z)
        dOm = (np.pi/180.)**2*params['dOm'] # Convert from sq degree to steradian
        dz = params['dz']
        Vi = (c/H)*(X**2)*dOm*dz
        Vi = Vi/(u.Mpc)**3
        # Add element to N matrix
        N[i,i] = float(Vi*ni)
        i+=1

    return np.nan_to_num(N)

def main(argv):
    # Read Config                                                                                           
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
    print fparams
    # Add special case for m_nu, need to use astropy.units.u.Quantity()

    # Save fiducials                                                                                        
    print "Calculating and saving fiducial cosmology..."
    fidN = clusterNum(fparams)
    print "Number of clusters: ",np.trace(fidN)
    #sys.exit("stop!")
    np.savetxt("output/"+fparams['exp_name']+'_'+fparams['hmf_model']+"_fN.csv",fidN,delimiter=",")

    # Make derivatives
    for paramName in paramList:

        h = stepSizes[paramName]

        if paramName in cosmo:
            print "Calculating forward difference for ", paramName
            pparams = fparams.copy()
            pparams['cosmo'][paramName] = fparams['cosmo'][paramName] + 0.5*h
            pN = clusterNum(pparams)
 
            print "Calculating backward difference for ", paramName
            mparams = fparams.copy()
            mparams['cosmo'][paramName] = fparams['cosmo'][paramName] - 0.5*h
            mN = clusterNum(mparams)
            
            dN = (pN-mN)/h
            
        else:
            print "Calculating forward difference for ", paramName
            pparams = fparams.copy()
            pparams[paramName] = fparams[paramName] + 0.5*h
            pN = clusterNum(pparams)
            
            print "Calculating backward difference for ", paramName
            mparams = fparams.copy()
            mparams[paramName] = fparams[paramName] - 0.5*h
            mN = clusterNum(mparams)
            
            dN = (pN-mN)/h
        np.savetxt("output/"+fparams['exp_name']+'_'+fparams['hmf_model']+"_"+paramName+"_dN.csv",dN,delimiter=",")

if (__name__ == "__main__"):
    main(sys.argv[1:])
