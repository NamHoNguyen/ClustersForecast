import sys
import numpy as np
import matlibplot.pyploy as plt
import math
import ConfigParser
from hmf.hmf import MassFunction as mf

def clusterNum(params):
    # Get Halo Mass Function
    halo = mf(Mmin=params['Mmin'],Mmax=params['Mmax'],dlog10m=params['dlog10m'],delta_c=params['delta_c'],delta_h=params['delta_h'],z=params['z'],hmf_model=params['hmf_model'])
    # Get Mass Selectrion Function
    P = (1./2.)*math.erfc( (np.log(params['Mth'])-np.log(halo.M)) / (np.sqrt(2)*params['sigmalnM']) )
    P =  np.nan_to_num(P)
    Ni = 0.
    for i in range(0,len(P)):
        Ni += halo.dlog10m*P[i]*halo.dndlog10m[i]
    return Ni*params['Vi']

def main(argv):
    print "Needs Ni, Sij, bi, b"
if (__name__ == "__main__"):
    main(sys.argv[1:])
