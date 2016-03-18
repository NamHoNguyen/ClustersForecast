import sys
import numpy as np
import matplotlib.pyplot as plt
import ConfigParser
from hmf import fitting_functions as ff

'''
#----- Check Table 2 -----
# defaults: delta_c = 1.686, delta_halo = 200
print "Table 2 - Tinker 2008"
for delta_halos in [200,300,400,600,800,1200,1600,2400,3200]:
    hmf = ff.Tinker08(nu2=nu2,delta_c=delta_c,delta_halo=delta_halos,z=z)
    print hmf.delta_halo,hmf.A,hmf.a,hmf.b,hmf.c

'''
def hmf_Tinker08(params):
    '''
    Get fsigma array for a given nu or sigma array at delta_halo and z
    '''
    if not('nu2' in params):
        params['nu2']=(params['delta_c']/10**(-params['neglogsigma']))**2
    hmf = ff.Tinker08(nu2=params['nu2'],delta_c=params['delta_c'],z=params['z'],delta_halo=params['delta_halo'])
    #print "log f(sigma) array for delta_halo = ",hmf.delta_halo," z = ",hmf.z
    #print np.log10(hmf.fsigma)
    return hmf

def main(argv):

    # Read general Config
    iniFile = 'input/halo_fit.ini'
    Config = ConfigParser.SafeConfigParser()
    Config.optionxform = str
    Config.read(iniFile)
    genparams={}
    paramList=[]
    for (key, val) in Config.items('Tinker08'):
        if ',' in val:
            start, end, step = val.split(',')
            genparams[key] = np.arange(float(start),float(end)+float(step),float(step))
        else:
            genparams[key] = float(val)
            paramList.append(key)

    print "----------------------------"
    #----- Check Table 2 -----
    # defaults: delta_c = 1.686, delta_halo = 200
    params = genparams.copy()
    print "Table 2 - Tinker 2008"
    for delta_halo in [200,300,400,600,800,1200,1600,2400,3200]:
        params['delta_halo'] = delta_halo
        hmf = hmf_Tinker08(params)
        print hmf.delta_halo,"\t",round(hmf.A,3),"\t",round(hmf.a,2),"\t",round(hmf.b,2),"\t",round(hmf.c,2)
    
    print "----------------------------"    
    # Produce 3 curves in fig 5, Tinker 2008
    print "Producing fig 5, Tinker 2008"
    params = genparams.copy()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for delta_halo in genparams['delta_halo']:
        params['delta_halo'] = delta_halo
        fsigma = hmf_Tinker08(params).fsigma
        if delta_halo == 200 or delta_halo == 800 or delta_halo == 3200:
            ax.plot(params['neglogsigma'],np.log10(fsigma),label='delta_halo='+str(delta_halo))
    ax.set_xlabel('log $1/\sigma$')
    ax.set_ylabel('log f($\sigma$)')
    ax.set_xlim([-0.65,0.5])
    ax.set_ylim([-3,0])
    ax.set_title('Halo Mass Function z = 0 (Tinker08 fit)')
    plt.legend()
    plt.grid()
    plt.show()


if (__name__ == "__main__"):
    main(sys.argv[1:])
