import sys
import numpy as np
import ConfigParser
import itertools

class haloFisher:

    def __init__(self,iniFile):

        self.Config = ConfigParser.SafeConfigParser()
        self.Config.read(iniFile)
       # self.Config.optionxform = str                                                                                                                                    

        print iniFile
        self.paramList = self.Config.get('general','paramList').split(',')
        self.numParams = len(self.paramList)

        self.totFisher = np.zeros((self.numParams,self.numParams))
         
    def calcFisher(self,verbose=True):

        Config = self.Config
 
        derivRoot = Config.get('general','derivRoot')

        fidN = np.loadtxt("output/"+derivRoot+"_fN.csv",delimiter=",")
        dN = {}
        for paramName in self.paramList:
            dN[paramName] = np.loadtxt("output/"+derivRoot+"_"+paramName+"_dN.csv",delimiter=",")
        
        # Calculate covariance matrix
        Cov = np.zeros([len(fidN),len(fidN)])
        for i in range(0,len(fidN)):
            for j in range(0,len(fidN)):
                # Add sample variance
                S = 0.
                Cov[i,j] = fidN[i,j] + S
        InvCov = np.linalg.inv(Cov)

        paramCombs = itertools.combinations_with_replacement(self.paramList,2)
        Fisher = np.zeros((self.numParams,self.numParams))

        for param1,param2 in paramCombs:
            if verbose: print "Parameter combination : ", param1,param2
            u = self.paramList.index(param1)
            v = self.paramList.index(param2)
            Fij = 0.
            for i in range(0,len(fidN)):
                for j in range(0,len(fidN)):
                    Fij += dN[param1][i,i]*InvCov[i,j]*dN[param2][j,j]
            Fisher[u,v] = Fij
            Fisher[v,u] = Fij
        
        self.totFisher += Fisher

        if verbose:
            print self.totFisher
            print np.linalg.det(self.totFisher)
            for param in self.paramList:
                print param, " ," ,self.margSigma(param)


        return self.totFisher

    def margSigma(self,param):
        i = self.paramList.index(param)
        return np.sqrt(np.linalg.inv(self.totFisher)[i,i])

def main(argv):

    try:
        iniFile = argv[0]
    except:
        iniFile = "input/halo_fisher.ini"

    F = haloFisher(iniFile)
    print "Calculating Fisher matrix..."
    F.calcFisher(verbose = True)

if (__name__ == "__main__"):
    main(sys.argv[1:])

