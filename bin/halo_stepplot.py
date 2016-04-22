import sys
import numpy as np
import matplotlib.pyplot as plt
import ConfigParser
from glob import glob

def getNvsParam(paramName,root,zrange):
    fnames=glob("output/step/"+root+"_"+paramName+"_*.csv")
    #print len(fnames)
    data={}
    for fname in fnames:
        f = np.loadtxt(fname,delimiter=",")
        string = fname.split('_')
        paramVal = string[len(string)-1].split('.')
        paramVal.pop(len(paramVal)-1)
        paramVal = float('.'.join(paramVal))
        N = []
        for i in range(len(f)):
            N.append(f[i,i])
        data[paramVal] = N
    for i in range(len(zrange)):
        z = zrange[i]
        plotdata = []
        for paramVal in data:
            plotdata.append([paramVal,data[paramVal][i]])
        plotdata = np.sort(np.array(plotdata),axis=0)
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        figName = "N vs. "+paramName+" (z="+str(z)+")"
        ax.set_title(figName, fontsize=20)
        ax.plot(plotdata[:,0],plotdata[:,1],'*')
        plt.savefig("output/"+root+"_"+figName+".png",format='png')
        #plt.show()
        plt.close()
        print "Saved figure ", figName

def main(argv):
    # Read Config
    iniFile = "input/halo_stepplot.ini"
    Config = ConfigParser.SafeConfigParser()
    Config.optionxform = str
    Config.read(iniFile)
    paramList = Config.get('general','paramList').split(',')
    root = Config.get('general','root')
    start,end,step = [float(x) for x in Config.get('general','zrange').split(',')]
    zrange = np.arange(start+step,end+step,step)
    print paramList,root,zrange
    for paramName in paramList:
        getNvsParam(paramName,root,zrange)
            

if (__name__ == '__main__'):
    main(sys.argv[1:])
