# Some part of this code has been taken from CosmoHammer

from copy import copy
import numpy as np
from scipy import stats
import logging
import matplotlib.pyplot as plt




class Params(object):
    """
    A key-value parameter store preserving the order of the params passed to the intializer.

    Examples::
    
        $ params = Params(("key1", [1,2,3]),
                          ("key2", [1,2,3]))
    
        $ print(params.keys)
        > ['key1', 'key2']
        
        $ print(params.key1)
        > [1, 2, 3]
        
        $ params[:,0] = 0
        
        $ print(params.values)
        > [[0 2 3]
           [0 2 3]]
           
        $ print(params[:,1])
        > [2 2]
        
    """
    
    def __init__(self, *args):

        values = []
        self._keys = []
        for k,v in args:
            if k in self._keys:
                raise KeyError("Duplicated key '%s'"%k)
            
            self.__dict__[k] = v
            self._keys.append(k)
            values.append(v)
        self._values = np.array(values)
 
    def write_paramnames_ranges_file(param_min_arr, param_max_arr, param_name_arr,param_latex_name_arr,file_path,file_root):
         file_paramnames = file_path +'/'+file_root +'.paramnames'
         file_ranges = file_path +'/'+file_root + '.ranges'
         #file_paramnames = open(file_path +'/'+file_root +'.paramnames', "w")
         #file_ranges = open(file_path +'/'+file_root + '.ranges', "w")
         ndim = len(param_min_arr)
         with open(file_ranges, 'w') as f:
             for i in range(ndim): f.write(param_name_arr[i] + ' ' + str(param_min_arr[i]) + ' ' + str(param_max_arr[i]) + '\n')
             f.write('tau'+ ' ' + str(0.05) + ' ' + str(0.1) + '\n')
         with open(file_paramnames, 'w') as f:
             for i in range(ndim): f.write(param_name_arr[i] + ' ' + param_latex_name_arr[i] + '\n')  
             f.write('tau' + ' '+ '\\tau' +  '\n' )

      
    def __getitem__(self, slice):
        return self.values[slice]
    
    def __setitem__(self, slice, value):
        self.values[slice] = value
    
    def __str__(self):
        return ",".join(("%s=%s"%(k,v) for k,v in zip(self.keys, self.values)))
    
    @property
    def keys(self):
        return copy(self._keys)

    @property
    def values(self):
        return self._values
    
    def get(self, key):
        return self.__dict__[key]
    
    def copy(self):
        return Params(*zip(self.keys, self.values))

def setspline_sigma(rho_0, sigma_sourav): #sigma_sourav is a method in NORMALIZE
	len_R=100
	logM=np.linspace(0.2,20.0,int(20.0/0.2))
	Rcube=3.0*(10**logM)/(4.0*np.pi*rho_0)
	R=Rcube**(1/3.0)
	logsig=np.zeros(len_R)
	vfunc = np.vectorize(sigma_sourav)
	logsig=np.log10(vfunc(R))
	return logM,logsig



def getLogger():
	return logging.getLogger(__name__)




def Bin_data(X,Y, statistics=None):  #for data binning
	X_unique=np.unique(X)
	#print(X_unique)
	tolX=np.abs(np.min(np.diff(X_unique)) / 2.0)
	nbinX=len(X_unique)
	Xbins_edges=np.zeros(len(X_unique)+1)
	Xbins_edges[0]=X_unique[0]-tolX
	Xbins_edges[-1]=X_unique[-1]+tolX
	Xbins_edges[1:nbinX]=(X_unique[1:]+X_unique[:-1])/2.0

	bin_Y, Xbins_edges, ignore=stats.binned_statistic(X, Y, statistic=statistics, bins=Xbins_edges )

	bin_width= abs(Xbins_edges[0]-Xbins_edges[1])
	bin_center = Xbins_edges[1:] - bin_width/2
	return bin_center, bin_Y/bin_width


def simple_plotting(X, Y):
	plt.plot(X,Y)
	plt.show()

