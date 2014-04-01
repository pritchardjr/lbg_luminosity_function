#
# Code for fisher matrix based on Cl i.e. 2D
#
#
#
#

import numpy as np
import fisher

class FisherCL(fisher.Fisher):
    """ Fisher class for 2D Cl derived calculations """
    def __init__(self,paramdict={'oml':0.7,'hubble':0.7,'omm':0.3,'omk':0.0,'Ombhh':0.0225,'nscal':0.95,'Ascal':25.0}):
        fisher.Fisher.__init__(self,paramdict)
        self._fishertype="FisherCL"

    def calcFisherMatrix(self):
        #self._fisher=np.identity(self._nparam)
        raise Exception("calcFisherMatrix not defined")
