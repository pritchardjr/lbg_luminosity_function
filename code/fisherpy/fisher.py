#
# Base class for Fisher matrix calculations
#
# Structure is as follows:
# 1) specify parameters in paramdict
#

import numpy as np
import scipy.linalg

#################################################################
# Fisher class definition
#################################################################

class Fisher:
    """ Base Fisher class for core manipulation """
    def __init__(self,paramdict={}):
        self._nparam=len(paramdict)

        self._fisher=None
        self._ifisher=None

        self._fishertype="Fisher"

        # Its important that matrix elements and parameters stay consistent
        # so assign each a unique ID
        self._paramdict={}
        for (i,key) in enumerate(paramdict.keys()):
            self._paramdict[key]={'value':paramdict[key],'id':i}

    def __repr__(self):
        return "%s with %s" % (self._fishertype,self._paramdict.keys())

############ variable access

    def getNParam(self):
        return self._nparam

    def getParam(self):
        return self._paramdict

    def fisherMatrix(self):
        return self._fisher

    def inverseFisherMatrix(self):
        if self._fisher is not None and self._ifisher is None:
            self.calcInverseFisher()

        return self._ifisher

########### methods

    def calcInverseFisher(self):
        if self._fisher is None:
            print "fisher matrix not calculated"
            return None
        
        try:
            self._ifisher=scipy.linalg.inv(self._fisher)
        except:
            print "problem inverting fisher"
            return None

        return self._ifisher

    def getErrors(self):
        """ Return Fisher errors """
        if self._ifisher is not None:
            errors=np.sqrt(np.diagonal(self._ifisher))
            
        errordict={}
        for (i,key) in enumerate(self._paramdict.keys()):
            errordict[key]=errors[self._paramdict[key]['id']]
            errordict[key]['value']=self._paramdict[key]['value']
        return errordict

    def calcFisherMatrix(self):
        #self._fisher=np.identity(self._nparam)
        raise Exception("calcFisherMatrix not defined")

    

##################################################################
# Helper functions
##################################################################
        
#################################################################
# Fisher matrix calculation
#################################################################
    def calculateFisher(self,tagdict):
        """ calculate Fisher matrix for parameters in tagdict """
        pass

