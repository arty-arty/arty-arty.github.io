import os
import sys

###BLOCKING CONSOLE OUTPUT
def redirect_stdout():
    sys.stdout.flush() # <--- important when redirecting to files
    newstdout = os.dup(1)
    devnull = os.open(os.devnull, os.O_WRONLY)
    os.dup2(devnull, 1)
    os.close(devnull)
    sys.stdout = os.fdopen(newstdout, 'w')

import ctypes
from ctypes import cdll
from ctypes import *

class MuMuPy:
    def __init__(self, pathToBinary = "/results/CrossSection.so", redirect = False):
        self.pathToBinary = pathToBinary
        if redirect: redirect_stdout()        
        ###IMPORTING THE COMPILED FORTRAN LIBRARY
        try:
            self.libc = cdll.LoadLibrary("/results/CrossSection.so")
            self.libc.wtcrs_.restype = ctypes.c_double
            self.libc.wtrcrs_.restype = ctypes.c_double
        except:
            self.libc = cdll.LoadLibrary(pathToBinary)
            self.libc.wtcrs_.restype = ctypes.c_double
            self.libc.wtrcrs_.restype = ctypes.c_double
            
    def calcCrossSection(self, iq1, iq2, Z, n, l, m):
        '''Calculate the total Cross-Section of dimuonium transitions from state (n, l, m) to some other states.
        The dimuonium scatters in an electric field of a nucleus with charge Z.

        Keyword arguments:
        iq1 (int) -- Quantization axis: value 0 - along transfered momentum; value 1 - along beam direction
        iq2 (int) -- Computation method: 1, 2, 3 
        Z (float) -- Nucleus electric charge
        n   (int) -- Dimuonium initial state main quantum number: n 
        l   (int) -- Dimuonium initial state angular quantum number: l
        m   (int) -- Dimuonium initial state angular momentum projection: m
        '''

        ###CHECKING ALLOWED RANGE
        error = True
        if (0 <= l <= n - 1) and (abs(m) <= l):  
            error = False

        ###CREATING THE PARAMETERS 
        iq1 = c_int(iq1)
        iq2 = c_int(iq2)

        Z = c_double(Z)
        n = c_int(n)
        l = c_int(l)
        m = c_int(m)
      
        ###PASSING THE PARAMETERS BY REFERENCE 
        args = [iq1, iq2, Z, n, l, m]
        refs = [ctypes.byref(arg) for arg in args]

        ###CALLING FORTRAN CODE
        if not error:
            result = self.libc.wtcrs_(*refs)
        if error:
            result = 0
        return result




    def calcTransportCrossSection(self, iq1, iq2, Z, n1, l1, m1, n2, l2, m2):
            '''Calculate the  Cross-Section of dimuonium transitions from state (n1, l1, m1) to (n2, l2, m2).
            The dimuonium scatters in an electric field of a nucleus with charge Z.

            Keyword arguments:
            iq1 (int) -- Quantization axis: value 0 - along transfered momentum; value 1 - along beam direction
            iq2 (int) -- Computation method: 1, 2, 3 
            Z (float) -- Nucleus electric charge
            n   (int) -- Dimuonium initial state main quantum number: n 
            l   (int) -- Dimuonium initial state angular quantum number: l
            m   (int) -- Dimuonium initial state angular momentum projection: m
            '''

            ###CHECKING ALLOWED RANGE
            error = True
            if (0 <= l1 <= n1 - 1) and (abs(m1) <= l1):  
                error = False

            if (0 <= l2 <= n2 - 1) and (abs(m2) <= l2):  
                error = False

            ###CREATING THE PARAMETERS 
            iq1 = c_int(iq1)
            iq2 = c_int(iq2)

            Z = c_double(Z)

            n1 = c_int(n1)
            l1 = c_int(l1)
            m1 = c_int(m1)

            n2 = c_int(n2)
            l2 = c_int(l2)
            m2 = c_int(m2)
        
            ###PASSING THE PARAMETERS BY REFERENCE 
            args = [iq1, iq2, Z, n1, l1, m1, n2, l2, m2]
            refs = [ctypes.byref(arg) for arg in args]

            ###CALLING FORTRAN CODE
            if not error:
                result = self.libc.wtrcrs_(*refs)
            if error:
                result = 0
            return result