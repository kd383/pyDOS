#!/usr/bin/env python
"""
This is a test function for eig_rand1 module

11/29/2016 - Created
11/29/2016 - Change path import

"""

import sys
from os.path import dirname,realpath
sys.path.insert(0,dirname(realpath(__file__))[:-11])
from pyDOS import *
import scipy.io as sio

def test_moment_comp(A,n,Z=0,x=0):
	if isinstance(Z,int):
		Z = np.sign(nr.randn(n,20))
	if isinstance(x,int):
		x = nr.randn(10,1)

	L = matrix_normalize(A)
	Lfun = mfunc_normalize(A)

	[c, cs] = moments_cheb_dos(L,5534,Z,20)
	[cl, csl] = moments_cheb_ldos(L,5534,Z,20)
	[cfun, csfun]= moments_cheb_dos(Lfun,5534,Z,20)
	[cfunl, csfunl] = moments_cheb_ldos(Lfun,5534,Z,20)
	cd = moments_delta(x,20)
	sio.savemat('test_mmt_cmp.mat',{'Z':Z,'c':c,'cs':cs,'cl':cl,'csl':csl,'cfun':cfun,
							'csfun':csfun,'cfunl':cfunl,'csfunl':csfunl,'x':x,'cd':cd})

if __name__ == '__main__':
	A = load_graph('erdos02-cc','../data/')
	test_moment_comp(A,A.shape[0])
