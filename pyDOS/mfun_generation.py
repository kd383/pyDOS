#!usr/bin/env python
"""
This module is a colletion of functions that generates matrix functions

8/1/2016 - Created
8/16/2016 - Commented
"""

import numpy as np
import scipy.sparse as ss

def mfunc_laplacian(W):
	"""
	Convert a weighted adjacency matrix function into a Laplacian operator
	
	Args:
		W: weighted adjacency matrix/operator

	Output:
		Lfun: graph Laplacian operator
	"""

	# Create a function handle if given a matrix
	if callable(W):
		Wfun = W
	else:
		if isinstance(W,np.ndarray):
			W = ss.csr_matrix(W)
		n = W.shape[0]
		Wfun = lambda x: W*x

	e = np.ones(n)
	d = ss.csr_matrix(Wfun(e)).transpose()

	# L= D-A
	Lfun = lambda x: d.multiply(x)-Wfun(x)
	return Lfun

def mfunc_modularity(W):
	"""
	Convert a weighted adjacency matrix function into a modularity operator
	
	Args:
		W: weighted adjacency matrix/operator

	Output:
		Bfun: graph modularity operator
	"""

	# Create a function handle if given a matrix
	if callable(W):
		Wfun = W
	else:
		if isinstance(W,np.ndarray):
			W = ss.csr_matrix(W)
		n = W.shape[0]
		Wfun = lambda x: W*x

	e = np.ones(n)
	d = ss.csr_matrix(Wfun(e))
	s = d.sum()

	# B = A-d*d'/s
	Bfun = lambda x: Wfun(x)-(d.transpose()/s)*(d*x)
	return Bfun

def mfunc_normalize(W,mode='s'):
	"""
	Convert a weighted adjacency matrix function into a normalized adjacency 
	operator

	Args:
		W: weighted adjacency matrix/operator
		mode: string indicating the style of normalization;
			's': Symmetric scaling by the degree (default)
			'r': Normalize to row-stochastic
			'c': Normalize to col-stochastic

	Output:
		Nfun: normalized adjacency operator
	"""

	# Create a function handle if given a matrix
	if callable(W):
		Wfun = W
	else:
		if isinstance(W,np.ndarray):
			W = ss.csr_matrix(W)
		n = W.shape[0]
		Wfun = lambda x: W*x

	# Calculate the degree of nodes
	e = np.ones(n)
	d = ss.csr_matrix(Wfun(e)).transpose()

	# Normalize by the desire style
	if mode =='s':
		d.data = 1/np.sqrt(d.data)
		Nfun = lambda x: d.multiply(Wfun(d.multiply(x)))
	elif mode == 'r':
		d.data = 1/d.data
		Nfun = lambda x: d.multiply(Wfun(x))
	elif mode == 'c':
		d.data = 1/d.data
		Nfun = lambda x: Wfun(d.multiply(x))
	else:
		raise ValueError('Unknown mode!')

	return Nfun

def mfunc_slaplacian(W):
	"""
	Convert a weighted adjacency matrix function into a signless Laplacian 
	operator
	
	Args:
		W: weighted adjacency matrix/operator

	Output:
		Lfun: graph signless Laplacian operator
	"""

	# Create a function handle if given a matrix
	if callable(W):
		Wfun = W
	else:
		if isinstance(W,np.ndarray):
			W=ss.csr_matrix(W)
		n = W.shape[0]
		Wfun = lambda x: W*x

	e = np.ones(n)
	d = ss.csr_matrix(Wfun(e)).transpose()

	# L = D+A
	Lfun = lambda x: d.multiply(x)+Wfun(x)
	return Lfun

def mfunc_cheb_poly(coeff,W):
	"""
	Construct the operator to evaluate the Chebyshev series for given matrix W
	and a vector of coefficients.

	Args:
		coeff: coefficients of Chebyshev polynomials
		W: The matrix or function input for Chebyshev polynomials

	Output:
		pWfun: Chebyshev series evaluation function
	"""

	# Create a function handle if given a matrix
	if callable(W):
		Wfun = W
	else:
		if isinstance(W,np.ndarray):
			W=ss.csr_matrix(W)
		Wfun = lambda x: W*x

	# Apply coefficients of Chebyshev series
	coeff = np.asarray(coeff)
	if coeff.shape == () or len(coeff) == 1:
		pWfun = lambda x: coeff*Wfun(x)
	else:
		pWfun = lambda x: mfunc_cheb_poly_apply(coeff,Wfun,x)

	return pWfun
	

def mfunc_cheb_poly_apply(coeff,Wfun,x):
	"""
	Auxiliary function for mfunc_cheb_poly_apply
	See description above
	"""

	Pm = x
	P0 = Wfun(x)
	y = coeff[0]*Pm

	# Apply recurrecne relation of Chebyshev polynomials
	for j in range(1,len(coeff)-1):
		y = coeff[j]*P0
		Pn = 2*Wfun(P0)-Pm
		Pm = P0
		P0 = Pn

	y = y+coeff[-1]*P0
	return y
	