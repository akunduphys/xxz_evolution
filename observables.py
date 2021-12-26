from __future__ import division
import h5py
import numpy as np

def loader(N,J1,J2):
		datafile = h5py.File('NNN_model_'+str(N)+'_J1_'+str(J1)+'_J2_'+str(J2)+'_data.h5', 'a')
		print "Loading..."+str(N)		
		#Anisotropy_12_0.5_1.5_201data
		eigenvalues = datafile['eigenvalues'][:]
		eigenvectors = datafile['eigenvectors'][:]
		eigenvectors = eigenvectors.T
		datafile.flush()
		datafile.close()
		return eigenvalues,eigenvectors

# flip between two bit values (i,j) in a . ^=XOR
def flip(a,i,j):
	
	return a^(2**i+2**j)

#Test i th bit of a. Returns 0,1
def btest(a,i):
	
	return (a&(1<<i))>>i #1<<i=2**i BiAND and then right shift ith position


#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

def Splus():
	return np.array([[ 0.,  1.],[ 0.,  0.]])

def Sminus():
	return np.array([[ 0.,  0.],[ 1.,  0.]])
	

def MatrixSp(i_arr,hdim):
	"""
	#returns the S+ matrix in direct sum basis, the S+ is at ith site array
	#i_arr = [1,N/2,N] #numbering from 1-N
	"""
	i_arr = np.asarray(i_arr)
	i_arr-=1
	Splus_arr = np.eye(hdim,hdim)
	for i in i_arr:
		Splus_arr[2*i:2*i+2,2*i:2*i+2] = Splus()
		
	return Splus_arr	
	
def MatrixSm(i_arr,hdim):
	"""
	#returns the S- matrix in direct sum basis, the S- is at ith site array
	#i_arr = [1,N/2,N] #numbering from 1-N
	"""
	i_arr = np.asarray(i_arr)
	i_arr-=1
	Sminus_arr = np.eye(hdim,hdim)
	for i in i_arr:
		Sminus_arr[2*i:2*i+2,2*i:2*i+2] = Sminus()
		
	return Sminus_arr

#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------


def opSp(i,state):
	"""
	Sp(i):
	Flips ith position of state if down. returns 0 if already up
	"""
	if(btest(state,i)==1):return 0 #null vector
	return state^(2**i)


def ElementSpBinaryBasis(i,m,n):
	"""
	returns <Xm|Sp(i)|Xn> of operator
	Xm and Xn are (0,1) basis
	"""
	#if (opSp(i,n) == m):
	#	return 1
	if(btest(n,i) == 0):
		flip = n^(2**i)
		if(flip == m):
			return 1
		
	return 0
	
def MatrixSpConfigBasis(i,hdim):	
	"""
	returns the Matrix Sp in Configuration Basis
	i<=N-1
	"""
	r_vec = np.arange(0,hdim)
	sp_matrix = np.zeros((hdim,hdim))
	for m in r_vec:
		for n in r_vec:
			sp_matrix[m,n] = ElementSpBinaryBasis(i,m,n)
	return sp_matrix


def MatrixSpEnergyBasis(i,hdim,evals,evec):	
	"""
	returns the Matrix Sp in Energy Basis
	i<=N-1
	Sp_energy = U.Sp.U+
	"""
	
	r_vec = np.arange(0,hdim)
	sp_matrix_bb = MatrixSpConfigBasis(i,hdim)#MatrixSp([i],hdim)#
	sp_matrix_eb = np.zeros((hdim,hdim))	
	#for m in r_vec:
	#	for n in r_vec:			
	#		sp_matrix_eb[m,n] = np.dot(evec[m], np.dot(sp_matrix_bb,evec[n].T))	
	sp_matrix_eb = np.dot(evec,np.dot(sp_matrix_bb,evec.T))
	return sp_matrix_eb	
	
def opSm(i,state):
	"""
	Sm(i):
	Flips ith position of state if up. returns 0 if already down
	"""
	if(btest(state,i)==0):return 0 #null vector
	return state^(2**i)


def opSz(i,state):
	"""
	test if the given ith bit is up or down and returns 1,0 respectively
	"""
	if(btest(state,i)==0):return -1 
	else: return 1
	
def ElementSmBinaryBasis(i,m,n):
	"""
	returns <Xm|Sm(i)|Xn> of operator
	Xm and Xn are (0,1) basis
	"""
#	if (opSm(i,n) == m):
#		return 1
	if(btest(n,i) == 1):
		flip = n^(2**i)
		if(flip == m):
			return 1
		
	return 0

def ElementSzBinaryBasis(i,m,n):
	"""
	returns <Xm|Sz(i)|Xn> of operator
	Xm and Xn are (0,1) basis
	"""
#	if (opSm(i,n) == m):
#		return 1
	if(m==n):
		if(btest(n,i) == 1):		
				return 1
		else:
			return -1
	return 0
	
def MatrixSzConfigBasis(i,hdim):	
	"""
	returns the Matrix Sz in Configuration Basis
	"""
	r_vec = np.arange(0,hdim)
	sz_matrix = np.zeros((hdim,hdim))
	
	for m in r_vec:
		for n in r_vec:
			sz_matrix[m,n] = ElementSzBinaryBasis(i,m,n)
	return sz_matrix


	
def MatrixSmConfigBasis(i,hdim):	
	"""
	returns the Matrix Sp in Configuration Basis
	"""
	r_vec = np.arange(0,hdim)
	sm_matrix = np.zeros((hdim,hdim))
	for m in r_vec:
		for n in r_vec:
			sm_matrix[m,n] = ElementSmBinaryBasis(i,m,n)
	return sm_matrix


def MatrixSmEnergyBasis(i,hdim,evals,evec):	
	r_vec = np.arange(0,hdim)
	sm_matrix_bb =MatrixSmConfigBasis(i,hdim)# MatrixSm([i],hdim)#
	sm_matrix_eb = np.zeros((hdim,hdim))
	#evec = evec.T
	#for m in r_vec:
	#	for n in r_vec:			
	#		sm_matrix_eb[m,n] = np.dot(evec[m], np.dot(sm_matrix_bb,evec[n].T))	
	sm_matrix_eb = np.dot(evec,np.dot(sm_matrix_bb,evec.T))
	return sm_matrix_eb	

def MatrixSzEnergyBasis(i,hdim,evals,evec):	
	r_vec = np.arange(0,hdim)
	sz_matrix_bb =MatrixSzConfigBasis(i,hdim)# MatrixSm([i],hdim)#
	sz_matrix_eb = np.zeros((hdim,hdim))
	#evec = evec.T
	#for m in r_vec:
	#	for n in r_vec:			
	#		sm_matrix_eb[m,n] = np.dot(evec[m], np.dot(sm_matrix_bb,evec[n].T))	
	sz_matrix_eb = np.dot(evec,np.dot(sz_matrix_bb,evec.T))
	return sz_matrix_eb	
	

def InProd(a,b):
	return np.dot(a.T,b)	

