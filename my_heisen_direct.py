from __future__ import division
import h5py
import numpy as np


# flip between two bit values (i,j) in a . ^=XOR
def flip(a,i,j):
	
	return a^(2**i+2**j)

#Test i th bit of a. Returns 0,1
def btest(a,i):
	
	return (a&(1<<i))>>i #1<<i=2**i BiAND and then right shift ith position


# construct Hamiltonian 
def Hamiltonian(delta):
	Jx=1
	Jz=1
	i=j=np.int32(0)
	a=b=i
	H=np.zeros((hdim,hdim))
	
	for a in xrange(0,hdim):		
		for i in xrange(0,N):
			j=(i+1)%N	#periodic boundary condition
			if (btest(a,i)==btest(a,j)):				
				H[a,a]=H[a,a] + Jz*delta*0.25
				
			else:				
				H[a,a]=H[a,a]-Jz*delta*0.25
				b=flip(a,i,j)				
				H[a,b]=Jx*0.5
	return H

# construct Hamiltonian in  next neareast neighbour
def Hamiltonian_NNN(J1,J2):	

	i=j=np.int32(0)
	a=b=i
	H=np.zeros((hdim,hdim))
	
	for a in xrange(0,hdim):		
		for i in xrange(0,N):
			j=(i+1)%N	#periodic boundary condition
			if (btest(a,i)==btest(a,j)):				
				H[a,a]=H[a,a] + J1*0.25
				
			else:				
				H[a,a]=H[a,a]-J1*0.25
				b=flip(a,i,j)				
				H[a,b]=J1*0.5
	
	for a in xrange(0,hdim):		
		for i in xrange(0,N):
			j=(i+2)%N	#periodic boundary condition
			if (btest(a,i)==btest(a,j)):				
				H[a,a]=H[a,a] + J2*0.25
				
			else:				
				H[a,a]=H[a,a]-J2*0.25
				b=flip(a,i,j)				
				H[a,b]=J2*0.5
	return H
				
# Solve for eigenvalues and eigenvectors
def solve(H):
	
	eigenvalues,eigenvectors=np.linalg.eigh(H)
	return eigenvalues,eigenvectors


# store
def store(eigenvalues,eigenvectors,J1,J2,N):

		datafile = h5py.File('NNN_model_'+str(N)+'_J1_'+str(J1)+'_J2_'+str(J2)+'_data.h5', 'a')
		datafile['eigenvalues'] = eigenvalues
		datafile['eigenvectors'] = (eigenvectors)
		datafile.close()
		print "Successfully saved eigenvalues and eigenvectors"


N=8  # no_of_sites
hdim=2**N
J1 = 1.0
J2 = 0.0
H = Hamiltonian_NNN(J1,J2)
evals,evec = solve(H)
store(evals,evec,J1,J2,N)
