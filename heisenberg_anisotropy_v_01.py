# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 21:09:14 2013

@author: aritrak
"""

from __future__ import division
import h5py
import numpy as np

N=12  # no_of_sites
hdim=2**N

# flip between two bit values (i,j) in a . ^=XOR
def flip(a,i,j):
	
	return a^(2**i+2**j)

#Test i th bit of a. Returns 0,1
def btest(a,i):
	
	return (a&(1<<i))>>i #1<<i=2**i BiAND and then right shift ith position

#find out fixed magnetisation basis
def makebasis(nu): #nu total number of up spins; in case: N/2 for m=0

	fm=0
	fma=[]
	for a in xrange(0,hdim):
		nus=0
		for i in xrange(0,N):
			if(btest(a,i)):nus=nus+1 # 1 is up
		if(nus==nu):
			fm=fm+1
			fma.append(a)
	
	return fm,fma

def findstate(sa,fm,fma): 
# find the position of state sa from the list of fma
# fm is the length of fma the reduced basis set
	bmin,bmax=0,fm-1
	#index=bisect_left(sa,fm)
	while 1:
		b=bmin+(bmax-bmin)//2
		if(sa<fma[b]):
			bmax=b-1
		elif (sa>fma[b]):
			bmin=b+1
		else:return b
	
	
 
	
# construct Hamiltonian in magnetized basis
def Hamiltonian(fm,fma,delta):
	Jx=1
	Jz=1
	i=j=np.int32(0)
	a=b=i
	H=np.zeros((fm,fm))
	
	for a in xrange(0,fm):
		sa=fma[a]
		for i in xrange(0,N):
			j=(i+1)%N	#periodic boundary condition
			if (btest(sa,i)==btest(sa,j)):
				H[a,a]=H[a,a]+Jz*delta*0.25
			else:
				H[a,a]=H[a,a]-Jz*delta*0.25
				sb=flip(sa,i,j)
				b=findstate(sb,fm,fma)
				H[a,b]=Jx*0.5
	return H
				
	

# Solve for eigenvalues and eigenvectors
def solve(H):
	
	eigenvalues,eigenvectors=np.linalg.eig(H)
	return eigenvalues,eigenvectors

# store
def store(eigenvalues,eigenvectors,delta,sites):

		datafile = h5py.File('Anisotropy_'+sites+'data.h5', 'a')
		datafile['eigenvalues_'+delta] = eigenvalues
		datafile['eigenvectors_'+delta] = (eigenvectors)
		datafile.close()
		print "Successfully saved eigenvalues and eigenvectors"

	
fm,fma=makebasis(N/2.)	
#print fm,fma
def anisotropic_solve():
	delarr=[1]
	for delta in delarr:# change spacing from .5 fpr more points and change in operations too
		print "delta: ",delta,"  constructing hamiltonian"
		H=Hamiltonian(fm,fma,delta)
		print " Solving Hamiltonian "
		eva,eve=solve(H)
		
		store(eva,eve,str(delta),str(N))
		print "eigenval/vec: stored ..."
