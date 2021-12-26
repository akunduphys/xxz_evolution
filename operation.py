"""
1.Find a hamiltonian such that the <a|A|a> - <a+1|A|a+1> is smooth fucntion of ordered eigenvalues
but the off diagoonal elements are comareabel to diagonal elements. See if ETH holds

2. Check if autocorrelation commutator satisfies ETH
"""

from __future__ import division
import numpy as np
from observables import *
import matplotlib.pyplot as plt

	
def CommutatorAvg2(i,t,a,evals,evec):
	"""
	returns avg of commutator of <a|[S+_i(t),S-i_(0)]|a>
	"""

	#r_vec = np.arange(0,hdim)
	sp_eb = MatrixSpEnergyBasis(i,hdim,evals,evec)
	sm_eb = MatrixSmEnergyBasis(i,hdim,evals,evec)
	
	comm_avg = 0 + 0j
	for k in xrange(hdim):
		for j in xrange(hdim):
			for n in xrange(hdim):
				temp = np.exp((evals[n]-evals[k])*t*1j)*sp_eb[n,k]*InProd(a,evec[j])
				temp -= np.exp((evals[j]-evals[n])*t*1j)*sp_eb[j,n]*InProd(a,evec[k])
				comm_avg += sm_eb[k,j]*InProd(a,evec[n])*temp
	return comm_avg
	
	
def CommutatorAvg(i,t,psi,evals,evec):
	"""
	|psi> = (0.5,1,0.....2^N)
	|a> = \sum_m psi_m*|m>
	is a vector giving coeffecients of the corrosponding energy basis
	returns avg of commutator of <|[S+_i(t),S-i_(0)]|a>
	"""
	#r_vec = np.arange(0,hdim)
	sp_eb = MatrixSpEnergyBasis(i,hdim,evals,evec)
	sm_eb = MatrixSmEnergyBasis(i,hdim,evals,evec)
	comm_avg = 0 + 0j
	for m in xrange(hdim):
		for n in xrange(hdim):
			for k in xrange(hdim):
				temp = np.exp((evals[m]-evals[n])*t*1j)*sp_eb[m,n]*psi[k]				
				comm_avg += (psi[m]*sm_eb[n,k] - psi[n]*sm_eb[k,m])*temp
				
	return comm_avg


def CommutatorMatrix(i,t,evals,evec):
	"""
	returns a matrix of [U Sp Ud , Sm] in energy basis
	"""
	sp = MatrixSpEnergyBasis(1,hdim,evals,evec)
	sm = MatrixSmEnergyBasis(1,hdim,evals,evec)
	
	U = np.diag(np.exp(-evals*t*1j))
	Ud = np.diag(np.exp(evals*t*1j))
	
	spt = np.dot(np.dot(U,sp),Ud)
	
	commute = np.dot(spt,sm) - np.dot(sm,spt) #this is in energy basis
	return commute

def CommutatorMatrixSq(i,t,evals,evec):
	"""
	returns a matrix of [U Sp Ud , Sm][U Sp Ud , Sm]d in energy basis
	"""
	sp = MatrixSpEnergyBasis(1,hdim,evals,evec)
	sm = MatrixSmEnergyBasis(1,hdim,evals,evec)
	
	U = np.diag(np.exp(-evals*t*1j))
	Ud = np.diag(np.exp(evals*t*1j))
	
	spt = np.dot(np.dot(U,sp),Ud)
	
	commute = np.dot(spt,sm) - np.dot(sm,spt) #this is in energy basis
	commuted = np.conj(commute).T
	
	commutesq = np.dot(commute,commuted)
	return commutesq
	

def CommutatorAvg3(i,t,psi,evals,evec):
	
	commMatrix = CommutatorMatrix(i,t,evals,evec)
	avgComm = np.dot(psi.T,np.dot(commMatrix,psi))
	return avgComm
	
	
	
	
def ThermalAvgCommutator(i,t,T,evals,evec):
	"""
	compute \sum e^(-beta En) <n|[]|n> /\sum e^(-beta E_n)
	"""
	beta = 1.0/T
	thavg =0
	norm = 0
	commMatrix = CommutatorMatrix(i,t,evals,evec)
	d = commMatrix.diagonal()
	thermal_weights = np.exp(-beta*evals)
	norm = np.sum(thermal_weights)
	thavg  = np.dot(thermal_weights,d)/norm
	
	#for n in xrange(hdim):
		#thavg += np.exp(-beta*evals[n])*commMatrix[n,n]
		#norm += np.exp(-beta*evals[n])	
	#thavg = thavg/norm
	return thavg	
	
	
def ThermalAvgCommutatorSq(i,t,T,evals,evec):
	"""
	compute \sum e^(-beta En) <n|[]|n> /\sum e^(-beta E_n)
	"""
	beta = 1.0/T
	thavg =0
	norm = 0
	commMatrix = CommutatorMatrixSq(i,t,evals,evec)
	
	d = commMatrix.diagonal()
	thermal_weights = np.exp(-beta*evals)
	norm = np.sum(thermal_weights)
	thavg  = np.dot(thermal_weights,d)/norm
	
	
	#for n in xrange(hdim):
		#thavg += np.exp(-beta*evals[n])*commMatrix[n,n]
		#norm += np.exp(-beta*evals[n])
	
	#thavg = thavg/norm
	return thavg	
	
	

		
N = 6
J1 = 1.0
J2 = 0.0
T  = 10
evals_arr,evec_arr = loader(N,J1,J2)
hdim = 2**N
#check_commutator()

time = np.arange(0,500,0.5)
comm_avg_ts = time*(1.0+0j)
comm_avg_ts2 = time*(1.0+0j)

#a is list of strenghts of different energy eigenbasis
#a = np.random.random(hdim)
#a = np.ones(hdim)
#a = init_thermal_state(hdim,1,evals_arr)
#a = evec_arr[0]#np.ones(hdim)
#print a

#statevec = np.dot(evec,a)
#testCommAvgT0(1,a,evals_arr,evec_arr)
#check_commutator_t0(1,a,evals_arr*1,evec_arr*1)
#check_commutator(1,evals_arr*1,evec_arr*1)



for idx,t in enumerate(time):
	print t
	
	comm_avg_ts2[idx]=ThermalAvgCommutatorSq(1,t,T,evals_arr,evec_arr) #for sums of energy states
	#comm_avg_ts[idx]=ThermalAvgCommutator(1,t,T,evals_arr,evec_arr) #for sums of energy states
	
	#comm_avg_ts[idx]=CommutatorAvg3(1,t,a,evals_arr,evec_arr) #for sums of energy states
	#comm_avg_ts2[idx]=CommutatorAvg(1,t,a,evals_arr,evec_arr)#for direct states

plt.plot(time,np.real(comm_avg_ts2))
#plt.plot(time,np.abs(comm_avg_ts))
plt.show()

	
def init_thermal_state(hdim,T,evals):
	beta = 1/T
	a = np.zeros(hdim)
	
	for i in xrange(hdim):
		a[i]  = np.exp(-beta*evals[i])
	
	norm = np.sum(a)
	a = a/norm
	return a	
	

def check_commutator_t0(i,a,evals,evec):
	"""
	check if time zero commutator works
	"""
	#commutator of [S+,S-] = 2Sz
#	print evals,evec
	state = np.zeros(hdim)
	for j in xrange(hdim):
		state += a[j]*evec_arr[j]
	#print "state",state
	
	#cm1 = CommutatorAvg2(i,0,state,evals,evec) # in energy basis
	cm1 = CommutatorAvg3(i,0,a,evals,evec) # in energy basis
	cm2 = CommutatorAvg(i,0,a,evals,evec) # in energy basis	
		
	sz_eb = MatrixSzEnergyBasis(i,hdim,evals,evec)
	
	#szexp = np.dot(state,np.dot(sz_eb,state.T))
	szexp = np.dot(a,np.dot(sz_eb,a))
	
	#sz_trans = np.dot(evec,np.dot(cm1,evec.T))
	print "commutator f:state\n ",cm1,"\n f:a \n",cm2
	
	print "sz exp\n ",szexp



	
def testCommAvgT0(i,a,evals,evec):
	"""
	
	"""
	state = np.zeros(hdim)
	for j in xrange(hdim):
		state += a[j]*evec_arr[j]

	sp_eb = MatrixSpEnergyBasis(i,hdim,evals,evec)
	sm_eb = MatrixSmEnergyBasis(i,hdim,evals,evec)
	commute = np.dot(sp_eb,sm_eb) - np.dot(sm_eb,sp_eb) #this is in energy basis
	
	avgc = np.dot(a,np.dot(commute,a))
	
	sz_eb = MatrixSzEnergyBasis(i,hdim,evals,evec)	
	szexp = np.dot(a,np.dot(sz_eb,a))
	
	print avgc,"\n",szexp
	
	

def check_commutator(i,evals,evec):
	"""
	check if commutator [s+,s-] = sz in energy badid.
	trasnform sz to original basis by Ud sz U = sz_con
	"""
	
	sp = MatrixSpEnergyBasis(1,hdim,evals,evec)
	sm = MatrixSmEnergyBasis(1,hdim,evals,evec)
	commute = np.dot(sp,sm) - np.dot(sm,sp) #this is in energy basis
	print "sm\n",sm
	print "Energy basis [s+,s-]=\n ",commute
	#transfrom back
	sz_eb = MatrixSzEnergyBasis(i,hdim,evals,evec)
	print "Matrix Sz in energy basis constructed \n",sz_eb
	
	# transform back is U+ S U
	sz_trans = np.dot(evec.T,np.dot(commute,evec))
	print "transformed back to config basis sz\n",sz_trans		
		

