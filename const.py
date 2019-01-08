#! /usr/bin/python
import numpy as np
um = np.float(1e-6);
mm = np.float(1e-3)
cm = np.float(1e-2);
fs = np.float(1e-15);
ns = np.float(1e-9);
nm = np.float(1e-9);
km = np.float(1e3);

pi = 3.1415926;
me = 9.1e-31;
lambda_c = 2.42e-12
c = np.float(3e8);
qe = 1.6e-19;
epsilon0 = 8.854e-12;
mu0 = 4.0*pi*1e-7; '''unit N/A^2'''
h_planck = 6.626069e-34;   '''unit J \cdot s''' 
hbar = h_planck/2.0/pi;
alpha_f = qe**2.0/hbar/c/4.0/pi/epsilon0; '''\sim 1/137 the fine structure constant'''
r_e = qe**2.0/4.0/pi/epsilon0/me/c**2;
Es = alpha_f*qe/4.0/pi/epsilon0/r_e**2.0;
Bs = Es/c
Mev = 1e6*qe;
ev = qe;
kb = 1.3806488e-23

def Get_omega(lambd=1*um):
	return 2*pi*c/lambd

def Get_T0(lambd=1*um):
	return(lambd/c)

def Get_norme(omega,a0 = 1):
	return a0*me*omega*c/qe;

def Get_nc(omega):
	return omega**2*me*epsilon0/qe**2;

#default setting 
#linestyle 
MYLinecolor = ['r','g','b','c','y','k']
MYLinestyle = ['o-','1-','v-','1-','s-','p-']
MYColormap =  ['jet','viridis','RdBu','bwr','seismic']
import sdf
import sdfread as sr
import drawfig as df
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import imp
import os 
import sys
