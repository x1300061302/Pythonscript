#! /usr/bin/python
um = 1e-6;
pi = 3.1415926;
me = 9.1e-31;
c = 3e8;
qe = 1.6e-19;
epsilon0 = 8.854e-12;
mu0 = 4.0*pi*1e-7; '''unit N/A^2'''
h_planck = 6.626069e-34;   '''unit J \cdot s''' 
hbar = h_planck/2.0/pi;
alpha = qe**2.0/hbar/c/4.0/pi/epsilon0; '''\sim 1/137 the fine structure constant'''
r_e = qe**2.0/4.0/pi/epsilon0/me/c**2;
E_s = alpha*qe/4.0/pi/epsilon0/r_e**2.0;
Mev = 1e6*qe;

def omega(lambda=1*um):
	return 2*pi*c/lambda

def T0(lambda=1*um):
	return(lambda/c)

def norme(a0=1,omega):
	return a0*me*omega*c/qe;

def nc(omega):
	return omega**2*me*epsilon0/qe**2;

