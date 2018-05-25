#! /usr/bin/python
import sdfread as sr
import drawfig as df
import numpy as np
from const import*

def draw_specturm(allfilename): 
	for filename in allfilename:
		dgam = sr.Get_particle_variable(filename,('Gamma'),'electron');
		gam = dgam.data;
		px = dpx.data
		aad = np.array((ad1,ad2));
		name = ('gamma','px')
		import drawfig as df
		df.draw_spectrum(aad,name,numl=2)

	print('ok')
	
def angle_distribution(filename):
	"""draw angle distribution of filename"""
	dpy,dpz,dgam = sr.Get_particle_variable(filename,['Py','Pz','Gamma'],'electron');
	py = dpy.data/me/c;
	pz = dpz.data/me/c;
	gam = dgam.data;
	theta = np.arctan2(py,pz)
	import drawfig as df
	df.draw_angle_distribution3d(theta,gamma);

def draw_field(filename,fieldname):
	grid = sr.Get_field_variable(filename,'Grid_mid');
	extent = grid.extents
	a1 = sr.Get_field_variable(filename,fieldname);
	df.draw_field_snapshot(a1.data,extent=extent,plain='xy',index=96,label=('x','y','z','Ex06'),figname='Ex06',Display=0)

def cal_circular_frequency(filename):
	[dpy,dpz,dgam,dgrid] = sr.Get_particle_variable(filename,['Py','Pz','Gamma','Grid'],'electron');
	py = dpy.data;
	pz = dpz.data;
	gam = dgam.data;
	y = dgrid.data[1]
	z = dgrid.data[2]
	r = y**2+z**2
	omega = (py*z + pz*y)/r**2/gam/me
	return omega 

def cal_laser_frequency(filename):
	'''laser info lambda = 1um omega = 2*np.pi/T0,vph=c'''
	[dpx,dgam]= sr.Get_particle_variable(filename,['Px','Gamma'],'electron');
	px = dpx.data;
	gam = dgam.data
	lambda_l = 1*um
	omega_l0 = 2*np.pi*c/lambda_l;
	vx = px/gam/me;
	omega_le = omega_l0*(1 - vx/c)
	return omega_le 



'''main'''



'''main'''

