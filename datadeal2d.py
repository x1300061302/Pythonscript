import sdfread as sr
import numpy as np
import drawfig as df
import sdf
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from const import *

###########Part A: is common functions which need some specifical parameter 
''' in this part functions should use sdfread to get data and combine them
	and then draw the figure with drawfig.py'''
def print_spectrum(filename,prefix ='',varname = 'Gamma',species='electron'): 
	'''a function for print spectrum of electron 
	in filename like 'pxxxx.sdf' and with prefix like 16T0
	in fact one can choose 
	varname as px,py,pz to draw the histogram, or different species 
	like species = 'electron' 'proton' '''
	a = sdf.read(filename)
	dgam = sr.Get_particle_variable(a,varname,species);
	gam = dgam.data;
	df.draw_spectrum(gam,label =('x','y',varname+prefix),figname = varname+prefix)
	print('ok')

def draw_field(filename,varname,prefix=''):
	'''filename = 'fxxxx.sdf' 
	varname = Ex,Ey or averaged and EkBar,Density
	prefix is like 16T0'''
	
	a = sdf.read(filename)
	extent =np.array(sr.Get_extent(a))/um;
	
	var = sr.Get_field_variable(a,varname);
	df.draw_field_snapshot(var.data.T,\
			extent=extent,\
			xylim=([0,100],[-20,20]),\
			label=('x','y',varname),\
			figname=varname+prefix)
def print_compared_spectrum(dirnames,filename):
	'''this function aims to compare different spectrum in different cases
	dirnames is a list or tuple of the names of directory
	filename is the pxxxx.sdf in these directory with out prefix'''
	gam = []
	for dirname in dirnames:
		a = sdf.read(dirname+'/'+filename)
		dgam = sr.Get_particle_variable(a,'Gamma','electron');
		gam.append(dgam.data)

	prefix = filename[3:5]
	df.draw_spectrum_nline(gam,dataname = dirnames,label = ('$\gamma$','N','$\gamma-N$'),weight=0,figname = 'gam'+prefix,numl = len(dirnames))

#####################
	
#####################Part B functions combination to set parameter 
	'''in this part one can set some specifical parameter through 
	the functions in Part B for those functions in A ,in generally in part B,
	you just need some operators for files'''
def print_field():
	'''print field_variable
	   default is Ex_averaged,Ey,EkBar_electron,Density_electorn				'''
	for files in sr.Get_file('f'):
		try:
			draw_field(files,'Ex_averaged',prefix=files[2:5]);
		except:
			print('Wrong in',files,'with variable Ex_averaged')
		try:
			draw_field(files,'Ey',prefix=files[2:5]);
		except:
			print('Wrong in',files,'with variable Ey')
		try:
			draw_field(files,'EkBar_electron',prefix=files[2:5]);
		except:
			print('Wrong in',files,'with variable EkBar')
		try:
			draw_field(files,'Density_electron',prefix=files[2:5]);
		except:
			print('Wrong in',files,'with variable Density')

def print_particle_spectrum():
	'''just print the electron's spectrum'''
	for files in sr.Get_file('p'):
		print(files)
		print_spectrum(files,prefix = files[2:5]);
	#	print('Wrong in',files,'with gam')

def print_particle_spectrum_compare():

	dirnames = ('alpha05','alpha1','alpha0')
	ff = sr.Get_file('p',dirnames[0]);
	print(ff)
	for files in sr.Get_file('p',dirnames[0]):
		print(files)
		print_compared_spectrum(dirnames,files)
	
	#	print('Wrong in',files,'with gam')

###############Part C: some special operate
''' in Part C some special operators for draw typical figures to anaylize the physics
	You MUST intepret the functions in this part to  make the function clear'''
def draw_partx_ex():
	'''this function aims to draw particles' position 'scatter figure'and the 1-D field 
	,for example, in order to analyze if particle is accelerated by ex, I'll draw the 
	higher energy particles in x with ex distribution'''
	import sdf 
	pd = sr.Get_file(prefix = 'p')
	fd = sr.Get_file(prefix = 'f')
	print(pd)
	print(fd)
	for i in range(0,len(pd)):
		#high energy electron x-pos
		a = sdf.read(pd[i]) 
		dgrid = sr.Get_particle_variable(a,'Grid','electron')
		grid = dgrid.data
		x = grid[0]
		print(type(x))
		dgam = sr.Get_particle_variable(a,'Gamma','electron')
		max_gam = np.max(dgam.data);
		print(max_gam,'choose',max_gam*2/3)
		gam = dgam.data[dgam.data > max_gam*2/3];
		xx = x[dgam.data > max_gam*2/3]/um;
		######## field ex 
		a = sdf.read(fd[i])
		exx = sr.Get_field_variable(a,'Ex_averaged')
		print(exx)
		ex = exx.data
		nx,ny = ex.shape #3-D test should be deleted
		print(nx,ny)
		extent = sr.Get_extent(a) #2-D
		ex = ex[:,ny/2];
		axx = np.linspace(extent[0],extent[1],nx)/um
		df.plot_line(axx,ex/np.max(ex),('x','ex','ex-x'),'ex',savefig = 0)
		plt.scatter(xx,gam/np.max(gam),c = 'b',norm = 1,edgecolors = 'none')
		plt.savefig('scatxex',dpi=300,facecolor='none',edgecolor='b')
		




###############main procedure 
draw_partx_ex()
