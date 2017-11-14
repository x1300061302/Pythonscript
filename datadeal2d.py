import sdfread as sr
import numpy as np
import drawfig as df
import sdf
import os
from const import *


def print_spectrum(filename,prefix =''): 
	a = sdf.read(filename)
	dgam = sr.Get_particle_variable(a,'Gamma','electron');
	gam = dgam.data;
	df.draw_spectrum(gam,label =('x','y','gam'+prefix),figname = 'gam'+prefix)
	print('ok')

def draw_field(filename,varname,prefix=''):
	a = sdf.read(filename)
	extent =np.array(sr.Get_extent(a))/um;
	
	var = sr.Get_field_variable(a,varname);
	df.draw_field_snapshot(var.data.T,\
			extent=extent,\
			xylim=([0,100],[-20,20]),\
			label=('x','y',varname),\
			figname=varname+prefix)

def print_compared_spectrum(dirnames,filename):
	gam = []
	for dirname in dirnames:
		a = sdf.read(dirname+'/'+filename)
		dgam = sr.Get_particle_variable(a,'Gamma','electron');
		gam.append(dgam.data)

	prefix = filename[3:5]
	df.draw_spectrum_nline(gam,dataname = dirnames,label = ('$\gamma$','N','$\gamma-N$'),weight=0,figname = 'gam'+prefix,numl = len(dirnames))

#####################
	

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


###############main procedure 

#mkdir 
#try:
#	os.mkdir('figure')
#except:
#	print("figure has been existed")
#print_field()
#print_particle()
print_particle_spectrum_compare();
print('Completed')
#os.system('mv *.png figure/')

