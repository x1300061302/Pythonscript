import sdfread as sr
import numpy as np
import drawfig as df
import sdf
import os
from const import *


def draw_specturm(filename,prefix =''): 
	a = sdf.read(filename)
	dgam = sr.Get_particle_variable(a,('Gamma'),'electron');
	print(dgam)
	gam = dgam.data;
	name = ('gamma')
	df.draw_spectrum(gam,name,figname = 'gam'+prefix,numl=1)
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

def print_particle():
	for files in sr.Get_file('p'):
		try:
			draw_spectrum(files,prefix = files[2:5]);
		except:
			print('Wrong in',files,'with gam')



###############main procedure 

#mkdir 
try:
	os.mkdir('figure')
except:
	print("figure has been existed")
#print_field()

print_particle()

try:
	os.system('mv *.png figure/')
except:
	print("Wrong mv")

