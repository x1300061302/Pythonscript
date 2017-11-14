#!/usr/bin/python
import sdf 

def Get_particle_variable(sdffile,var,species):
	'''def Get_particle_variable(sdffile,var,species):
	sdffile e.g = 'pxxxx.sdf'
	var = [Px,Py,Gamma,Grid] etc.
	species is defined in input like electron trackele proton 
	also for BlockPointVariable'''
	import re

	rege_var = '';nofind = 1;
	''' only Grid '''
	if var == 'Grid':
		rege_var = r'^'+var + '_Particles' + '\w*'+species+'$' 
	else:
		rege_var = r'^Particles_'+var+'\w*'+species+'$'
	dic = sdffile.__dict__
	for key in dic: 
		if re.match(rege_var,key):
			data = dic[key];
			nofind = 0;
			break;

	if (nofind):
		print('Can Not Find '+var+' of '+species)
		print(dic.keys())
	return data


def Get_field_variable(sdffile,var):
	
	'''def Get_field_variable(sdf.read,var):
	var = 1 should be Ex,Ey,Ex_averaged ...
	or Grid_Grid_mid
	or Density_electron ... 
	or Ekbar_electron...'''
	import re
	rege_var ='';nofind = 1;data = 0;
	
	rege_var = r'^' + '\w*'+var+'$' 

	dic = sdffile.__dict__
	for key in dic: 
		if re.match(rege_var,key):
			data = dic[key];
			nofind = 0;
			break;
	if (nofind):
		print('*******Can Not Find '+var+' in '+'***********')
		print(dic.keys())
		print('******************************************************')
	return data
def Get_extent(sdffile):
	import numpy as np;
	grid = Get_field_variable(sdffile,'Grid');
	index = np.array([0,2,1,3]);
	extent = np.array(grid.extents)
	myextent = extent[index]
	return myextent
def Get_file(prefix):
	import os 
	import re
	filename = [];
	regu_var = r'^'+prefix+'\w+.sdf$';
	cwd = os.getcwd()
	for files in sorted(os.listdir(cwd)):
		if re.match(regu_var,files):
			filename.append(files)
	return filename

