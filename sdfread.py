#!/usr/bin/python
import sdf 
import re
import numpy as np
import os 
PARTVARNAME_RC = ['Px','Py','Pz','Gamma','Grid']
SPECIES_RC = 'electron'
FIELDNAME_RC = ['Ex_averaged','Ey','Bz_averaged','Density_electron','EkBar_electron']

def Get_partvar_npy(varname=PARTVARNAME_RC,species=SPECIES_RC):
	fps = Get_file('p') 
	for i in range(0,len(fps)):
		a = sdf.read(fps[i])
		for var in varname:
			data = Get_particle_variable(a,var,SPECIES_RC)
			data = data.data;
			print(var)
			np.save(var+fps[i][3:5]+'.npy',data)

def Get_fieldvar_npy(varname=FIELDNAME_RC,species=SPECIES_RC):
	ffs = Get_file('f') 
	for i in range(1,len(ffs)):
		a = sdf.read(ffs[i])
		for var in varname:
			data = Get_field_variable(a,var)
			print(var)
			data = data.data;
			np.save(var+str(i)+'.npy',data)

def Get_particle_variable(sdffile,var,species):
	'''def Get_particle_variable(sdffile,var,species):
	sdffile e.g = 'pxxxx.sdf'
	var = [Px,Py,Gamma,Grid] etc.
	species is defined in input like electron trackele proton 
	also for BlockPointVariable'''
	data = -1

	rege_var = '';nofind = 1;
	# only Grid 
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
	return data.data


def Get_field_variable(sdffile,var):
	
	'''def Get_field_variable(sdf.read,var):
	var = 1 should be Ex,Ey,Ex_averaged ...
	or Grid_Grid_mid
	or Density_electron ... 
	or Ekbar_electron...'''
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

	return data.data

def Get_extent(sdffile):
	'''return grid.extents which is a list of [xmin,ymin,xmax,ymax] for 2D
	 or xmin,ymin,zmin,xmax,ymax,zmax for 3D but I just write the version for 2D'''
	grid = sdffile.__dict__['Grid_Grid_mid'];
	print(grid)

	dims = len(grid.dims)
	if (dims == 1):
		index = np.array([0,1])
		extent = np.array(grid.extents)
		myextent = extent[index]
	if (dims == 2):
		index = np.array([0,2,1,3]);
		extent = np.array(grid.extents)
		myextent = extent[index]
	if (dims == 3):
		index = np.array([0,3,1,4,2,5]);
		extent = np.array(grid.extents)
		myextent = extent[index]

	return myextent

def Get_file(prefix,dirc=''):
	''' prefix = 'p' or 'f' 
	dirc is the name directory of the sdf file
	return sdf filename without dirc added'''
	filename = [];
	regu_var = r'^'+prefix+'\d+.sdf$';
	cwd = os.getcwd()
	for files in sorted(os.listdir(cwd+'/'+dirc)):
		if re.match(regu_var,files):
			filename.append(files)
	return filename

##test region
#Get_partvar_npy()



