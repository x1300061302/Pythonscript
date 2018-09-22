#!/usr/bin/python
import sdf 
import re
import numpy as np
import os 
PARTVARNAME_RC = ['Px','Py','Pz','Gamma','Grid']
SPECIES_RC = 'electron'
FIELDNAME_RC = ['Ex_averaged','Ey','Bz_averaged','Density_electron','EkBar_electron']

def Get_index(data,region):
#     print(len(grids[0]))
    index = np.array(data[0] < region[0][1])*np.array(data[0] >region[0][0])*np.array(data[1]<region[1][1])*np.array(data[1]>region[1][0])
    if (len(data[0][index]) == 0):
        print("Warning the length of index = 0,\n Make sure unit of data is right")# error attention
    return index

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
def Get_particle_theta(sdffile,species,dim):
    if (dim == 2):
        py = Get_particle_variable(sdffile,'Py',species);
        px = Get_particle_variable(sdffile,'Px',species);
        theta = np.arctan2(py,px);
        #make sure theta is in (0,2*pi)
    elif (dim == 3):
        py = Get_particle_variable(sdffile,'Py',species);
        pz = Get_particle_variable(sdffile,'Pz',species);
        px = Get_particle_variable(sdffile,'Px',species);
        theta = np.arctan2(np.sqrt(py**2+pz**2),px);

    index = (theta < 0);
    theta[index] = theta[index] + 2*np.pi;
    return theta
        

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

def Get_laser_en(sdffile):
    data = sdffile.__dict__['Absorption_Total_Laser_Energy_Injected__J_'];
    data = data.data;
    return data;

###############--------------------############some operation:
def Get_hist_var(var,bins=500,normed=False):
    hd, ad = np.histogram(var, bins=bins, normed=normed)
    hd = hd/(ad[1]-ad[0])
    axx = 0.5*(ad[1:]+ad[:-1]);
    return hd,axx
def Get_hist2d_var(x,y,bins=[100,200],normed = False):
    h_xy, xedge, yedge = np.histogram2d(x, y, bins=bins)
    xax = np.linspace(np.min(x), np.max(x), bins[0])
    yax = np.linspace(np.min(y), np.max(y), bins[1])
    xx, yy = np.meshgrid(xax, yax)
    return h_xy,xx,yy

    



##test region
#Get_partvar_npy()



