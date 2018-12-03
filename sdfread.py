#!/usr/bin/python
import sdf 
import re
import numpy as np
import os 

PARTVARNAME_RC = ['Px','Py','Pz','Gamma','Grid']
SPECIES_RC = 'electron'
FIELDNAME_RC = ['Ex_averaged','Ey','Bz_averaged','Density_electron','EkBar_electron']


# def generate_key(str_list):
#     varname = str_list[0]
#     if (varname == 'Number_Density'):
#         return 'Derived+''_'.join(str_list)

def Get_EMTensor(a,averaged = True,dim = 2):
    c = 3e8
    '''
    Get F_{ik} = [[0 Ex Ey Ez],
                  [-Ex,0,-Hz,Hy],
                  [-Ey,Hz,0,-Hx],
                  [-Ez,-Hy,Hx,0]]
         F_{ik} = [[F_ik[0] F_ik[1]  ....F_ik[15]  
    '''
    if (averaged):
        Bx = a.Magnetic_Field_Bx_averaged.data
        By = a.Magnetic_Field_By_averaged.data
        Bz = a.Magnetic_Field_Bz_averaged.data
        Ex = a.Electric_Field_Ex_averaged.data
        Ey = a.Electric_Field_Ey_averaged.data
        Ez = a.Electric_Field_Ez_averaged.data
    else:
        Bx = a.Magnetic_Field_Bx.data
        By = a.Magnetic_Field_By.data
        Bz = a.Magnetic_Field_Bz.data
        Ex = a.Electric_Field_Ex.data
        Ey = a.Electric_Field_Ey.data
        Ez = a.Electric_Field_Ez.data
    if (dim == 2):
        nx,ny =  Bx.shape;
    F_ik = np.zeros([nx,ny,16]);
    for i in range(0,nx):
        for j in range(0,ny):
            F_ik[i,j,0] = 0;
            F_ik[i,j,1] = Ex[i,j];
            F_ik[i,j,2] = Ey[i,j];
            F_ik[i,j,3] = Ez[i,j];
            F_ik[i,j,4] = -Ex[i,j];
            F_ik[i,j,5] = 0;
            F_ik[i,j,6] = -c*Bz[i,j];
            F_ik[i,j,7] = c*By[i,j];
            F_ik[i,j,8] = -Ey[i,j];
            F_ik[i,j,9] = c*Bz[i,j];
            F_ik[i,j,10] = 0;
            F_ik[i,j,11] = -c*Bx[i,j];
            F_ik[i,j,12] = -Ez[i,j];
            F_ik[i,j,13] = -c*By[i,j];
            F_ik[i,j,14] = c*Bx[i,j];
            F_ik[i,j,15] = 0;
    
    return F_ik
                 

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
def Get_max_var(key,prefix = '',dirc = '', ffs = []):
    ''' Input:
            key : a.__dict__key or a.key
            prefix and dirc is the file list information like xxx/f0010.sdf prefix = 'f',dirc = 'xxx'
        Output 
            return the maximum value of key in the list of ffs
            and the time at each time.
    '''
    if (len(ffs) == 0):
        ffs = Get_file(prefix = prefix,dirc = dirc);
    key = key.split('.')[1]
    max_var = []
    time = []
    for i in range(0,len(ffs)):
        a = sdf.read(dirc+ffs[i])
        max_var.append(np.max(a.__dict__[key].data));
        time.append(a.Header['time'])
    return max_var,time

def Get_Compare_var(sdflist, key,return_time = False):
    key = key.split('.')[1]
    vars = []
    time = []
    for sdffile in sdflist:
        a = sdf.read(sdffile); 
        try:
            time.append(a.Header['time'])
            vars.append(a.__dict__[key].data)
        except:
            print('Wrong Can not Read')
    if (return_time): 
        return vars,time
    else:
        return vars
    

def Get_field_energy(sdffile,key='all'):
    '''
    just for MR case;
    '''
    #Const 
    c = 3e8
    epsilon0 = 8.854e-12
    a = sdf.read(sdffile);
    dx = a.Grid_Grid.data[0][1]-a.Grid_Grid.data[0][0]
    dy = a.Grid_Grid.data[1][1]-a.Grid_Grid.data[1][0]
    if (key == 'all'):
        ez = a.Electric_Field_Ez.data;
        bx = a.Magnetic_Field_Bx.data;
        by = a.Magnetic_Field_By.data;
        bz = a.Magnetic_Field_Bz.data;
    else:
        ez = a.Electric_Field_Ez_averaged.data;
        bx = a.Magnetic_Field_Bx_averaged.data;
        by = a.Magnetic_Field_By_averaged.data;
        bz = a.Magnetic_Field_Bz_averaged.data;
    nx,ny = bx.shape
    s ={'bx':0,'by':0,'bz':0,'ex':0,'ey':0,'ez':0}
    for ix in range(0,nx):
        for iy in range(0,ny):
            s['bx'] = s['bx'] + c**2 * bx[ix,iy]**2
            s['by'] = s['by'] + c**2 * by[ix,iy]**2
            s['bz'] = s['bz'] + c**2 * bz[ix,iy]**2
#             sb = sb + sbx + sby + sbz;
            s['ez'] = s['ez'] + ez[ix,iy]**2
    for key in s.keys():
        s[key] = s[key] * 0.5*epsilon0*dx*dy
#     sb = 0.5*epsilon0*sb*dx*dy
#     se = 0.5*epsilon0*se*dx*dy
    return s
#     grids = a.

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
def Get_time(prefix=''):
    ffs = Get_file(prefix = prefix);
    time = []
    for i in range(0,len(ffs)):
        a = sdf.read(ffs[i]);
        time.append(a.Header['time']);
    return np.array(time) 
def Get_curlB(sdffile):
    a = sdf.read(sdffile)
    try:
        Bx = a.Magnetic_Field_Bx;
        By = a.Magnetic_Field_By;
        dx = Bx.grid.data[0][1] - Bx.grid.data[0][0];
        dy = Bx.grid.data[1][1] - Bx.grid.data[1][0]
    except:
        print('Wrong, can not read Magnetic Field')
    curlB = curl(Ax = Bx,Ay = By,dx = dx,dy = dy,order = 2);
    return curlB
def curl(Ax,Ay,dx = 1.0,dy = 1.0,Az=0,order=2):
    nx,ny = Ax.shape
    curlAz = np.zeros(Ax.shape)
    cx = 1.0/(dx) 
    cy = 1.0/(dy) 
    for ix in range(1,nx):
        for iy in range(1,ny):
            curlAz[ix,iy] = cx*(Ay[ix  , iy] - Ay[ix-1, iy])-cy*(Ax[ix  , iy]   - Ax[ix  , iy-1])
    return curlAz
def Get_energy(prefix='',dirc = ''):
    ffs = Get_file(prefix = prefix,dirc = dirc);
    TFE = []
    TPE = []
    for i in range(0,len(ffs)):
        a = sdf.read(dirc + ffs[i]);
        TFE.append(a.Total_Field_Energy_in_Simulation__J_.data)
        TPE.append(a.Total_Particle_Energy_in_Simulation__J_.data)
    return np.array(TFE),np.array(TPE)

###############--------------------############some operation:
def Get_hist_var(var,weights = 0, bins=500,normed=False):
    if (type(weights) == int):
        weights = np.ones(var.shape)+weights;
    hd, ad = np.histogram(var, weights = weights,bins=bins, normed=normed)
    hd = hd/(ad[1]-ad[0])
    axx = 0.5*(ad[1:]+ad[:-1]);
    return hd,axx
def Get_hist2d_var(x,y,bins=[100,200],normed = False):
    h_xy, xedge, yedge = np.histogram2d(x, y, bins=bins)
    xax = np.linspace(np.min(x), np.max(x), bins[0])
    yax = np.linspace(np.min(y), np.max(y), bins[1])
    xx, yy = np.meshgrid(xax, yax)
    return h_xy,xx,yy

#############-------------Tracer script------------##########
def tp_pos(ID,fts=0):
    if (type(fts) == np.int):
        fts = Get_file('tr');
    #length of pos   
    l_pos = len(fts);
    l_id = len(ID);
    print(l_pos,l_id)
        
#     print('a',ID[0])
    pos_list = np.zeros([l_pos,l_id]);
    print(pos_list.shape)
    for i in range(0,l_pos):
  
        a = sdf.read(fts[i]);
        ID_list = Get_particle_variable(sdffile=a,species='tra_ele',var='ID')
#         print(i,len(ID_list))
#         print(len(ID_list))
        for j in range(0,l_id):
            pos = np.where(ID_list == ID[j])[0][0]
#             print(i,j,pos)
            pos_list[i,j] = pos
            
            
    return pos_list

def tp_var(pos,varname,species='tra_ele',fts=0,dim = 2):
    if (type(fts) == np.int):
        fts = Get_file('tr');
    l_pos = len(pos);
# print(pos[0,:])
    tvar = []
    for i in range(0,l_pos):
        a = sdf.read(fts[i]);
        var = Get_particle_variable(sdffile=a,var=varname,species=species);
        int_pos = int(pos[i,:])
#         print(int_pos)
#         print(varname)
        if (varname == 'Grid'):
#               print(len(var))
#             print('a')
            tx= var[0][int_pos];
            ty= var[1][int_pos];
            if (dim == 3):
                tz= var[2][int_pos];
                tvar.append([tx,ty,tz]);
            tvar.append([tx,ty]);
        else:
            tvar.append(var[int_pos]);
        
    return tvar


##test region
#Get_partvar_npy()



