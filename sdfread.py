#!/usr/bin/python
import sdf 
import re
import numpy as np
import os 
import gc
### const Part
c = 3e8
qe = 1e-16


PARTVARNAME_RC = ['Px','Py','Pz','Gamma','Grid']
SPECIES_RC = 'electron'
FIELDNAME_RC = ['Ex_averaged','Ey','Bz_averaged','Density_electron','EkBar_electron']


# def generate_key(str_list):
#     varname = str_list[0]
#     if (varname == 'Number_Density'):
#         return 'Derived+''_'.join(str_list)
def get_FWHM(data,dx=1.0,peak=None):
    base = data[0]
    data = data-base
    if peak == None:
        p = np.max(data);
    else:
        p = peak;
    p_2 = p/2;
    xx = np.linspace(0,len(data),len(data))
    pos = np.where(data > p_2)[0]
    p1 = pos[0];
    pos2 = np.zeros(len(pos));
    for i in range(len(pos)):
        if i == 0:
            pos2[i] = 1
        else:
            pos2[i] = pos[i] - pos[i-1]
    try:
        p2 = np.where(pos2 > 1)[0][0]
    except:
        p2 = 0
    FWHM = p2 - p1
    FWHM =FWHM*dx
    x1 = p1*dx
    x2 = p2*dx
    return FWHM,x1,x2,p/2+base
def get_mvp(sdffile,avg = True):
    a = sdf.read(sdffile);
    try:
        xx = a.Grid_Grid_mid.data[0]
        yy = a.Grid_Grid_mid.data[1]
        dx = xx[1] - xx[0]
        dy = yy[1] - yy[0]
    except:
        print('Grid wrong')
        dx = 1.0
        dy = 1.0
    if avg:
        bx = a.Magnetic_Field_Bx_averaged.data;
        by = a.Magnetic_Field_By_averaged.data;
    else:
        bx = a.Magnetic_Field_Bx.data;
        by = a.Magnetic_Field_By.data;
    mvp = np.zeros(bx.shape);
    Nx,Ny = mvp.shape
    integral_y = 0.0
    for iy in range(Ny):
        integral_y = integral_y + bx[0,iy]
        integral_x = 0.0
        for ix in range(Nx):
            integral_x = integral_x - by[ix,iy]
            mvp[ix,iy] = integral_x*dx + integral_y*dy 
    return np.array(mvp),xx,yy

def Calculate_Flux(a,energy_range,direction = None,species = 'subset_pho_photon',csp_limit = None):
    '''
    Input:a,energy_range
    direction = None
    species = 'subset_pho_photon'
    Output:
    sum_Flux
    sum_energy_flux
    '''
    sum_Flux = 0;
    px = Get_particle_variable(sdffile=a,var='Px',species=species)
    py = Get_particle_variable(sdffile=a,var='Py',species=species)
    pz = Get_particle_variable(sdffile=a,var='Pz',species=species)
    ekp2 = Get_particle_variable(sdffile = a, var = 'Ek',species = species)
    weisp = Get_particle_variable(sdffile = a, var = 'Weight',species = species)
    #1 energy_range
    index = (ekp2> energy_range[0])*(ekp2 < energy_range[1])
    if (direction !=None):
        cos_theta = (px*direction[0] + py*direction[1] + pz*direction[2])/np.sqrt(px**2 + py**2 + pz**2)
        Omega = 2*np.pi*(1 -  cos_theta);
        index2 = Omega < 4*np.pi/100;
        index = index*index2
    if (csp_limit != None ):
        csp = Get_particle_variable(sdffile = a,var='Compton_Scatter_Times',species = species);
        index3 = csp > csp_limit
        index = index*index3
    sum_Flux = sum(weisp[index])
    sum_Energy_flux = sum(ekp2[index]*weisp[index])
    return sum_Flux,sum_Energy_flux
def makedirs(dirc):
    try:
        os.makedirs(dirc)
    except:
        print('exist')
    return dirc + '/'
def Track2D(ID,ffs,species):
    nn = len(ffs)
    xpos = np.zeros(nn)
    ens = np.zeros(nn)
    ypos = np.zeros(nn)
    tts = np.zeros(nn)
    for i in range(0,nn):
        a = sdf.read(ffs[i]);
        grids = Get_particle_variable(sdffile=a,var='Grid',species=species)
        en = Get_particle_variable(sdffile = a,var = 'Gamma',species=species)
        ids = Get_particle_variable(sdffile = a,var = 'ID',species = species)
        xx = grids[0]
        yy = grids[1]
        pos = np.where(ids == ID)[0]
        xpos[i] = xx[pos]
        ypos[i] = yy[pos]
        ens[i] = en[pos]
        tts[i] = a.Header['time']
    try:
        return xpos,ypos,ens,tts
    finally:
        del xpos
        del ypos
        del ens
        del tts
        del grids
        del en
        del ids
        del xx
        del yy
        del pos
        del a

def get_times(files,ran = None):
    '''
    Input: files of sdf file
           ran = range(0,len(files),1)
    Output: times: np.array(len(files))
    '''
    times = []
    if (ran == None):
        ran = range(0,len(files))
    for i in ran:
        a = sdf.read(files[i]);
        times.append(a.Header['time'])
    return np.array(times)




def smooth_data(data,sigma = 2):
    from scipy import ndimage
    sd = ndimage.gaussian_filter(data,sigma = sigma)
    return sd

def Get_EMTensor(a,averaged = True,dim = 2):
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
	return data

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
def Get_max_var(key,prefix = '',dirc = '', ffs = [],region=0):
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
        a = sdf.read(ffs[i])
        if (type(region) == np.int):
            region[1],region[3] = a.__dict__[key].data.shape();
            region[0],region[2] = [0,0];
        max_var.append(np.max(a.__dict__[key].data[region[0]:region[1],region[2]:region[3]]));
        time.append(a.Header['time'])
    return np.array(max_var),np.array(time)

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
        return np.array(vars),np.array(time)
    else:
        return np.array(vars)
    

def Get_field_energy(sdffile,key='all'):
    '''
    just for MR case;
    '''
    #const 
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

def Get_file(prefix,dirc='',suffix = '.sdf'):
    ''' prefix = 'p' or 'f' 
    dirc is the name directory of the sdf file
    return sdf filename without dirc added'''
    filename = [];
    regu_var = r'^'+prefix+'\d+'+suffix+'$';
    for files in sorted(os.listdir(dirc)):
        if re.match(regu_var,files):
            filename.append(dirc+files)
    if (len(filename) == 0):
        print('No files');
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
def Get_curlB(sdffile,dx = 1.0,dy = 1.0):
    a = sdf.read(sdffile)
    try:
        Bx = a.Magnetic_Field_Bx.data;
        By = a.Magnetic_Field_By.data;
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
def Get_optical_stokes(prefix = '',dirc = '',files = None):
    '''Get Total_Optical_I1 and Total_optical_I2 and U
    input: prefix and dirc
           ret means if return time and files
    output: toi1,toi2,tou,times
    '''
    ffs = Get_file(prefix = prefix,dirc = dirc);
    toi1 = []
    toi2 = []
    tou = []
    time = []
    if (files !=None):
        ffs = files
    for i in range(0,len(ffs)):
        a = sdf.read(ffs[i]);
        time.append(a.Header['time'])
        toi1.append(a.Total_Optical_Intensity_I1_in_Simulation.data)
        toi2.append(a.Total_Optical_Intensity_I2_in_Simulation.data)
        tou.append(a.Total_Optical_Intensity_U_in_Simulation.data)
    return np.array(toi1),np.array(toi2),np.array(tou),np.array(time)

def Get_energy(prefix='',dirc = '',ret = False,files = None):
    '''Get TFE and TPE
    input: prefix and dirc
           ret means if return time
    output: TFE TPE
    
    '''
    ffs = Get_file(prefix = prefix,dirc = dirc);
    TFE = []
    TPE = []
    time = []
    if (files !=None):
        ffs = files
    for i in range(0,len(ffs)):
        a = sdf.read(ffs[i]);
        time.append(a.Header['time'])
        TFE.append(a.Total_Field_Energy_in_Simulation__J_.data)
        TPE.append(a.Total_Particle_Energy_in_Simulation__J_.data)
    if (ret):
        return np.array(TFE),np.array(TPE),np.array(time)
    else:
        return np.array(TFE),np.array(TPE)

###############--------------------############some operation:
def Get_hist_var(data,weights = 0, bins=500,normed=False):
    '''
    Input: data,weights  
    Output: hd,axx
    '''
    if (type(weights) == int):
        weights = np.ones(data.shape)+weights;
    hd, ad = np.histogram(data, weights = weights,bins=bins, normed=normed)
    hd = hd/(ad[1]-ad[0])
    axx = 0.5*(ad[1:]+ad[:-1]);
    return axx,hd
def Get_hist2d_var(x,y,bins=[100,200],normed = False,weights = None):
    '''
    Input:
        x,y,
        bins = [100,200]
        normed = False
    Output:
        xx,yy,h_xy
    '''
    h_xy, xedge, yedge = np.histogram2d(x, y, bins=bins,weights = weights)
    xax = np.linspace(np.min(x), np.max(x), bins[0])
    yax = np.linspace(np.min(y), np.max(y), bins[1])
    xx, yy = np.meshgrid(xax, yax)
    return xx,yy,h_xy

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

def read_const(dirc):
    info = {}
    try:
        f = open(dirc+'/const.status','r')
    except:
        print('can not find const.status')
        return;
    lines = f.readlines()
    for line in lines:
        exp = line.split()
        info[exp[0]] = float(exp[2])
    f.close()
    return info
    

##########Class Simu_info()-------#################

class simu_info():
    '''
    This Class is aimed to get the information of the Simulation.
    And normalized parameter in Simulations.
    '''
    def __init__(self,dirc='',name = ''):
        self.name = name;
        self.axis={};
        self.map={};
        self.const={'di':1,'T0':1,'B0':1,'E0':1,'J0':1};
        self.Nx,self.Ny = [1,1]
        self.read_input(dirc); #read nx,ny
        self.read_const(dirc = dirc)
        try:
            self.axis['x'] = np.linspace(self.const['xmin']/self.const['di'],self.const['xmax']/self.const['di'],self.Nx);
            self.axis['y'] = np.linspace(self.const['ymin']/self.const['di'],self.const['ymax']/self.const['di'],self.Ny);
            self.map['xx'],self.map['yy'] = np.meshgrid(self.axis['x'],self.axis['y']);  
        except:
            print('No parameters')
        
    def read_input(self,dirc):
        '''
        Input:  directory of input.deck 
        Output:
        nx,ny 
        '''
        try:
            f = open(dirc + '/input.deck','r')
        except:
            print('can not find input.deck')
            return;
        lines = f.readlines()
        for line in lines:
            exp = line.split()
            if (not exp):#if empty
                continue;
            if (exp[0] in ['nx','ny']): 
                self.const[exp[0]] = int(exp[2])
        f.close()
        self.Nx =  self.const['nx']
        self.Ny =  self.const['ny']

    def read_const(self,dirc):
        try:
            f = open(dirc+'/const.status','r')
        except:
            print('can not find const.status')
            return;
        lines = f.readlines()
        for line in lines:
            exp = line.split()
            self.const[exp[0]] = float(exp[2])
        self.const['E0'] = self.const['B0']*c
        self.const['J0'] = self.const['drift_V']*self.const['ne']*qe

        f.close()
    def get_extent(self,sdffile):
        a = sdf.read(sdffile);
        extent = Get_extent(a);
        return np.array(extent)

    
    
#####------Class mrc ----------######  
def get_juuj(a):
    '''
    Input:
    a:sdffile
    Output:
    the Jv+Jv term for calculating the Ohm
    '''
    j = {}
    u = {}
    juuj={}
    j['x'] = a.Current_Jx.data
    j['y'] = a.Current_Jy.data
    j['z'] = a.Current_Jz.data
    u['x'] = a.Derived_Velocity_Ux.data
    u['y'] = a.Derived_Velocity_Uy.data
    u['z'] = a.Derived_Velocity_Uz.data
    for dir1 in ['x','y','z']:
        for dir2 in ['x','y','z']:
            juuj[dir1+dir2] = j[dir1]*u[dir2]+u[dir1]*j[dir2];
    return juuj

def cross(a,b):
    c = np.zeros(3)
    c[2] = a[0]*b[1] - a[1]*b[0]
    c[1] = a[1]*b[2] - a[2]*b[1]
    c[0] = a[2]*b[1] - a[1]*b[2]
    return c 
def calc_cross(Ax,Ay,Bx,By):
    Nx,Ny = Ax.shape
    C = np.zeros([Nx,Ny])
    for ix in range(Nx):
        for iy in range(Ny):
            a = [Ax[ix,iy],Ay[ix,iy],0]
            b = [Bx[ix,iy],By[ix,iy],0]
            C[ix,iy] = cross(a,b)[2]
    return C

    
def calc_mean_var(data,x,y, region):
    ''' Calculate the mean value in the specific region
    Input: data 
           x,y,region with same unit 
    Output: meanvalue
    '''
#     var = a.Electric_Field_Ez_averaged;
#     x = ez.grid_mid.data[0]
#     y = ez.grid_mid.data[1]
    # print(len(x))
    mv = 0.0

    dx = x[1]-x[0]
    dy = y[1]-y[0]
    xmin,xmax,ymin,ymax = region;
    index_x = np.array(x < xmax)*np.array(x > xmin); 
    index_y = np.array(y < ymax)*np.array(y > ymin); 
    cell_x0 = int((xmin - x[0])/dx)
    cell_x1 = int((xmax - x[0])/dx)

    cell_y0 = int((ymin - y[0])/dy)
    cell_y1 = int((ymax - y[0])/dy)

    w = 0.0
    for ix in range(cell_x0,cell_x1):
        for iy in range(cell_y0,cell_y1):
            mv = mv + data[ix,iy]
            w = w +1
    return (mv/w)

def calc_mean_var_line(data,direction,cell,delta_grid = 1):
    '''
    Input:data = np.array(Nx,Ny)
    direction = 'x' is to calculate the mean value in y direction
    direction = 'y' is to calculate the mean value in x direction
    delta_grid = default =  1
    cell = centre cell
    Output:data_mean_line around [cell - delta_grid, cell + delta_grid]
    
    '''
    Nx,Ny = data.shape
    if (direction == 'x'):
        data_mean_line = np.zeros(Nx)
        for iy in range(cell-delta_grid,cell+delta_grid+1):
            data_mean_line[:] = data_mean_line[:] + data[:,iy]
    if (direction == 'y'):
        data_mean_line = np.zeros(Ny)
        for ix in range(cell-delta_grid,cell+delta_grid+1):
            data_mean_line[:] = data_mean_line[:] + data[ix,:]
    return data_mean_line/(delta_grid*2+1)

def calc_dot_2d(Ax,Ay,region,dx,dy,delta_nx,delta_ny):
        nx,ny = Ax.shape;
        gx = np.zeros(Ax.shape);
        gy = np.zeros(Ax.shape);
        for i in range(region[0],region[1]):
            for j in range(region[2],region[3]):
                if (i == 0):
                    gx[i,j] = (Ax[i+delta_nxx,j]-Ax[i,j])/delta_nx/dx;
                elif (i == nx-1):
                    gx[i,j] = (Ax[i,j]-Ax[i-delta_nx,j])/delta_nx/dx;
                else:
                    gx[i,j] = (Ax[i+delta_nx,j]-Ax[i-delta_nx,j])/2/delta_nx/dx;
                if (j == 0):
                    gy[i,j] = (Ay[i,j+delta_ny]-Ay[i,j])/delta_ny/dy;
                elif (j == ny-1):
                    gy[i,j] = (Ay[i,j]-Ay[i,j-delta_ny])/delta_ny/dy;
                else:  
                    gy[i,j] = (Ay[i,j+delta_ny]-Ay[i,j-delta_ny])/2/delta_ny/dy;
        return gx+gy
            
def calc_gradient(data,dx,dy=0,region = 0,dim = 1,delta_n = 1,delta_ny=1):
    '''
    Input data is a scalar variation 2-D or 1-D
    Output data is a vector 2-D or 1-D
    '''
    if (dim == 1):
        nx = len(data);
        gx = np.zeros(nx);
        if (type(region) == np.int):
            cell_x0 = 0;
            cell_x1 = nx
        else:
            cell_x0,cell_x1 = region

        for i in range(cell_x0,cell_x1):
            #raise TypeError()
            if (i == 0):
                gx[i] = (data[i+delta_n]-data[i])/delta_n/dx;
            elif (i == nx-1):
                gx[i] = (data[i]-data[i-delta_n])/delta_n/dx;
            else:
                gx[i] = (data[i+delta_n]-data[i-delta_n])/2/delta_n/dx;
        return gx
    if (dim == 2):
        nx,ny = data.shape;

        if (type(region) == np.int):
            region =[0,nx,0,ny]
        else:
            cell_x0,cell_x1,cell_y0,cell_y1 = region

        gx = np.zeros(data.shape);
        gy = np.zeros(data.shape);
        for i in range(region[0],region[1]):
            for j in range(region[2],region[3]):
                if (i == 0):
                    gx[i,j] = (data[i+delta_n,j]-data[i,j])/delta_n/dx;
                elif (i == nx-1):
                    gx[i,j] = (data[i,j]-data[i-delta_n,j])/delta_n/dx;
                else:
                    gx[i,j] = (data[i+delta_n,j]-data[i-delta_n,j])/2/delta_n/dx;
                if (j == 0):
                    gy[i,j] = (data[i,j+delta_ny]-data[i,j])/delta_ny/dy;
                elif (j == ny-1):
                    gy[i,j] = (data[i,j]-data[i,j-delta_ny])/delta_ny/dy;
                else:  
                    gy[i,j] = (data[i,j+delta_ny]-data[i,j-delta_ny])/2/delta_ny/dy;
        return gx,gy
            
def Magnetic_Flux(sdfname,region,xx=0,yy=0):
        ''' Input:
            sdfname: 0000.sdf 
            region = [xmin,xmax,ymin,ymax] of your integral region
            Output:
            Fx =int By*dx
            Fy =int Bx*dy 
        '''
        a = sdf.read(sdfname)
        Bx = a.Magnetic_Field_Bx_averaged;
        By = a.Magnetic_Field_By_averaged;
        nx,ny = Bx.data.shape
        if (type(xx) !=np.int):
            x = xx;
            y = yy;
        else:
            extents = Bx.grid_mid.extents;
            x = np.linspace(extents[0],extents[2],nx)
            y = np.linspace(extents[1],extents[3],ny)
        if (nx != len(x)):
            print("The shape nx,ny is",nx,ny,', but the length x,y is ',len(x),len(y))
            return 0;
        #return cell_x and cell_y 
        dx = x[1] - x[0];
        dy = y[1] - y[0];
        #print('dx=',dx,'dy=',dy)
        xmin,xmax,ymin,ymax = region;
        index_x = np.array(x < xmax)*np.array(x > xmin); 
        #print(xmax,xmin)
        index_y = np.array(y < ymax)*np.array(y > ymin); 
        cell_x0 = int((xmin - x[0])/dx)
        cell_y0 = int((ymin - y[0])/dy)
        #print(cell_x0,cell_y0)
        #Flux_y = int^{y_max}_{y_min} Bx dy
        Bx2 = Bx.data[cell_x0,:];
        By2 = By.data[:,cell_y0];
        Flux_y = np.trapz(Bx2[index_y],y[index_y])
        Flux_x = np.trapz(By2[index_x],x[index_x])
        return Flux_x,Flux_y    

def Find_p(sdffile, region,index = 128):
    a = sdf.read(sdffile)
    bx = a.Magnetic_Field_Bx_averaged.data;
    by = a.Magnetic_Field_By_averaged.data;
    abs_b = np.sqrt(bx**2 + by**2)
    data = abs_b[:,index];
    pos = np.where(data[region[0]:region[1]] == np.min(data[region[0]:region[1]]))[0]
    return pos + region[0]
def get_slope(Xi,Yi,p0,ran = None,logscale = None):
    '''Input: 
        Xi,Yi
       Output:
       p0 
    '''
    if (logscale == 'x'):
        Xi = np.log10(Xi)
    if (logscale == 'y'):
        Yi = np.log10(Yi)
    if (logscale == 'xy'):
        Xi = np.log10(Xi)
        Yi = np.log10(Yi)
    if (ran != None):
        index = np.array(Xi > ran[0])*np.array(Xi < ran[1])
        Xi = Xi[index]
        Yi = Yi[index]
        print(len(Xi),len(Yi))
    from scipy.optimize import leastsq
    def func(p,x):
        k,b=p
        return k*x+b

    def error(p,x,y):
        return func(p,x)-y #x、y都是列表，故返回值也是个列表
    #初始值
    p0=p0
    Para=leastsq(error,p0,args=(Xi,Yi)) #把error函数中除了p以外的参数打包到args中
    k,b=Para[0]
    print("k=",k,'\n',"b=",b)
    return k,b
    
def get_mag_flux(region, prefix = '',dirc = '',files = None):
    '''Input:
        region: in real unit like [xmin,xmax,ymin,ymax]
        prefix
        dirc
       Output:
       return fxs,fys,time(datatype: np.array)
    '''

    ffs = Get_file(prefix = prefix,dirc = dirc);
    time = []
    fxs = []
    fys = []
    if files != None:
        ffs = files
    for i in range(1,len(ffs),1):
        time.append(sdf.read(ffs[i]).Header['time']);
        if (type(region) == np.ndarray):
            fx,fy = Magnetic_Flux(ffs[i],region = region);
        elif(type(region) == list):
            fx,fy = Magnetic_Flux(ffs[i],region = region[i]);
        fxs.append(fx);
        fys.append(fy);
    return np.array(fxs),np.array(fys),np.array(time)
def Get_De(sdffile):
    a = sdf.read(sdffile)
    Jx = a.Current_Jx_averaged.data;
    Jy = a.Current_Jy_averaged.data;
    Jz = a.Current_Jz_averaged.data
    Bx = a.Magnetic_Field_Bx_averaged.data
    By = a.Magnetic_Field_By_averaged.data
    Bz = a.Magnetic_Field_Bz_averaged.data
    Ex = a.Electric_Field_Ex_averaged.data
    Ey = a.Electric_Field_Ey_averaged.data
    Ez = a.Electric_Field_Ez_averaged.data
    Nume = a.Derived_Number_Density_electron.data;

    nx,ny = Nume.shape
    print(nx,ny)
    Vx = Jx/Nume/qe
    Vy = Jy/Nume/qe
    Vz = Jz/Nume/qe

    Dex = Jx*(Ex+Vy*Bz-Vz*By)
    Dey = Jy*(Ey+Vz*Bx-Vx*Bz)
    Dez =  Jz*(Ez + Vx*By - Vy*Bx)

    VdotE = Vx*Ex + Vy*Ey + Vz*Vz
    De = Dex+ Dey + Dez
    return De
##test region
#Get_partvar_npy()

########People's model 
def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if(x.ndim != 1):
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if (x.size < window_len):
        raise ValueError("Input vector needs to be bigger than window size.")


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y



#########demonstrate how to use it

