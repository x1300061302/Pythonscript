import sdfread as sr
import numpy as np
import drawfig as df
import sdf
import os
import matplotlib
import re 
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from const import *
linecolorstyle = ['ro-','g1-','bv-','c1-','ys-','kp-']

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
	dweights = sr.Get_particle_variable(a,'Weight',species);
	gam = dgam.data;
	weights = dweights.data;
	df.draw_spectrum(gam,label =('x','y',varname+prefix),weights = weights, figname = varname+prefix)
	print('ok')

def print_vad(filename,prefix):
	a = sdf.read(filename)
	varname = ('Gamma','Py','Pz','Grid')
	species = 'electron'
	gam = sr.Get_particle_variable(a,varname[0],species);
	gam = gam.data;
	py = sr.Get_particle_variable(a,varname[1],species);
	print(py)
	pz = sr.Get_particle_variable(a,varname[2],species);
	py = py.data;
	pz = pz.data;
	####choose condition 
	grid = sr.Get_particle_variable(a,varname[3],species);	
	y = grid.data[1]
	cr =np.array(y/um>-5)*np.array(y/um<5);
	theta = np.arctan2(py,pz)
	
	df.draw_angle_distribution3d(T = theta[cr],\
							  R = gam[cr],\
							  rlim= 0 ,\
							  figname = 'v_gam_ad'+prefix);

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
	weights = []
	for dirname in dirnames:
		a = sdf.read(dirname+'/'+filename)
		dgam = sr.Get_particle_variable(a,'Gamma','electron');
		gam.append(dgam.data)
		dweight = sr.Get_particle_variable(a,'Weight','electron');
		weights.append(dweight.data)

	prefix = filename[3:5]
	df.draw_spectrum_nline(gam,dataname = dirnames,label = ('$\gamma$','N','$\gamma-N$'),weights=weights,figname = 'gam'+prefix,numl = len(dirnames))

#####################
	
#####################Part B functions combination to set parameter 
	'''in this part one can set some specifical parameter through 
	the functions in Part B for those functions in A ,in generally in part B,
	you just need some operators for files'''
def print_field(varname = 0):
	'''print field_variable
	   default is Ex_averaged,Ey,EkBar_electron,Density_electorn and Bz_averaged for 2-D			'''
	if (type(varname) == np.int):
		varname = []
		varname.append('Ex_averaged')
		varname.append('Ey')
		varname.append('EkBar_electron')
		varname.append('Density_electron')
		varname.append('Bz_averaged')
	else:
		for files in sr.Get_file('f'):
			if (files =='f0000.sdf'):
				pass
			else:
				for var in varname:
					draw_field(files,var,prefix=files[2:5]);
		
def print_particle_spectrum():
	'''just print the electron's spectrum'''
	for files in sr.Get_file('p'):
		print(files)
		print_spectrum(files,prefix = files[2:5]);
	#	print('Wrong in',files,'with gam')

def print_particle_spectrum_compare():

	dirnames = ('a20n01w5','a20n03w5','a20n05w5','a20n08w5','a20n1w5')
	ff = sr.Get_file('p',dirnames[0]);
	for files in ff:
		print_compared_spectrum(dirnames,files)
	
	#	print('Wrong in',files,'with gam')

###############Part C: some special operate
''' in Part C some special operators for draw typical figures to anaylize the physics
	You MUST intepret the functions in this part to  make the function clear'''
def IDrecord():
	ff = 'p0016.sdf'
	a = sdf.read(ff)
	gam = sr.Get_particle_variable(a,'Gamma','electron')
	gam = gam.data
	ID =sr.Get_particle_variable(a,'ID','electron')
	ID = ID.data;
	grid =sr.Get_particle_variable(a,'Grid','electron')
	gdx = grid.data[0]
	index = np.array(gam > np.max(gam)*2/3)*np.array()
	IDre = ID[gam > np.max(gam)*2/3];
	return IDre

def trackparticle_ex_max(IDre):
	fp = sr.Get_file('p')
	ff = sr.Get_file('f')
	xx = np.zeros([10,len(ff)])
	tt = np.linspace(0,len(ff)*5,len(ff)-1)
	xemax = np.zeros([1,len(ff)])

	for i in range(0,len(ff)):
		if re.match(r'^w0000.sdf$',fp[i]):
			pass
		else:
			a = sdf.read(fp[i])
			ID =sr.Get_particle_variable(a,'ID','electron')
			ID =ID.data;
			grid = sr.Get_particle_variable(a,'Grid','electron')
			x = grid.data[0]
	#record the high energy position 
			for idd in range(0,10):
				pos = np.where(ID == IDre[idd])
				xx[idd,i] = x[pos[0][0]]
	
	#ex max position
		if re.match(r'^w0000.sdf$',ff[i]):
			pass
		else:
			a= sdf.read(ff[i])
			ex = sr.Get_field_variable(a,'Ex_averaged');
			grid = sr.Get_field_variable(a,'Grid_Grid_mid');
			ex = ex.data
			nx,ny = grid.dims 
			ex = ex[:,ny//2]
			pos_emax = np.where(ex == np.max(ex));
			xemax[0,i] = grid.data[0][pos_emax[0][0]]		
		
	print(xemax[0,:])
#	print(xx[0,:])
	
#draw figure
	fig = plt.figure()	
	ax = fig.add_subplot(111)
	plt.plot(tt,xemax[0,1:]/um,'r-',linewidth =1.5,label='xemax')
	plt.plot(tt,xx[0,1:]/um,'y-',linewidth =1.5,label='part1')
	plt.plot(tt,xx[1,1:]/um,'g-',linewidth =1.5,label='part2')
	plt.plot(tt,xx[2,1:]/um,'k-',linewidth =1.5,label='part3')
	plt.legend()
	plt.xlabel('t/T0')
	plt.ylabel('x/um')
	plt.savefig('trackpart_exmax'+'.png',dpi=300,facecolor='none',edgecolor='b')


def enterdir(dirname):
	import os 
	os.chdir(dirname)
	return os.getcwd()

def outdir():
	import os
	os.chdir('..')
	return os.getcwd()
def ex_max_evolution():
	ex_max = []
	for file in sr.Get_file('f'):
		print(file)
		a = sdf.read(file)
		dex = sr.Get_field_variable(a,'Ex_averaged')
		ex_max.append(np.max(dex.data))
	return ex_max
def print_xemax():
	dataname = ('n01','n03','n05','n08','n1')
	xemax = []
	for name in dataname:
		xemax.append(np.load(name+'xemax.npy')[0])
	#dirnames = ('.')
	dt = 5; #print frequency 
	tt = dt*np.arange(0,len(xemax[0]))
	#draw_nline
	fig = plt.figure(figsize=[10,8])
	ax = fig.add_subplot(111)
	v = np.zeros(len(dataname))
	for i in range(0,len(dataname)):
		v[i] = (xemax[i][16]-xemax[i][10])/um/30;
		ax.plot(tt,xemax[i][0:25]/um,linecolorstyle[i],label = dataname[i]+' v='+str(v[i])[0:5])
	
	#figureset 
	plt.xlabel('t/T_0')
	plt.ylabel('x/um')
	plt.title('$pos_{exmax}-t$')
	plt.legend()
	plt.show()
	fig.savefig('xemax.png',dpi=300,facecolor='w',edgecolor='b')
	print(v)
	print('Completed')

def print_ex_max_evolutions():
	dataname = ('n01','n03','n05','n08','n1')
	ex_max = []
	for name in dataname:
		ex_max.append(np.load(name+'ex_t.npy'))
	#dirnames = ('.')
	dt = 5; #print frequency 
	tt = dt*np.arange(0,len(ex_max[0]))
	print(len(tt))
	print(len(ex_max[0]))
	#draw_nline
	fig = plt.figure(figsize=[10,8])
	ax = fig.add_subplot(111)

	for i in range(0,len(dataname)):
		ax.plot(tt,ex_max[i][0:25],label = dataname[i])

	plt.legend()
	fig.savefig('ex_t.png',dpi=300,facecolor='w',edgecolor='b')
	print('Completed')

def printdirs_partx_ex():
	dirnames = ('a20n1w5','a20n01w5','a20n03w5','a20n05w5','a20n08w5')
	maker
	for dirname in dirnames:
		en_fail = 0
		try:
			nowdir = enterdir(dirname)
		except:
			en_fail = 1
			print('no',dirname)
	
		draw_partx_ex()
	
		if (en_fail==0):
			nowdir = outdir()
	
		print(nowdir)
	
def draw_partx_ex(beg =0):
	'''this function aims to draw particles' position 'scatter figure'and the 1-D field 
	,for example, in order to analyze if particle is accelerated by ex, I'll draw the 
	higher energy particles in x with ex distribution'''
	import sdf 
	pd = sr.Get_file(prefix = 'p')
	fd = sr.Get_file(prefix = 'f')
	print(pd)
	print(fd)
	for i in range(beg,len(pd)):
		#high energy electron x-pos
		a = sdf.read(pd[i]) 
		dgrid = sr.Get_particle_variable(a,'Grid','electron')
		grid = dgrid.data
		x = grid[0]
		y = grid[1]
		print(type(x))
		dgam = sr.Get_particle_variable(a,'Gamma','electron')
		max_gam = np.max(dgam.data);
		gam = dgam.data;

		x_id = np.array(y/um > -5.0)*np.array(y/um < 5.0)\
				  *np.array(gam > np.max(gam)*2/3);
		xx = x[x_id]/um;
		ggam = gam[x_id];
		####distribution of high energy particles in xdirection
		h_xx,a_xx=np.histogram(xx,bins=4800,normed = False) 
		h_xx = h_xx/np.max(h_xx)
		######## field ex,ey 
		a = sdf.read(fd[i])
		exx = sr.Get_field_variable(a,'Ex_averaged')
		ex = exx.data
		nx,ny = ex.shape #3-D test should be deleted
		eyy = sr.Get_field_variable(a,'Ey')
		ey = eyy.data
		#####ekbar density 
		ekbar = sr.Get_field_variable(a,'Ekbar_electron')
		ekbar = ekbar.data;
		ekbar = ekbar[:,ny//2]

		####charge_density
		dnume = sr.Get_field_variable(a,'Density_electron')
		dnump = sr.Get_field_variable(a,'Density_proton')
		nume = dnume.data;
		nump = dnump.data;
		rho_q = -qe*(nume[:,ny//2] - nump[:,ny//2])

		extent = sr.Get_extent(a) #2-D
		ex = ex[:,ny//2];
		axx = np.linspace(extent[0],extent[1],nx)/um
		#ex 
		#plt.plot(axx,ex/np.max(ex),'r-',linewith = 1.5,label='ex')
		#plt.plot(axx,ey/np.max(ey),'g-',linewidth =1.5,label='ey')
		#plt.scatter(xx,ggam/np.max(gam),c = 'b',norm = 1,edgecolors = 'none',label='')
		#plt.plot(0.5*(a_xx[1:]+a_xx[:-1]),h_xx,'b-',linewidth =1.5,label='hist_x')
		plt.plot(axx,rho_q/np.max(abs(rho_q)),'y-',linewidth =1.5,label='rho_q')
		plt.plot(axx,ekbar/np.max(ekbar),'k--',linewidth = 1.0,label='ekbar')
		plt.legend()
		plt.savefig('ex-partx'+str(i)+'.png',dpi=300,facecolor='none',edgecolor='b')
		plt.close()
def print_partgam_r(prefix = ''):
	'''aims to draw the histogram figure of particle in gam and r '''
	files = sr.Get_file(prefix = 'p')
	for ff in files:
		a = sdf.read(ff)
		Gam = sr.Get_particle_variable(a,'Gamma','electron')
		Gam = Gam.data;
		Grid = sr.Get_particle_variable(a,'Grid','electron')
		Grid = Grid.data;
		y = Grid[1];
		Gamid = np.array(Gam > 1)#np.max(Gam)/3*2)
		df.draw_histogram2d(y[Gamid]/um,Gam[Gamid],\
							label=('y','$\gamma$','dist_y_$\gamma$'),\
							bins = [300,400],\
							#xylim = ([-10,10],[np.min(Gam[Gamid]),np.max(Gam[Gamid])]),\
							figname = 'dist_y_gam'+prefix,\
							display = 1,\
							savefig = 1)
def plot_trackpart():
	dataname = 'n05'
	track = np.load(dataname+'.npy')	
	fig = plt.figure(figsize=[10,8])
	ax = fig.add_subplot(111)
	tt = 5*np.linspace(0,120,25)
	print(track.shape)
	print(len(track[0]))
	print(len(tt))
	for i in range(0,10):
		ax.plot(tt,track[i][:]/um,label='part'+str(i))
	
	#figureset 
	plt.xlabel('t/T_0')
	plt.ylabel('x/um')
	plt.title('$pos_{trackele}-t$')
	plt.legend()
	plt.show()
	fig.savefig(dataname+'.png',dpi=300,facecolor='w',edgecolor='b')
	print('Completed')
def check_region_gamma_value(species='electron'):
	a = sdf.read('p0016.sdf')
	xemax = np.load('xemax.npy')
	xemax = xemax[0][16]
	gam = sr.Get_particle_variable(a,'Gamma',species=species)
	gam = gam.data;
	grid = sr.Get_particle_variable(a,'Grid',species=species)
	gdx = grid.data[0]
	gdy = grid.data[1]
	index = np.array(gdy > -5*um)*np.array(gdy<5*um)
	index = index*np.array(gdx > xemax-0.5*um)*np.array(gdx<xemax+0.5*um)
	mgam = np.mean(gam[index])
	print(mgam)
	np.save('mgam.npy',mgam)
##########Part D dealwith npy file
def get_npy():
	try:
		os.mkdir('data')
	except:
		pass
	sr.Get_partvar_npy()
	sr.Get_fieldvar_npy()
	os.system('mv *.npy data/')
###############Main procedure 

#df.draw_spectrum(gamma*0.511,label=saxis)

