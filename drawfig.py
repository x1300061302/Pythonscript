#!/usr/bin/python

def draw_angle_distribution3d(T,R,rlim = 0,figname='fig',binT=360,binR=1000):
	'''T,R is N-array 
	rlim = 0 means default else rlim = [r_min,r_max]
	draw angle distribution histogram'''
	try:
		import matplotlib
		matplotlib.use('Agg')
	except:
		print('use matplotlib.use before import plt')        
	import matplotlib.cm as cm
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D
	import numpy as np 
	'''histogram'''
	H,Tedge,Redge = np.histogram2d(T,R,bins=[binT,binR])
	log10H = np.log10(H)
	print("get histogram")
	'''change T,R to meshgrid'''
	rad = np.linspace(R.min(),R.max(),binR)
	the = np.linspace(0,2*np.pi,binT)
	r,t = np.meshgrid(rad,the)
	fig = plt.figure()
	ax = plt.subplot(projection = 'polar')
	pcmesh = ax.pcolormesh(t,r,log10H,vmin = 0,vmax = log10H.max())
	if (rlim != 0):
		ax.set_rlim(rlim)
	cbar = fig.colorbar(pcmesh)
	ffigname = '.'.join([figname,'png']);
	fig.savefig(ffigname,dpi = 300,facecolor='w',edgecolor='b')
	print(figname,'has been printed');

def draw_spectrum(data,dataname,weight=0,figname='fig',numl=1,logx=0,display=0):
	'''data is 2-D x numl array
	dataname is numl length array ['','',''] 
	figname default = fig
	make sure numl is given'''
	if (display==0):
		try:
			import matplotlib
			matplotlib.use('Agg')
		except:
			print('use matplotlib.use before import plt')
	
	import matplotlib.pyplot as plt
	import numpy as np 
	line = []
	colorline = ('red','black','blue')
	fig = plt.figure(figsize=(8,4))
	ax = plt.subplot(111);
	if weight == 0: 
		print('No Weight hist')
		for i in range(0,numl):
			hd,ad = np.histogram(data[i][:],bins = 500,normed = False)
			if logx==1:
				ad = np.log10(ad)
			pl = plt.semilogy(.5*(ad[1:]+ad[:-1]),hd,color=colorline[i],label=dataname[i],linewidth=2) 
			line.append(pl);
			print("has printed",dataname[i])
	else:
		print('Weight hist')
		for i in range(0,numl):
			hd,ad = np.histogram(data[i][:],bins = 500,normed = False,weight=weight[i][:])
			if logx==1:
				ad = np.log10(ad)
			plt.semilogy(.5*(ad[1:]+ad[:-1]),hd,color=colorline[i],label=dataname[i],linewidth=2) 
			print("has printed",dataname[i])
			
	plt.legend(tuple(line),dataname,loc='upper right')
	ffigname = '.'.join([figname,'png']);
	fig.savefig(ffigname,dpi = 300,facecolor='w',edgecolor='b')
	if (display):
		plt.show()
	return fig,line 

def draw_field_snapshot(data,extent,label,numtick=5,Display=0,figname='fig'):
	'''data should be 2-D array (nx,ny) 
	extent should be ([x_min,x_max],[y_min,y_max])
	label = ['xlabel/Unit','ylabel/Unit','title'] 
	Display dicide whether to show
	figname default = title'''
	if (Display==0):
		try:
			import matplotlib
			matplotlib.use('Agg')
		except:
			print('use matplotlib.use before import plt')        	
	import matplotlib.cm as cm
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D
	import numpy as np 

	fig = plt.figure(figsize = (6,6));
	ax = fig.add_subplot(111)

	nx,ny = data.shape
	aspect = ny/nx*0.8

	gci = ax.imshow(data,extent = extent,origin='lower',cmap='jet',\
			 vmax=data.max(),vmin=data.min(),interpolation='spline36')

	#set x_ticklabel 
	axx = ax.xaxis;
	xticklocs = np.linspace(0,ny,5)
	xtick = np.linspace(extent[0],extent[1],numtick)
	xticks = [str(k) for k in xtick]
	print(xticks)
	axx.set_ticks(xticklocs)
	axx.set_ticklabels(xticks)
	#set y_ticklabel 
	ayy = ax.yaxis;
	yticklocs = np.linspace(0,nx,5)
	ayy.set_ticks(yticklocs)
	ytick = np.linspace(extent[2],extent[3],numtick)
	yticks = [str(k) for k in ytick]
	ayy.set_ticklabels(yticks)
	#label 
	plt.xlabel(label[0])
	plt.ylabel(label[1])
	plt.title(label[2])
	ax.set_ylim([0,5.0])
	ax.set_xlim([0,5.0])
	#colorbar
	plt.colorbar(gci);

	if (Display):
		plt.show();
	else:
		fig.savefig(figname,dpi=300,facecolor='w',edgecolor='b');	
	
		

