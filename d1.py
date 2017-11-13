import sdfread as sr
import numpy as np
import drawfig as df
import sdf
from const import *

def draw_field(filename,varname):
	a = sdf.read(filename)
	extent =np.array(sr.Get_extent(a))/um;
	print(extent)
	var = sr.Get_field_variable(a,varname);

	df.draw_field_snapshot(var.data.T,\
			extent=extent,\
			xylim=([0,100],[-20,20]),\
			label=('x','y',varname),\
			figname=varname)
	


#main procedure 
draw_field('f0041.sdf','Ex_averaged');
draw_field('f0041.sdf','Ey')
