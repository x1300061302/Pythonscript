#!/usr/bin/python
import sdf 

def Get_particle_variable(deckname,varcombine,species):
	'''def Get_particle_variable(deckname,var,species):
	deckname e.g = 'pxxxx.sdf'
	var = [Px,Py,Gamma,Grid] etc.
	species is defined in input like electron trackele proton 
	also for BlockPointVariable'''
	import re
	data = [];
	for var in varcombine:
		rege_var = '';nofind = 1;
		''' only Grid '''
		if var == 'Grid':
			rege_var = r'^'+var + '_Particles' + '\w*'+species+'$' 
		else:
			rege_var = r'^Particles_'+var+'\w*'+species+'$'
		sdffile = sdf.read(deckname)
		dic = sdffile.__dict__
		for key in dic: 
			if re.match(rege_var,key):
				data.append(dic[key]);
				nofind = 0;
				break;
		
		if (nofind):
			print('Can Not Find '+var+' of '+species+' in '+deckname)
			print(dic.keys())
	return data


def Get_field_variable(deckname,var):
	
	'''def Get_field_variable(deckname,var):
	var = 1 should be Ex,Ey,Ex_averaged ...
	or Grid_Grid_mid
	or Density_electron ... 
	or Ekbar_electron...'''
	import re
	rege_var ='';nofind = 1;data = 0;
	
	rege_var = r'^' + '\w*'+var+'$' 

	sdffile = sdf.read(deckname)
	dic = sdffile.__dict__
	for key in dic: 
		if re.match(rege_var,key):
			data = dic[key];
			nofind = 0;
			break;
	if (nofind):
		print('*******Can Not Find '+var+' in '+deckname+'***********')
		print(dic.keys())
		print('******************************************************')
	return data

'''test'''
