import sys
import sdf 
import matplotlib.pyplot as plt
import os 
import re 

a = r'^tp\d{4}.sdf$'

for root,dirs,files in os.walk(os.getcwd()):
	print(files)
'''above read file'''

for filename in files:
	if (re.match(a,filename)):
		ff = sdf.read(filename)
		'''read sdf'''
		keys = ff.__dict__.keys()
		print(keys)
		print('Next--------------');
		dataname = r'^\w+_Gamma_\w+$'
		for key in keys: 
			if (re.match(dataname,key)):
				print(key)
				data = ff.__dict__[key] 
				






