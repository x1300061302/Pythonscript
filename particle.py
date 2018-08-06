#This is the class of particle
import numpy as np  
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
g = 1e-3
mol = 6.02e23
A = 1e-10
ps = 1e-12
kb =1.38e-23
dim = 3
class particle:
	pos = np.zeros(3)
	mass = 1 
	v = np.zeros(3)
	a = np.zeros(3)
	ID = 0 
	def __init__(self,mass=1,pos=[0,0,0],v =[0,0,0],a=[0,0,0],ID = 0):
		self.mass = mass;
		self.pos = pos;
		self.v = v
		self.a = a 
		self.id = ID
	def print_info(self):
		print("mass =",self.mass);
		print('pos =',self.pos);
		print('v = ',self.v);
		print('a = ',self.a);
		print('id =',self.id);
	def normv(self):
		return np.sqrt(self.v[0]**2+self.v[1]**2+self.v[2]**2)

class particle_list:
	Number = 0
	part_list=[]
	def __init__(self):
		self.part_list = []
		self.Number = 0
	def add_particle(self,particle = particle):
		self.Number = self.Number + 1
		self.part_list.append(particle)
	def print_info(self):
		print('Num=',self.Number);
		for i in range(0,self.Number):
			self.part_list[i].print_info()
	def update_part_info(self,n,v,pos):
		self.part_list[n].pos = pos;
		self.part_list[n].v = v;
	def temperature(self):
		#all
		T = 0;
		for part in self.part_list:
			T += 2/3*1/2*part.mass*(g/mol)*part.normv()**2*(A/ps)**2
		T = T/kb/self.Number
		#print('Temp calculated')
		return T
	def draw(self,ax):
		#ax = fig.add_subplot(111,projection='3d')
		for i in range(0,self.Number):
			part = self.part_list[i];
			ax.scatter(part.pos[0],part.pos[1],part.pos[2])
