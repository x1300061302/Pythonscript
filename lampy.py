#!/anaconda3/bin/python
from const import *
from particle import particle,particle_list


#unit = metal 
#mass = grams/mole
f1 = 'testFe2.lammpstrj'

times = 0; #dump times
fh = open(f1,'r')
lines = fh.readlines()
def readbf(fh):
	'''to read the basis information'''
	i = 0;
	while(i < 10):
		line = lines[i];
		data = line.split()
		if ((data[0] == 'ITEM:') and (data[1] == 'TIMESTEP')):
			i = i + 1;
			times = int(lines[i]);
		if ((data[0] == 'ITEM:') and (data[1] == 'NUMBER')):
			i = i + 1;
			Num = int(lines[i]);
		if ((data[0] =='ITEM:') and (data[1] == 'ATOMS')):
			begin_line = i + 1;
		i = i + 1 
	Times = int(len(lines)/(Num+begin_line));
	return Num,begin_line,Times

Num, begin_line,Times = readbf(fh)
print(readbf(fh))
mass = 1
p1 = particle_list() #center
p2 = particle_list() #wall
pt1 =[]
pt2 =[]
T1 =[]
T2 =[]
for j in range(0,Times):
	for i in range(begin_line,begin_line + Num):
		line = lines[i];
		data = line.split()
		vx = float(data[5]);
		vy = float(data[6]);
		vz = float(data[7]);
		id = float(data[0]);
		typeid = data[1]
		x = float(data[2]);
		y = float(data[3]);
		z = float(data[4]);
		a = particle(mass = 56,pos=[x,y,z],v=[vx,vy,vz],ID=id)
		if ((typeid == '2')):
		#	p1.append([vx,vy,vz])   # id = 2 
			p2.add_particle(particle = a);
		elif((typeid == '1')):
			p1.add_particle(particle = a);
	pt1.append(p1)
	pt2.append(p2)
	T1.append(p1.temperature())
	T2.append(p2.temperature())

np.save('pt1.npy',pt1)
np.save('pt2.npy',pt1)
np.save('T1.npy',T1)
np.save('T2.npy',T2)

print(T1,T2)
fh.close()
