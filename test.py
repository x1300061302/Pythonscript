import numpy as np

a = np.linspace(0,1000,1000)
b = np.linspace(0,1000,1000)
c = np.linspace(0,1000,1000)

d = a > 900;
c = a < 990;
m = d*c; 
m2 = np.array(a>90)*np.array(a<900)
print(m2)
print(a[m2])
print('****')
