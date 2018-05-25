import sdfread as sr
a = sr.Get_field_variable('f0006.sdf','Ex')

ex = a.data
ex_xy = ex[:,:,96]
print(ex_xy.shape)
ex_yz = ex[200,:,:]
print(ex_yz.shape)
ex_xz = ex[:,96,:]
print(ex_xz.shape)

import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(2,1,1)
print(dir(ax))
ax.imshow(ex_xy.T);
ax.set_aspect(960/192)

plt.xlabel('x')
plt.ylabel('y')

plt.show()
