#### Calculate Ohm component
import matplotlib
matplotlib.use('Agg')
from const import *
from scipy import ndimage

# cell_y0 = Ny//2d
# imp.reload(sr)
dirc='qed-2-nop/'
sdffile = '0060.sdf'
# Nx,Ny = [500,500]
Nx,Ny = [1000,1000]

mratio = 1
betau = 1 + 1/mratio;
betam = 1 - 1/mratio;


def calc_djzdt(sdf1,sdf2):
    a = sdf.read(sdf1)
    time1 = a.Header['time']
    fv1 = a.Derived_Four_Velocity_Uz_electron.data;
    gfv1 = ndimage.gaussian_filter(fv1,sigma=1.0)
    nmm1 = a.Derived_Number_Density.data/2;
    
    a = sdf.read(sdf2)
    fv0 = a.Derived_Four_Velocity_Uz_electron.data;
    gfv0 = ndimage.gaussian_filter(fv0,sigma=1.0)
    
    nmm0 = a.Derived_Number_Density.data/2;
    time0 = a.Header['time']

    djzdt = -(gfv1-gfv0)/(time1-time0)*me/qe
    return djzdt

cell_y0 = Ny//4
cell_x0 = Nx//2
dg =  8
region = [20,Nx-20,20,Ny-20]
q = df.quick_draw(dirc = dirc,sdffile=sdffile,Nx=Nx,Ny=Ny)
s = q.s;
a = q.a;
va = q.s.const['B0']/np.sqrt(q.s.const['ne']*s.const['mi']*mu0)
# va = c
dxx = 2*s.const['xmax']/Nx;
dyy = dxx;
num_e = a.Derived_Number_Density_electron;
num = a.Derived_Number_Density;
num_p = a.Derived_Number_Density_positron;

bx = a.Magnetic_Field_Bx
by = a.Magnetic_Field_By
bz = a.Magnetic_Field_Bz
try:
    bx_avg = a.Magnetic_Field_Bx_averaged
    by_avg = a.Magnetic_Field_By_averaged
    bz_avg = a.Magnetic_Field_Bz_averaged
    convention_ez = a.Derived_Convection_Ez_averaged
    ez_avg = a.Electric_Field_Ez_averaged
except:
    print('---Has no averaged ')

jz = a.Current_Jz
jx = a.Current_Jx
jy = a.Current_Jy
jx_e = a.Derived_Jx_electron
jy_e = a.Derived_Jy_electron
jz_e = a.Derived_Jz_electron
jz_p = a.Derived_Jz_positron
temp = a.Derived_Temperature_electron

# print(ne.data.shape)
###uu term 
uui = {}
uui['xz'] = a.Derived_Velocity_Ux_positron.data*a.Derived_Four_Velocity_Uz_positron.data
uui['yz'] = a.Derived_Velocity_Uy_positron.data*a.Derived_Four_Velocity_Uz_positron.data

uue = {}
uue['xz'] = a.Derived_Velocity_Ux_electron.data*a.Derived_Four_Velocity_Uz_electron.data
uue['yz'] = a.Derived_Velocity_Uy_electron.data*a.Derived_Four_Velocity_Uz_electron.data
#############-draw-------------
## dg = 5
curldg = 5
### Ez
ez_avg_ml = sr.calc_mean_var_line(data = ez_avg.data,cell=cell_y0,delta_grid=dg,direction='x')
####convection term vxb
convention_ez_ml = sr.calc_mean_var_line(cell=cell_y0,delta_grid=dg,data=convention_ez.data,direction='x') #/2 for num_e
#####non-diagnal term
Pre = {}
Pre['xz'] = a.Derived_Pressure_Pxz_electron.data   
Pre['yz'] = a.Derived_Pressure_Pyz_electron.data  
####non-diagnal term for i 
Pri = {}
Pri['xz'] = a.Derived_Pressure_Pxz_positron.data;
Pri['yz'] = a.Derived_Pressure_Pyz_positron.data;
###calculate the derivation
dPre = sr.calc_dot_2d(Ax=Pre['xz'],Ay = Pre['yz'],dx = dxx,dy = dyy,delta_nx=curldg,delta_ny=curldg,region=region)*1/(betau*qe*(num.data/2))
g_dPre = ndimage.gaussian_filter(dPre,sigma=5.0)
dPre_ml = sr.calc_mean_var_line(data = g_dPre,cell=cell_y0,delta_grid=dg,direction='x')

dPri = sr.calc_dot_2d(Ax=Pri['xz'],Ay = Pri['yz'],dx = dxx,dy = dyy,delta_nx=curldg,delta_ny=curldg,region=region)*1/(betau*qe*(num.data/2))
g_dPri = ndimage.gaussian_filter(dPri,sigma=5.0)
dPri_ml = sr.calc_mean_var_line(data = g_dPri,cell=cell_y0,delta_grid=dg,direction='x')
Pre_term = dPri_ml/mratio - dPre_ml
# plt.plot(dPre_ml)
#### juuj
juuj = sr.get_juuj(a)
juuj_ez = sr.calc_dot_2d(Ax = juuj['xz'],Ay = juuj['yz'],dx = dxx,dy = dyy,delta_nx=curldg,delta_ny=curldg,region=region)*me/qe**2/(num.data/2)
g_juuj = ndimage.gaussian_filter(juuj_ez,sigma=2.0)
juuj_ez_ml = sr.calc_mean_var_line(data = g_juuj,cell=cell_y0,delta_grid=dg,direction='x')
####hall effect
jxb = betam/betau*(jx.data*by_avg.data - jy.data*bx_avg.data)/(num.data/2)/qe #/总的数密度
g_jxb = ndimage.gaussian_filter(jxb,sigma = 2.0)

# jxb = sr.calc_cross(Ax=jx.data,Ay = jy.data,Bx = bx.data,By=by.data)/num.data/qe
jxb_ml = sr.calc_mean_var_line(data = g_jxb,cell=cell_y0,delta_grid=dg,direction='x')
djzdt = calc_djzdt(dirc + '0070.sdf',dirc+'0060.sdf');
g_djzdt = ndimage.gaussian_filter(djzdt,sigma=2.0)
djzdt_ml = sr.calc_mean_var_line(data = g_djzdt,cell=cell_y0,delta_grid=dg,direction='x')
# djdt = a.Current_Jz_derivated.data
# gdjdt = ndimage.gaussian_filter(djdt,sigma = 2.0);
# djdt_ml = sr.calc_mean_var_line(data = gdjdt,cell=cell_y0,delta_grid=dg,direction='x')
###uui term 
uudg = 5
numm =  (num_e.data + num_p.data)/2
nme = num_e.data 
nmi = num_p.data
uui_ez = sr.calc_dot_2d(Ax = nmi*uui['xz'],\
                       Ay = nmi*uui['yz'],\
                       dx=dxx,dy=dyy,\
                       delta_nx = curldg,delta_ny=curldg,\
                       region=region)*me/(num.data/2)/qe
####
uue_ez = sr.calc_dot_2d(Ax = nme*uue['xz'],\
                       Ay = nme*uue['yz'],\
                       dx=dxx,dy=dyy,\
                       delta_nx = curldg,delta_ny=curldg,\
                       region=region)*me/(num.data/2)/qe
uui_uue = ndimage.gaussian_filter(uui_ez-uue_ez,sigma = 4.0)
uu_ml = 1/betau*sr.calc_mean_var_line(cell=cell_y0,data=uui_uue,delta_grid=uudg,direction='x')
# inertial_term = uu_ml
#For relativistic term v\dot nabla u
uz = a.Derived_Four_Velocity_Uz_electron.data;
guz = ndimage.gaussian_filter(uz,sigma=2.0)
gx,gy = sr.calc_gradient(data = guz,delta_n = curldg,delta_ny=curldg,dim=2,dx = dxx,dy =dyy,region=region)
vx = a.Derived_Velocity_Ux_electron.data;
vy = a.Derived_Velocity_Uy_electron.data;
vdu = -(vx*gx+vy*gy)*me/qe/2

gvdu = ndimage.gaussian_filter(vdu,sigma=4.0)
gvdu_ml = sr.calc_mean_var_line(cell=cell_y0,data=gvdu,delta_grid=uudg,direction='x')
inertial_term = gvdu_ml

#----------------X direction ----------- #########
xx = s.axis['x']
fig,ax = df.Create_Figure([6,4])
ax.plot(xx,ez_avg_ml/s.const['E0'])
ax.plot(xx,convention_ez_ml/s.const['E0'])
ax.plot(xx,Pre_term/s.const['E0'])
ax.plot(xx,djzdt_ml/s.const['E0'],label = 'djz/dt')
ax.plot(xx,inertial_term/s.const['E0'],label = 'inertial term')

ez_cal = convention_ez_ml + Pre_term  + inertial_term + djzdt_ml 
ax.plot(xx,ez_cal/s.const['E0'],label = r'ez_Ohm_calc')

df.Axis_set(ax = ax,axesname=['x/di','$E/B_0c$','Each Component of Ez in X direction'])
plt.legend(['ez_avg','convection_ez','Pre','djzdt','Inertial term','ez_calc'])

plt.savefig(dirc + 'Ohm-X.eps',dpi = 200)

####------------Y direction:
cell_x0 = Nx//2
dg = 20
ez_avg_mly = sr.calc_mean_var_line(data = ez_avg.data,cell=cell_x0,delta_grid=dg,direction='y')
convention_ez_mly = sr.calc_mean_var_line(cell=cell_x0,delta_grid=dg,data=convention_ez.data,direction='y') #/2 for num_e

# dPre = sr.calc_dot_2d(Ax=Pre['xz'],Ay = Pre['yz'],dx = dxx,dy = dyy,delta_nx=curldg,delta_ny=curldg,region=region)*1/(betau*qe*(num.data/2))
# g_dPre = ndimage.gaussian_filter(dPre,sigma=2)
dPre_mly = sr.calc_mean_var_line(data = g_dPre,cell=cell_x0,delta_grid=dg,direction='y')

# dPri = sr.calc_dot_2d(Ax=Pri['xz'],Ay = Pri['yz'],dx = dxx,dy = dyy,delta_nx=curldg,delta_ny=curldg,region=region)*1/(betau*qe*(num.data/2))
# g_dPri = ndimage.gaussian_filter(dPri,sigma=2)
dPri_mly = sr.calc_mean_var_line(data = g_dPri,cell=cell_x0,delta_grid=dg,direction='y')
Pre_y = dPri_mly - dPre_mly

uu_mly = 1/betau*sr.calc_mean_var_line(cell=cell_x0,data=uui_uue,delta_grid=uudg,direction='y')

djzdt_mly = sr.calc_mean_var_line(data = g_djzdt,cell=cell_x0,delta_grid=dg,direction='y')

gvdu_mly = sr.calc_mean_var_line(cell=cell_x0,data=gvdu,delta_grid=uudg,direction='y')

fig,ax = df.Create_Figure([6,4])
ax.plot(xx,ez_avg_mly/s.const['E0'])
ax.plot(xx,convention_ez_mly/s.const['E0'])
ax.plot(xx,Pre_y/s.const['E0'])
# plt.plot(uu_mly)
ax.plot(xx,djzdt_mly/s.const['E0'])
ax.plot(xx,gvdu_mly/s.const['E0'])
ax.plot(xx,(convention_ez_mly+Pre_y+djzdt_mly+gvdu_mly)/s.const['E0'])
df.Axis_set(ax = ax, axesname=[r'$y/d_i$','ez','Each Component of Ez in Y direction'])
plt.legend(['ez_avg','convection_ez','Pre','djzdt','Inertial term','ez_calc'])
plt.savefig(dirc + 'Ohm-Y.eps',dpi = 200)



