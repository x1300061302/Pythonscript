#This is the script for Magnetic Reconnection 
#In 2D-simulation
import scipy.integrate as integrate
import numpy as np 
from const import *

class Simu_info():
    '''
    This Class is aimed to get the information of the Simulation.
    '''
    def __init__(self,name):
        self.name = name;
        self.T0 = 0;
        self.di = 0
        self.ne = 0
        self.omega_pe = 0
        self.omega_ce = 0
        self.B0 = 0;
        self.nx = 0;
        self.ny = 0;
        self.deck_read()
        self.box_info('0000.sdf')
    def deck_read(self,deck_name='const.status'):
        f = open(deck_name,'r')
        lines = f.readlines()
        for line in lines:
            exp = line.split()
#             print(exp)
            if (exp[0] == 'T0'):
                self.T0 = float(exp[2]); 
            if (exp[0] == 'di'):
                self.di = float(exp[2]);
            if (exp[0] == 'ne'):
                self.ne = float(exp[2]);
            if (exp[0] == 'omega_pe'):
                self.omega_pe = float(exp[2]);
            if (exp[0] == 'B0'):
                self.B0 = float(exp[2]);
            if (exp[0] == 'omega_ce'):
                self.omega_ce = float(exp[2]);
            if (exp[0] == 'nx'):
                self.nx = float(exp[2]);
            if (exp[0] == 'ny'):
                self.ny = float(exp[2]);
#             print(exp[0])
        f.close()
        
    def box_info(self,sdfname):      
        a = sdf.read(sdfname);
        extent = sr.Get_extent(a)
        self.xmin,self.xmax,self.ymin,self.ymax = np.array(extent)/self.di
#         print(xmin,xmax,ymin,ymax)
        if (self.nx == 0):
            self.nx = len(a.Grid_Grid_mid.data[0])
            self.ny = len(a.Grid_Grid_mid.data[1])
        self.x = np.linspace(self.xmin,self.xmax,self.nx)
        self.y = np.linspace(self.ymin,self.ymax,self.ny)
        self.xx,self.yy = np.meshgrid(self.x,self.y)

class quick_draw():
    
    def __init__(self,sdfname):
        self.name = sdfname;
        self.a = sdf.read(sdfname);
        
    
    def draw_Ekbar(ax):
        Ekbar = 
        ax,gci = df.draw_field_snapshot(ax=ax,\
                                                data=Ekbar.T/1e-12,\
                                                extent=np.array(extent)/di,\
                                                cmap='jet',\
                                        caxis=[0,1]
                                       )
    #                                             )
        df.Axis_set(ax,\
                        axesname=['$x/d_i$','$y/d_i$','$Ekbar$ at '+ str((i)*dT)+'$T_0$'],\
                        xylims = xylims,\
                       )
        acs = df.Colorbar_set(ax,gci)
    #     plt.show()
        return acs

    def draw_photon_Ekbar(ax):
        ax,gci = df.draw_field_snapshot(ax=ax,\
                                                data=pho_Ekbar.T/1e-12,\
                                                extent=np.array(extent)/di,\
                                                cmap='jet',\
                                        caxis=caxis
                                       )
    #                                             )
        df.Axis_set(ax,\
                        axesname=['$x/d_i$','$y/d_i$','$Ekbar$ at '+ str((i)*dT)+'$T_0$'],\
                        xylims = xylims,\
                       )
        df.Colorbar_set(ax,gci)
    #     plt.show()


    def draw_Bx(ax):
        ax,gci = df.draw_field_snapshot(ax=ax,\
                                                data=Bx.T,\
                                                extent=np.array(extent)/di,\
                                                cmap='seismic',\
                                       )
    #                                             )
        df.Axis_set(ax,\
                        axesname=['$x/d_i$','$y/d_i$','$B_x$ at '+ str((i)*dT)+'$T_0$'],\
    #         #             xylims = xylims,\
                       )
        df.Colorbar_set(ax,gci)
        plt.show()
        #plt.savefig('bz'+str(i)+'.eps')


        #contour
        # fig,ax = df.Create_Figure()
        # ax.contour(Bx.T)

        #Jx
    def draw_Jx(ax):
        ax,gci = df.draw_field_snapshot(ax=ax,\
                                                data=Jx.T/Jz0,\
                                                extent=np.array(extent)/di,\
                                                cmap='seismic',\
                                            caxis=[-1,1]
                                           )
        df.Axis_set(ax,\
                        axesname=['$x/d_i$','$y/d_i$','$J_x$ at '+ str((i)*dT)+'$T_0$'],\
                          xylims = xylims,\
                       )
        df.Colorbar_set(ax,gci)
        #Jy
    def draw_Jy(ax):
        ax,gci = df.draw_field_snapshot(ax=ax,\
                                                data=Jy.T/Jz0,\
                                                extent=np.array(extent)/di,\
                                                cmap='seismic',\
                                            caxis=[-1,1]
                                           )
        df.Axis_set(ax,\
                        axesname=['$x/d_i$','$y/d_i$','$J_y$ at '+ str((i)*dT)+'$T_0$'],\
                          xylims = xylims,\
                       )
        df.Colorbar_set(ax,gci)



        #Jz
    def draw_Jz(ax):
        ax,gci = df.draw_field_snapshot(ax=ax,\
                                                data=Jz.T/Jz0,\
                                                extent=np.array(extent)/di,\
                                                cmap='seismic',\
                                            caxis=[-1,1]
                                           )
        df.Axis_set(ax,\
                        axesname=['$x/d_i$','$y/d_i$','$J_z$ at '+ str((i)*dT)+'$T_0$'],\
                          xylims = xylims,\
                       )
        df.Colorbar_set(ax,gci)

        ##Bz
    def draw_Bz(ax):
        ax,gci = df.draw_field_snapshot(ax=ax,\
                                                data=Bz.T,\
                                                extent=np.array(extent)/di,\
                                                cmap='seismic')
        df.Axis_set(ax,\
                        axesname=['$x/d_i$','$y/d_i$','$B_z$ at '+ str((i)*dT)+'$T_0$'],\
                        xylims = xylims,\
                       )
        df.Colorbar_set(ax,gci)
        ###By
    def draw_By(ax):
        ax,gci = df.draw_field_snapshot(ax=ax,\
                                                data=By.T,\
                                                extent=np.array(extent)/di,\
                                                cmap='seismic')
        df.Axis_set(ax,\
                        axesname=['$x/d_i$','$y/d_i$','$B_y$ at '+ str((i)*dT)+'$T_0$'],\
                        xylims = xylims,\
                       )
        df.Colorbar_set(ax,gci)
        #Ez 
    def draw_Ez(ax):
        ax,gci = df.draw_field_snapshot(ax=ax,\
                                                data=Ez.T/Ez0,\
                                                extent=np.array(extent)/di,\
                                                cmap='seismic',\
                                                caxis=[-10,10])
        df.Axis_set(ax,\
                        axesname=['$x/d_i$','$y/d_i$','$E_z$ at '+ str((i)*dT)+'$T_0$'],\
                        xylims = xylims,\
                       )
        df.Colorbar_set(ax,gci)



        ####Nume
    def draw_Nume(ax):
        ax,gci = df.draw_field_snapshot(ax=ax,\
                                                data=Nume.T/ne,\
                                                extent=np.array(extent)/di,\
                                                cmap='jet',\
                                                caxis = [0,3],\
                                                )
        df.Axis_set(ax,\
                        axesname=['$x/d_i$','$y/d_i$','$Num_e$ at '+ str((i)*dT)+'$T_0$'],\
                        xylims = xylims,\
                       )
        axc = df.Colorbar_set(ax,gci)
            # axc.set_clim([0,1.5])
        #     plt.show()
    def draw_pho_Nume(ax):
        ax,gci = df.draw_field_snapshot(ax=ax,\
                                                data=Nume_pho.T/ne,\
                                                extent=np.array(extent)/di,\
                                                cmap='jet',\
                                                caxis = caxis,\
                                                )
        df.Axis_set(ax,\
                        axesname=['$x/d_i$','$y/d_i$','$Num_pho_e$ at '+ str((i)*dT)+'$T_0$'],\
                        xylims = xylims,\
                       )
        axc = df.Colorbar_set(ax,gci)
            # axc.set_clim([0,1.5])
        #     plt.show()

    def draw_streamplot(ax):
        ax.streamplot(x,y,Bx.T,By.T,\
                          density = 1,\
                          linewidth = 3.0,\
                          cmap = 'jet')
    def draw_contour(ax):
        C = ax.contour(xx,yy,np.sqrt(Bx.T**2+By.T**2)/B0,\
                           extent = np.array(extent)/di,\
                           linewidths = 3.0,\
                           cmap='Oranges_r')
        df.Axis_set(ax,axesname=['$x/d_l$','$y/d_l$','Field Line at '+ str(i*dT)+'$T_0$'],\
                    xylims = xylims)
        plt.clabel(C,fontsize=20)

    def draw_De(ax,xylims):
        ax,gci = df.draw_field_snapshot(ax=ax,\
                                                data=De.T,\
                                                extent=np.array(extent)/di,\
                                                cmap='jet',\
                                                caxis = [-1e41,1e41],\
                                                )
        ax.streamplot(x,y,Bx.T,By.T,density=2,linewidth = 3.0)
        xylims = xylims
        df.Axis_set(ax,\
                        axesname=['$x/d_i$','$y/d_i$',r'$D_e$ at '+ str((i)*dT)+'$T_0$'],\
                        xylims = xylims,\
                        legend = False
                       )
        axc = df.Colorbar_set(ax,gci)
    def draw_Temp(ax):
        ax,gci = df.draw_field_snapshot(ax=ax,\
                                                data=Temp.T/1e12,\
                                                extent=np.array(extent)/di,\
                                                cmap='jet',\
                                                caxis = [0,5],\
                                                )
        #     ax.streamplot(x,y,Bx.T,By.T)
        #     xylims = xylims
        df.Axis_set(ax,\
                        axesname=['$x/d_i$','$y/d_i$',r'$Temp$ at '+ str((i)*dT)+'$T_0$'],\
                        xylims = xylims,\
                        legend = False
                       )
        axc = df.Colorbar_set(ax,gci)

        axc.set_clim([0,5])
    def Get_De(a):
        Jx = a.Current_Jx_averaged.data;
        Jy = a.Current_Jy_averaged.data;
        Jz = a.Current_Jz_averaged.data
        Bx = a.Magnetic_Field_Bx_averaged.data
        By = a.Magnetic_Field_By_averaged.data
        Bz = a.Magnetic_Field_Bz_averaged.data
        Ex = a.Electric_Field_Ex_averaged.data
        Ey = a.Electric_Field_Ey_averaged.data
        Ez = a.Electric_Field_Ez_averaged.data
        Nume = a.Derived_Number_Density_electron.data;

        nx,ny = Nume.shape
    #     print(nx,ny)
        extent = sr.Get_extent(a)
        xmin,xmax,ymin,ymax = np.array(extent)/di
    #     print(xmin,xmax,ymin,ymax)
        x = np.linspace(xmin,xmax,nx)
        y = np.linspace(ymin,ymax,ny)
        xx,yy = np.meshgrid(x,y)
        Vx = Jx/Nume/qe
        Vy = Jy/Nume/qe
        Vz = Jz/Nume/qe

        Dex = Jx*(Ex+Vy*Bz-Vz*By)
        Dey = Jy*(Ey+Vz*Bx-Vx*Bz)
        Dez =  Jz*(Ez + Vx*By - Vy*Bx)

        VdotE = Vx*Ex + Vy*Ey + Vz*Vz
        De = Dex+ Dey + Dez
        return De





class MR_Calc(): 
    ''' This is Magnetic Reconnection Calculation Class aimed to calculate the physics parameters and draw 
    some figure quickly
    
    '''
    def __init__(self):
        print('MR_calc')
    def print_sdfname(self):
        print(self.sdfname)
    def box(self):
        '''This function is aimed to get the box information 2d 
        include: 
        xmin,xmax,ymin,ymax
        extent
        
        '''
        a = sdf.read(self.sdfname)
        
        extent = sr.Get_extent(a)
        
        return extent
        
        
        
    def Magnetic_Flux(self,x,y,Bx,By,region):
        '''B should be 2-D arrayddddd
           x,y is the box setup   
           region = [xmin,xmax,ymin,ymax] of your integral region
        '''
        nx,ny = Bx.shape;
        if (nx != len(x)):
            print("The shape nx,ny is",nx,ny,', but the length x,y is ',len(x),len(y))
            return 0;
        #return cell_x and cell_y 
        dx = x[1] - x[0];
        dy = y[1] - y[0];
        xmin,xmax,ymin,ymax = region;
        index_x = np.array(x < xmax)*np.array(x > xmin); 
        index_y = np.array(y < ymax)*np.array(y > ymin); 
        cell_x0 = int((xmin - x[0])/dx)
        cell_y0 = int((ymin - y[0])/dy)
        #print(cell_x0,cell_y0)
        #Flux_y = int^{y_max}_{y_min} Bx dy
        Bx2 = Bx[cell_x0,:];
        By2 = By[:,cell_y0];
        Flux_y = np.trapz(Bx2[index_y],y[index_y])
        Flux_x = np.trapz(By2[index_x],x[index_x])
        return Flux_x,Flux_y        
    
    
    


    


def calc_gradient(dx,dy,data):
    nx,ny = data.shape;
    gx = np.zeros(data.shape);
    gy = np.zeros(data.shape);
    for i in range(0,nx):
        for j in range(0,ny):
            if (i == 0):
                gx[i,j] = (data[i+1,j]-data[i,j])/dx;
            elif (i == nx-1):
                gx[i,j] = (data[i,j]-data[i-1,j])/dx;
            else:
                gx[i,j] = (data[i+1,j]-data[i-1,j])/2/dx;
            if (j == 0):
                gy[i,j] = (data[i,j+1]-data[i,j])/dy;
            elif (j == ny-1):
                gy[i,j] = (data[i,j]-data[i,j-1])/dy;
            else:  
                gy[i,j] = (data[i,j+1]-data[i,j-1])/2/dy;
    return gx,gy
    
def Magnetic_Flux(x,y,Bx,By,region):
        '''B should be 2-D arrayddddd
           x,y is the box setup   
           region = [xmin,xmax,ymin,ymax] of your integral region
        '''
        nx,ny = Bx.shape;
        if (nx != len(x)):
            print("The shape nx,ny is",nx,ny,', but the length x,y is ',len(x),len(y))
            return 0;
        #return cell_x and cell_y 
        dx = x[1] - x[0];
        dy = y[1] - y[0];
        xmin,xmax,ymin,ymax = region;
        index_x = np.array(x < xmax)*np.array(x > xmin); 
        index_y = np.array(y < ymax)*np.array(y > ymin); 
        cell_x0 = int((xmin - x[0])/dx)
        cell_y0 = int((ymin - y[0])/dy)
        #print(cell_x0,cell_y0)
        #Flux_y = int^{y_max}_{y_min} Bx dy
        Bx2 = Bx[cell_x0,:];
        By2 = By[:,cell_y0];
        Flux_y = np.trapz(Bx2[index_y],y[index_y])
        Flux_x = np.trapz(By2[index_x],x[index_x])
        return Flux_x,Flux_y    
    
def test():	
    #test
    x = np.linspace(-10,10,1000);   
    y = np.linspace(-10,10,1000);
    region = [0,5,0,5]
    xx,yy = np.meshgrid(x,y)
    Bx = np.tanh(yy/5);
    By = np.zeros(Bx.shape)
    print(Magnetic_Flux(x,y,Bx,By,region))
	
