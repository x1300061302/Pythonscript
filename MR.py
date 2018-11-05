#This is the script for Magnetic Reconnection 
#In 2D-simulation
import scipy.integrate as integrate
import numpy as np 
from const import *

class Simu_info():
    '''
    This Class is aimed to get the information of the Simulation.
    And normalized parameter in Simulations.
    '''
    def __init__(self,deckfile,name = ''):
        self.name = name;
        self.T0 = 0;
        self.di = 0
        self.ne = 0
        self.omega_pe = 0
        self.omega_ce = 0
        self.B0 = 0;
        self.E0 = 0;
        self.J0 = 0;
        self.nx = 0;
        self.ny = 0;
        self.dt = 20;
        self.drift_V = 0;
        print('Begin Read Deck')
        self.deck_read(deck_name = deckfile)
        
        #self.box_info('0000.sdf')]
    def deck_read(self,deck_name):
        f = open(deck_name,'r')
        print('Open file',deck_name);
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
            if (exp[0] == 'temp'):
                self.temp = float(exp[2]);
            if (exp[0] == 'drift_V'):
                self.drift_V = float(exp[2]);
            if (exp[0] == 'kd'):
                self.kd = float(exp[2]);
            ##
            self.E0 = self.B0*c
            self.J0 = self.drift_V*self.ne*qe
                
#             print(exp[0])
        f.close()
        print('Close File')
    def set_dT(self,dt):
        self.dt = dt;
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


class quick_draw(object):
    '''
    I want to achieve quick_draw figure in this class as 
    Q = quick_draw()
    Q.draw_Ekbar()
    Q.draw_photon_Ekbar()
    like this
    '''
    def __init__(self, sdffile,name = '',deckfile='const.status'):
        if name == '':
            self.name = str(sdffile.Header['time']);
        else:
            self.name = name;
        self.a = sdffile;
        self.S = Simu_info(deckfile);  
        self.extent = np.array(sr.Get_extent(self.a))/self.S.di;
        self.para={
            'norm':1,
            'caxis':[0,1],
            'cmap':'jet',
            'xylims':[[self.extent[0],self.extent[1]],[self.extent[2],self.extent[3]]],
            'save':False,
            'axesname': [r'$x/d_i$',r'$y/d_i$',r'$title$'+self.name]
        }

    def get_S(self):
        return self.S;
#    @staticmethod
#     def get_norm(self,key=''):
    def get_para(self):
        return self.para;
    def set_para(self,para):
#         print(para)
        for key in para.keys():
            self.para[key] = para[key];
#             print(key)
#             print(para[key]);
#         return self.para;
    def draw(self,ax,key,para):
        #setting 
        #judge whether the key exist in a. and get data.
        #judge the parameter setting.
        # if default  
        # or not default set_para(para)



        if (type(para) != np.int):
            self.set_para(para);
        if (key == 'Derived_Number_Density_electron'):
            para = {
                    'norm':self.S.ne,
                    'caxis':[0,3],
                    'cmap':'jet',
                    'xylims':[[self.extent[0],self.extent[1]],[self.extent[2],self.extent[3]]],
                    'save':False,
                    'axesname': [r'$x/d_i$',r'$y/d_i$',r'$N_e$'+self.name]
                }
#             key = 
            
        if (key == 'Derived_EkBar_electron'):           
            para = {
                    'norm':self.S.ne,
                    'caxis':[0,3],
                    'cmap':'jet',
                    'xylims':[[self.extent[0],self.extent[1]],[self.extent[2],self.extent[3]]],
                    'save':False,
                    'axesname': [r'$x/d_i$',r'$y/d_i$',r'$N_e$'+self.name]
                }
            
#         self.set_para(para)
        if (type(ax) == np.int):
            fig,ax = df.Create_Figure()
            
        var = self.a.__dict__[key].data;  
#         di = self.S.di;
        ax,gci = df.draw_field_snapshot(ax=ax,\
                                        data=var.T/self.para['norm'],\
                                        extent=self.extent,\
                                        cmap=self.para['cmap'],\
                                        caxis=self.para['caxis']
                                       )
    #                                             )
#         axesname[2] = axesname[2] + self.name;
        df.Axis_set(ax,\
                        axesname=self.para['axesname'],\
                        xylims = self.para['xylims'],\
                       )
        acs = df.Colorbar_set(ax,gci)
        if (self.para['save']):
            plt.savefig(self.para['axesname'][2]+'.png',dpi = 200);
        else:
            plt.show()


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
	
