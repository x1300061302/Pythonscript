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
    def __init__(self,deckfile='const.status',Nx=0,Ny=0,Nz=0,name = ''):
        self.name = name;
        self.Const={};
        self.axis={};
        self.map={}
        print('Begin Read Deck')
        self.deck_read(deck_name = deckfile)
        self.Nx = Nx;
        self.Ny = Ny;
        self.axis['x'] = np.linspace(self.Const['xmin']/self.Const['di'],self.Const['xmax']/self.Const['di'],self.Nx);
        self.axis['y'] = np.linspace(self.Const['ymin']/self.Const['di'],self.Const['ymax']/self.Const['di'],self.Ny);
        self.map['xx'],self.map['yy'] = np.meshgrid(self.axis['x'],self.axis['y']);  
    def deck_read(self,deck_name):
        f = open(deck_name,'r')
        print('Open file',deck_name);
        lines = f.readlines()
        for line in lines:
            exp = line.split()
#             print(exp)
            self.Const[exp[0]] = float(exp[2])
            ##
        self.Const['E0'] = self.Const['B0']*c
        self.Const['J0'] = self.Const['drift_V']*self.Const['ne']*qe
                
#             print(exp[0])
        f.close()
        print('Close File')
    def get_extent(self,sdffile):
        a = sdf.read(sdffile);
        extent = sr.Get_extent(a);
        return np.array(extent)

class quick_draw(object):
    '''
    I want to achieve quick_draw figure in this class as 
    Q = quick_draw()
    Q.draw_Ekbar()
    Q.draw_photon_Ekbar()
    like this
    '''
    def __init__(self, sdffile,name = '',deckfile='const.status'):
        a = sdf.read(sdffile)
        if name == '':
            self.name = str(a.Header['time']);
        else:
            self.name = name;
        self.a = a;
        self.S = Simu_info(deckfile,\
                           Nx=len(a.Grid_Grid.data[0])-1,\
                           Ny=len(a.Grid_Grid.data[1])-1,\
                          );  
        self.extent = np.array(sr.Get_extent(self.a))/self.S.Const['di'];
        self.para={
            'norm':1,
            'caxis':0,
            'cmap':'jet',
            'xylims':[[self.extent[0],self.extent[1]],[self.extent[2],self.extent[3]]],
            'save':False,
            'axesname': [r'$x/d_i$',r'$y/d_i$',r'$title$'+str(a.Header['time'])],
            'density': 1,
            'linewidth': 2.0
        }
        self.default_para = self.para;
    def get_S(self):
        return self.S;
#    @staticmethod
#     def get_norm(self,key=''):
    def get_para(self):
        return self.para;
    def set_para(self,para):
        self.para.update(para)
#         return self.para;
    def draw_MagneticLine(self,ax=0,para={}):
        self.para.update(self.default_para)
        if (ax == 0):
            fig,ax = df.Create_Figure();
        self.set_para(para);
        Bx = self.a.Magnetic_Field_Bx_averaged.data;
        By = self.a.Magnetic_Field_By_averaged.data;
        ax.streamplot(self.S.axis['x'],self.S.axis['y'],Bx.T,By.T,\
                      density = self.para['density'],\
                      linewidth = self.para['linewidth'],\
                      cmap = self.para['cmap'])
    def draw_spectrum(self,ax,key, para={}, weight = 0):
        key = key.split('.')[1]
        self.para.update(self.default_para)
        if (ax == 0):
            fig,ax = df.Create_Figure();
        var = self.a.__dict__[key]; 
        speciesname = var.name.split('/')[3];# species
        print(np.min(var.data),np.max(var.data))
        self.set_para(para={'axesname':['E/Mev','dN/dE',speciesname+str(self.a.Header['time'])],\
                            'xylims':0});
        if (type(weight) == np.int):
            df.draw_spectrum(ax=ax,data = var.data);
#         else:
#             keyw =  var.name.split('/');
#             weight = self.a.dict__[key]
        df.Axis_set(ax,\
                        axesname=self.para['axesname'],\
                        xylims = self.para['xylims'],\
                       )
        
    def draw(self,ax,key,para={}):
        '''
        Input: 
            ax is the axis of the figure 
            key is the key value of sdf file, like key = 'Derived_Number_Density_electron'
            para is the dictionary paramater setting which includes:
                'norm':1,
                'caxis':0,
                'cmap':'jet',
                'xylims':[[self.extent[0],self.extent[1]],[self.extent[2],self.extent[3]]],
                'save':False,
                'axesname': [r'$x/d_i$',r'$y/d_i$',r'$title$'+self.name]
        Output:
            Return Fig,ax 
            Figure
        '''
                
        #setting 
        #judge whether the key exist in a. and get data.
        #judge the parameter setting.
        # if default  
        # or not default set_para(para)
#        if (key == 'Derived_Number_Density_electron'):
#            para = {
#                    'norm':self.S.Const['ne'],
#                    'caxis':[0,3],
#                    'cmap':'jet',
#                    'xylims':[[self.extent[0],self.extent[1]],[self.extent[2],self.extent[3]]],
#                    'save':False,
#                    'axesname': [r'$x/d_i$',r'$y/d_i$',r'$N_e$'+self.name]
#                }
##             key = 
#            
#        if (key == 'Derived_EkBar_electron'):           
#            para = {
#                    'norm':self.S.Const['ne'],
#                    'caxis':[0,3],
#                    'cmap':'jet',
#                    'xylims':[[self.extent[0],self.extent[1]],[self.extent[2],self.extent[3]]],
#                    'save':False,
#                    'axesname': [r'$x/d_i$',r'$y/d_i$',r'$N_e$'+self.name]
#                }
#            
#         self.set_para(para)
        key = key.split('.')[1]
        if (type(ax) == np.int):
            fig,ax = df.Create_Figure()
            
        var = self.a.__dict__[key].data;  
        vmin = np.min(var);
        vmax = np.max(var);
        #self.para.update(self.default_para)
        self.para['caxis'] = 0; # autosetting vmin - vmax;
        
        if (type(para) != np.int):
            self.set_para(para);
#         di = self.S.di;
        print(self.para['caxis']);
        ax,gci,cb = df.draw_field_snapshot(ax=ax,\
                                        data=var.T/self.para['norm'],\
                                        extent=self.extent,\
                                        cmap=self.para['cmap'],\
                                        caxis=self.para['caxis'],\
                                       )
    #                                             )
        df.Axis_set(ax,\
                        axesname=self.para['axesname'],\
                        xylims = self.para['xylims'],\
                       )
        #autoseting colorcaxis:
        if (self.para['save']):
            plt.savefig(self.para['axesname'][2]+'.png',dpi = 200);

        return ax 


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
	
