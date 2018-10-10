#This is the script for Magnetic Reconnection 
#In 2D-simulation
import scipy.integrate as integrate
import numpy as np 
from const import *

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
	
