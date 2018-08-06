#This is the script for Magnetic Reconnection 
#In 2D-simulation
import scipy.integrate as integrate
import numpy as np 

def Magnetic_Flux(x,y,Bx,By,region):
    '''xx,yy,B should be 2-D array
       region = [xmin,xmax,ymin,ymax]
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
	
