#!/usr/bin/python
from const import *
import numpy as np
colorline = MYLinecolor
figw, figh = [10,10]
ix,iy = [1,1]
iwidth, iheight= [8,8]


def TextPos(xylims,scale = [0.1,0.9],islog = False):
    if (islog):
        Dy = np.log(xylims[1][1]/xylims[1][0])
        Dx = xylims[0][1] - xylims[0][0]
        posx = xylims[0][0]+ scale[0]*Dx;
        posy = xylims[1][0]+np.exp(scale[1]*(Dy))
    else:
        Dy = xylims[1][1] - xylims[1][0]
        Dx = xylims[0][1] - xylims[0][0]
        posx = xylims[0][0]+ scale[0]*Dx;
        posy = xylims[1][0]+ scale[1]*Dy;
    return posx,posy
def gifmake(gifname='gif.gif', duration=0.1, beg=0, end=0, prefix='.png'):
    import matplotlib.pyplot as plt
    import imageio
    import os
    images = []
    filenames = sorted((fn for fn in os.listdir('.') if fn.endswith(prefix)))
    if (end == 0):
        end = len(filenames)
    for i in range(beg, end):
        filename = filenames[i]
        images.append(imageio.imread(filename))
        imageio.mimsave(gifname, images, duration=duration)
# ----------Create Figure---------$$$$$$$$$$$$

def Create_Figure2(n,m,fw=10,fh=10,sp=2,sw=8,sh=8,ratio=1):
    '''
    default fw = 10, fh = 10, sp = 2, sw = sh = 8
    '''
    import matplotlib.pyplot as plt
    spw = sp/fw;
    sph = sp/fh;
    rsw = sw/fw;
    rsh = sh/fh;
    ax_list = []
    #zoom
    fw = fw*ratio;
    fh = fh*ratio;
    
    fig = plt.figure(figsize= [fw,fh])
    for i in range(0,n):
        for j in range(0,m):
            ax = fig.add_axes([spw*(j+1)+rsw*j,sph*(i+1)+rsh*i,rsw,rsh]) #横排从左往右,#竖排从下往上
            ax_list.append(ax);       
    return ax_list
#     plt.show()

def Create_Figure(figsize=[8,8], x=1, y=1, n=1, polar=False):
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=figsize)
    if (polar):
        ax = fig.add_subplot(x, y, n,projection='polar')
    else:
        ax = fig.add_subplot(x, y, n)
    return fig,ax


# ------------Figure/axis setting------$$$$$$$$$
def Axis_set(ax,axesname=['x','y',''],fs=20.0,xticklabel=0,xtickrange=0,yticklabal=0,ytickrange=0,grid=False,legend=False,xylims = 0,ax_style='d',lw = 2.0,showtick = True, ticklength = 10):
    import matplotlib.pyplot as plt
    if (axesname == 'sp'):   #Quick set 
        axesname = ['E/MeV','dN/dE',''];
        
    if (ax_style == 'd'): 
        if (type(xylims) != np.int):
            ax.set_xlim(xylims[0]);
            ax.set_ylim(xylims[1]);
        # x-axis
        ax.set_xlabel(axesname[0],fontsize=fs); 
#         ax.set_xticks(xticklabel); 
        #y-axis
        ax.set_ylabel(axesname[1],fontsize=fs);
#         ax.set_yticks(fontsize=fs);
	#title
        ax.set_title(axesname[2],fontsize=fs);
        ax.grid(grid);
        if (legend):
            plt.legend(fontsize =fs);
	
        ax.spines['bottom'].set_linewidth(lw)
        ax.spines['left'].set_linewidth(lw)
        ax.spines['right'].set_linewidth(lw)
        ax.spines['top'].set_linewidth(lw)
        ax.tick_params(which = 'both',width = lw,colors = 'black',labelsize=fs)
        ax.tick_params(direction = 'in')
        ax.tick_params(length  = ticklength)
        if (showtick):
            ax.tick_params(top = True,right = True)
		
    elif (ax_style == 'p'):
        if (type(xylims) != np.int):
            ax.set_rlim(xylim[1])
        # x-axis
        ax.set_xlabel(axesname[0],fontsize=fs); 
        ax.set_xticks(fontsize=fs);
        #y-axis
        ax.set_ylabel(axesname[1],fontsize=fs);
        ax.set_yticks(fontsize=fs);
        #title
        ax.set_title(axesname[2],fontsize=fs);
        plt.grid(grid);
        ax.spines['polar'].set_linewidth(lw)
        if (legend):
            plt.legend(fontsize =fs,shadow=True);
    else:
        print('Can not understand ax style'+ax_style);

def Colorbar_set(ax,gci,pos='right',size='3%',pad=0.1):
    import matplotlib.cm as cm
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import matplotlib.pyplot as plt
    divider = make_axes_locatable(ax)
    cax = divider.append_axes(pos, size=size, pad=pad)
    cbar = plt.colorbar(gci, cax)
    #cbar.set_ticks() and other operation
    return cbar

def Legend_outside(ax,loc,ncol = 1,bbox=(1.0,1.0),fs=20):
    box = ax.get_position();
    ax.set_position([box.x0,box.y0,box.width,box.height]);
    ax.legend(loc = loc,ncol = ncol, bbox_to_anchor=bbox,fontsize=fs)
    

#***********Draw----------------$$$$$$$$$$$$$    
def draw_histogram2d(ax, x, y, bins=[300, 300], xylim=0, caxis=0, fontsize=20,cmap = 'jet',log10h =True):
    if (ax == 0):
        fig,ax = Create_Figure();
    '''to draw 2-D histogram figure 
       handle = ax
       data = [x,y]
       bins =[nx,ny]
       xylim = [[xlim],[ylim]] means the region to show
       caxis = region of colobar
       log10h = log10(h)
       '''
    import matplotlib.pyplot as plt
 
    h_xy, xedge, yedge = np.histogram2d(x, y, bins=bins)
    if (log10h):
        log10h = np.log10(h_xy)
    else:
        log10h = h_xy

    if (type(caxis) != np.int):
        vmin = caxis[0]
        vmax = caxis[1]
    else:
        vmin = np.min(log10h)
        vmax = np.max(log10h)
    print('get histogram')
    # axis
    xax = np.linspace(np.min(x), np.max(x), bins[0])
    yax = np.linspace(np.min(y), np.max(y), bins[1])
    xx, yy = np.meshgrid(xax, yax)
    # draw
    gci = ax.pcolormesh(xx, yy, log10h.T, cmap=cmap,
                        vmin=vmin, vmax=vmax)
    return ax,gci
    
    # colorbar

def draw_angle_distribution3d(ax, T, R, tlim=0, rlim=0, binT=360, binR=1000, caxis=0,log10=True,cmap='jet'):
    if (ax == 0):
        fig,ax = Create_Figure(polar=True);
    ''' ax should be polarization
    T,R is N-array 
    rlim = 0 means default else rlim = [r_min,r_max]
    draw angle distribution histogram'''
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    #histogram
    H, Tedge, Redge = np.histogram2d(T, R, bins=[binT, binR])
    if (log10):
        log10H = np.log10(H)
    else:
        log10H = H;
    print("get histogram")
    #change T,R to meshgrid
    rad = np.linspace(R.min(), R.max(), binR)
    if(type(tlim) != np.int):
        the = np.linspace(tlim[0], tlim[1], binT)
    else:
        the = np.linspace(0, 2*np.pi, binT)
    r, t = np.meshgrid(rad, the)
    if (type(caxis) != np.int):
        vmin = caxis[0]
        vmax = caxis[1]
    else:
        vmin = 0
        vmax = log10H.max()
    pcmesh = ax.pcolormesh(t, r, log10H, cmap=cmap, vmin=vmin, vmax=vmax)
    if (type(rlim) != np.int):
        ax.set_rlim(rlim)
    else:
        ax.set_rlim(np.min(R), np.max(R))
    cbar = plt.colorbar(pcmesh)
    return ax,pcmesh

def draw_spectrum(data, ax = 0, label='', weights=0, grid=True,gethd=0,lw=3.0,cl='',factor=1.0):
    ''' Input: data 
               ax = 0 default figure and axis;
               label = '' 
               weights = 0 
               grid =True
               lw means the linewidth of the curve
               cl means the color of the line 
               factor means the factor * hd 
        Output:
           return ax  
    '''
    if (ax == 0):
        fig,ax = Create_Figure();
    import matplotlib.pyplot as plt
    if (cl == ''):
        cl = MYLinecolor[0];
    if (type(weights) == np.int):
        print('No Weights hist')
        hd, ad = np.histogram(data, bins=500, normed=False)
        hd = hd/(ad[1]-ad[0])
        axx = 0.5*(ad[1:]+ad[:-1]);
        pl = plt.semilogy(axx, factor*hd,\
                          color=cl,label=label, linewidth=lw)
    else:
        print('Weights hist')
        hd, ad = np.histogram(data, bins=500, normed=False, weights=weights)
        hd = hd/(ad[1]-ad[0])
        axx = 0.5*(ad[1:]+ad[:-1]);
        pl = plt.semilogy(axx, factor*hd,\
                          color=cl, label=label, linewidth=lw)
    return ax


def draw_spectrum_nline(data, ax = 0, label=0, weights=0, lw=3.0, cl=MYLinecolor,ls = '-', bins=500, normed=False):
    ''' Input:
            data is l x N array
            label is the legend name
            weights = 0 is default meaning no weight
            lw is the width of the line
            cl is the color of the line 
            and ls is the linestyle of the line
            figname default = fig
            make sure numl is given'''
    if (ax == 0):
        fig,ax = Create_Figure();
    numl = len(data)
    if (type(label)==np.int):
        label = []
        for i in range(0,numl):
            label.append(str(i))
    if (display == 0):
        try:
            import matplotlib
            matplotlib.use('Agg')
        except:
            print('use matplotlib.use before import plt')
    import matplotlib.pyplot as plt
    if (type(weights) == np.int):
        print('No Weights hist')
        for i in range(0, numl):
            hd, ad = np.histogram(data[i][:], bins=500, normed=False)
            plsy = plt.semilogy(.5*(ad[1:]+ad[:-1]), hd,
                                color=cl[i],
                                linestyle = ls,
                                label=label[i],
                                linewidth=lw)
    else:
        print('Weights hist')
        for i in range(0, numl):
            hd, ad = np.histogram(
                data[i][:], bins=bins, normed=normed, weights=weights[i][:])
            plsy = plt.semilogy(.5*(ad[1:]+ad[:-1]), hd,
                                color=cl[i], 
                                linestyle = ls,
                                label=label[i], linewidth=lw)
    return ax


def draw_field_snapshot(data, extent, ax = 0, caxis=0,cmap='jet',bar_set=True):
    '''data should be 2-D array (nx,ny) 
    extent should be ([x_min,x_max],[y_min,y_max])
    xylim = ([xmin,xmax],[ymin,ymax])
    label = ['xlabel/Unit','ylabel/Unit','title'] 
    Display dicide whether to show
    figname default = title'''
    if (ax == 0):
        fig,ax = Create_Figure();
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from mpl_toolkits.mplot3d import Axes3D
    nx, ny = data.shape
    if (type(caxis) != np.int):
        vmin = caxis[0]
        vmax = caxis[1]
    else:
        vmin = data.min()
        vmax = data.max()
    gci = ax.imshow(data, extent=extent, origin='lower', cmap=cmap,
                    vmax=vmax, vmin=vmin, interpolation='spline36')
    if (bar_set):
        cb = Colorbar_set(ax,gci);
        return ax,gci,cb;
    
    # colorbar
    return ax,gci;


def plot(xx, data, ax = 0,label='', lw=2.0):
    if (ax == 0):
        fig,ax = Create_Figure();
    ax.plot(xx, data, 'r-', label=label, linewidth=lw)
    return ax


def draw_nline(ax, xx, data, label, cl=MYLinecolor, lw=2.0):
    ''' xx = nlxlength np.array datatype 
    '''
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from mpl_toolkits.mplot3d import Axes3D
    if (ax == 0):
        fig,ax = Create_Figure();
    nl = len(data)
    for i in range(0, nl):
        ax.plot(xx, data[i][:], cl[i], label=label[i], linewidth=lw)
    return ax


#**************main test
