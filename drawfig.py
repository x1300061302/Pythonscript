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

def Create_Figure2(n,m,fw=10,fh=10,sp=2,sw=8,sh=8):
    '''
    default fw = 10, fh = 10, sp = 2, sw = sh = 8
    '''
    import matplotlib.pyplot as plt
    spw = sp/fw;
    sph = sp/fh;
    rsw = sw/fw;
    rsh = sh/fh;
    ax_list = []
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
    if (ax_style == 'd'): 
        if (type(xylims) != np.int):
            plt.xlim(xylims[0]);
            plt.ylim(xylims[1]);
        # x-axis
        plt.xlabel(axesname[0],fontsize=fs); 
        plt.xticks(fontsize=fs);
        #y-axis
        plt.ylabel(axesname[1],fontsize=fs);
        plt.yticks(fontsize=fs);
	#title
        plt.title(axesname[2],fontsize=fs);
        plt.grid(grid);
        if (legend):
            plt.legend(fontsize =fs);
	
        ax.spines['bottom'].set_linewidth(lw)
        ax.spines['left'].set_linewidth(lw)
        ax.spines['right'].set_linewidth(lw)
        ax.spines['top'].set_linewidth(lw)
        ax.tick_params(which = 'both',width = lw,colors = 'black')
        ax.tick_params(direction = 'in')
        ax.tick_params(length  = ticklength)
        if (showtick):
            ax.tick_params(top = True,right = True)
		
    elif (ax_style == 'p'):
        if (type(xylims) != np.int):
            ax.set_rlim(xylim[1])
        # x-axis
        plt.xlabel(axesname[0],fontsize=fs); 
        plt.xticks(fontsize=fs);
        #y-axis
        plt.ylabel(axesname[1],fontsize=fs);
        plt.yticks(fontsize=fs);
        #title
        plt.title(axesname[2],fontsize=fs);
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

def Legend_outside(ax,loc,bbox=(1.35,1.0),fs=20):
    box = ax.get_position();
    ax.set_position([box.x0,box.y0,box.width,box.height]);
    ax.legend(loc = loc,bbox_to_anchor=bbox,fontsize=fs)
    

#***********Draw----------------$$$$$$$$$$$$$    
def draw_histogram2d(ax, x, y, bins=[100, 200], xylim=0, caxis=0, fontsize=20,cmap = 'jet'):
    '''to draw 2-D histogram figure 
       handle = ax
       data = [x,y]
       bins =[nx,ny]
       xylim = [[xlim],[ylim]] means the region to show
       caxis = region of colobar
       '''
    import matplotlib.pyplot as plt
 
    h_xy, xedge, yedge = np.histogram2d(x, y, bins=bins)
    log10h = np.log10(h_xy)

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


def draw_spectrum(ax, data, label, weights=0, logx=0, grid=True,gethd=0):
    ''' ax should be input
    data is 1-D  array
    label = name of legend 
    '''
    import matplotlib.pyplot as plt
    if (type(weights) == np.int):
        print('No Weights hist')
        hd, ad = np.histogram(data, bins=500, normed=False)
        hd = hd/(ad[1]-ad[0])
        axx = 0.5*(ad[1:]+ad[:-1]);
        if logx == 1:
            ad = np.log10(ad)
        pl = plt.semilogy(axx, hd,\
                          color=MYLinecolor[0],label=label, linewidth=2)
    else:
        print('Weights hist')
        hd, ad = np.histogram(data, bins=500, normed=False, weights=weights)
        hd = hd/(ad[1]-ad[0])
        axx = 0.5*(ad[1:]+ad[:-1]);
        if logx == 1:
            ad = np.log10(ad)
        pl = plt.semilogy(axx, hd,\
                          color=MYLinecolor[0], label=label, linewidth=2)
    return hd,axx


def draw_spectrum_nline(ax,data, label=0, weights=0, lw=2.0, logx=0, cl=MYLinecolor, bins=500, normed=False):
    '''data is 2-D x numl array
    dataname is numl length array ['','',''] 
    figname default = fig
    make sure numl is given'''
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
            if logx == 1:
                ad = np.log10(ad)
            plsy = plt.semilogy(.5*(ad[1:]+ad[:-1]), hd,
                                color=cl[i],
                                label=label[i],
                                linewidth=lw)
    else:
        print('Weights hist')
        for i in range(0, numl):
            hd, ad = np.histogram(
                data[i][:], bins=bins, normed=normed, weights=weights[i][:])
            if logx == 1:
                ad = np.log10(ad)
            plsy = plt.semilogy(.5*(ad[1:]+ad[:-1]), hd,
                                color=cl[i], label=label[i], linewidth=lw)


def draw_field_snapshot(ax, data, extent, caxis=0,cmap='jet'):
    '''data should be 2-D array (nx,ny) 
    extent should be ([x_min,x_max],[y_min,y_max])
    xylim = ([xmin,xmax],[ymin,ymax])
    label = ['xlabel/Unit','ylabel/Unit','title'] 
    Display dicide whether to show
    figname default = title'''
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
    # colorbar
    return ax,gci;


def draw_line(ax, xx, data, label, lw=2.0):
    ax.plot(xx, data, 'r-', label=label, linewidth=lw)


def draw_nline(xx, data, label, cl=MYLinecolor, lw=2.0):
    ''' xx = nlxlength np.array datatype 
    '''
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from mpl_toolkits.mplot3d import Axes3D
    nl = len(data)
    for i in range(0, nl):
        ax.plot(xx, data[i][:], cl[i], label=label[i], linewidth=lw)


#**************main test
