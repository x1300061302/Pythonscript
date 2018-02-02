#!/usr/bin/python
from const import *
import numpy as np
colorline = MYLinecolor


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


def Create_Figure(figsize=[15, 10], x=1, y=1, n=1, polar='False'):
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=[15, 10])
    ax = fig.add_subplot(x, y, n)
    return ax


# ------------Figure/axis setting------$$$$$$$$$
def Figure_setting():
    # x-axis
    import matplotlib.pyplot as plt
    plt.xlabel('x');
    # ***********Draw----------------$$$$$$$$$$$$$


def draw_histogram2d(ax, x, y, label, bins=[100, 200], xylim=0, caxis=0, fontsize=20):
    '''to draw 2-D histogram figure 
       handle = ax
       data = [x,y]
       label = [xlabel,ylabel,title]
       bins =[nx,ny]
       xylim = [[xlim],[ylim]] means the region to show
       caxis = region of colobar
       '''
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    h_xy, xedge, yedge = np.histogram2d(x, y, bins=bins)
    log10h = np.log10(h_xy)
    print('get histogram')
    # axis
    xax = np.linspace(np.min(x), np.max(x), bins[0])
    yax = np.linspace(np.min(y), np.max(y), bins[1])
    xx, yy = np.meshgrid(xax, yax)
    # draw
    gci = ax.pcolormesh(xx, yy, log10h.T, cmap='jet',
                        vmin=0, vmax=np.max(log10h))
    # colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='3%', pad=0.1)
    cbar = plt.colorbar(gci, cax)


def draw_angle_distribution3d(ax, T, R, label=0, tlim=0, rlim=0, binT=360, binR=1000, caxis=0):
    ''' ax should be polarization
    T,R is N-array 
    rlim = 0 means default else rlim = [r_min,r_max]
    draw angle distribution histogram'''
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    #histogram
    H, Tedge, Redge = np.histogram2d(T, R, bins=[binT, binR])
    log10H = np.log10(H)
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
    pcmesh = ax.pcolormesh(t, r, log10H, cmap='jet', vmin=vmin, vmax=vmax)
    if (type(rlim) != np.int):
        ax.set_rlim(rlim)
    else:
        ax.set_rlim(np.min(R), np.max(R))
    cbar = fig.colorbar(pcmesh)


def draw_spectrum(ax, data, label, Linecolor='b', weights=0, logx=0, grid=True):
    ''' ax should be input
    data is 1-D  array
    label = name of legend 
    '''
    import matplotlib.pyplot as plt
    if (type(weights) == np.int):
        print('No Weights hist')
        hd, ad = np.histogram(data, bins=500, normed=False)
        hd = hd/(ad[1]-ad[0])
        if logx == 1:
            ad = np.log10(ad)
        pl = plt.semilogy(.5*(ad[1:]+ad[:-1]), hd,
                          color=color, label=label, linewidth=2)
    else:
        print('Weights hist')
        hd, ad = np.histogram(data, bins=500, normed=False, weights=weights)
        hd = hd/(ad[1]-ad[0])
        if logx == 1:
            ad = np.log10(ad)
        pl = plt.semilogy(.5*(ad[1:]+ad[:-1]), hd,
                          color=color, label=label, linewidth=2)


def draw_spectrum_nline(data, dataname, label, weights=0, lw=2.0, logx=0, cl=MYLinecolor, bins=500, normed=False):
    '''data is 2-D x numl array
    dataname is numl length array ['','',''] 
    figname default = fig
    make sure numl is given'''
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
        for i in range(0, len(dataname)):
            hd, ad = np.histogram(
                data[i][:], bins=bins, normed=normed, weights=weights[i][:])
            if logx == 1:
                ad = np.log10(ad)
            plsy = plt.semilogy(.5*(ad[1:]+ad[:-1]), hd,
                                color=colorline[i], label=label[i], linewidth=lw)


def draw_field_snapshot(ax, data, extent, caxis=0):
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
    gci = ax.imshow(data, extent=extent, origin='lower', cmap='jet',
                    vmax=vmax, vmin=vmin, interpolation='spline36')
    # colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='3%', pad=0.1)
    cbar = plt.colorbar(gci, cax)


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
