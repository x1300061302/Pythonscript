#!/usr/bin/python
from const import *
import numpy as np
colorline = MYLinecolor; 

def gifmake(gifname = 'gif.gif',duration = 0.1,beg =0,end=0,prefix='.png'):
    import matplotlib.pyplot as plt
    import imageio,os
    images=[]
    filenames=sorted((fn for fn in os.listdir('.') if fn.endswith(prefix)))
    if (end == 0 ):
        end = len(filenames);
    for i in range(beg,end):
        filename = filenames[i];
        images.append(imageio.imread(filename))
        imageio.mimsave(gifname, images,duration=duration)


def draw_histogram2d(x,y, label, bins =[100,200],xylim=0, caxis=0, figname='fig', display=0,  savefig=1):
    '''to draw 2-D histogram figure '''

    try:
        import matplotlib
        matplotlib.use('Agg')
    except:
        print('use matplotlib.use before import plt')
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
    # figure
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    gci = ax.pcolormesh(xx,yy, log10h.T, cmap='jet',
                           vmin=0, vmax=np.max(log10h))
    # axis
    if (type(xylim) == np.int):
        pass 
    else:
        ax.set_xlim(xylim[0])
        ax.set_ylim(xylim[1])
    # label
    plt.xlabel(label[0],fontsize=20)
    plt.ylabel(label[1],fontsize=20)
    plt.title(label[2],fontsize=20)
    plt.xticks(fontsize=20);
    plt.yticks(fontsize=20);
    # colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='3%', pad=0.1)
    cbar = plt.colorbar(gci, cax)

    ffigname = '.'.join([figname, 'png'])

    if (display):
        plt.show()
    if (savefig):
        fig.savefig(ffigname, dpi=300, facecolor='none', edgecolor='b')
        plt.close()


def draw_angle_distribution3d(T, R, label = 0,tlim = 0,rlim=0, figname='fig', binT=360, binR=1000,caxis=0,Display=1):
    '''T,R is N-array 
    rlim = 0 means default else rlim = [r_min,r_max]
    draw angle distribution histogram'''
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    '''histogram'''
    H, Tedge, Redge = np.histogram2d(T, R, bins=[binT, binR])
    log10H = np.log10(H)
    print("get histogram")
    '''change T,R to meshgrid'''
    rad = np.linspace(R.min(), R.max(), binR)
    if(type(tlim)!=np.int):
	    the = np.linspace(tlim[0],tlim[1],binT)
    else:
        the = np.linspace(0, 2*np.pi, binT)
    r, t = np.meshgrid(rad, the)
    fig = plt.figure(figsize=(10, 10))
    ax = plt.subplot(projection='polar')
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
	    ax.set_rlim(np.min(R),np.max(R));
    cbar = fig.colorbar(pcmesh)
    ffigname = '.'.join([figname, 'png'])
    plt.show()
    if (type(label)  != np.int):
        plt.title(label[2]);
    else:
        plt.title(figname)
    fig.savefig(ffigname, dpi=300, facecolor='none', edgecolor='b')
    print(ffigname, 'has been printed')
    plt.close()


def draw_spectrum(data, label, color='black', weights=0, figname='fig', logx=0, display=0,grid = True):
    '''data is 1-D  array
    label = ['xlabel','ylabel','title'] 
    figname default = fig
    '''
    if (display == 0):
        try:
            import matplotlib
            matplotlib.use('Agg')
        except:
            print('use matplotlib.use before import plt')

    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(20, 10))
    ax = plt.subplot(111)
    if (type(weights) == np.int):
        print('No Weights hist')
        hd, ad = np.histogram(data, bins=500, normed=False)
        hd = hd/(ad[1]-ad[0])
        if logx == 1:
            ad = np.log10(ad)
        pl = plt.semilogy(.5*(ad[1:]+ad[:-1]), hd,
                          color=color, label=label[2], linewidth=2)
    else:
        print('Weights hist')
        hd, ad = np.histogram(data, bins=500, normed=False, weights=weights)
        hd = hd/(ad[1]-ad[0])
        if logx == 1:
            ad = np.log10(ad)
        pl = plt.semilogy(.5*(ad[1:]+ad[:-1]), hd,
                          color=color, label=label[2], linewidth=2)

    plt.xlabel(label[0],fontsize=20)
    plt.ylabel(label[1],fontsize=20)
    plt.title(label[2],fontsize = 20)
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.legend(fontsize=20)
    plt.grid(grid);

    ffigname = '.'.join([figname, 'png'])
    fig.savefig(ffigname, dpi=300, facecolor='none', edgecolor='b')
    print(ffigname, 'has been printed')
    if (display):
        plt.show()
    return fig, pl


def draw_spectrum_nline(data, dataname, label, xylim = 0,weights=0, figname='fig', numl=1, logx=0, display=0,grid = True):
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
    fig = plt.figure(figsize=(20, 10))
    ax = plt.subplot(111)
    if (type(weights) == np.int):
        print('No Weights hist')
        for i in range(0, numl):
            hd, ad = np.histogram(data[i][:], bins=500, normed=False)
            if logx == 1:
                ad = np.log10(ad)
            plsy = plt.semilogy(.5*(ad[1:]+ad[:-1]), hd,
                                color=colorline[i], label=dataname[i], linewidth=2)
    else:
        print('Weights hist')
        for i in range(0, numl):
            hd, ad = np.histogram(
                data[i][:], bins=500, normed=False, weights=weights[i][:])
            if logx == 1:
                ad = np.log10(ad)
            plsy = plt.semilogy(.5*(ad[1:]+ad[:-1]), hd,
                                color=colorline[i], label=dataname[i], linewidth=2)

    if (type(xylim) != np.int):
        plt.xlim(xylim[0]);
        plt.ylim(xylim[1]);

    plt.xlabel(label[0],fontsize=20)
    plt.ylabel(label[1],fontsize=20)
    plt.title(label[2],fontsize = 20)
    plt.xticks(fontsize = 20)
    plt.legend(fontsize=20)
    plt.grid(grid);

    ffigname = '.'.join((figname, 'png'))
    fig.savefig(ffigname, dpi=300, facecolor='none', edgecolor='b')
    print("has printed", ffigname)

    if (display):
        plt.show()
    return ax 


def draw_field_snapshot(data, extent, label, caxis=0, xylim=0, Display=0, figname='fig'):
    '''data should be 2-D array (nx,ny) 
    extent should be ([x_min,x_max],[y_min,y_max])
    xylim = ([xmin,xmax],[ymin,ymax])
    label = ['xlabel/Unit','ylabel/Unit','title'] 
    Display dicide whether to show
    figname default = title'''
    if (Display == 0):
        try:
            import matplotlib
            matplotlib.use('Agg')
        except:
            print('use matplotlib.use before import plt')
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(111)

    nx, ny = data.shape
    if (type(caxis) != np.int):
        vmin = caxis[0]
        vmax = caxis[1]
    else:
        vmin = data.min()
        vmax = data.max()

    gci = ax.imshow(data, extent=extent, origin='lower', cmap='jet',\
                    vmax=vmax, vmin=vmin, interpolation='spline36')

    # set axis
    if (type(xylim) == np.int):
        pass
    else:
        ax.set_xlim(xylim[0])
        ax.set_ylim(xylim[1])
    plt.xlabel(label[0])
    plt.ylabel(label[1])
    plt.title(label[2])
    # colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='3%', pad=0.1)
    cbar = plt.colorbar(gci, cax)

    if (Display):
        plt.show()
    else:
        fig.savefig(figname, dpi=300, facecolor='none', edgecolor='b')
        plt.close()


def plot_line(xx, data, label, figname, savefig=1):
    try:
        import matplotlib
        matplotlib.use('Agg')
    except:
        print('use matplotlib.use before import plt')

    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111)
    ax.plot(xx, data, 'r-', label=label[2], linewidth=2.0)
    # figure setting
    plt.xlabel(label[0])
    plt.ylabel(label[1])
    plt.title(label[2])

    plt.legend()


def plot_nline(xx, data, label, display=0, figname='fig'):
    ''' xx = nlxlength np.array datatype 
    '''
    if (display):
        pass
    else:
        try:
            import matplotlib
            matplotlib.use('Agg')
        except:
            print('use matplotlib.use before import plt')
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from mpl_toolkits.mplot3d import Axes3D
    nl = len(data)
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111)
    for i in range(0, nl):
        ax.plot(xx, data[i][:], colorline[i], label=label[2], linewidth=2.0)

    # figure setting
    plt.xlabel(label[0])
    plt.ylabel(label[1])
    plt.title(label[2])

    plt.legend()
    if (display):
        pass
    else:
        fig.savefig(figname+'.png', dpi=300, facecolor='none', edgecolor='b')

#**************main test
