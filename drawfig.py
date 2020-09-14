#!/usr/bin/python
from const import *
import numpy as np
colorline = MYLinecolor
font = {'family':'serif'}
figw, figh = [10,10]
ix,iy = [1,1]
iwidth, iheight= [8,8]


def draw_dist_fn(data,ax = 0,axis_set = True,label = ''):
    x = data.grid.data[0];
    y = data.data
    if (ax == 0):
        fig,ax = Create_Figure()
    ax.plot(x/Mev,y,label = label)
    if (axis_set): 
        Axis_set(ax = ax,axesname = 'sp',logscale = 'xy')
        Line_set(ax = ax )
    return ax

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

def Create_Artist_Figure(para={'figsize':(10,10)}):
    '''Input:

    Output:fig,ax
    '''
    import mpl_toolkits.axisartist as axisartist
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize = para['figsize'])
    ax = axisartist.Subplot(fig,111)
    fig.add_axes(ax)
    ax.axis[:].set_visible(False)
    ax.axis["x"] = ax.new_floating_axis(0,0)
    ax.axis["x"].set_axisline_style("->",size = 1.0)
    ax.axis["y"] = ax.new_floating_axis(1,0)
    ax.axis["y"].set_axisline_style("-|>",size = 1.0)
    ax.axis["x"].set_axis_direction("bottom")
    ax.axis["y"].set_axis_direction("left")
    return fig,ax
def Create_Figure2(h,w,aw=5,ah=5,sp_w =1,sp_h = 1,p0 = [0.5,0.5]):
    '''
    Input:
    h x w , h is row and w is column
    default aw = 10, ah = 10, sp_w(space in column) = 2,sp_h(space in row) = 2, sw = sh = 8
    p0 is the initial position(rate) default = [2,2]
    Output:
    fig
    ax_list
    '''

    import matplotlib.pyplot as plt
    fw = aw*w + (w - 1)*sp_w + p0[0]*2
    fh = ah*h + (h - 1)*sp_h + p0[1]*2
    rx0 = p0[0]/fw
    ry0 = p0[1]/fh
    spw = sp_w/fw
    sph = sp_h/fh

    #calculate the space for figure
    rsw = aw/fw
    rsh = ah/fh

    ax_list = []
    
    fig = plt.figure(figsize= [fw,fh])
    for i in range(0,h):
        for j in range(0,w):
            ax = fig.add_axes([rx0+spw*(j)+rsw*j,ry0+sph*(i)+rsh*i,rsw,rsh]) #横排从左往右,#竖排从下往上
            ax_list.append(ax);       
    return fig,ax_list
#     plt.show()

def Create_Figure(axes =[6,4],p0 = None ,projection=None):
    '''
    Input:figsize is auto setted by axes + 2*p0
          axes is the size of axes default is [6,4]
          p0 is the initial position of axes default is 0.1,0.1 as ratio
          projection is to decide the polar axes, projection = None,'polar'
    Output:fig,ax
    '''
    import matplotlib.pyplot as plt
    #modification:
    if (p0 == None):
       p0 = np.zeros(2);
       p0[0] = 0.2
       p0[1] = 0.2

    figsize = [axes[0]/(1 - 2*p0[0]) , axes[1]/(1 - 2*p0[1])]
    fig = plt.figure(figsize=figsize)
    if (projection == None):
        ax = fig.add_axes([p0[0],p0[1],axes[0]/figsize[0],axes[1]/figsize[1]])
    else:
        #ax = fig.add_subplot(x, y, n)
        ax = fig.add_subplot(1,1,1,projection = projection)
    return fig,ax


# ------------Figure/axis setting------$$$$$$$$$
def Line_set(ax,linewidths = None,linestyles = None,labels = None,markers =None):
    '''Input: ax: is the axis you draw line
              linewidth:single value of a list 
              linestyle:
              label:
    '''
    lines = ax.get_lines()
    for i in range(0,len(lines)):
        line = lines[i]

        #linestyles
        if (linestyles == None):
            linestyle = '-';
        else:
            linestyle = linestyles[i];
        #linewidths
        if (linewidths == None):
            linewidth = 3.0
        else:
            linewidth = linewidths[i]

        #linemarkers
        if (markers == None):
            marker = None 
        else:
            marker = markers[i]

        if (labels != None):
            line.set_label(labels[i])
            
        line.set_linewidth(linewidth)
        line.set_linestyle(linestyle)
        line.set_marker(marker)
#         line.set_label(label[i])
    
    
    
def Axis_set(ax,axesname=[r'$x/d_i$',r'$y/d_i$',''],fs=35.0,xticklabel=0,xtickrange=0,yticklabal=0,ytickrange=0,grid=False,legend=False,xylims = None,ax_style='d',lw = 2.0,showtick = False, ticklength = 10,set_line = False,logscale = '',family = 'sans-serif'):
    import matplotlib.pyplot as plt
    plt.rc({'family':family})
    if (axesname == 'sp'):   #Quick set 
        axesname = ['E/MeV','dN/dE',''];
    elif (axesname == 'cb'): # for colorbar
        axesname = ['','','']
        
    if (ax_style == 'd'): 
        if (xylims != None):
            try:
                ax.set_xlim(xylims[0],xylims[1]);
            except:
                print('No xlims')
            try:
                ax.set_ylim(xylims[2],xylims[3]);
            except:
                print('No ylims')
        # x-axis
        ax.set_xlabel(axesname[0],fontsize=fs+3); 
        #ax.set_xticks(xticklabel); 
        #y-axis
        ax.set_ylabel(axesname[1],fontsize=fs+3);
        #ax.set_yticks(fontsize=fs);
	#title
        ax.set_title(axesname[2],fontsize=fs-3);
        ax.grid(grid);
        ax.minorticks_on()
        if (legend):
            plt.legend(fontsize =fs);
	
        ax.spines['bottom'].set_linewidth(lw)
        ax.spines['left'].set_linewidth(lw)
        ax.spines['right'].set_linewidth(lw)
        ax.spines['top'].set_linewidth(lw)
        #ax.tick_params(which = 'both',width = lw,colors = 'black',labelsize=fs,direction = 'in',length = ticklength)
        ax.tick_params(axis = 'x',labelsize = fs*0.8)
        ax.tick_params(axis = 'y',labelsize = fs*0.8)
        ax.tick_params(which = 'minor',direction = 'in',length = ticklength/2,width = lw/2)
        ax.tick_params(which = 'major',direction = 'in',length = ticklength,width=lw/2)
        if (showtick):
            ax.tick_params(which = 'minor',top = True,right = True)
            ax.tick_params(which = 'major',top = True,right = True)
        if (set_line):  
            Line_set(ax = ax)
        if (logscale =='x'):
            ax.set_xscale('log')
        elif (logscale =='y'): 
            ax.set_yscale('log')
        elif (logscale =='xy'):
            ax.set_yscale('log')
            ax.set_xscale('log')
		
    elif (ax_style == 'p'):
        if (type(xylims) != np.int):
            ax.set_rlim(xylim[1])
        # x-axis
        ax.set_xlabel(axesname[0],fontsize=fs+3); 
        ax.set_xticks(fontsize=fs);
        #y-axis
        ax.set_ylabel(axesname[1],fontsize=fs+3);
        ax.set_yticks(fontsize=fs);
        #title
        ax.set_title(axesname[2],fontsize=fs-3);
        plt.grid(grid);
        ax.spines['polar'].set_linewidth(lw)
        if (legend):
            plt.legend(fontsize =fs,shadow=True);
    elif (ax_style == '3d'):
        ax.set_xlabel('X',fontsize = fs)
        ax.set_ylabel('Y',fontsize = fs)
        ax.set_zlabel('Z',fontsize = fs)
        ax.spines['bottom'].set_linewidth(lw)
        ax.spines['left'].set_linewidth(lw)
        ax.spines['right'].set_linewidth(lw)
        ax.spines['top'].set_linewidth(lw)
        #ax.tick_params(which = 'both',width = lw,colors = 'black',labelsize=fs,direction = 'in',length = ticklength)
        ax.tick_params(axis = 'x',labelsize = fs*0.8)
        ax.tick_params(axis = 'y',labelsize = fs*0.8)
        ax.tick_params(axis = 'z',labelsize = fs*0.8)
        
    else:
        print('Can not understand ax style'+ax_style);

def Colorbar_set(ax,gci,pos='right',size='3%',pad=0.1,ori = 'vertical'):
    '''
    ori = vertical or horizantal
    '''
    import matplotlib.cm as cm
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import matplotlib.pyplot as plt
    divider = make_axes_locatable(ax)
    cax = divider.append_axes(pos, size=size, pad=pad)
    cbar = plt.colorbar(gci, cax,orientation = ori)
    #cbar.set_ticks() and other operation
    return cbar

def Colorbar_set2(fig,ax,gci,size = 0.05,sp = 0.01, position = 'right' , ori='vertical'):
    '''
    Input: fig, ax ,gci size,sp,position,ori
    Output:
    cbar,cax

    '''
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    pos = ax.get_position()
    #for position
    if (position == 'right'):
        bbox = [pos.x1 + sp,pos.y0,size*pos.width,pos.height]
    elif (position == 'left'):
        bbox = [pos.x0-size*pos.width-sp,pos.y0,size*pos.width,pos.height]
    elif (position =='top'):
        bbox = [pos.x0,pos.y1+sp,pos.width,pos.height*size]
    elif (position =='bottom'):
        bbox = [pos.x0,pos.y0-size-sp,pos.width,pos.height*size]
    
    cax = fig.add_axes(bbox)

    cbar = plt.colorbar(gci,cax,orientation= ori)

    #for ticks 
    if (position == 'right'):
        pass;
    elif (position == 'left'):
        cax.yaxis.tick_left()
    elif (position =='top'):
        cax.xaxis.tick_top()

    return cbar,cax

def Legend_outside(ax,loc,ncol = 1,bbox=(1.0,1.0),fs=20):
    box = ax.get_position();
    ax.set_position([box.x0,box.y0,box.width,box.height]);
    ax.legend(loc = loc,ncol = ncol, bbox_to_anchor=bbox,fontsize=fs)
    

#***********Draw----------------$$$$$$$$$$$$$    
def draw_histogram2d(x, y, ax = 0, bins=[500, 500], xylim=0, caxis=0, fontsize=20,cmap = 'jet',log10h =True,typefigure = 'imshow',weights = None,fig = None):
    '''to draw 2-D histogram figure 
        Input:
        x,y = 1-D np.ndarray
        handle = ax
        bins =[nx,ny]
        xylim = [[xlim],[ylim]] means the region to show
        caxis = region of colobar
        log10h = log10(h)
        Output:
        ax,gci 
        
    
    '''
    import matplotlib.pyplot as plt

    if (ax == 0):
        fig,ax = Create_Figure();
 
    h_xy, xedge, yedge = np.histogram2d(x, y, bins=bins,weights = weights)
    if (log10h):
        log10h = np.log10(h_xy)
    else:
        log10h = h_xy

    if (type(caxis) == type([])):
        vmin = caxis[0]
        vmax = caxis[1]
    elif (caxis == 'sym_bar'):
        vmin = np.min(np.min(log10h))
        vmax = np.max(np.max(log10h))
        value = np.max(abs(vmin),abs(vmax));
        vmin = -value
        vmax = value
    elif (caxis == 0):
        vmin = None
        vmax = None
    print('get histogram')
    # axis
    xax = np.linspace(np.min(x), np.max(x), bins[0])
    yax = np.linspace(np.min(y), np.max(y), bins[1])
    xx, yy = np.meshgrid(xax, yax)
    # draw
    #gci = ax.pcolormesh(xx, yy, log10h.T, cmap=cmap,
    #                    vmin=vmin, vmax=vmax)
    if (typefigure == 'imshow'):
        gci = ax.imshow(log10h.T,origin = 'lower',extent = [xedge[0],xedge[-1],yedge[0],yedge[-1]],aspect = 'auto',cmap = cmap,vmin = vmin,vmax = vmax)
    elif (typefigure == 'contourf'):
        gci = ax.contourf(xx,yy,log10h.T) 
        
    cb = Colorbar_set2(fig,ax,gci);
    return ax,gci,cb
    
    # colorbar

def draw_angle_distribution3d(ax, T, R, weights = 0,tlim=0, rlim=0, binT=360, binR=1000, caxis=0, log10=True,cmap='jet'):
    ''' 
    Input: ax
           T = Theta  1-D np.ndarray or list
           R = Radius 1-D np.ndarray or list 
           weights = Weight 1-D np.ndarray or list
           tlim,rlim,binT,binR,caxis,log10,cmap
    Output:
           ax,pcmesh
    Note:
    make sure ax style is polarization
    '''
    if (ax == 0):
        fig,ax = Create_Figure(polar=True);
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    if (type(weights) == np.int):
        weights = np.ones(T.shape)
    #histogram
    H, Tedge, Redge = np.histogram2d(T, R, weights=weights, bins=[binT, binR])
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

def draw_spectrum(data, ax = 0, label='', weights=0, grid=True,gethd=0,lw=3.0,cl=None,factor=1.0,normed = False,smooth_data = False,bins = 500,axis_set = True,logdata = False,plot = True):

    ''' Input: data 
               ax = 0 default figure and axis;
               label = '' 
               weights = 0 
               grid =True
               lw means the linewidth of the curve
               cl means the color of the line 
               factor means the factor * hd 
        Output:
           return ax,x,y
    '''
    if (ax == 0 and plot):
        fig,ax = Create_Figure();
    import matplotlib.pyplot as plt
    if (logdata):
        data = np.log10(data);
    if (type(weights) == np.int):
        #print('No Weights hist')
        hd, ad = np.histogram(data, bins=bins, normed=normed)
    else:
        #print('Weights hist')
        hd, ad = np.histogram(data, bins=bins, normed=normed, weights=weights)
    hd = hd/(ad[1]-ad[0]);
    axx = 0.5*(ad[1:]+ad[:-1]);
    if (smooth_data):
        hd = sr.smooth_data(hd,sigma = 2)
    if (logdata):
        axx = 10**axx
    if (plot):
        ax.plot(axx, factor*hd,\
            color=cl, label=label, linewidth=lw)
    #df.Axis_set(ax = ax,axesname = 'sp',logscale = 'xy')
    if (plot and axis_set):
        Axis_set(ax = ax, axesname = ['E','dN/dE',''],logscale = 'xy')
    return ax,axx,factor*hd
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
    numl = len(data)
    if (ax == 0):
        fig,ax = Create_Figure();
    if (type(label)==np.int):
        label = []
        for i in range(0,numl):
            label.append(str(i))
    import matplotlib.pyplot as plt
    if (type(weights) == np.int):
        print('No Weights hist')
        for i in range(0, numl):
            hd, ad = np.histogram(data[i][:], bins=500, normed=normed)
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


def draw_field_snapshot(data, extent=None, ax = 0, caxis=None,cmap='jet',bar_set=True,axis = True,name = None,sym = True):
    '''Input:
    data should be 2-D array (nx,ny) 
    extent should be ([x_min,x_max],[y_min,y_max])
    xylim = ([xmin,xmax],[ymin,ymax])
    label = ['xlabel/Unit','ylabel/Unit','title'] 
    Display dicide whether to show
    figname default = title
       Output:
        or return ax,gci 
        return ax,gci,cb;
    
    '''
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from mpl_toolkits.mplot3d import Axes3D
    #For python it is strange for nx,ny 
    ny,nx = data.shape
    if (caxis != None):
        vmin = caxis[0]
        vmax = caxis[1]
    elif (sym):
        vmin = -max(abs(data.min()),abs(data.max()))
        vmax = -vmin
    else:
        vmin = data.min()
        vmax = data.max()
    if (extent == None):
        extent = [0,nx,0,ny]

    if (ax == 0):
        ly = extent[3]-extent[2]
        lx = extent[1]-extent[0]
        ratio = ly/lx
        fig,ax = Create_Figure(axes = [8.0,8.0*ratio]);

    gci = ax.imshow(data, extent=extent, origin='lower', cmap=cmap,
                    vmax=vmax, vmin=vmin, interpolation='spline36')

    if (axis):
        Axis_set(ax = ax)
    if(name != None):
        Axis_set(ax = ax, axesname = [r'$X/d_i$',r'$Y/d_i$',name])

    if (bar_set):
        cb = Colorbar_set(ax,gci);
        return ax,gci,cb;

    
    # colorbar
    return ax,gci;


def plot(data, xx = None,ax = 0,label='', lw=2.0):
    if (ax == 0):
        fig,ax = Create_Figure();
    if (xx == None):
        xx = np.arange(len(data))
    ax.plot(xx,data,'r-', label=label, linewidth=lw)
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

def draw_MagneticLine(bx,by,x,y,ax=0,avg = True,di = 1,mypara = {}):
    '''
    Default Para:
    density: 1
    lw: 1.5
    arrowstyle : ->
    
    '''
    #default para
    para = {'density':1,'lw':1.5,'arrowstyle':'->','color':'r','minlength':0.5,'maxlength':10.0}
    para.update(mypara)
    #print(para)
    
    if (ax == 0):
        fig,ax = Create_Figure();
    sp = ax.streamplot(x,y,bx.T,by.T,\
                   color = para['color'],\
                   density = para['density'],\
                   linewidth = para['lw'],\
                   arrowstyle = para['arrowstyle'],\
                   minlength = para['minlength'],\
                   maxlength = para['maxlength'],\
                   )
    return ax,sp
##########Quick Draw Part:
class quick_draw(object):
    import matplotlib.pyplot as plt
    '''
    I want to achieve quick_draw figure in this class as 
    Q = quick_draw()
    Q.draw_Ekbar()
    Q.draw_photon_Ekbar()
    like this
    '''
    def __init__(self,prefix='',dirc = '',name = 'sdf',info = {}):
        self.dirc = dirc
        self.name = name
        self.files = sr.Get_file(prefix = prefix,dirc = dirc)
        try:
            self.s = sr.read_const(dirc);  
        except:
            print('No const.deck')

        self.para={
                'x_axis':np.linspace(0,1,100),
                'y_axis':np.linspace(0,1,100),
            'figsize':[10,5],
            'norm':1,
            'caxis':0,
            'cmap':'jet',
            'xylims':0,
            'extent':[0,1,0,1],
            'save':False,
            'axesname': [r'$x/d_i$',r'$y/d_i$',r'$title$'],
            'density': 1,
            'linewidth': 2.0,
            'bins':[500,500]
        }
        #update extent
        self.s.update(info)
        try:
            x1,x2 = [self.s['xmin']/self.s['di'],self.s['xmax']/self.s['di']]
            y1,y2 = [self.s['ymin']/self.s['di'],self.s['ymax']/self.s['di']]
            extent = [x1,x2,y1,y2]
            xx = np.linspace(x1,x2,self.s['nx'])
            yy = np.linspace(x1,x2,self.s['ny'])
            para = {'extent':extent,\
                    'xylims':[[x1,x2],[y1,y2]],\
                    'x_axis':xx,\
                    'y_axis':yy,\
                    }
            self.para.update(para)
        except:
            print('no xmin,xmax,ymin,ymax or di information in const.deck')
        
        self.default_para = self.para;
    def get_s(self):
        return self.s;
    def update_s(self,info):
        self.s.update(info)
#    @staticmethod
#     def get_norm(self,key=''):
    def get_para(self):
        return self.para;
    def set_para(self,para):
        self.para.update(para)
#         return self.para;
    def qdraw_MagneticLine(self,ax=0,n = 0, para={},avg = True):
        a = sdf.read(self.files[n])
        if (ax == 0):
            fig,ax = Create_Figure();
        if (avg):
            bx = sr.Get_field_variable(a,var = 'Bx_averaged')
            by = sr.Get_field_variable(a,var = 'By_averaged')
        else:
            bx = sr.Get_field_variable(a,var = 'Bx')
            by = sr.Get_field_variable(a,var = 'By')
        draw_MagneticLine(bx=bx.data,by=by.data,x=self.para['x_axis'],y=self.para['y_axis'],ax=ax,avg = True,di = self.s['di'])
        Axis_set(ax,\
                        axesname=self.para['axesname'],\
                        xylims = self.para['xylims'],\
                       )
        return ax
    def qdraw_phasespace(self,x,y,ax = 0,para={}):
        self.set_para(para);
        ax,gci = draw_histogram2d(x = x.data/self.para['norm'],y = y.data/self.para['norm'],ax = ax,bins = self.para['bins'],caxis = self.para['caxis'])
        xname = x.name.split('/')[1]
        yname = y.name.split('/')[1]
        df.Axis_set(ax,axesname = [xname,yname,'Phase Figure'],xylims = self.para['xylims'])
        df.Colorbar_set(ax,gci)
        #self.set_para(para={'axesname':[])

    def qdraw_spectrum(self,ax,key, para={}, weight = 0):
        key = key.split('.')[1]
        self.para.update(self.default_para)
        if (ax == 0):
            fig,ax = Create_Figure();
        var = self.a.__dict__[key]; 
        speciesname = var.name.split('/')[3];# species
        print(np.min(var.data),np.max(var.data))
        self.set_para(para={'axesname':['E/Mev','dN/dE',speciesname+str(self.a.Header['time']/self.s.const['T0'])[0:5]+'T0'],\
                            'xylims':self.para['xylims']});
        if (type(weight) == np.int):
            draw_spectrum(ax=ax,data = var.data/Mev);
#         else:
#             keyw =  var.name.split('/');
#             weight = self.a.dict__[key]
        Axis_set(ax,\
                        axesname=self.para['axesname'],\
                        xylims = self.para['xylims'],\
                       )

    def savefig(self,name=None,dpi = 200):
        import matplotlib.pyplot as plt
        '''
        Input:
        name is the name of the saved file
        and dpi = 200 is default
        '''
        if (name == None):
            name = self.dirc+self.para['axesname'][2]+'.png'; 
        plt.savefig(name,dpi = dpi)
        
    def qdraw(self,ax,n,key,para={}):
        import matplotlib.pyplot as plt
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
            Return fig,ax,cb
        '''
        ##########data 
        
        a = sdf.read(self.files[n])
        var = sr.Get_field_variable(a,var = key)
        print(var)
        #############
        self.para['caxis'] = 0; # autosetting vmin - vmax;
        self.para['axesname'][2] = key + ' '+str(a.Header['time']/self.s['T0'])[0:5] + 'T0'

        vmin = np.min(var);
        vmax = np.max(var);

        data = var.data
        #update draw para
        self.para.update(para)
        
        if(self.para['norm'] == 0):
            self.para['norm'] == np.max(data)


        ax,gci,cb = draw_field_snapshot(ax=ax,\
                                 data=data.T/self.para['norm'],\
                                 extent=self.para['extent'],\
                                 cmap=self.para['cmap'],\
                                 caxis=self.para['caxis'],\
                                       )
    #                                             )
        Axis_set(ax,\
                        axesname=self.para['axesname'],\
                        xylims = self.para['xylims'],\
                       )
        #autoseting colorcaxis:
        if (self.para['save']):
            savefig(name = self.dirc+self.para['axesname'][2]+'.png',dpi = 200);

        self.para.update(self.default_para)
        return ax 
#**************main test
