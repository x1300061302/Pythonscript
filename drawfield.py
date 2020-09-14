#draw numberdensity line
#yhbatch -N 1 ./pxiey.sh drawfield.py Jz dirc prefix
import matplotlib
matplotlib.use('Agg')
from const import *
from scipy import ndimage
import sys,getopt
#from mpi4py import MPI

def draw_field(dirc,varname,norm,ffs,i,log_data = False,caxis = None,cmap='jet'):
    a = sdf.read(ffs[i])
    a1 = sdf.read(ffs[i-1])
    var = sr.Get_field_variable(a,varname)  
    if (log_data):
        data = np.nan_to_num(np.log(var.data/norm  + 1e-40))
    else:
        data = var.data/norm
    name = varname + ' at ' + str(a.Header['time']/info['T0'])[0:6] + 'T0'
    ax = df.draw_field_snapshot(ax = 0,data = data.T,caxis = caxis, cmap = cmap); 
    df.Axis_set(ax = ax[0],axesname = ['X/di','Y/di',name],xylims = xylims)
    plt.savefig(savedir+varname + str(i) + '.png', dpi = 300)
    plt.close('all');

#Some examples
try:
    opts,args = getopt.getopt(sys.argv[1:],'hd:p:v:n:c:lm:x:',['help','dirc=','prefix=','varname=','norm=','caxis=','log_data','cmap=','xylims='])
except:
    getopt.GetoptError

for o, a in opts:
    if o in ("-h", "--help"):
        print('example: python drawfield.py -d dirc -p prefix -v varname -n norm -x xylims([-1,1,-1,1]) -c caxis([-1,1]) -l logdata --cmap=')
    if o in ("-d", "--dirc"):
        dirc = a
    if o in ("-p", "--prefix"):
        prefix = a
    if o in ("-v", "--varname"):
        varname = a
    if o in ("-norm", "--norm"):
        norm = float(a)
    else:
        norm = 1.0
    if o in ("-c", "--caxis"):
        caxis = [float(cc) for cc in a[1:-1].split(',')]
    else:
        caxis = [-1.0,1.0]
    if o in ("-l", "--logdata"):
        log_data  = True
    else:
        log_data = False;
    if o in ("-m", "--cmap"):
        cmap = a
    else:
        cmap = 'jet'
    if o in ("-x", "--cmap"):
        xylims = [float(cc) for cc in a[1:-1].split(',')]
    else:
        xylims = [-1,1,-1,1]

savedir  = sr.makedirs(dirc + varname)

#	
#comm = MPI.COMM_WORLD
#rank = comm.Get_rank()
#size = comm.Get_size()
info = sr.read_const(dirc)
##norm = info['drift_V']*info['ne']*qe
##norm = info['B0']*c*info['drift_V']*info['ne']*qe
##norm  = c*info['B0']
##norm = me*c**2
##norm  = 1e6
##info['ne']
ffs = sr.Get_file(prefix = prefix,dirc=dirc)
info = sr.read_const(dirc)
print(ffs)
draw_field(dirc,varname,norm,ffs,0,log_data, caxis,cmap)
#for i in range(len(ffs)//size+1):
#    if (rank == 0):
#        print(i)
#    if (i*size + rank == 0):
#         continue;
#    if (i*size + rank < len(ffs)):
#        draw_field(dirc,varname,info['ne'],ffs,i*size + rank,log_data_list[key],caxis_list[key])
