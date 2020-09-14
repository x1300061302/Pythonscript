from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np
from const import *
import drawfig as df
from scipy import integrate

def diff(time,xx):
    dx = np.zeros(len(xx));
    for i in range(1,len(xx)):
        dx[i-1] = (xx[i]- xx[i-1])/((time[i]-time[i-1])*omega0)
    dx[len(xx)-1] = dx[len(xx)-2]
    return dx
def generate_x(r0, v0):
    x = np.zeros(7)
    x[0:3] = r0
    gam = np.sqrt(1 / (1 - np.sqrt(np.sum(v0**2)) / c))
    x[3:6] = gam * v0
    x[6] = gam
    return x

def motion_equation(x, t):
    part_q = -qe
    lamb = 1 * um
    omega0 = 2 * np.pi * c / lamb
    phi0 = np.pi/2
    a0 = 10
    vph = c
    duration = 30*T0;
    # above ------------parameter set up
    Ex, Ey, Ez, Bx, By, Bz = laser(
        a0, t, x, omega0, 'LP', dura = duration, \
        phi0=phi0, vph=vph, \
        eff=True,alpha = alpha)
    # right - hand
    Dx = np.zeros(7)
    Dx[0] = x[3] / x[6]  # dx/dt = gamma vx /gamma
    Dx[1] = x[4] / x[6]
    Dx[2] = x[5] / x[6]
    Dx[3] = part_q / me * (Ex + (x[4] / x[6] * Bz - x[5] / x[6] * By))
    Dx[4] = part_q / me * (Ey - (x[3] / x[6] * Bz - x[5] / x[6] * Bx))
    Dx[5] = part_q / me * (Ez + (x[3] / x[6] * By - x[4] / x[6] * Bx))
    Dx[6] = part_q / me / c**2 * \
        (x[3] / x[6] * Ex + x[4] / x[6] * Ey + x[5] / x[6] * Ez)
    return Dx


def laser(a0, t, xx, omega0, format, dura, phi0=0, vph=c, eff=False,alpha = [0,0,0]):
    #
    f = lambda x : 1
    #np.sin(np.pi*x/dura)**2
    x = xx[0]
    Phi = omega0*(t-x/vph);
    if((Phi > 0) and (Phi < dura/T0*2*np.pi )):
        E0 = a0 * me * omega0 * c / qe
        B0 = E0 / vph
    else:
        E0 = 0
        B0 = 0
    if(eff == True):
        ef = cal_ef(a0, t, xx,omega= omega0,alpha=alpha)
    else:
        ef = [0, 0, 0, 0, 0, 0]
    if (format == 'CP_right'):
        
        Ex = 0 + ef[0]
        Ey = E0 * f(t-x/c)*np.cos(omega0 * (t - x / vph) + phi0) + ef[1]
        Ez = E0 * f(t-x/c)*np.sin(omega0 * (t - x / vph) + phi0) + ef[2]
        Bx = 0 + ef[3]
        By = -B0 * f(t-x/c)*np.sin(omega0 * (t - x / vph) + phi0) + ef[4]
        Bz = 0 + B0 * f(t-x/c)*np.cos(omega0 * (t - x / vph) + phi0) + ef[5]
    if (format == 'LP'):
        Ex = 0 + ef[0]
        Ey = E0 * f(t-x/c)*np.cos(omega0 * (t - x / vph) + phi0) + ef[1]
        Ez = 0 + ef[2]
        Bx = 0 + ef[3]
        By = 0 + ef[4]
        Bz = B0 * f(t-x/c)*np.cos(omega0 * (t - x / vph) + phi0) + ef[5]
    return Ex, Ey, Ez, Bx, By, Bz


def cal_ef(a0, t, xx,omega,alpha):
    r0 = r00;
    r = np.sqrt(xx[1]**2 + xx[2] ** 2)
    E0 = me * omega * c / qe
    B0 = E0/c
    Bm = alpha[0]*B0;
    #print(Bm)

    Ex = alpha[2]*E0
    Ey = 0
    Ez = 0 
    Bx = alpha[1]*B0 # axial magnetic field 
    By = Bm*xx[2]/r0*np.exp(r**2/r0**2+1/2)
   
    Bz = -Bm*xx[1]/r0*np.exp(r**2/r0**2+1/2)
    ef = [Ex,Ey,Ez,Bx,By,Bz]
    return ef


def main():
    global T0
    global omega0
    r0 = np.array([0, 0, 0])
    v0 = np.array([0.0*c, 0, 0])
    x0 = generate_x(r0, v0)
    time = np.linspace(time_start, time_end, N_step)
    #print(x0)
    xx = odeint(motion_equation, x0, time)
    #draw_figure(xx, time)
    return xx,time


def draw_figure(xx, time):
    fig = plt.figure(figsize=[20,20])
    ax = fig.add_subplot(311)
    ax.plot(time/T0, xx[:,0]/um,'b-',label = 'x-t');
    line1,label1 = ax.get_legend_handles_labels()
    df.Axis_set(ax,axesname=['','$x/um$','$ex = -0.03E_0, bx = 0.05$'],legend=False);
    ax2 = ax.twinx()
    ax2.plot(time/T0, xx[:,1]/um,'r-',label = 'y-t');
    line2,label2 = ax2.get_legend_handles_labels()
    df.Axis_set(ax,axesname=['t/T0','$y/um$',''],legend=False);
    plt.legend(line1+line2,label1+label2,fontsize=20)

    #
    vx = xx[:,3]/xx[:,6]
    ax = fig.add_subplot(312)
    #ax.plot(time/T0, xx[:, 3] / xx[:, 6]/c,label = 'vx-t')
    dvx = diff(time=time,xx=vx/c);
    ax.plot(time/T0,dvx/omega0,'r-',label = 'dvx')
    plt.ylim()
    df.Axis_set(ax,axesname=['t/T0','$v_x/c$',''],legend=True);
  
    ax = fig.add_subplot(313)
    ax.plot(time/T0,xx[:,6],label ='gamma-t')
    df.Axis_set(ax,axesname=['t/T0','$\gamma$',''],legend=True);
    plt.show()
def angle_distribution(xx,time):
    #data deal 
    
    pp = np.sqrt(xx[:,4]**2 + xx[:,5]**2)
    theta = np.arctan(pp/xx[:,3]);
    fig = plt.figure(figsize = [20,5]);
    ax = fig.add_subplot(111)
    ax.plot(time/T0,theta,'r-',label = '$\Theta-t$');
    df.Axis_set(ax,axesname=['','$\Theta$','$ex = -0.03E_0, bx = 0.05$'],legend=False);
    plt.show()
def draw_frequency(xx,time):
    r0 = r00
    omega0 = 2*np.pi*c/(1*um);
    vx = xx[:,3]/xx[:,6];
    E0 = me * omega0 * c / qe
    B0 = E0/c
    Bm = alpha[0]*B0;
    r = np.sqrt(xx[:,1]**2+xx[:,2]**2)
    #r = xx[:,1];
    omega_L = omega0*(1 - vx/c);
    omega_beta = np.sqrt(qe*vx/xx[:,6]/me*Bm/r0*(1-r**2/r0**2)*np.exp(-r**2/2/r0**2+1/2));
    omega_beta2 = np.sqrt(qe*vx/xx[:,6]/me*Bm/r0);
    Phi_L = integrate.cumtrapz(omega_L[5000:],time[5000:]);
    Phi_beta = integrate.cumtrapz(omega_beta[5000:],time[5000:]);
    Phi_beta2 = integrate.cumtrapz(omega_beta2[5000:],time[5000:]);
    fig = plt.figure(figsize=[10,10])
    ax = fig.add_subplot(211);
    
    ax.plot(time[5000:-1]/T0,Phi_L,label ='Phi_L')
    ax.plot(time[5000:-1]/T0,Phi_beta,label = 'Phi_{beta}');
    ax.plot(time[5000:-1]/T0,Phi_beta2,label = 'Phi_{beta2}');
    ax.plot(time[5000:-1]/T0,Phi_L-Phi_beta,label='L-beta');
    ax.plot(time[5000:-1]/T0,xx[:,6][5000:-1]/30,label = 'gamma/30')
    plt.legend()
    
    ax2 = fig.add_subplot(212);
    ax2.plot(time[5000:-1]/T0,omega_L[5000:-1],label = 'L');
    ax2.plot(time[5000:-1]/T0,omega_beta[5000:-1],label ='beta');
    ax2.plot(time[5000:-1]/T0,omega_beta2[5000:-1],label ='beta2');
    
    
    plt.legend()
    plt.show()

#init
r00 = 2*um    
T0 = 1 * um / c
omega0 = 2*np.pi/T0;
time_start, time_end, N_step = [0.0,400*T0,10000]
# main-------------------------------
print('begin')
#alpha1= np.linspace(-5,0,500);
#gam = []
#for i in range(0,len(alpha1)):
#    print(i)
#    alpha =[10**(alpha1[i]),0,0.1]
#    xx,time = main()
#    gam.append(np.max(xx[:,6]))
alpha=[3,0,0.3]
xx,time = main()
draw_frequency(xx,time)
draw_figure(xx,time)
angle_distribution(xx,time)
print('end')