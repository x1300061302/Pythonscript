##This script is aimed to dealwith the data of 2d Magnetic reconnection.

def draw_Ekbar(ax,norm):
    ax,gci = df.draw_field_snapshot(ax=ax,\
                                            data=Ekbar.T/norm,\
                                            extent=np.array(extent)/di,\
                                            cmap='jet',\
                                    caxis=[0,1]
                                   )
#                                             )
    df.Axis_set(ax,\
                    axesname=['$x/d_i$','$y/d_i$','$Ekbar$ at '+ str((i)*dT)+'$T_0$'],\
                    xylims = xylims,\
                   )
    acs = df.Colorbar_set(ax,gci)
#     plt.show()
    return acs

def draw_photon_Ekbar(ax):
    ax,gci = df.draw_field_snapshot(ax=ax,\
                                            data=pho_Ekbar.T/1e-12,\
                                            extent=np.array(extent)/di,\
                                            cmap='jet',\
                                    caxis=caxis
                                   )
#                                             )
    df.Axis_set(ax,\
                    axesname=['$x/d_i$','$y/d_i$','$Ekbar$ at '+ str((i)*dT)+'$T_0$'],\
                    xylims = xylims,\
                   )
    df.Colorbar_set(ax,gci)
#     plt.show()


def draw_Bx(ax):
    ax,gci = df.draw_field_snapshot(ax=ax,\
                                            data=Bx.T,\
                                            extent=np.array(extent)/di,\
                                            cmap='seismic',\
                                   )
#                                             )
    df.Axis_set(ax,\
                    axesname=['$x/d_i$','$y/d_i$','$B_x$ at '+ str((i)*dT)+'$T_0$'],\
                    xylims = xylims,\
                   )
    df.Colorbar_set(ax,gci)
    plt.show()
    #plt.savefig('bz'+str(i)+'.eps')


    #contour
    # fig,ax = df.Create_Figure()
    # ax.contour(Bx.T)

    #Jx
def draw_Jx(ax):
    ax,gci = df.draw_field_snapshot(ax=ax,\
                                            data=Jx.T/Jz0,\
                                            extent=np.array(extent)/di,\
                                            cmap='seismic',\
                                        caxis=[-1,1]
                                       )
    df.Axis_set(ax,\
                    axesname=['$x/d_i$','$y/d_i$','$J_x$ at '+ str((i)*dT)+'$T_0$'],\
                      xylims = xylims,\
                   )
    df.Colorbar_set(ax,gci)
    #Jy
def draw_Jy(ax):
    ax,gci = df.draw_field_snapshot(ax=ax,\
                                            data=Jy.T/Jz0,\
                                            extent=np.array(extent)/di,\
                                            cmap='seismic',\
                                        caxis=[-1,1]
                                       )
    df.Axis_set(ax,\
                    axesname=['$x/d_i$','$y/d_i$','$J_y$ at '+ str((i)*dT)+'$T_0$'],\
                      xylims = xylims,\
                   )
    df.Colorbar_set(ax,gci)



    #Jz
def draw_Jz(ax):
    ax,gci = df.draw_field_snapshot(ax=ax,\
                                            data=Jz.T/Jz0,\
                                            extent=np.array(extent)/di,\
                                            cmap='seismic',\
                                        caxis=[-1,1]
                                       )
    df.Axis_set(ax,\
                    axesname=['$x/d_i$','$y/d_i$','$J_z$ at '+ str((i)*dT)+'$T_0$'],\
                      xylims = xylims,\
                   )
    df.Colorbar_set(ax,gci)

    ##Bz
def draw_Bz(ax):
    ax,gci = df.draw_field_snapshot(ax=ax,\
                                            data=Bz.T,\
                                            extent=np.array(extent)/di,\
                                            cmap='seismic')
    df.Axis_set(ax,\
                    axesname=['$x/d_i$','$y/d_i$','$B_z$ at '+ str((i)*dT)+'$T_0$'],\
                    xylims = xylims,\
                   )
    df.Colorbar_set(ax,gci)
    ###By
def draw_By(ax):
    ax,gci = df.draw_field_snapshot(ax=ax,\
                                            data=By.T,\
                                            extent=np.array(extent)/di,\
                                            cmap='seismic')
    df.Axis_set(ax,\
                    axesname=['$x/d_i$','$y/d_i$','$B_y$ at '+ str((i)*dT)+'$T_0$'],\
                    xylims = xylims,\
                   )
    df.Colorbar_set(ax,gci)
    #Ez 
def draw_Ez(ax):
    ax,gci = df.draw_field_snapshot(ax=ax,\
                                            data=Ez.T/Ez0,\
                                            extent=np.array(extent)/di,\
                                            cmap='seismic',\
                                            caxis=[-10,10])
    df.Axis_set(ax,\
                    axesname=['$x/d_i$','$y/d_i$','$E_z$ at '+ str((i)*dT)+'$T_0$'],\
                    xylims = xylims,\
                   )
    df.Colorbar_set(ax,gci)



    ####Nume
def draw_Nume(ax):
    ax,gci = df.draw_field_snapshot(ax=ax,\
                                            data=Nume.T/ne,\
                                            extent=np.array(extent)/di,\
                                            cmap='jet',\
                                            caxis = [0,10],\
                                            )
    df.Axis_set(ax,\
                    axesname=['$x/d_i$','$y/d_i$','$Num_e$ at '+ str((i)*dT)+'$T_0$'],\
                    xylims = xylims,\
                   )
    axc = df.Colorbar_set(ax,gci)
        # axc.set_clim([0,1.5])
    #     plt.show()
def draw_pho_Nume(ax):
    ax,gci = df.draw_field_snapshot(ax=ax,\
                                            data=Nume_pho.T/ne,\
                                            extent=np.array(extent)/di,\
                                            cmap='jet',\
                                            caxis = caxis,\
                                            )
    df.Axis_set(ax,\
                    axesname=['$x/d_i$','$y/d_i$','$Num_pho_e$ at '+ str((i)*dT)+'$T_0$'],\
                    xylims = xylims,\
                   )
    axc = df.Colorbar_set(ax,gci)
        # axc.set_clim([0,1.5])
    #     plt.show()

def draw_streamplot(ax):
    ax.streamplot(x,y,Bx.T,By.T,\
                      density = 1,\
                      linewidth = 3.0,\
                      cmap = 'jet')
def draw_contour(ax):
    C = ax.contour(xx,yy,np.sqrt(Bx.T**2+By.T**2)/B0,\
                       extent = np.array(extent)/di,\
                       linewidths = 3.0,\
                       cmap='Oranges_r')
    df.Axis_set(ax,axesname=['$x/d_l$','$y/d_l$','Field Line at '+ str(i*dT)+'$T_0$'],\
                xylims = xylims)
    plt.clabel(C,fontsize=20)

def draw_De(ax,xylims):
    ax,gci = df.draw_field_snapshot(ax=ax,\
                                            data=De.T,\
                                            extent=np.array(extent)/di,\
                                            cmap='jet',\
                                            caxis = [-1e41,1e41],\
                                            )
    ax.streamplot(x,y,Bx.T,By.T,density=2,linewidth = 3.0)
    xylims = xylims
    df.Axis_set(ax,\
                    axesname=['$x/d_i$','$y/d_i$',r'$D_e$ at '+ str((i)*dT)+'$T_0$'],\
                    xylims = xylims,\
                    legend = False
                   )
    axc = df.Colorbar_set(ax,gci)
def draw_Temp(ax):
    ax,gci = df.draw_field_snapshot(ax=ax,\
                                            data=Temp.T/1e12,\
                                            extent=np.array(extent)/di,\
                                            cmap='jet',\
                                            caxis = [0,5],\
                                            )
    #     ax.streamplot(x,y,Bx.T,By.T)
    #     xylims = xylims
    df.Axis_set(ax,\
                    axesname=['$x/d_i$','$y/d_i$',r'$Temp$ at '+ str((i)*dT)+'$T_0$'],\
                    xylims = xylims,\
                    legend = False
                   )
    axc = df.Colorbar_set(ax,gci)

    axc.set_clim([0,5])
def Get_De(a):
    Jx = a.Current_Jx_averaged.data;
    Jy = a.Current_Jy_averaged.data;
    Jz = a.Current_Jz_averaged.data
    Bx = a.Magnetic_Field_Bx_averaged.data
    By = a.Magnetic_Field_By_averaged.data
    Bz = a.Magnetic_Field_Bz_averaged.data
    Ex = a.Electric_Field_Ex_averaged.data
    Ey = a.Electric_Field_Ey_averaged.data
    Ez = a.Electric_Field_Ez_averaged.data
    Nume = a.Derived_Number_Density_electron.data;

    nx,ny = Nume.shape
#     print(nx,ny)
    extent = sr.Get_extent(a)
    xmin,xmax,ymin,ymax = np.array(extent)/di
#     print(xmin,xmax,ymin,ymax)
    x = np.linspace(xmin,xmax,nx)
    y = np.linspace(ymin,ymax,ny)
    xx,yy = np.meshgrid(x,y)
    Vx = Jx/Nume/qe
    Vy = Jy/Nume/qe
    Vz = Jz/Nume/qe

    Dex = Jx*(Ex+Vy*Bz-Vz*By)
    Dey = Jy*(Ey+Vz*Bx-Vx*Bz)
    Dez =  Jz*(Ez + Vx*By - Vy*Bx)

    VdotE = Vx*Ex + Vy*Ey + Vz*Vz
    De = Dex+ Dey + Dez
    return De
    