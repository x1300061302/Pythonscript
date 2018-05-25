#import sdf
#import numpy as np
#import sys
#import matplotlib.cm as cm
#from constants import constants as const
#from mpl_toolkits.axes_grid1 import make_axes_locatable
class sdf_class(object):
    '''
    This class contanin some functions to visulize sdf data file.
    '''
#initialization
    def __init__(self):
        pass
#line plot parameter
    def line_set(self):
        '''
        This function stored some predefined line plot parameters.
        '''
        col_sty = ['r-','g--','b-.','c:','k-','y-','m-']
        #col_sty = ['r-','b-.','c:','k-','y-','m-']
        return col_sty
#get file name list
    def get_list(self,filename):
        '''
        This function is used to get file name list.
        parameters:
        filename----filenames,integer or integer list.
        '''
        import sys
        if(type(filename) == type(1)):
            return [filename]
        elif(type(filename) == type([1,2])):
            return filename
        else:
            sys.exit(0)
#get data from a sdf file
    def get_data(self,filename,prefix='1'):
        '''
        This function is used to get the data dictionary from a sdf file.
        parameters:
        filename----a sdf file name is needed, defualt:0001.sdf
        prefix------file name prefix, default:1
        '''
        import sdf
        #construct file name
        if(prefix != 'None'):
            sdfname = prefix+str(filename).zfill(4)+'.sdf'
        else:
            sdfname = str(filename).zfill(4)+'.sdf'
        sdfdata = sdf.read(sdfname)
        data_dict = sdfdata.__dict__
        return data_dict
#use a loop to get the exact physical field
    def get_field(self,field='bx',keys=["Magnetic_Field_Bx"]):
        '''
        This function is used to get the exact physical field name in keys.
        parameters:
        field-------the field name. default:'bx'
        keys--------the field names list, default:"Magnetic_Field_Bx"
        '''
        import sys
        ifornot = False
        for eachone in keys:
            if (field.lower() in eachone.lower()):
                field = eachone
                ifornot = True
                break
        if(ifornot == False):
            print "There is not such a field name %s."%(field)
            sys.exit(0)
        return field
#use a dictionary to chose the normalization constant
    def get_constant(self,field='Magnetic_Field_Bx'):
        '''
        This function is used to chose the normalization constant according to the field.
        parameters:
        field-------the field name. default:'Magnetic_Field_Bx'
        '''
        from constants import bubble_mr as const
        normal = {"magnetic":const.B0,"electric":const.E0,\
                  "current":const.J0,"density":const.N0,\
                  "temperature":const.T0,"ekbar":const.T0,\
                  "xz_u":const.V0,"xz_t":const.T0,"axis":const.D0,\
                  "derived_j":const.J0}
        normal_s = normal.keys()
        for eachone in normal_s:
            if(eachone in field.lower()):
                factor = normal[eachone]
                break
        return factor
#if magnitude is True, this function will return a vector's module
    def get_module(self,data_dict,field='magnetic',gf=False,g_field=150):
        '''
        This function is used to calculate a vector's module.
        parameters:
        data_dict---a data dictionary read from a sdf file. 
        field-------a vector field, default:magnetic
        gf----------guide field, default:False
        g_field-----guide field, default:150
        '''
        import numpy as np
        mode = 0
        const = 0
        s = data_dict.keys()
        for eachone in s:
            if((field.lower() in eachone.lower()) and \
              ('rank' not in eachone.lower())):
                const = const + 1
                data_field = data_dict[eachone]
                array = data_field.data
                if((gf == True) and (eachone == 'Magnetic_Field_Bz')):
                    array = array - g_field
                mode = mode + array*array
        array = np.sqrt(mode)
        if(const < 3):
            print "Warning! There is not enough componens in this vector!"
        return np.transpose(array)
#get line index
    def get_index(self,data_dict,field='Magnetic_Field_Bx',axes='x'):
        '''
        This function is used to get line data from an array.
        parameters:
        data_dict---data dictionary from a sdf file.
        field-------physical field, default:'Magnetic_Field_Bx'
        axes--------axes,'x' or 'y', default:'x'
        '''
        import numpy as np
        field_data = data_dict[field]
        data = field_data.data
        array = np.transpose(data)
        shape = array.shape
        if (axes == 'x'):
            print "There are %d rows in this field array, please chose proper one."%(shape[0])
        else:
            print "There are %d cloumns in this field array, please chose proper one."%(shape[1])
        index = raw_input("Please input a proper integer number:")
        index = int(index)
        return index
#get line array
    def get_line(self,data_dict,field='Magnetic_Field_Bx',axes='x',index=1):
        '''
        This function is used to get line data from an array.
        parameters:
        data_dict---data dictionary from a sdf file.
        field-------physical field, default:'Magnetic_Field_Bx'
        axes--------axes,'x' or 'y', default:'x'
        index-------array index, default:0
        '''
        import numpy as np
        field_data = data_dict[field]
        data = field_data.data
        array = np.transpose(data)
        if (axes == 'x'):
            vector = array[index-1,:]
        else:
            vector = array[:,index-1]
        return vector
#get module line
    def get_module_line(self,array,axes='x',index=0):
        '''
        This function is  used to get a line from an array.
        parameters:
        array-------input array
        axes--------axes,'x' or 'y', default:'x'
        index-------array index, default:0
        '''
        if (axes == 'x'):
            vector = array[index-1,:]
        else:
            vector = array[:,index-1]
        return vector
#get extent
    def get_extent(self,data_dict):
        '''
        This function is used to get array extent.
        parameters:
        data_dict---data dictionary read from sdf file.
        '''
        from constants import bubble_mr as const
        #import numpy as np
        grid = data_dict["Grid_Grid_mid"]
        #use this method, sometimes contour core dumped
        #grid = data_grid.data
        #x = grid[0]/const.di
        #y = grid[1]/const.di
        #extent = [np.min(x),np.max(x),np.min(y),np.max(y)]
        grid_extent = grid.extents
        xmin = grid_extent[0]/const.D0
        xmax = grid_extent[2]/const.D0
        ymin = grid_extent[1]/const.D0
        ymax = grid_extent[3]/const.D0
        extent = [xmin,xmax,ymin,ymax]
        return extent
#get array
    def get_array(self,data_dict,field='Magnetic_Field_Bx',gf=False,g_field=150):
        '''
        This function is used to get array frim a sdf dictionary.
        parameters:
        data_dict---data dictionary.
        field-------field, default:'Magnetic_Field_Bx'
        gf----------guide field, default:False
        g_field-----guide field, default:150
        '''
        import numpy as np
        data_array = data_dict[field]
        array = data_array.data
        array = np.transpose(array)
        if((gf == True) and (field == 'Magnetic_Field_Bz')):
            array = array - g_field
        return array
#get cordinate
    def get_axis(self,data_dict,axes='x'):
        '''
        This function is used to get 'x', 'y' cordinate.
        parameters:
        data_dict---data dictionary.
        axes--------axes, 'x' or 'y', default:'x'
        '''
        import numpy as np
        extent = sdf_class.get_extent(self,data_dict)
        grid = data_dict["Grid_Grid_mid"]
        dimens = grid.dims
        if(axes == 'x'):
            axis = np.linspace(extent[0],extent[1],dimens[0])
        else:
            axis = np.linspace(extent[2],extent[3],dimens[1])
        return axis
#line integrate
    def line_integrate(self,namelist,field='bx',axes='y',magnitude=False,semi=True,\
                       max_or_min=False,prefix='1'):
        '''
        Thismagnitude---in integrate a vector's module, set True, default:False
 function is used to integrate field along a line,'x' or 'y'.
        parameters:
        namelist----sdf name list.
        field-------physical field to be integrated, default:'bx'
        axes--------axes, 'x' or 'y' ,default:'y'.
        magnitude---in integrate a vector's module, set True, default:False
        semi--------if integrate semi-axis, set True, default:True
        max_or_min--when semi is True, if find extreme min, set False, else set True.
        prefix------file name prefix, default:1
        '''
        import numpy as np
        from constants import bubble_mr as const
        #use sample sdfget parameters
        data_dict = sdf_class.get_data(self,namelist[0],prefix=prefix)
        keys = data_dict.keys()
        if(magnitude == False):
            field = sdf_class.get_field(self,field=field,keys=keys)
        else:
            field_d = field
            field = sdf_class.get_field(self,field=field,keys=keys)
        data = np.array((data_dict[field]).data)
        dimen = data.shape
        a = dimen[1]
        b = dimen[0]
        axis = sdf_class.get_axis(self,data_dict,axes=axes)
        dx = axis[1]-axis[0]
        constant = sdf_class.get_constant(self,field=field)
        n = len(namelist)
        integrate = np.zeros(n)
        #use a loop to integrate filed
        for i in range(n):
            data_dict = sdf_class.get_data(self,namelist[i],prefix=prefix)
            if(magnitude == False):
                array = sdf_class.get_array(self,data_dict,field=field)/constant
            else:
                array = sdf_class.get_module(self,data_dict,field=field_d)/constant
            if(axes == 'x'):
                vector = (array[a/2-1,:] + array[a/2,:] + array[a/2-2,:] + array[a/2+1,:])/4.0
            else:
                vector = (array[:,b/2-1] + array[:,b/2] + array[:,b/2-2] + array[:,b/2+1])/4.0
            if(semi == True):
                if(i <= 26):
                    index = sdf_class.get_local_extreme(self,vector,max_or_min=False)
                else:
                    index = len(vector)/2
                sub_vector = vector[index:]
                #ilength = len(sub_vector)
                integrate[i] = abs(np.sum(sub_vector)*dx)
                #integrate[i] = index
            else:
                integrate[i] = abs(np.sum(vector)*dx/(const.D0))
        return integrate 
#find local min index
    def get_local_extreme(self,vector,max_or_min = False):
        '''
        This function is used to find a vector's local extreme value.
        parameters:
        vector------input dector
        max_or_min--if find extreme min, set False, else set True
        '''
        import numpy as np
        n = len(vector)
        #find start index
        vector_d = vector
        sub1 = np.argmax(vector_d)
        sub2 = np.argmin(vector_d)
        if(sub1 < sub2):
            index = sub1 + np.argmin(abs(vector[sub1:sub2]))
        else:
            index = sub2 + np.argmin(abs(vector[sub2:sub1]))
        return index
        #while(True):
        #   base = vector[index]
        #    left = vector[index-1]
        #    right = vector[index+1]
        #    if(max_or_min == False):
        #        if((left < base) and (base < right)):
        #            index = index-1
        #        elif(left > base) and (base > right):
        #            index = index+1
        #        else:
        #            return index
        #            break
        #    else:
        #        if((left < base) and (base < right)):
        #            index = index+1
        #        elif(left > base) and (base > right):
        #            index = index-1
        #        else:
        #            return index
        #            break
#find current sheet
    def get_file(self,namelist,field='jz',axes='y',magnitude=False,find_max=True,prefix='1'):
        '''
        This function is used to find max value in which file, and return the index.
        parameters:
        namelist----sdf name list.
        field-------physical field to be integrated, default:'jz'.
        axes--------axes, 'x' or 'y' ,default:'y'.
        magnitude---in integrate a vector's module, set True, default:False
        find_max----if to find max file, set True, default:True
        prefix------file name prefix, default:1
        '''
        import numpy as np
        #use sample sdf to get parameters
        data_dict = sdf_class.get_data(self,namelist[0],prefix=prefix)
        keys = data_dict.keys()
        if(magnitude == False):
            field = sdf_class.get_field(self,field=field,keys=keys)
        else:
            field_d = field
            field = sdf_class.get_field(self,field=field,keys=keys)
        data = np.array((data_dict[field]).data)
        dimen = data.shape
        a = dimen[1]
        b = dimen[0]
        constant = sdf_class.get_constant(self,field=field)
        n = len(namelist)
        info = np.zeros((3,n),np.float)
        #use a loop to open each file and find max value an it's location.
        for i in range(n):
            data_dict = sdf_class.get_data(self,namelist[i],prefix=prefix)
            if(magnitude == False):
                array = sdf_class.get_array(self,data_dict,field=field)/constant
            else:
                array = sdf_class.get_module(self,data_dict,field=field_d)/constant
            if(axes == 'x'):
                vector = (array[a/2-1,:] + array[a/2,:])/2.0
            else:
                vector = (array[:,b/2-1] + array[:,b/2])/2.0
            max_value = np.max(vector)
            index = np.argmax(vector)
            #save into info
            info[0,i] = namelist[i]
            info[1,i] = max_value
            info[2,i] = index
        #find what we want
        if(find_max == True):
            index = np.argmax(info[1,:])    
            return (info[0,index],info[2,index])  
        else:
            return info
#find FWHM
    def get_fwhm(self,filenumber=0,field='Current_Jz',axes='x',magnitude=False,index=256,peak=1.0,\
                 constant=1.0,prefix='1'):
        '''
        This function is used to find FWHM.
        parameters:
        filenumber--sdf file.
        field-------physical field to be integrated, default:'Current_Jz'.
        axes--------axes, 'x' or 'y' ,default:'x'.
        magnitude---in integrate a vector's module, set True, default:False.
        index-------line array index, default=256.
        peak--------max value, default=1.0.
        constant----normalization constant.
        prefix------file name prefix, default:1
        '''
        import numpy as np
        data_dict = sdf_class.get_data(self,filenumber,prefix=prefix)
        if(magnitude == True):
            array = sdf_class.get_module(self,data_dict,field=field)/constant
        else:
            array = sdf_class.getdata_array = data_dict[field]
        array = data_array.data

        shape = array.shape
        if(axes == 'x'):
            vector = array[index,:]
            mid = shape[1]/2
        else:
            vector = array[:,index]
            mid = shape[0]
        #then to find FWHM
        vector_abs = np.abs(vector-peak/2.0)
        sub1 = np.argmin(vector_abs)
        vector_abs[sub1] = peak/2.0
        step = 0
        while(True):
            step = step + 1
            sub2 = np.argmin(vector_abs)
            vector_abs[sub2] = peak/2.0
            if(abs(sub2-sub1) > (sub1-mid)):
                break
            else:
                continue
            if(step > 2*mid):
                print 'Find Nothing!'
                break
        if(sub1 < sub2):
            return (sub1,sub2)
        else:
            return (sub2,sub1)
#calculate electron dissipation region.
    def cal_dissipation(self,filenumber,charge=[-1,1,-1,1,-1,1],prefix='1'):
        '''
        This function is used to calculate electron dissipation region.
        parameters:
        filenumber--sdf file name.
        charge------species electric charge, default:[-1,1,1,1].
        prefix------file name prefix, default:1
        '''
        import numpy as np
        from constants import bubble_mr as const
        data_dict = sdf_class.get_data(self,filenumber,prefix=prefix)
        #get constant
        j0 = const.J0
        e0 = const.E0
        b0 = const.B0
        v0 = const.V0
        n0 = const.N0
        qe = const.qe
        c = const.c
        #read current
        jx = sdf_class.get_array(self,data_dict,field='Current_Jx')
        jy = sdf_class.get_array(self,data_dict,field='Current_Jy')
        jz = sdf_class.get_array(self,data_dict,field='Current_Jz')
        #read electric field
        ex = sdf_class.get_array(self,data_dict,field='Electric_Field_Ex_averaged')
        ey = sdf_class.get_array(self,data_dict,field='Electric_Field_Ey_averaged')
        ez = sdf_class.get_array(self,data_dict,field='Electric_Field_Ez_averaged')
        #read magnetic field
        bx = sdf_class.get_array(self,data_dict,field='Magnetic_Field_Bx')
        by = sdf_class.get_array(self,data_dict,field='Magnetic_Field_By')
        bz = sdf_class.get_array(self,data_dict,field='Magnetic_Field_Bz')
        #read electron number density
        #n = len(charge)
        rho_ele1 = sdf_class.get_array(self,data_dict,field='Derived_Number_Density_ele1')
        rho_ele2 = sdf_class.get_array(self,data_dict,field='Derived_Number_Density_ele2')
        rho_ele3 = sdf_class.get_array(self,data_dict,field='Derived_Number_Density_ele3')
        rho_pro1 = sdf_class.get_array(self,data_dict,field='Derived_Number_Density_pro1')
        rho_pro2 = sdf_class.get_array(self,data_dict,field='Derived_Number_Density_pro2')
        rho_pro3 = sdf_class.get_array(self,data_dict,field='Derived_Number_Density_pro3')
        rho = charge[3]*rho_pro1 + charge[4]*rho_pro2 + charge[5]*rho_pro3 + charge[0]*rho_ele1+\
              charge[1]*rho_ele2 + charge[2]*rho_ele3
        #read velocity
        vx = sdf_class.get_array(self,data_dict,field='Derived_xz_ux_averaged')
        vy = sdf_class.get_array(self,data_dict,field='Derived_xz_uy_averaged')
        vz = sdf_class.get_array(self,data_dict,field='Derived_xz_uz_averaged')
        #nie = sdf_class.get_array(self,data_dict,field='Derived_Number_Density_'+species[i])/sdf_class.get_array(self,data_dict,field='Derived_Number_Density_ele')
        #calculate gamma
        v_module = np.sqrt(vx*vx + vy*vy + vz*vz)/c
        gamma = np.sqrt(1/(1-v_module*v_module))
        #calculate dissipation terms
        #j*e
        je = gamma*(jx*ex + jy*ey + jz*ez)
        #j*v*b
        jvb = gamma*(jx*(vy*bz-vz*by) + jy*(vz*bx-vx*bz) + jz*(vx*by-vy*bx))
        #rho*v*e
        rhove = gamma*(qe*rho*(vx*ex + vy*ey + vz*ez))
        #dissipation scaler
        d = je + jvb - rhove
        constant = j0*b0*v0
        return (d/constant,(je+jvb)/constant,je/constant,jvb/constant,-rhove/constant)
#general ohm theory
    def ohm_theory(self,filenumber,axes='x',cut=[512,1024],prefix='1'):
        '''
        This funvtion is used to calculate ohm theory.
        parameters:
        filenumber--sdf file number.
        axes--------axes, 'x' or 'y', default:'x'
        prefix------file name prefix, default:1
        '''
        import numpy as np
        import string
        from constants import bubble_mr as const
        s = 'xy'
        s_index = s.find(axes)
        new_axes = s[1-s_index]
        file_info = sdf_class.get_file(self,namelist=[filenumber],field='jz',axes=new_axes,\
                                       find_max=False)
        index = int(file_info[2])
        a = cut[0]
        b = cut[1]
        data_dict = sdf_class.get_data(self,filenumber,prefix=prefix)
        #read data
        #read magnetic field
        ez = sdf_class.get_array(self,data_dict,field='Electric_Field_Ez_averaged')
        bx = sdf_class.get_array(self,data_dict,field='Magnetic_Field_Bx')
        by = sdf_class.get_array(self,data_dict,field='Magnetic_Field_By')
        jx = sdf_class.get_array(self,data_dict,field='Current_Jx')
        jy = sdf_class.get_array(self,data_dict,field='Current_Jy')
        vx = sdf_class.get_array(self,data_dict,field='Derived_xz_ux_averaged')
        vy = sdf_class.get_array(self,data_dict,field='Derived_xz_uy_averaged')
        rho = sdf_class.get_array(self,data_dict,field='Derived_Number_Density_ele')
        dpxz = sdf_class.get_array(self,data_dict,field='Derived_xz_DxTxz_averaged_ele')
        dpyz = sdf_class.get_array(self,data_dict,field='Derived_xz_DyTyz_averaged_ele')
        if(axes == 'x'):
            ez = (ez[index-1,a:b] + ez[index,a:b] + ez[index+1,a:b])/3.0
            bx = (bx[index-1,a:b] + bx[index,a:b] + bx[index+1,a:b])/3.0
            by = (by[index-1,a:b] + by[index,a:b] + by[index+1,a:b])/3.0
            jx = (jx[index-1,a:b] + jx[index,a:b] + jx[index+1,a:b])/3.0
            jy = (jy[index-1,a:b] + jy[index,a:b] + jy[index+1,a:b])/3.0
            vx = (vx[index-1,a:b] + vx[index,a:b] + vx[index+1,a:b])/3.0
            vy = (vy[index-1,a:b] + vy[index,a:b] + vy[index+1,a:b])/3.0
            rho = (rho[index-1,a:b] + rho[index,a:b] + rho[index+1,a:b])/3.0
            dpxz = (dpxz[index-1,a:b] + dpxz[index,a:b] + dpxz[index+1,a:b])/3.0
            dpyz = (dpyz[index-1,a:b] + dpyz[index,a:b] + dpyz[index+1,a:b])/3.0
        else:
            pass
        e0 = const.E0
        ez = ez/e0
        vb = -(vx*by - vy*bx)/e0
        pxz = -dpxz/rho/const.qe/e0
        pyz = -dpyz/rho/const.qe/e0
        jb = (jx*by - jy*bx)/rho/const.qe/e0
        array = [ez,vb,pxz,pyz,jb]
        #average if necessary
        n = len(ez)
        new_array = array
        for i in range(3):
            for j in range(5):
                for k in range(n-2):
                    new_array[j][k+1] = (array[j][k] + array[j][k+1] + array[j][k+2])/3.0
        return array
#reconnection rate
    def reconnection_rate(self,namelist,field='Electric_Field_Ez_averaged',axes='y',magnitude=False,\
                          semiwidth=2,prefix='1'):
        '''
        This function is used to calculate reconnection rate.
        parameters:
        namelist----sdf name list.
        field-------physical field to be integrated, default:'Electric_Field_Ez_averaged'.
        axes--------axes, 'x' or 'y' ,default:'y'.
        magnitude---in integrate a vector's module, set True, default:False.
        semiwidth---average semi-width, default:2.
        prefix------file name prefix, default:1
        '''
        import numpy as np
        from constants import bubble_mr as const
        #use sample sdf to get parameters
        data_dict = sdf_class.get_data(self,namelist[0],prefix=prefix)
        keys = data_dict.keys()
        if(magnitude == False):
            field = sdf_class.get_field(self,field=field,keys=keys)
        else:
            field_d = field
            field = sdf_class.get_field(self,field=field,keys=keys)
        data = np.array((data_dict[field]).data)
        dimen = data.shape
        a = dimen[1]
        b = dimen[0]
        constant = sdf_class.get_constant(self,field=field)
        n = len(namelist)
        rate = np.zeros(n,dtype=np.float)
        #determine row or cloumn
        if(axes == 'x'):
            index = a/2
        else:
            index = b/2
        for i in range(n):
            data_dict = sdf_class.get_data(self,namelist[i],prefix=prefix)
            #line_1 = sdf_class.get_line(self,data_dict,field=field,axes=axes,index=index-1)
            #line_2 = sdf_class.get_line(self,data_dict,field=field,axes=axes,index=index)
            #line_3 = sdf_class.get_line(self,data_dict,field=field,axes=axes,index=index-2)
            #line_4 = sdf_class.get_line(self,data_dict,field=field,axes=axes,index=index+1)
            #line = (line_1 + line_2 + line_3 + line_4)/4.0
            ez = sdf_class.get_array(self,data_dict,field=field)
            #calculate convect ez:v*b
            bx = sdf_class.get_array(self,data_dict,field='Magnetic_Field_Bx')
            by = sdf_class.get_array(self,data_dict,field='Magnetic_Field_By')
            vx = sdf_class.get_array(self,data_dict,field='Derived_xz_ux_averaged')
            vy = sdf_class.get_array(self,data_dict,field='Derived_xz_uy_averaged')
            ez = ez + (vx*by - vy*bx)
            if(axes == 'y'):
                part_ez = ez[:,index-semiwidth:index+semiwidth-1]
                line_ez = np.sum(part_ez,axis=1)/(2.0*semiwidth)
            else:
                part_ez = ez[index-semiwidth:index+semiwidth-1,:]
                line_ez = np.sum(part_ez,axis=0)/(2.0*semiwidth)
            sub_max = np.argmax(line_ez)
            sub_ez = line_ez[sub_max-semiwidth:sub_max+semiwidth-1]
            rate[i] = np.sum(sub_ez)/(len(sub_ez)*1.0)/constant
        return rate
#general ohm theory
    def dissipation(self,filenumber,prefix='1'):
        '''
        This funvtion is used to calculate ohm theory.
        parameters:
        filenumber--sdf file number.
        prefix------file name prefix, default:1
        '''
        import numpy as np
        import string
        from constants import bubble_mr as const
        data_dict = sdf_class.get_data(self,filenumber,prefix=prefix)
        #read data
        #read magnetic field
        ez = sdf_class.get_array(self,data_dict,field='Electric_Field_Ez_averaged')
        bx = sdf_class.get_array(self,data_dict,field='Magnetic_Field_Bx')
        by = sdf_class.get_array(self,data_dict,field='Magnetic_Field_By')
        vx = sdf_class.get_array(self,data_dict,field='Derived_xz_ux_averaged')
        vy = sdf_class.get_array(self,data_dict,field='Derived_xz_uy_averaged')
        vb = vx*by - vy*bx
        e0 = const.E0
        array = [ez/e0,vb/e0]
        return array
#u*b
    def cal_ub(self,filenumber,component=3,prefix='1'):
        '''
        This function is used to calculate u*b vector
        parameters:
        filenumber--sdf file number.
        component---x, y, z = 1,2,3 respectively, default:3.
        prefix------file name prefix, default:1
        '''
        import numpy as np
        from constants import bubble_mr as const
        data_dict = sdf_class.get_data(self,filenumber,prefix=prefix)
        bx = sdf_class.get_array(self,data_dict,field='Magnetic_Field_Bx')
        by = sdf_class.get_array(self,data_dict,field='Magnetic_Field_By')
        bz = sdf_class.get_array(self,data_dict,field='Magnetic_Field_Bz')
        vx = sdf_class.get_array(self,data_dict,field='Derived_xz_ux_averaged')
        vy = sdf_class.get_array(self,data_dict,field='Derived_xz_uy_averaged')
        vz = sdf_class.get_array(self,data_dict,field='Derived_xz_uz_averaged')
        ex = sdf_class.get_array(self,data_dict,field='Electric_Field_Ex_averaged')
        ey = sdf_class.get_array(self,data_dict,field='Electric_Field_Ey_averaged')
        ez = sdf_class.get_array(self,data_dict,field='Electric_Field_Ez_averaged')
        ub1 = (vy*bz - vz*by)
        ub2 = (vz*bx - vx*bz)
        ub3 = (vx*by - vy*bx)
        ub = [ex+ub1,ey+ub2,ez+ub3]
        return ub[component-1]/const.E0
#charge density
    def charge_density(self,filenumber,species=3,charge=[1,1,1],prefix='1'):
        '''
        This function is used to calculate charge density.
        parameters:
        filenumber--sdf file number.
        species-----pro species.
        charge------electric charge for each species.
        prefix------file name prefix, default:1
        '''
        data_dict = sdf_class.get_data(self,filenumber,prefix=prefix)
        density_ele = sdf_class.get_array(self,data_dict,field='Derived_Number_Density_ele')
        charge_density = -1*density_ele
        for i in range(species):
            field = 'Derived_Number_Density_pro' + str(i+1)
            array = sdf_class.get_array(self,data_dict,field=field)
            charge_density += charge[i]*array
        return charge_density
#calculate electron dissipation region.
    def cal_dissipation_s(self,filenumber,prefix='1'):
        '''
        This function is used to calculate dissipation scaler for different species.
        parameters:
        filenumber--sdf file name.
        prefix------file name prefix, default:1
        '''
        import numpy as np
        from constants import bubble_mr as const
        data_dict = sdf_class.get_data(self,filenumber,prefix=prefix)
        #get constant
        j0 = const.J0
        e0 = const.E0
        b0 = const.B0
        v0 = const.V0
        n0 = const.N0
        qe = const.qe
        c = const.c
        #read current
        jx = sdf_class.get_array(self,data_dict,field='Current_Jx')
        jy = sdf_class.get_array(self,data_dict,field='Current_Jy')
        jz = sdf_class.get_array(self,data_dict,field='Current_Jz')
        #read electric field
        ex = sdf_class.get_array(self,data_dict,field='Electric_Field_Ex_averaged')
        ey = sdf_class.get_array(self,data_dict,field='Electric_Field_Ey_averaged')
        ez = sdf_class.get_array(self,data_dict,field='Electric_Field_Ez_averaged')
        #read magnetic field
        bx = sdf_class.get_array(self,data_dict,field='Magnetic_Field_Bx')
        by = sdf_class.get_array(self,data_dict,field='Magnetic_Field_By')
        bz = sdf_class.get_array(self,data_dict,field='Magnetic_Field_Bz')
        species = ['ele','pro1','pro2','pro3']
        n = len(species)
        #read electron number density
        density_s = []
        rho = 0
        for i in range(n):
            density = sdf_class.get_array(self,data_dict,field='Derived_Number_Density_'+species[i])
            density_s = density_s + [density]
            rho = rho + density_s[i]
        rho = rho - 2*density_s[0]
        #read velocity and calculate D
        D = []
        total_pro = 0
        for i in range(n):
            vx = sdf_class.get_array(self,data_dict,field='Derived_xz_ux_averaged_'+species[i])
            vy = sdf_class.get_array(self,data_dict,field='Derived_xz_uy_averaged_'+species[i])
            vz = sdf_class.get_array(self,data_dict,field='Derived_xz_uz_averaged_'+species[i])
            nie = density_s[i]/density_s[0]
            #calculate gamma
            v_module = np.sqrt(vx*vx + vy*vy + vz*vz)/c
            gamma = np.sqrt(1/(1-v_module*v_module))
            #calculate dissipation terms
            #j*e
            je = nie*gamma*(jx*ex + jy*ey + jz*ez)
            #j*v*b
            jvb = nie*gamma*(jx*(vy*bz-vz*by) + jy*(vz*bx-vx*bz) + jz*(vx*by-vy*bx))
            #rho*v*e
            rhove = nie*gamma*(qe*rho*(vx*ex + vy*ey + vz*ez))
            #dissipation scaler
            d = je + jvb - rhove
            constant = j0*b0*v0
            D = D + [d/constant]
            if(i > 0):
                total_pro = total_pro + d/constant
        D = D + [total_pro]
        return D
#energy spectrum
    def cal_enspe(self,filenumber,species='ele1',info='momentum',ndomain=500,prefix='2',\
                  mass_ratio=1,g_max=1.05,average=False):
        '''
        This function is used to calculate energy spectrum.
        parameters:
        filenumber--sdf file name.
        species-----pro species, default:ele1
        info--------particle information, default:momentum.
        ndomain-----axis step, default:500.
        prefix------file name prefix, default:2
        mass_ratio--mass ratio to electron, default:1
        g_max-------max gamma, default:1.05
        average-----if average, default:False
        '''
        from constants import bubble_mr as const
        import numpy as np
        data_dict = sdf_class.get_data(self,filenumber,prefix=prefix)
        #read data
        if(info == 'momentum'):
            m = mass_ratio*const.me
            p0 = m*const.c
            px = sdf_class.get_array(self,data_dict,field='Particles_Px_subset_'+species[:3]+'_'+species)
            py = sdf_class.get_array(self,data_dict,field='Particles_Py_subset_'+species[:3]+'_'+species) 
            pz = sdf_class.get_array(self,data_dict,field='Particles_Pz_subset_'+species[:3]+'_'+species)
            p = np.sqrt(px*px + py*py + pz*pz)
            gamma = np.sqrt(1 + (p/p0)**2)
        if(info == 'gamma'):
            gamma = sdf_class.get_array(self,data_dict,field='Particles_Gamma_subset_'+species[:3]+'_'+species)
        if(info == 'energy'):
            m = mass_ratio*const.me
            e0 = m*const.c*const.c
            energy = sdf_class.get_array(self,data_dict,field='Particles_Energy_subset_'+species[:3]+'_'+species)
            gamma = energy/e0
        #calculate particle number
        weight = sdf_class.get_array(self,data_dict,field='Particles_Weight_subset_'+species[:3]+'_'+species)
        #dg = (max(gamma) - min(gamma))/(ndomain*1.0)
        #g_test = max(gamma)
        #if(g_test < g_max):
        #    g_max = g_test
        #n_step = int((max(gamma) - min(gamma))/ndomain)
        dg = np.linspace(min(gamma),max(gamma),ndomain+1)
        dn = np.zeros(ndomain+1,np.float)
        dgx = dg[1] - dg[0]
        lg = len(gamma)
        ld = len(dg)
        particles = 0
#        for i in range(ndomain):
#            nparticle = 0
#            for j in range(lg):
#                if((gamma[j] > dg[i]) and (gamma[j] <= dg[i+1])):
#                       nparticle += weight[j]
#                       particles += 1
#            dn[i] = nparticle
        for j in range(lg):
            for i in range(ndomain):
                if((gamma[j] > dg[i]) and (gamma[j] <= dg[i+1])):
                    dn[i] += weight[j]
                    break
        #average, if True
        if(average == True):
            for each in range(1,ndomain-1):
                dn[each] = (dn[each-1] + dn[each] + dn[each+1])/3.0
        total_n = sum(dn)
        for each in range(ndomain+1):
            dn[each] = dn[each]/float(total_n)/dgx
            dg[each] = dg[each] - 1
        return (dg,dn,total_n)
#read 3d array
    def get_array_3d(self,data_dict,field='Magnetic_Field_Bx',axis='x',index=0):
        '''
        This function is used to get array frim a sdf dictionary.
        parameters:
        data_dict---data dictionary.
        field-------field, default:'Magnetic_Field_Bx'
        axis--------slice axis, default:'x'
        index-------slice index, default:0
        '''
        import numpy as np
        data_array = data_dict[field]
        array = data_array.data
        if(axis == 'x'):
            resu = array[index,:,:]
        elif(axis == 'y'):
            resu = array[:,index,:]
        else:
            resu = array[:,:,index]
        return np.transpose(resu)
#get 3d extent
    def get_extent_3d(self,data_dict,axis='x'):
        '''
        This function is used to get array extent.
        parameters:
        data_dict---data dictionary read from sdf file.
        axis--------slice axis, default:'x'
        '''
        from constants import laser_mr as const
        #import numpy as np
        grid = data_dict["Grid_Grid_mid"]
        grid_extent = grid.extents
        xmin = grid_extent[0]/const.la
        xmax = grid_extent[3]/const.la
        ymin = grid_extent[1]/const.la
        ymax = grid_extent[4]/const.la
        zmin = grid_extent[2]/const.la
        zmax = grid_extent[5]/const.la
        if(axis == 'x'):
            extent = [ymin,ymax,zmin,zmax]
        elif(axis == 'y'):
            extent = [xmin,xmax,zmin,zmax]
        else:
            extent = [xmin,xmax,ymin,ymax]
        return extent
#use a dictionary to chose the normalization constant
    def get_constant_3d(self,field='Magnetic_Field_Bx'):
        '''
        This function is used to chose the normalization constant according to the field.
        parameters:
        field-------the field name. default:'Magnetic_Field_Bx'
        '''
        from constants import laser_mr as const
        normal = {"magnetic":const.B0,"electric":const.E0,\
                  "current":const.J0,"axis":const.la,\
                  "derived_j":const.J0,"density":const.nc}
        normal_s = normal.keys()
        for eachone in normal_s:
            if(eachone in field.lower()):
                factor = normal[eachone]
                break
        return factor
#calculate magnetic field tensor force
    def tensor_force(self,filenumber,prefix='1',component=2):
        '''
        This function is used to calculate magnetic field tensor force.
        parameters:
        filenumber--sdf file name.
        prefix------file name prefix, default:1
        component---vector component, default:2(y)
        '''
        from constants import bubble_mr as const
        import numpy as np
        data_dict = sdf_class.get_data(self,filenumber,prefix=prefix)
        bx = sdf_class.get_array(self,data_dict,field='Magnetic_Field_Bx')
        by = sdf_class.get_array(self,data_dict,field='Magnetic_Field_By')
        bz = sdf_class.get_array(self,data_dict,field='Magnetic_Field_Bz')
        b = [bx,by,bz]
        #get dx
        grid = data_dict['Grid_Grid_mid']
        extent = grid.extents
        [m,n] = bx.shape
        dx = (extent[2] - extent[0])/(1.0*n)
        force = np.zeros([m,n])
        #for loop to calculate force
        for i in range(1,m-2):
            for j in range(1,n-2):
                force[i,j] = b[0][i,j]*(b[component-1][i,j+1] - b[component-1][i,j-1])/2.0/dx + \
                             b[1][i,j]*(b[component-1][i+1,j] - b[component-1][i-1,j])/2.0/dx
        force[0,:] = force[1,:]
        force[m-1,:] = force[m-2,:]
        force = force/(const.B0*const.B0/const.di)
        return force
#calculate jb
    def cal_jb(self,filenumber,component=3,prefix='1'):
        '''
        This function is used to calculate u*b vector
        parameters:
        filenumber--sdf file number.
        component---x, y, z = 1,2,3 respectively, default:3.
        prefix------file name prefix, default:1
        '''
        import numpy as np
        from constants import bubble_mr as const
        data_dict = sdf_class.get_data(self,filenumber,prefix=prefix)
        bx = sdf_class.get_array(self,data_dict,field='Magnetic_Field_Bx')
        by = sdf_class.get_array(self,data_dict,field='Magnetic_Field_By')
        bz = sdf_class.get_array(self,data_dict,field='Magnetic_Field_Bz')
        jx = sdf_class.get_array(self,data_dict,field='Current_Jx')
        jy = sdf_class.get_array(self,data_dict,field='Current_Jy')
        jz = sdf_class.get_array(self,data_dict,field='Current_Jz')
        rho_ele1 = sdf_class.get_array(self,data_dict,field='Derived_Number_Density_ele1')
        rho_ele2 = sdf_class.get_array(self,data_dict,field='Derived_Number_Density_ele2')
        rho_ele3 = sdf_class.get_array(self,data_dict,field='Derived_Number_Density_ele3')
        rho = rho_ele1 + rho_ele2 + rho_ele3
        jb1 = (jy*bz - jz*by)
        jb2 = (jz*bx - jx*bz)
        jb3 = (jx*by - jy*bx)
        jb = [jb1,jb2,jb3]
        return jb[component-1]/rho/const.qe/const.E0
#particle information
    def get_parinfo(self,namelist,field='Px',species='ele',prefix='3',component=1,\
                    id_index=1,mass_ratio=1,time_factor=1.0):
        '''
        This function is used to get particle information.
        parameters:
        namelist----sdf file list.
        field-------physical field, default:Px. #particles_px_subset_tracer_p_tracer_ele.
        species-----species, default:ele
        prefix------file name prefix, default:3
        component---axis, default:1(x).
        id_index----particle id.
        mass_ratio--mass ratio, default:1
        time_factor-time factor, an output sdf file represent time step, default:1.
        '''
        import numpy as np
        from constants import bubble_mr as const
        n = len(namelist)
        array = np.zeros(n,np.float) 
        vz = np.zeros(n,np.float) 
        for i in range(n):
            filenumber = namelist[i]
            data_dict = sdf_class.get_data(self,filenumber,prefix=prefix)
            #get ID
            ID = data_dict['Particles_ID_subset_tracer_p_tracer_'+species]
            ID = ID.data
            len_ID = len(ID)
            if((len_ID < id_index) and (i == 0)):
                print 'Too large id_index!'
                break
            else: 
                if((species == 'pro') and (i == 0)):
                    id_index += len_ID
                for j in range(len_ID):
                    if(id_index == ID[j]):
                        index = j
                if(field == 'Grid'):
                    if(component < 3):
                        grid = data_dict['Grid_Particles_subset_tracer_p_tracer_'+species]
                        grid = grid.data
                        #print id_index,'    ',index
                        array[i] = grid[component-1][index]/const.di
                    else:
                        pz = data_dict['Particles_Pz_subset_tracer_p_tracer_'+species]
                        g = data_dict['Particles_Gamma_subset_tracer_p_tracer_'+species]
                        pz = pz.data
                        g = g.data
                        vz[i] = pz[index]/g[index]/const.me
                        if(i > 0):
                            array[i] = array[i-1] + ((vz[i] + vz[i-1])/2.0/const.omi*float(time_factor))/const.di
                else:
                    data = data_dict['Particles_'+field+'_subset_tracer_p_tracer_'+species]
                    data = data.data
                    #print id_index,'    ',index
                    array[i] = data[index]/float(const.Pe0*mass_ratio)
        return array

