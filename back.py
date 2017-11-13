
	if (dims == 3):
		nx,ny,nz = data.shape;
		if (plain == 'xy'): 
			ddata = data[:,:,index];
			ddata = ddata.T;
			new_extent = extent[0:4];
			xlabel = label[0];
			ylabel = label[1];
			myaspect = nx/ny;
		if (plain == 'yz'):
			ddata = data[index,:,:].reshape(ny,nz);
			new_extent = extent[2:6];
			xlabel = label[1];
			ylabel = label[2];
		if (plain == 'xz'):
			ddata = data[:,index,:].reshape(nx,nz);
			ddata = ddata.T;
			new_extent = (extent[0],extent[1],extent[4],extent[5]);
			aspect = nz/nx;
			xlabel = label[0];
			ylabel = label[2];
	else:
		nx,ny = data.shape;
		ddata = ddata.T;
		new_extent = extent;
		aspect = nx/ny;
		xlabel = label[0];
		ylabel = label[1];
	title = label[3];	
