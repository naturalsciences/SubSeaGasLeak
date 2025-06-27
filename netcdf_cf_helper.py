from netCDF4 import Dataset
import numpy as np

class netcdf_file:

    def __init__(self, path, title, lons, lats):
        """
        Object for creating a standard netcdf cf file.

        params
        -----
        path : the path of the file
        title : the title of the file
        lons : an array with all the longitude for the data
        lats : an array with all the latitude for the data
        """
        self.nc_file = Dataset(path, 'w', format='NETCDF4_CLASSIC')

        self.nc_file.Conventions = "CF-1.8"
        self.nc_file.title = title
        self.nc_file.source = "OSERIT"
        self.nc_file.institution = 'RBINS'
        self.create_dimension('lon', "longitude", "longitude", "degrees_east", lons)
        self.create_dimension('lat', "latitude", "latitude", "degrees_north", lats)


    def add_time(self, times):
        """
        Add the time dimension to the file, times is an array with
        all the times since 1970-1-1 00:00:00
        """
        self.create_dimension('time', 'time', 'time', "seconds since 1970-01-01 00:00:00", times)

    def create_dimension(self, name, std_name, lg_name, units, data):
        """
        Create a dimension and the variable linked to it

        params
        -----
        name is the dimension/variable name
        std_name is the standard_name of the variable
        lg_name is the long_name of the variable
        units is the units of the variable
        data is the array to be put as dimension
        """
        self.nc_file.createDimension(name,len(data))
        self.nc_file.createVariable(name,np.float64, (name,))
        var = self.nc_file.variables[name]
        var.standard_name = std_name
        var.long_name = lg_name
        var.units = units
        self.nc_file.variables[name][:] = data


    def add_info_var(self,name, std_name, lg_name, units, valid_range = [0,1000000000000]):
        """
        Add the attribute to a variable

        params
        -----
        name : the variable name
        std_name : the standard_name of the variable
        lg_name : the long_name of the variable
        units : the units of the variable
        valid_range : the valid_range, by default [0,1000000000000]
        """
        var = self.nc_file.variables[name]
        var.standard_name = std_name
        var.long_name = lg_name
        var.units = units
        var.valid_range = valid_range

    def add_variable_lon_lat(self, name, data, fill_value = 0):
        """
        Add a variable 2d to the file
        name is the variable name
        data is the array(lon,lat) to put in the netcdf file
        fill_value is the fill value, by default 0
        """
        self.nc_file.createVariable(name,np.float64, ('lat','lon'), fill_value = fill_value)
        self.nc_file.variables[name][:,:] = np.transpose(data, (1,0))


    def add_variable_lon_lat_time(self, name, data, fill_value = 0):
        """
        Add a variable 3d with the the time to the file (lon, lat, time)

        params
        -----
        name : the variable name
        data : the array(lon,lat) to put in the netcdf file
        fill_value : the fill value, by default 0
        """
        self.nc_file.createVariable(name,np.float64, ('time','lat','lon'), fill_value = fill_value)
        self.nc_file.variables[name][:,:,:] = np.transpose(data, (2,1,0))

    def add_variable_lon_lat_x(self, name, data, dim, fill_value = 0):
        """
        Add a variable 3d with the a arbitrary dimension to the file (lon, lat dim)

        params
        -----
        name : the variable name
        data : the array(lon,lat) to put in the netcdf file
        dim : the arbitrary dimension
        fill_value : the fill value, by default 0
        """
        self.nc_file.createVariable(name,np.float64, (dim,'lat','lon'), fill_value = fill_value)
        self.nc_file.variables[name][:,:,:] = np.transpose(data, (2,1,0))


    def close(self):
        self.nc_file.close()

def combine_files(paths_in, path_out,title):
    """
    Combine netcdf files in the list paths_in together in the "path_out"
    support only lon lat time and concentration as commun dimension
    """
    datasets = []
    for path in paths_in:
        datasets.append(Dataset(path))

    lons_id = np.full((len(datasets)), -1) #-1 means this variable is unused
    lats_id = np.full((len(datasets)), -1)
    times_id = np.full((len(datasets)), -1)
    concs_id = np.full((len(datasets)), -1)

    lons = []
    lats = []
    times = []
    concs = []

    for index, dataset in enumerate(datasets):
        try:
            find_dim(dataset,index,'lon',lons, lons_id)
        except IndexError as e:
            pass
        try:
            find_dim(dataset,index,'lat',lats, lats_id)
        except IndexError as e:
            pass
        try:
            find_dim(dataset,index,'time',times, times_id)
        except IndexError as e:
            pass
        try:
            find_dim(dataset,index,'concentration_thr',concs, concs_id)
        except IndexError as e:
            pass

    nc_file = Dataset(path_out, 'w', format='NETCDF4_CLASSIC')

    nc_file.Conventions = "CF-1.8"
    nc_file.title = title
    #nc_file.source = "OSERIT"
    nc_file.institution = 'RBINS'

    for i in range(len(lons)):
        create_dim(nc_file, f'lon{i}', "longitude", "longitude", "degrees_east", lons[i])
    for i in range(len(lats)):
        create_dim(nc_file, f'lat{i}', "latitude", "latitude", "degrees_north", lats[i])
    for i in range(len(times)):
        create_dim(nc_file, f'time{i}', 'time', 'time', "seconds since 1970-01-01 00:00:00", times[i])
    for i in range(len(concs)):
        create_dim(nc_file, f'concentration_thr{i}', "concentration_threshold", "concentration threshold ppmv", "1e-6", concs[i])

    for i_d, dataset in enumerate(datasets):
        for i_var, var in enumerate(dataset.variables):
            if str(var) not in ['lon', 'lat', 'time', 'concentration_thr']:
                dim = list(dataset[var].dimensions)
                for j in range(len(dim)):
                    if dim[j] == 'lat':
                        dim[j] = f'lat{lats_id[i_d]}'
                    elif dim[j] == 'lon':
                        dim[j] = f'lon{lons_id[i_d]}'
                    elif dim[j] == 'time':
                        dim[j] = f'time{times_id[i_d]}'
                    elif dim[j] == 'concentration_thr':
                        dim[j] = f'concentration_thr{concs_id[i_d]}'

                nc_file.createVariable(var,np.float64,tuple(dim), fill_value = dataset[var]._FillValue)

                if len(dim) == 1:
                    nc_file.variables[var][:] = dataset[var][:]
                elif len(dim) == 2:
                    nc_file.variables[var][:,:] = dataset[var][:,:]
                elif len(dim) == 3:
                    nc_file.variables[var][:,:,:] = dataset[var][:,:,:]
                elif len(dim) == 4:
                    nc_file.variables[var][:,:,:,:] = dataset[var][:,:,:,:]
                nc_file.variables[var].standard_name = dataset[var].standard_name
                nc_file.variables[var].long_name = dataset[var].long_name
                nc_file.variables[var].units = dataset[var].units
                nc_file.variables[var].valid_range = dataset[var].valid_range
        dataset.close()
    nc_file.close()

def find_dim(dataset, index, name_dim, list_val, list_id):
    """
    Find the dimension of name_dim in dataset at the index, ad the values at list_val
    if not already defined and at list_id to find where it is
    """
    data = dataset[name_dim][:]
    i = -1
    find = False
    for i in range(len(list_val)):
        if len(list_val[i]) == len(data):
            if (list_val[i] == data).all():
                find = True
                break

    if not find:
        list_val.append(data)
        list_id[index] = i+1
    else:
        list_id[index] = i

def create_dim(file, name, std_name, lg_name, units, data):
    """
    Create a dimension and the variable linked to it

    params
    -----
    name is the dimension/variable name
    std_name is the standard_name of the variable
    lg_name is the long_name of the variable
    units is the units of the variable
    data is the array to be put as dimension
    """
    file.createDimension(name,len(data))
    file.createVariable(name,np.float64, (name,))
    var = file.variables[name]
    var.standard_name = std_name
    var.long_name = lg_name
    var.units = units
    file.variables[name][:] = data

