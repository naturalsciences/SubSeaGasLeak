import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import netcdf_cf_helper as nc_h
import math
from shapely.geometry import Polygon, box
import os
import zipfile

class Particle:
    def __init__(self):
        self.x =[]
        self.y = []
        self.z = []
        self.radius = []
        self.height = []
        self.number_bubble = []
        self.mole_content = []
        self.timesteps = []

class Release:
    def __init__(self, x, y, radius, quantity, time):
        self.x = x
        self.y = y
        self. radius = radius
        self.quantity = quantity
        self.time = time

def rectangle_area(X, Y):
    # For a rectangle with coordinates in clockwise or counterclockwise order, Shoelace formula
    n = len(X)
    area = 0
    
    for i in range(n):
        j = (i + 1) % n
        area += X[i] * Y[j]
        area -= Y[i] * X[j]
    
    return abs(area) / 2

def calculate_grid_coverage(rect_x, rect_y, grid_x, grid_y):
    """
    Calculate the fraction of rectangle area in each grid cell.
    
    Parameters:
    - rect_x: Vector of x-coordinates for rectangle vertices
    - rect_y: Vector of y-coordinates for rectangle vertices
    - grid_x: Vector of x-coordinates defining grid cell boundaries
    - grid_y: Vector of y-coordinates defining grid cell boundaries
    
    Returns:
    - 2D numpy array with area fractions for each grid cell
    - Total sum of fractions will be 1.0
    """
    # Validate input
    if len(rect_x) != 4 or len(rect_y) != 4:
        raise ValueError("Rectangle must have exactly 4 points")
    
    # Create polygon from rectangle points
    rect_points = list(zip(rect_x, rect_y))
    rect_polygon = Polygon(rect_points)
    
    # Calculate total rectangle area
    total_rect_area = rect_polygon.area
    
    # Create grid coverage array
    coverage = np.zeros((len(grid_x)-1, len(grid_y)-1))
    
    # Iterate through each grid cell
    for x in range(len(grid_x)-1):
        for y in range(len(grid_y)-1):
            # Create a box representing the grid cell using exact grid coordinates
            grid_cell = box(grid_x[x], grid_y[y], grid_x[x+1], grid_y[y+1])
            
            # Calculate intersection area
            intersection = rect_polygon.intersection(grid_cell)
            
            # Store fraction of rectangle area in this cell
            coverage[x, y] = intersection.area / total_rect_area
    
    return coverage

def load_layers(path_layers):
    """
    Load layers in a csv to several numpy array
    """
    data_array = np.loadtxt(
            path_layers,
            skiprows=1,
            delimiter = ","
        )
    layer_id = np.unique(data_array[:,1])
    time_id = np.unique(data_array[:,0])
    #create layer matrices and populate them
    layer_matrix_x = np.zeros((len(layer_id), len(time_id), 4))
    layer_matrix_y = np.zeros((len(layer_id), len(time_id), 4))
    layer_matrix_quantity = np.zeros((len(layer_id), len(time_id), int((data_array.shape[1]-10)/2)))
    layer_release = np.zeros((len(time_id),int((data_array.shape[1]-10)/2)))
    # Ensure layer_release is a 2D array
    if layer_release.ndim == 1:
        layer_release = layer_release.reshape(len(time_id), -1)

    for i in range(data_array.shape[0]):
        id_time = (np.where(time_id == data_array[i,0]))[0][0]
        id_part = (np.where(layer_id == data_array[i,1]))[0][0]

        layer_matrix_x[id_part, id_time] = data_array[i,2:6]
        layer_matrix_y[id_part, id_time] = data_array[i,6:10]
        layer_matrix_quantity[id_part, id_time] = data_array[i,10:10+int((data_array.shape[1]-10)/2)]
        layer_release[id_time, :] = data_array[i,10+int((data_array.shape[1]-10)/2):data_array.shape[1]]

    return layer_matrix_x, layer_matrix_y, layer_matrix_quantity, layer_release, time_id

def load_depth_layers(file_path):
    """
    Load a JSON file and extract the depth_layers as a numpy array.
    
    Parameters:
    file_path (str): Path to the JSON file
    
    Returns:
    numpy.ndarray: Numpy array containing the depth_layers values
    """
    try:
        # Open and load the JSON file
        with open(file_path, 'r') as file:
            data = json.load(file)
        
        # Extract the depth_layers array
        depth_layers = data.get('depth_layers')
        
        if depth_layers is None:
            raise KeyError("The JSON file does not contain a 'depth_layers' key")
        
        # Convert to numpy array
        depth_array = np.array(depth_layers)
        return depth_array
    
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found")
        return None
    except json.JSONDecodeError:
        print(f"Error: File '{file_path}' is not valid JSON")
        return None
    except Exception as e:
        print(f"Error: {str(e)}")
        return None
    
def filter_strings_starting_with(strings, prefix):
    """
    only select strings with prefix in a list of string
    """
    return [s for s in strings if s.startswith(prefix)]

def remove_prefix_from_list(strings, prefix):
    """
    remove a prefix from a string in a list
    """
    return [s[len(prefix):] if s.startswith(prefix) else s for s in strings]

def load_particles(path_particle):
    with open(path_particle, 'r') as f:
        particle_data = json.load(f)

    array_tmtsp = particle_data['simulation']
    array_time = np.zeros(len(array_tmtsp))
    
    keys = list(array_tmtsp[0].keys())
    list_key_gas = filter_strings_starting_with(keys, "mol_bubble_part")
    name_gas = remove_prefix_from_list(list_key_gas, "mol_bubble_part")
    particles = []
    for i in range(len(array_tmtsp)):
        array_time[i] = array_tmtsp[i]['current_time']

        active_part = np.where(np.array(array_tmtsp[i]["active_part"]) == 1)

        for index_p in active_part[0]:
            id_part = array_tmtsp[i]["part_id"][index_p]
            while len(particles) < id_part+1:
                particles.append(Particle())

            list_of_mol = []
            for gas in list_key_gas:
                list_of_mol.append(array_tmtsp[i][gas][index_p])
            
            particles[id_part].x.append(array_tmtsp[i]["x_part"][index_p])
            particles[id_part].y.append(array_tmtsp[i]["y_part"][index_p])
            particles[id_part].z.append(array_tmtsp[i]["z_part"][index_p])
            particles[id_part].radius.append(array_tmtsp[i]["radius_part"][index_p])
            particles[id_part].height.append(array_tmtsp[i]["height_part"][index_p])
            particles[id_part].number_bubble.append(array_tmtsp[i]["number_bubble_part"][index_p])
            particles[id_part].mole_content.append(list_of_mol)
            particles[id_part].timesteps.append(array_time[i])
    return particles, name_gas, array_time

def load_releases(path_release):
    releases = []
    with open(path_release, 'r') as f:
        release_data = json.load(f)
        array_tmtsp = release_data['simulation']
        for tmstp in array_tmtsp:
            time = tmstp["current_time"]
            for release in tmstp["releases"]:
                releases.append(Release(release["x"], release["y"], release["radius_part"], release["quantity_released"], time))
    return releases

    
def save_quantity_layers(time, z_layers, quantity, title, units, path):
    plt.figure(figsize=(12, 6))


    # For flat shading, we need the coordinates to be one larger than the data dimensions
    time_centers = np.linspace(time[0], time[-1], quantity.shape[0] + 1)
    depth_centers = z_layers  # Already correct size
    time_mesh2, depth_mesh2 = np.meshgrid(time_centers, depth_centers)
    pcm = plt.pcolormesh(time_mesh2, depth_mesh2, quantity.T, 
                        cmap='viridis', shading='flat')

    # Add a colorbar
    cbar = plt.colorbar(pcm)
    cbar.set_label(units)

    # Set labels and title
    plt.xlabel('Time [s]')
    plt.ylabel('Depth [m]')
    plt.title(title)

    plt.tight_layout()
    plt.savefig(path)

def q_released_bubble(times, releases, id_g):
    amount_vola = np.zeros(len(times))
    for release in releases:
        id_time = np.searchsorted(times, release.time, side='left')-1 
        amount_vola[id_time] += release.quantity[id_g]

    return np.cumsum(amount_vola)

def q_released_layer(times, time_l, releases_l, id_g):
    amount_vola = np.zeros(len(times))
    for i in range(len(time_l)):
        id_time = np.searchsorted(times, time_l[i], side='left')-1 
        amount_vola[id_time] += releases_l[i,id_g]

    return np.cumsum(amount_vola)

def q_dis_layer(times, time_l, layer_q, id_g):
    amount_dis = np.zeros(len(times))
    for i in range(len(times)):
        for j in range(len(time_l)):
            if times[i] == time_l[j]:
                for k in range(len(layer_q[:,j,0])):
                    amount_dis[i] += layer_q[k,j,id_g]
                continue

    return amount_dis

def q_fresh(times, parts, id_g):
    amount_fresh = np.zeros(len(times))
    for part in parts:
        for i in range(len(part.timesteps)):
            found = False
            for j in range(len(times)):
                if part.timesteps[i] - times[j] < 1e-6:
                   amount_fresh[j] += part.mole_content[i][id_g] * part.number_bubble[i]
                   found = True
                   break
            if not found:
                print(part.timesteps[i], times[-1])
    return  amount_fresh

def compute_vola(layer_x_all, layer_y_all, layer_release_all, time_layers, releases, x_res, y_res, alias_factor, times, id_g):
    
    layer_release = layer_release_all[:, id_g]
    layer_x = layer_x_all[-1,:,:]
    layer_y = layer_y_all[-1,:,:]
    quantity_tr_bubble = np.zeros((len(times), x_res, y_res))
    quantity_tr_layer = np.zeros((len(times), x_res, y_res))
    added_factor = 1.05 # make the grid a bit bigger than the circle

    max_x = -1e10
    min_x = 1e10
    max_y = -1e10
    min_y = 1e10

    if len(releases) > 0:
        for release in releases:
            if release.x - release.radius * added_factor < min_x:
                min_x = release.x - release.radius * added_factor
            if release.x + release.radius * added_factor > max_x:
                max_x = release.x + release.radius * added_factor

            if release.y - release.radius * added_factor < min_y:
                min_y = release.y - release.radius * added_factor
            if release.y + release.radius * added_factor > max_y:
                max_y = release.y + release.radius * added_factor
    for i in range(len(layer_release)):
        if layer_release[i] > 0:
            if np.min(layer_x[i,:]) < min_x:
                min_x = np.min(layer_x[i,:])
            if np.min(layer_y[i,:]) < min_y:
                min_y = np.min(layer_y[i,:])
            if np.max(layer_x[i,:]) > max_x:
                max_x = np.max(layer_x[i,:])
            if np.max(layer_y[i,:]) > max_y:
                max_y = np.max(layer_y[i,:])

    dx = (max_x-min_x) / x_res
    dy = (max_y-min_y) / y_res

    dx_sub = dx/alias_factor
    dy_sub = dy/alias_factor

    xs = np.linspace(min_x + dx_sub/2, max_x - dx_sub/2, x_res * alias_factor)
    ys = np.linspace(min_y + dy_sub/2, max_y - dy_sub/2, y_res * alias_factor)


    is_on_circle = np.zeros((len(xs), len(ys)), dtype=bool)
    #adding particle
    for release in releases:
        id_time = np.searchsorted(times, release.time, side='left')-1
        is_on_circle[:] = False
        d_horiz = release.x - xs
        d_vert = release.y - ys
        # making arrays the good shape
        d_horiz = np.tile(d_horiz, (y_res * alias_factor, 1)).T
        d_vert = np.transpose(np.tile(d_vert, (x_res * alias_factor, 1)).T)

        distsq = (d_horiz)**2 + (d_vert)**2


        is_on_circle[distsq < release.radius**2] = True
        alias_factor_mask = np.zeros((x_res, y_res), dtype=int)
        for k in range(x_res):
            for l in range(y_res):
                alias_factor_mask[k, l] = np.sum(is_on_circle[k*alias_factor:k*alias_factor + alias_factor, l*alias_factor:l*alias_factor + alias_factor])
        alias_factor_mask = alias_factor_mask / alias_factor**2
        total = np.sum(alias_factor_mask)
        if total > 0:
            quantity_tr_bubble[id_time,:,:] += alias_factor_mask * release.quantity[id_g] / total

    #adding layers
    dx = (max_x-min_x) / x_res
    dy = (max_y-min_y) / y_res

    xs = np.linspace(min_x , max_x+dx , x_res+1)
    ys = np.linspace(min_y , max_y+dy , y_res+1)
    for i in range(len(time_layers)):
        
        if layer_release[i] == 0:
            continue

        id_time = np.searchsorted(times, time_layers[i], side='left')-1
        
        quantity_tr_layer[id_time,:,:] += layer_release[i] * calculate_grid_coverage(layer_x[i,:], layer_y[i,:], xs, ys)

    return quantity_tr_layer, quantity_tr_bubble, np.linspace(min_x, max_x, x_res), np.linspace(min_y, max_y, y_res)

def save_vola(quantity, xs, ys, lon, lat, name, unit, times, path):
    EARTH_RADIUS = 6371000 #[m]
    m_to_deg_lat = 360/(EARTH_RADIUS*math.pi*2)
    lat_radius = math.cos(lat*math.pi /180) * EARTH_RADIUS
    m_to_deg_lon = 360 / (lat_radius * math.pi * 2)

    lons = xs * m_to_deg_lon + lon
    lats = ys * m_to_deg_lat + lat

    nc_file_h = nc_h.netcdf_file(path,
                            f"bubble_rising_output_{name}_volatilization",
                            lons,
                            lats)
    nc_file_h.add_time(times)
    nc_file_h.add_variable_lon_lat_time(f"{name}_volatilization",np.transpose(quantity, (1,2,0)))
    nc_file_h.add_info_var(f"{name}_volatilization", f"{name}_volatilization", f"{name} volatilization", unit)
    nc_file_h.add_variable_lon_lat(f"{name}_volatilization_total",np.sum(quantity,0))
    nc_file_h.add_info_var(f"{name}_volatilization_total", f"{name}_volatilization_total", f"{name} volatilization total", unit)

    q_sqm = quantity/((xs[1]-xs[0])*(ys[1]-ys[0]))
    flow_rate_sqm = q_sqm / np.ediff1d(np.append(times, 2*times[-1]-times[-2])).reshape(-1, 1, 1)
    
    nc_file_h.add_variable_lon_lat_time(f"{name}_volatilization_rate",np.transpose(flow_rate_sqm, (1,2,0)))
    nc_file_h.add_info_var(f"{name}_volatilization_rate", f"{name}_volatilization_rate", f"{name} volatilization rate", f"{unit}/m²s")

def create_stackplot(times, vars, names, title, ylabel, path):
    plt.figure(figsize=(10, 6))


    colors = ['#1f77b4', '#1fb4af', '#2ca02c', '#d62728', '#9467bd']
    plt.figure(figsize=(10, 6))
    plt.stackplot(times, vars, labels=names, colors=colors)
    plt.xlabel('Time [s]')
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend(loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(path)

def run_postproc(path_metadata, path_particle, path_layers, path_release, id_gas_interest, x_res, y_res, alias_factor, mw, path_out, out_zip, path_timeseries, lon, lat, epoch_time):
    
    z_layers = load_depth_layers(path_metadata)
    dz = np.ediff1d(z_layers)
    layer_matrix_x, layer_matrix_y, layer_matrix_quantity, layer_release, time_id_layers = load_layers(path_layers)
    particles, name_gas_part, time_id_particles = load_particles(path_particle)
    releases = load_releases(path_release)

    #check time
    for time_l in time_id_layers:
        found = False
        for time_p in time_id_particles:
            if time_l == time_p:
                found = True

        if not found:
            raise Exception(f"No particle time for layer time {time_l}")


    nbr_layers = len(z_layers)-1
    nbr_tmstp = len(time_id_particles)

    quantity_dis = np.zeros((nbr_tmstp, nbr_layers))
    conc_dis = np.zeros((nbr_tmstp, nbr_layers))

    for i in range(nbr_tmstp):
        ids_layer = np.where(time_id_layers == time_id_particles[i])[0]
        if len(ids_layer) == 0:
            continue
        id_layer = ids_layer[0]
        for j in range(nbr_layers):
            if j >= layer_matrix_quantity.shape[0]:
                continue
            quantity_dis[i,j] = layer_matrix_quantity[j, id_layer, id_gas_interest]
            volume = rectangle_area(layer_matrix_x[j, id_layer], layer_matrix_y[j, id_layer])*dz[j]
            if volume > 0:
                conc_dis[i,j] = quantity_dis[i,j] / volume * 1000 #ppm

    quantity_dis_kg = quantity_dis * (mw[id_gas_interest]) #to kg
    #concentration

   
    save_quantity_layers(time_id_particles, z_layers, quantity_dis_kg,'Dissolved mass as a function of depth','kg',f'{path_timeseries}quantity_profile.png')
    save_quantity_layers(time_id_particles, z_layers, conc_dis,'Concentration as a function of depth','ppm',f'{path_timeseries}concentration_profile.png')


    quantity_tr_layer, quantity_tr_bubble, xs, ys = compute_vola(layer_matrix_x, layer_matrix_y, layer_release, time_id_layers, releases, x_res, y_res, alias_factor, time_id_particles, id_gas_interest)
    quantity_tr_layer = quantity_tr_layer * (mw[id_gas_interest]) #to kg
    quantity_tr_bubble = quantity_tr_bubble * (mw[id_gas_interest]) #to kg
    times = time_id_particles + epoch_time
    paths_vola = [path_out+"_vola_layer.nc",path_out+"_bubble_layer.nc", path_out+"_vola_total.nc"]
    save_vola(quantity_tr_layer, xs, ys, lon, lat, "layer", "kg", times, paths_vola[0])
    save_vola(quantity_tr_bubble, xs, ys, lon, lat, "bubble", "kg", times, paths_vola[1])
    save_vola(quantity_tr_layer+quantity_tr_bubble, xs, ys, lon, lat, "total", "kg", times, paths_vola[2])

    nc_h.combine_files(paths_vola, path_out+"_vola.nc","volatilization")
    for f in paths_vola:
        if os.path.isfile(f):
            os.remove(f)
    # plt.figure()
    # plt.plot(np.sum(quantity_tr_layer+quantity_tr_bubble,(1,2)))
    # plt.show()
    # partition ratio
    rel_b = q_released_bubble(time_id_particles, releases, id_gas_interest)
    rel_l = q_released_layer(time_id_particles, time_id_layers, layer_release, id_gas_interest)
    dis = q_dis_layer(time_id_particles, time_id_layers, layer_matrix_quantity, id_gas_interest)
    fresh = q_fresh(time_id_particles, particles, id_gas_interest)
    tot = rel_b + rel_l + dis + fresh
    create_stackplot(time_id_particles, [rel_b*mw[id_gas_interest],rel_l*mw[id_gas_interest],dis*mw[id_gas_interest],fresh*mw[id_gas_interest]], ["Atmosphere from bubble", "Atmosphere from volatilization", "Dissolved", "In bubble"], "Mass repartition", "kg",f'{path_timeseries}mass_balance.png')
    create_stackplot(time_id_particles, [rel_b/tot*100,rel_l/tot*100,dis/tot*100,fresh/tot*100], ["Atmosphere from bubble", "Atmosphere from volatilization", "Dissolved", "In bubble"], "Repartition", "%",f'{path_timeseries}repartition.png')

    path_to_files = [path_out+"_vola.nc",f'{path_timeseries}quantity_profile.png',f'{path_timeseries}concentration_profile.png',f'{path_timeseries}mass_balance.png',f'{path_timeseries}repartition.png']
    with zipfile.ZipFile(out_zip, 'w') as zipf:
        for file_path in path_to_files:
            # Add file to zip, preserving only the file name (not full path)
            zipf.write(file_path, arcname=os.path.basename(file_path))
