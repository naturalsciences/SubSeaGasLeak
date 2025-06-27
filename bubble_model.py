from enum import Enum
import math
import numpy as np
import os
import netcdf_cf_helper as nc_h

import json

import time

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import datetime
import csv

from shapely.geometry import Polygon, box

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


def to_ose_time(start_time, s_since_start):
    time = start_time + datetime.timedelta(seconds = s_since_start)
    return time.strftime("%Y-%m-%dT%H:%M:%S")

def molar_volume_Le_Bas(MW, organic = True):
    """
    Return the molar volume [m³/mol] from the molecular weight
    source : (HNS-MS)
    Parameters
    ----------
    MW : Molecular weight [kg/mol]
    organic : Chemical type, True if organic, the default is True.
    """

    if organic:
        return 4.9807  * (MW*1000)**0.6963 / (100**3)
    else:
        return 2.8047 * (MW*1000)**0.651 / (100**3)


def diffusion_coefficient(mol_vol, wat_viscosity = 10**-3):
    """
    Return the diffusion coefficent [m²/s]
    source : (HNS-MS)
    Parameters
    ----------
    mol_vol : Molar volume of the component [mol/m³]
    wat_viscosity : Dynamic viscosity of water [Pa s]
    """
    return (13.26*10**-5)/((wat_viscosity*1000)**1.14 * (mol_vol*100**3)**0.589)/(100*100)


def k_dis(Dc, De, wb, clean = True):
    """
    Return the mass transfert coefficient dissolution of a bubble [m/s]
    source : (HNS-MS)
    Parameters
    ----------
    Dc : diffusion coefficent[m²/s]
    De : equivalent diameter  [m]
    wb : terminal velocity [m/s]
    clean: if clean bubble: True, else: False
    """
    wb = abs(wb)
    if clean:
        n = 0.5
    else:
        n = 2/3

    if De < 0.5:
        K = 1.13*(Dc**n)*math.sqrt(wb/(0.45+0.2*De))
    elif De < 1.3:
        K = 6.5*Dc**n
    else:
        K = 6.94 * (De**-1/4) * (Dc**n)
    return K / 100#from cm to m



def area_sphere(V):
    """
    return the area of a sphere of a volume V
    """
    return 4*math.pi*((3*V/(4*math.pi))**(1/3))**2

def radius_sphere(V):
    """
    return the raduis of a sphere of a volume V
    """
    return (3*V/(4*math.pi))**(1/3)

def V_gaz_parfait(nbr_mol, T, P = 101325, R = 8.314):
    """
    return the volume (m³) from the number of mol
    Parameters
    ----------
    nbr_mol : number of mol [mol]
    T: Temperature[K]
    P: pressure [Pa]
    R: perfect gas constant [m³ Pa/(K mol)]
    """
    return nbr_mol * R * T / P

def radius_to_mol(radius, T, P = 101325, R = 8.314):
    """
    return the number of mole consiring the radius of a sphere
    Parameters
    ----------
    nbr_mol : number of mol [mol]
    T: Temperature[K]
    P: pressure [Pa]
    R: perfect gas constant [m³ Pa/(K mol)]
    radius: radius [m]
    """
    return P * (4/3 * math.pi * radius**3) / (R * T)


def flux(K, A, S, C):
    """
    Compute the flux [mol/s] of dissolution or volatilization
    Parameters
    ----------
    K: diffusion coefficient[m/s]
    A: area [m²]
    S: solubility [mol/m³]
    C: concentration in the water [mol/m³]
    """
    return K * A *  (S - C)


def intf_tens(rho_hydroc, rho_w, T, critical_T):
    """
    compute the interfacial tension of hydrocarbon [N/m]
    source: https://www.sciencedirect.com/topics/engineering/water-gas-interfacial-tension
    Parameters
    ----------
    rho_hydroc: density of the hydrocarbon[kg/m³]
    rho_w: density of water [kg/m³]
    T: Temperature
    critical_T: critical temperature of the hydrocarbon
    """
    return 0.111*(rho_w/1000-rho_hydroc/1000)**1.024*(T/critical_T)**-1.25

def nondimensional_henry(henry, temperature, R = 8.2e-5):
    """
    Return the nondimensional henry constant []
    source : (Lyman et al., 1990)

    params
    ------
    henry: henry constant [atm m³/mol]
    temperature: temperature [K]
    R: gaz constant 8.2e-5 [atm m³/mol]
    """
    return henry / (temperature * R)

def liquid_phase_exchange_coef(molar_mass, wind_speed):
    """
    Return the liquid phase exchange coefficent [cm/hr]
    source : (Lyman et al., 1990)

    params
    ------
    wind_speed: wind speed [m/s]
    molar_mass: molar mass [kg/mol]
    """
    molar_mass = molar_mass * 1000
    if(molar_mass < 65):
        return 20 * math.sqrt(44/molar_mass)
    elif wind_speed < 3:
        return 2.5
    elif wind_speed < 6:
        return 10
    elif wind_speed < 10:
        return 23
    else:
        return 50


def gas_phase_exchange_coef(molar_mass, wind_speed, current_speed):
    """
    Return the gas phase exchange coefficent [cm/hr]
    source : (Lyman et al., 1990)

    params
    ------
    wind_speed: wind speed [m/s]
    molar_mass: molar mass [kg/mol]
    current_speed: speed of the current [m/s]
    """
    molar_mass = molar_mass * 1000
    if(molar_mass < 65):
        return 3000 * math.sqrt(18/molar_mass)
    else:
        return 1137.5*(wind_speed+current_speed) * math.sqrt(18/molar_mass)


def mass_transfer_coefficient(nd_henry, gas_coef, liquid_coef):
    """
    Return the mass transfer coefficient [m/s]
    source : (Lyman et al., 1990)

    params
    ------
    nd_henry: wind speed []
    gas_coef: gas phase exchange coefficient [cm/hr]
    liquid_coef: liquid phase exchange coefficient [cm/hr]
    """
    a = 1 / (100 * 3600) #from cm/hr to m/s
    return a * (nd_henry * gas_coef * liquid_coef) / (nd_henry * gas_coef + liquid_coef)


def mass_flux_lyman(K, C, H, molar_mass, P = 0, atm_pressure = 101325):
    """
    Return the mass flux [kg/m² s]
    source : (Lyman et al., 1990)

    params
    ------
    K: mass transfer coefficent [m/s]
    C: concentration [kg/m³]
    H: henry constant [atm m³/mol]
    molar_mass: molar mass [kg/m³]
    P: vapor pressure (default is neglected at 0)[Pa]
    atm_pressure: atmospheric pressure, defaut is 101325 [Pa]
    """
    psh = P/(H * atm_pressure / molar_mass)
    return K * (C - psh)


def clift_terminal_velocity(rho_hns, T, critical_T, De, rho_w = 1000, g = 9.81, mu_w=0.001, satif=0.01, max_it=100):
    """
    Compute the terminal velocity
    source: Clift
    Parameters
    ----------
    rho_hns: hns density [kg/m³]
    T: temperature [K]
    critical_T: critical temperature of the HNS [K]
    De: bubble diameter [m]
    rho_w: water density [kg/m³]
    g: earth gravity [m/s²]
    mu_w: water viscosity [Pa s]
    satif: maximum variation between two iteration under wich iterative process stops []
    max_it: maximum number of iteration
    """
    intf_t = intf_tens(rho_hns, rho_w, T, critical_T)


    d_density = abs(rho_w - rho_hns)
    wd = 1
    sign = 1
    if rho_w < rho_hns:
        sign = -1
    if De < 0.001: #less than 1mm
        wd_old = 0
        while abs((wd - wd_old)/wd) > satif and max_it > 0:
            wd_old = wd
            Re = wd * rho_w * De / mu_w
            max_it -= 1
            if Re < 1: #stockes
                wd = g * De**2 * d_density / (18 * mu_w)
            elif Re > 750 and Re < 350000:
                wd = 1.73 * math.sqrt(g * De * d_density / rho_w)
            else: #CLIFT 5.2
                lg = math.log10(Re)
                if Re <= 20:
                    Cd = 24 / Re * (1 + 0.1315 * Re**(0.82 - 0.05 * lg))
                elif Re <= 260:
                    Cd = 24 / Re * (1 + 0.1935 * Re ** 0.6305)
                elif Re <= 1500:
                    Cd = 10 ** (1.6435 - 1.1242 * lg + 0.1558 * lg * lg)
                elif Re <= 400000:
                    Cd = 29.78 - 5.3 * lg
                elif Re <= 1000000:
                    Cd = 0.1 * lg - 0.49
                else:
                    Cd = 0.19 - 800000 / Re
                wd = math.sqrt((4 * d_density * g * De) / (3 * rho_w * Cd))

    else:
        Eo = abs(g * (rho_w - rho_hns) * De * De / intf_t)  #Eotvos
        Mo = abs(g * (rho_w - rho_hns) * mu_w / (rho_w * rho_w * intf_t**3))#Morton
        do_cap = True #If it will follow a spherical cap regime
        if De < 0.015 and Eo < 40 and Mo< 0.001:
            H = 4/3 * Eo * Mo**-0.149*(mu_w/0.009)**-0.14
            if H > 2:
                if H > 59.3:
                    J = 3.42 * H**0.441

                else:
                    J = 0.94*H**0.757
                wd = mu_w/(De*rho_hns)*Mo**-0.149*(J-0.857)#ellipsoidal
                do_cap = False
        if do_cap:
            wd = 0.711 * math.sqrt(g * d_density*De/rho_w)

    return wd * sign

def P_from_depth(depth, p_atm = 101325, rho_w = 1000, g=9.81):
    """
    Return the pressure (Pa) at the depth (m), the P_atm is the pressure of the atmosphere
    atmosphere (Pa) and rho_w (kg/m³) is the ambient water density
    """
    return depth*g*rho_w + p_atm


def velocity_hole(k, rho_hns, rho_water, hole_height, g=9.81):
    """
    Return the velocity (m/s) of hns exiting from a sunken vessel
    Parameters
    ----------
    k : parameters depending on simple or double breach (0.246 vs 0.703)
    rho_hns : density of hns [kg/m³]
    rho_water : density of ambient water [kg/m³]
    hole_height : height of the hole to the bottom of the tank [m]
    g : earth acceleration [m²/s]
    """
    return k * math.sqrt(2*g*hole_height*(rho_water-rho_hns)/rho_water)

def velocity_pipeline(mw_hns, pressure_pipe, pressure_w, temperature, R=8.314):
    """
    Estimation of the velocity (m/s) of the gaz exiting the pipeline
    Parameters
    ----------
    mw_hns : molar weigth of the hns [kg/mol]
    pressure_pipe : pressure in the pipeline [Pa]
    temperature : temperature of the water [K]
    pressure_w : pressure of the water at the depth [Pa]
    R : perfect gas constant [j/mol K]
    """
    return math.sqrt(2*R*temperature/(pressure_pipe*mw_hns)*(pressure_pipe-pressure_w))

def density_from_pressure(mw_hns, pressure_w, temperature, R=8.314):
    """
    Compute the density of a gas from its depth (using perfect gas law)
    Parameters
    ----------
    mw_hns : molar weigth of the hns [kg/mol]
    temperature : temperature of the water [K]
    pressure_w : pressure of the water at the depth [Pa]
    R : perfect gas constant [j/mol K]
    """
    return pressure_w * mw_hns /(R*temperature)

def radius_from_volume(volume, height):
    """
    Return the radius of a cylinder from its volume and heigth
    """
    return math.sqrt(volume/(height*math.pi))



def data_for_cylinder_along_z(center_x,center_y, center_z,radius,height_z):
    z = np.linspace(center_z-height_z/2, center_z+height_z/2, 2)
    theta = np.linspace(0, 2*np.pi, 20)
    theta_grid, z_grid=np.meshgrid(theta, z)
    x_grid = radius*np.cos(theta_grid) + center_x
    y_grid = radius*np.sin(theta_grid) + center_y
    return x_grid,y_grid,z_grid

def clamp(value, min_value = 0, max_value = 1):
    return min(max(value, min_value), max_value)


release_type = Enum('release_type',[
                                    '1_HOLE',
                                    '2_HOLE',
                                    'PIPELINE'
                                ])


def clean_folder(folder_path):
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)

        if os.path.isfile(file_path) or os.path.islink(file_path):
            # Remove the file
            os.unlink(file_path)

def rotate_point(point, angle):
    """
    Rotate a point around the origin
    """
    x, y = point
    x_new = x * math.cos(angle) - y * math.sin(angle)
    y_new = x * math.sin(angle) + y * math.cos(angle)
    return x_new, y_new

def estimate_timestep(dt_per_part, init_radius, exit_velocity):
    """
    Estimate the timestep based on 
    - to have a element not higher than large
    Parameters
    ---------
    dt_per_part : number of time steps per particle release
    init_radius : initial radius of the particle (breach radius)
    exit_velocity : velocity of the particle at the exit
    """
    max_height = init_radius /2
    dt_hv = max_height / exit_velocity / dt_per_part
    return dt_hv

def find_free(active_part):
    """
    Find the index of the first free particle. If none are return -1
    """
    for i in range(len(active_part)):
        if active_part[i] == 0:
            return i
    return -1


def run_model(
            TOTAL_DEPTH, WIND_SPEED, RHO_W, T_W, VISCO_W, VEL_W, CLEAN_W, P_ATM, AREA_BREACH, MASSIC_FLOW_RATE,
                OSE_LON, OSE_LAT, OSE_TIME, EVENT_DURATION, N_LAYERS, TMSTP_PER_PART, MAX_SHAPE_CHANGE, CARACTERISTIC_HEIGHT,
                KZ, KXY, NAME_HNS, AIR_CONC_HNS, INITIAL_FRACTION, MW_HNS, CRIT_T_HNS, H_HNS, MOLAR_VOLUME_HNS, DC_HNS, MU, SIGMA,
                plot_gif, out_files, init_parts = 10, GUI = None
           ):
    """
    Run the model
    
    parameters
    ---------
    TOTAL_DEPTH : depth of the release [m]
    WIND_SPEED : wind speed [m/s]
    RHO_W : density of ambient water [kg/m³]
    T_W : temperature of ambient water [K]
    VISCO_W : viscosity of ambient water [Pa s]
    VEL_W : velocity of ambient water (u, v, w) [m/s]
    CLEAN_W : True if the water is "clean" (no solutes)
    P_ATM : atmospheric pressure [Pa]
    AREA_BREACH : area of the breach [m²]
    MASSIC_FLOW_RATE : massic flow rate of the leak[kg/s]
    OSE_LON : longitude of the release [°]
    OSE_LAT : latitude of the release [°]
    OSE_TIME : time of the release (datetime)
    EVENT_DURATION : maximum duration of the simulation [s]
    N_LAYERS : number of layers
    TMSTP_PER_PART : number of timesteps per particle released
    MAX_SHAPE_CHANGE : maximum shape change of the particle at the surface before being considered in the atmosphere
    CARACTERISTIC_HEIGHT : heigth after which the regime change fully from jet to advection [m]
    KZ : vertical diffusivity [m2/s]
    KXY : horizontal diffusivity [m2/s]
    NAME_HNS : name of the hns (list)
    AIR_CONC_HNS : air concentration of the hns (list)
    INITIAL_FRACTION : initial fraction of the hns in the water (list)
    MW_HNS : molar weight of the hns (list) [kg/mol]
    CRIT_T_HNS : critical temperature of the hns (list) [K]
    H_HNS : henry constant of the hns (list) [mol/m³ Pa]
    MOLAR_VOLUME_HNS : molar volume of the hns (list) [m3/mol]
    DC_HNS : diffusion coefficient of the hns (list) [m2/s]
    MU : mean of the lognormal distribution of bubble size [Pa s]
    SIGMA : standard deviation of the lognormal distribution of bubble size [m]
    plot_gif : True if the gif is wanted
    out_files : path of the outputed files
    init_parts : initial number of particles in the array
    GUI : GUI if a gui is present, otherwise None
    """

    path_json_part = out_files + "_particles.json"
    path_csv_layers = out_files + "_layers.csv"
    path_metadata = out_files + "_metadata.json"
    path_release_event = out_files + "_part_release.json"

    if os.path.exists(path_json_part):
        os.remove(path_json_part)
    if os.path.exists(path_csv_layers):
        os.remove(path_csv_layers)
    if os.path.exists(path_metadata):
        os.remove(path_metadata)
    if os.path.exists(path_release_event):
        os.remove(path_release_event)



    NORM_VEL_W = math.sqrt(VEL_W[0]**2+VEL_W[1]**2+VEL_W[2]**2)
    ANGLE_W = math.atan(VEL_W[1]/VEL_W[0])

    EARTH_RADIUS = 6371000 #[m]
    m_to_deg_lat = 360/(EARTH_RADIUS*math.pi*2)
    lat_radius = math.cos(OSE_LAT*math.pi /180) * EARTH_RADIUS
    m_to_deg_lon = 360 / (lat_radius * math.pi * 2)

    # NBR_PART = int(EVENT_DURATION / DT / TMSTP_PER_PART)
    # NBR_TIMESTEP = int(EVENT_DURATION / DT)

    INITIAL_FRACTION = np.array(INITIAL_FRACTION)
    MW_HNS = np.array(MW_HNS)
    mw = sum(INITIAL_FRACTION * MW_HNS) / sum(INITIAL_FRACTION)

    #compute the water concentration for gas dissolved in high amount already
    water_conc_hns = [0]*len(NAME_HNS)#mol/m³
    for i in range(len(water_conc_hns)): #for all the chemical in equilibrium with the atmosphere
        water_conc_hns[i] = (AIR_CONC_HNS[i]/sum(AIR_CONC_HNS)) * H_HNS[i] * P_ATM

    array_name_csv = []
    for i in range(len(NAME_HNS)):
        if water_conc_hns[i] ==0:
            array_name_csv.append(NAME_HNS[i] + " [mol]")
    for i in range(len(NAME_HNS)):
        if water_conc_hns[i] ==0:
            array_name_csv.append(NAME_HNS[i] + " volatilized [mol]")
    # create the csv
    with open(path_csv_layers, 'w') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(['Current time [s]',
                             "id",
                             "x layers 0 [m]",
                             "x layers 1 [m]",
                             "x layers 2 [m]",
                             "x layers 3 [m]",
                             "y layers 0 [m]",
                             "y layers 1 [m]",
                             "y layers 2 [m]",
                             "y layers 3 [m]"] + array_name_csv)
    with open(path_json_part, "a") as f:
        f.write('{\n"simulation": [\n')
    with open(path_release_event, "a") as f:
        f.write('{\n"simulation": [\n')

    #cleaning
    clean_folder("tmp")

    
    init_density = density_from_pressure(mw, P_from_depth(TOTAL_DEPTH, p_atm = P_ATM, rho_w = RHO_W), T_W)

    flow_rate = MASSIC_FLOW_RATE / init_density

    exit_velocity = flow_rate / AREA_BREACH

    radius_init = math.sqrt(AREA_BREACH / math.pi)

    DT = estimate_timestep(TMSTP_PER_PART, radius_init, exit_velocity)
    
    #initial properties
    masse_gas_init_hns = MASSIC_FLOW_RATE * DT * TMSTP_PER_PART #initial mass of gas
    mol_gas_init_hns = masse_gas_init_hns / mw#init
    height_init = exit_velocity * DT * TMSTP_PER_PART #initial heigth of the particle
    volume_gas_init_hns = masse_gas_init_hns/mw*8.314 * T_W / P_from_depth(TOTAL_DEPTH, p_atm = P_ATM, rho_w = RHO_W)#initial volume of the particle
    ratio_init = radius_init / height_init #initial ratio

    #particle intialization

    initial_diameter = np.random.lognormal(MU, SIGMA, init_parts)
    initial_radius_3 = (initial_diameter/2)**3 # initial radius ^3
    number_bubble_part = volume_gas_init_hns / (4/3*math.pi*initial_radius_3)
    mol_bubble_part = np.zeros((init_parts,len(NAME_HNS)))
    for i in range(len(NAME_HNS)):
        mol_bubble_part[:,i] =masse_gas_init_hns /(mw*number_bubble_part)*INITIAL_FRACTION[i]

    part_init_d = -TOTAL_DEPTH+ height_init/2

    z_part = np.full((init_parts,), part_init_d)
    x_part = np.zeros_like(z_part)
    y_part = np.zeros_like(z_part)
    active_part = np.zeros_like(z_part)
    phase_part = np.zeros_like(z_part) #phase advection (1) or jet (0)
    part_id = np.full((init_parts,), -1)

    phase_change_t = np.full((init_parts,), 1e20)

    radius_part = np.full((init_parts,),radius_init)
    height_part = np.full((init_parts,),height_init)
    fully_inside_water_part = np.full((init_parts,),True) # have not touched surface yet

    regime_change_depth = part_init_d + CARACTERISTIC_HEIGHT

    # boxes for layers
    x_layers = np.zeros((N_LAYERS,4))
    y_layers = np.zeros((N_LAYERS,4))


    quantity_layers = np.zeros((N_LAYERS, len(NAME_HNS)))

    z_layers = np.linspace(-TOTAL_DEPTH, 0, N_LAYERS+1)

    mask_layers = np.full((N_LAYERS,), False)

    dz = np.diff(z_layers)

    l_part_in_layer = np.zeros((init_parts, N_LAYERS))

    activated_layers = np.zeros((N_LAYERS,))

    volume_init = radius_init * radius_init * math.pi * height_init

    volume_part = np.full((init_parts,), volume_init)

    # mol_bubble_part_all = np.zeros((NBR_TIMESTEP+1, NBR_PART,len(NAME_HNS)))
    # z_part_all = np.zeros((NBR_TIMESTEP+1, NBR_PART))
    # x_part_all = np.zeros_like(z_part_all)
    # y_part_all = np.zeros_like(z_part_all)
    # radius_part_all = np.zeros_like(z_part_all)
    # height_part_all = np.zeros_like(z_part_all)
    # active_part_all = np.zeros_like(z_part_all)


    # x_layers_all = np.zeros((NBR_TIMESTEP+1,N_LAYERS,4))
    # y_layers_all = np.zeros((NBR_TIMESTEP+1,N_LAYERS,4))

    vola_mol = np.zeros_like(H_HNS)

    # quantity_layers_all = np.zeros((NBR_TIMESTEP+1,N_LAYERS, len(NAME_HNS)))
    # volume_layer_all = np.zeros((NBR_TIMESTEP+1,N_LAYERS))

    release_events = [] # event of gas released in the atmosphere by bubble
    #volatilizations_mol = np.zeros((NBR_TIMESTEP+1,len(NAME_HNS))) # number of mol volatilized
    vola_since_last_out = np.zeros(len(NAME_HNS))

    rac_2 = math.sqrt(2)

    clift_terminal_velocity_vecto = np.vectorize(clift_terminal_velocity)
    k_dis_vecto = np.vectorize(k_dis)

    animate = True

    nbr_img = 1
    image_every = 0.25
    last_image = -1000
    last_output = -1000

    released_part = 0

    sqr2 = math.sqrt(2)

    bench_select_var = 0
    bench_elemnt_lagr = 0
    bench_jet_srfc = 0
    bench_set_layer = 0
    bench_transfer_layer = 0
    bench_diss = 0
    bench_vola = 0
    bench_img = 0
    bench_json = 0

    id_to_show = 0
    mol_dissolved = [0]
    mol_srfc = [0]
    mol_vola = [0]
    mol_fresh = [0]
    times = [0]

    number_past_tmstp = 10
    last_tmstps = np.zeros(number_past_tmstp)

    #only of the id_to_show
    # mol_dissolved = np.zeros(NBR_TIMESTEP+1)
    # mol_srfc = np.zeros(NBR_TIMESTEP+1)
    # mol_vola = np.zeros(NBR_TIMESTEP+1)

    current_time = 0
    simulation_duration = EVENT_DURATION
    dummy1_ar = np.array([1])
    dummy0_ar = np.array([0])
    dummyTrue_ar = np.array([True])
    dummy_large_arr = np.array([1e20])
    part_init_d_ar = np.array([part_init_d])
    radius_init_ar = np.array([radius_init])
    height_init_ar = np.array([height_init])
    dummy_layer_ar = np.zeros((N_LAYERS,))
    dummy_hns_ar = np.zeros((len(NAME_HNS),))
    dummy_volume_ar = np.array([volume_init])
    i = 0

    first_output = True
    first_release = True

    #write the metadata:
    metadata = {
        "depth_layers" :z_layers.tolist()
    }
    with open(path_metadata, "w") as f:
        json.dump(metadata, f, indent=4)



    while current_time <= simulation_duration: 
        time_st_tsmt = time.perf_counter()
        if GUI is not None:
            average_tmstp = np.mean(last_tmstps)
            s_pers_s = DT/average_tmstp
            if simulation_duration < 1e10:
                GUI.model_status.setText(f"Currently at {float(current_time):0.2f} s/{simulation_duration:0.2f} s of simulation. Currently {len(z_part)} elements. Simulation speed {s_pers_s:0.2f} s/s")  
            else:
                GUI.model_status.setText(f"Currently at {float(current_time):0.2f} s, highest element at {max(z_part):0.2f} m. Currently {len(z_part)} elements. Simulation speed {s_pers_s:0.2f} s/s")             
        times.append(current_time)
        # if i > 0:
        #     mol_dissolved[i+1] = mol_dissolved[i]
        #     mol_srfc[i+1] = mol_srfc[i]
        #     mol_vola[i+1] = mol_vola[i]


        start_bench = time.perf_counter()
        # mol_bubble_part_all[i] = mol_bubble_part
        # z_part_all[i] = z_part
        # x_part_all[i] = x_part
        # y_part_all[i] = y_part
        # radius_part_all[i] = radius_part
        # height_part_all[i] = height_part
        # x_layers_all[i] = x_layers
        # y_layers_all[i] = y_layers
        # active_part_all[i] = active_part
        # quantity_layers_all[i] = quantity_layers
        
        #double check because of the floored total number of particle
        if i // TMSTP_PER_PART > released_part:
            id_free = find_free(active_part)
            if id_free >= 0:
                part_id[id_free] = released_part
                active_part[id_free] = 1
                z_part[id_free] = part_init_d_ar
                x_part[id_free] = 0
                y_part[id_free] = 0
                phase_part[id_free] = 0
                phase_change_t[id_free] = 1e20
                radius_part[id_free] = radius_init
                height_part[id_free] = height_init
                fully_inside_water_part[id_free] = True
                l_part_in_layer[id_free] = dummy_layer_ar[:]
                volume_part[id_free] = volume_init
                initial_diameter_l = np.random.lognormal(MU, SIGMA, 1) 
                initial_radius_3_l = (initial_diameter_l/2)**3
                number_bubble_part[id_free] = volume_gas_init_hns / (4/3*math.pi*initial_radius_3_l)
                for j in range(len(NAME_HNS)):
                    mol_bubble_part[id_free, j] =masse_gas_init_hns /(mw*number_bubble_part[id_free])*INITIAL_FRACTION[j]
            else:
                active_part = np.concatenate((active_part,dummy1_ar))
                z_part = np.concatenate((z_part, part_init_d_ar))
                x_part = np.concatenate((x_part, dummy0_ar))
                y_part = np.concatenate((y_part, dummy0_ar))
                phase_part = np.concatenate((phase_part, dummy0_ar))

                phase_change_t = np.concatenate((phase_change_t, dummy_large_arr))

                radius_part = np.concatenate((radius_part, radius_init_ar))
                height_part = np.concatenate((height_part, height_init_ar))
                fully_inside_water_part = np.concatenate((fully_inside_water_part, dummyTrue_ar))
                l_part_in_layer = np.vstack((l_part_in_layer, dummy_layer_ar))
                volume_part = np.concatenate((volume_part, dummy_volume_ar))

                initial_diameter_l = np.random.lognormal(MU, SIGMA, 1)
                initial_radius_3_l = (initial_diameter_l/2)**3 # initial radius ^3
                number_bubble_part = np.concatenate((number_bubble_part,volume_gas_init_hns / (4/3*math.pi*initial_radius_3_l)))
                for j in range(len(NAME_HNS)):
                    dummy_hns_ar[j] =masse_gas_init_hns /(mw*number_bubble_part[-1])*INITIAL_FRACTION[j]
                mol_bubble_part = np.vstack((mol_bubble_part, dummy_hns_ar))
                part_id = np.concatenate((part_id, np.array([released_part])))
        
            released_part += 1
        mask_active = (active_part ==1)

        select_var_bench = time.perf_counter()
        bench_select_var += select_var_bench - start_bench
        tot_srfc = 0
        if np.sum(active_part) > 0:
            start_bench = time.perf_counter()

            phase_part[z_part > regime_change_depth] = 1

            phase_change_t[np.logical_and(phase_part == 1, phase_change_t == -1)] = i * DT

            #rising

            pressure_part = P_from_depth(-z_part[mask_active], p_atm = P_ATM, rho_w = RHO_W)
            
            # weigthing of the mw per the number of mol
            average_mw = np.sum(mol_bubble_part[mask_active] * MW_HNS, axis = 1) / np.sum(mol_bubble_part[mask_active], axis = 1)

            rho_part = density_from_pressure(average_mw, pressure_part, T_W, R=8.314)
            volume_bubble_part = V_gaz_parfait(np.sum(mol_bubble_part[mask_active], axis = 1), T_W, P = pressure_part)
            De_bubble_part = 2*(3*volume_bubble_part/(4*math.pi))**(1/3)
            ws_bubble_part = clift_terminal_velocity_vecto(rho_part, T_W, CRIT_T_HNS[0], De_bubble_part, rho_w = RHO_W, mu_w = VISCO_W)

            x_part[mask_active] += VEL_W[0] * DT
            y_part[mask_active] += VEL_W[1] * DT

            advection_phase_dist = (VEL_W[2] + ws_bubble_part) * DT

            active_jet = np.logical_and(mask_active, phase_part == 0)
            active_jet_active = phase_part[mask_active] == 0

            active_adv = np.logical_and(mask_active, phase_part == 1)
            active_adv_active = phase_part[mask_active] == 1

            interp_c = (z_part[active_jet]-regime_change_depth)/(part_init_d-regime_change_depth)

            z_part[active_jet] += interp_c * (exit_velocity * DT) + (1-interp_c) * advection_phase_dist[active_jet_active]

            z_part[active_adv] += advection_phase_dist[active_adv_active]

            bench_elemnt_lagr += time.perf_counter() - start_bench
            start_bench = time.perf_counter()

            #size increase jet

            mask_fully_underwater_jet = np.logical_and(np.logical_and(fully_inside_water_part, mask_active), phase_part == 0) #on all the particle
            mask_fully_underwater_jet_active = mask_fully_underwater_jet[mask_active] #only on the active particle
            F = 2*ws_bubble_part[mask_fully_underwater_jet_active]/(np.sqrt(9.81*(RHO_W-rho_part[mask_fully_underwater_jet_active])/RHO_W*radius_part[mask_fully_underwater_jet]))
            alpha = rac_2 * (0.057+0.554*np.sin(ANGLE_W)/F*F) / (1+5*NORM_VEL_W / ws_bubble_part[mask_fully_underwater_jet_active])

            Qs = 2*math.pi*radius_part[mask_fully_underwater_jet]*height_part[mask_fully_underwater_jet]*alpha*ws_bubble_part[mask_fully_underwater_jet_active]

            volume_part[mask_fully_underwater_jet] += Qs * DT
            height_part[mask_fully_underwater_jet] = (volume_part[mask_fully_underwater_jet]/(math.pi*ratio_init*ratio_init)) ** (1/3)
            radius_part[mask_fully_underwater_jet] = height_part[mask_fully_underwater_jet]*ratio_init

            #size increase adv, but only if the time since the phase change is long enough
            dt_since_phase = i * DT - phase_change_t
            mask_fully_underwater_adv = np.logical_and(np.logical_and(fully_inside_water_part, mask_active), dt_since_phase > 0) #on all the particle

            #cf wikipedia random walk for k_z
            height_part[mask_fully_underwater_adv] += sqr2 * np.sqrt(KZ/dt_since_phase[mask_fully_underwater_adv]) * DT
            radius_part[mask_fully_underwater_adv] = np.sqrt(radius_part[mask_fully_underwater_adv]*radius_part[mask_fully_underwater_adv]+4*KXY*DT/math.pi)

            #reaching surface
            surfc_reached = np.logical_and(active_part ==1,(height_part/2 + z_part)>0)
            surfc_reached_active = surfc_reached[mask_active]

            hi_tmp = height_part[surfc_reached]/2 - z_part[surfc_reached]
            fully_inside_water_part[surfc_reached] = False
            bubble_before = number_bubble_part[surfc_reached]
            number_bubble_part[surfc_reached] = number_bubble_part[surfc_reached] * (hi_tmp/height_part[surfc_reached])
            if np.sum(surfc_reached) > 0:
                #simulation shoudl stop 4x after the time a particle reached surface
                if simulation_duration > 4 * current_time:
                    simulation_duration = 4 * current_time
                for j, index_part in enumerate(np.where(surfc_reached)[0]):
                    #timestep, particle, mol gone in atmosphere
                    release_events.append((index_part,mol_bubble_part[index_part]*(bubble_before[j]-number_bubble_part[index_part])))
                    tot_srfc += (mol_bubble_part[index_part]*(bubble_before[j]-number_bubble_part[index_part]))[id_to_show]
                    #mol_srfc[i] += mol_bubble_part[index_part,id_to_show]*(bubble_before[j]-number_bubble_part[index_part])
                    # in_atm += mol_bubble_part[index_part,0]*(bubble_before[j]-number_bubble_part[index_part])
            volume_part[surfc_reached] = volume_part[surfc_reached] - volume_bubble_part[surfc_reached_active] * (1 - hi_tmp/height_part[surfc_reached])
            z_part[surfc_reached] = - height_part[surfc_reached]/2
            height_part[surfc_reached] = hi_tmp
            #to avoid zero height, causing dividing by 0
            idx = np.where(surfc_reached)[0]
            mask_no_zero_height = height_part[idx] != 0
            radius_part[idx[mask_no_zero_height]] = np.sqrt(volume_part[idx[mask_no_zero_height]]/(math.pi*height_part[idx[mask_no_zero_height]]))
            
            active_part[height_part < 0] = 0
            active_part[np.logical_and(surfc_reached, radius_part/height_part > ratio_init*MAX_SHAPE_CHANGE)] = 0

            bench_jet_srfc += time.perf_counter() - start_bench
            start_bench = time.perf_counter()
            

            
            #layer and dissolution
            x_layers[mask_layers,:] += VEL_W[0] * DT
            y_layers[mask_layers,:] += VEL_W[1] * DT

            l_part_in_layer[:,:] = 0

            #looking for the heigth of the particle at each level
            for j in np.where(active_part == 1)[0]:
                for k in range(N_LAYERS):
                    l_part_in_layer[j,k] = abs(clamp((z_layers[k+1]-z_part[j]-height_part[j]/2)/dz[k])-clamp((z_layers[k+1]-z_part[j]+height_part[j]/2)/dz[k]))

        activated_layers[np.where(np.sum(l_part_in_layer, axis=0)>0)] = 1
        activated_layers[np.where(np.sum(quantity_layers, axis=0)>0)] = 1

        mask_layers = (activated_layers == 1)

        start_bench = time.perf_counter()
        mol_srfc.append(tot_srfc+mol_srfc[-1])
        #making sure the layers encompass all the particles
        for j in np.where(np.sum(l_part_in_layer, axis=0)>0)[0]:
            for k in np.where(l_part_in_layer[:,j]>0)[0]:
                # the point are rotate, to have the rectangle aligned with the x and y axis, to simplify the calculation
                center_rot = rotate_point((x_part[k],y_part[k]), ANGLE_W)

                x_low, y_low = rotate_point((x_layers[j,0],y_layers[j,0]), ANGLE_W)
                x_high, y_high = rotate_point((x_layers[j,2],y_layers[j,2]), ANGLE_W)

                x_low = min(x_low, center_rot[0]-radius_part[k])
                x_high = max(x_high, center_rot[0]+radius_part[k])
                y_low = min(y_low, center_rot[1]-radius_part[k])
                y_high = max(y_high, center_rot[1]+radius_part[k])


                (x_layers[j,0],y_layers[j,0]) = rotate_point((x_low, y_low), -ANGLE_W)
                (x_layers[j,1],y_layers[j,1]) = rotate_point((x_high, y_low), -ANGLE_W)
                (x_layers[j,2],y_layers[j,2]) = rotate_point((x_high, y_high), -ANGLE_W)
                (x_layers[j,3],y_layers[j,3]) = rotate_point((x_low, y_high), -ANGLE_W)

        #making shure the layers are at least as big as the layers under
        for j in range(N_LAYERS-1):
            x_low, y_low = rotate_point((x_layers[j+1,0],y_layers[j+1,0]), ANGLE_W)
            x_high, y_high = rotate_point((x_layers[j+1,2],y_layers[j+1,2]), ANGLE_W)

            x_low_below, y_low_below = rotate_point((x_layers[j,0],y_layers[j,0]), ANGLE_W)
            x_high_below, y_high_below = rotate_point((x_layers[j,2],y_layers[j,2]), ANGLE_W)

            x_low = min(x_low, x_low_below)
            x_high = max(x_high, x_high_below)
            y_low = min(y_low, y_low_below)
            y_high = max(y_high, y_high_below)

            (x_layers[j+1,0],y_layers[j+1,0]) = rotate_point((x_low, y_low), -ANGLE_W)
            (x_layers[j+1,1],y_layers[j+1,1]) = rotate_point((x_high, y_low), -ANGLE_W)
            (x_layers[j+1,2],y_layers[j+1,2]) = rotate_point((x_high, y_high), -ANGLE_W)
            (x_layers[j+1,3],y_layers[j+1,3]) = rotate_point((x_low, y_high), -ANGLE_W)


        bench_set_layer += time.perf_counter() - start_bench
        start_bench = time.perf_counter()
        # transfer between layers

        width_layers = np.sqrt((x_layers[:,0] - x_layers[:,1]) * (x_layers[:,0] - x_layers[:,1]) + (y_layers[:,0] - y_layers[:,1]) * (y_layers[:,0] - y_layers[:,1]))
        depth_layers = np.sqrt((x_layers[:,1] - x_layers[:,2]) * (x_layers[:,1] - x_layers[:,2]) + (y_layers[:,1] - y_layers[:,2]) * (y_layers[:,1] - y_layers[:,2]))
        area_layers = width_layers * depth_layers
        
        volume_layers = area_layers * dz
        #volume_layer_all[i] = volume_layers
        conc_layers = np.where(volume_layers[:,np.newaxis] > 0, quantity_layers / volume_layers[:,np.newaxis], [0]*len(MW_HNS))

        #difference between each layer and the one above
        diff_layer = np.diff(conc_layers, axis=0)
        flux_diff = KZ * area_layers[:-1,np.newaxis] * diff_layer / dz[:-1,np.newaxis]
        # conc_layers[0,:] += KZ * (diff_layer[0,:]/dz[0]/dz[0]) * DT
        # conc_layers[1:-1,:] += KZ * (diff_layer[1:,:]/dz[1:-1,np.newaxis]-diff_layer[:-1,:]/dz[1:-1,np.newaxis]) /dz[1:-1,np.newaxis] * DT
        # conc_layers[-1,:] += KZ * (-diff_layer[-1,:]/dz[0]/dz[0]) * DT

        # quantity_layers_tmp = conc_layers * volume_layers[:,np.newaxis]

        # #correction to keep mass balance
        # total_q_diss = np.sum(quantity_layers, axis=0)
        # total_q_diss_tmp = np.sum(quantity_layers_tmp, axis=0)
        
        # ratio = total_q_diss/total_q_diss_tmp
        # ratio = np.where(total_q_diss_tmp == 0, 0, ratio)
        # quantity_layers = quantity_layers_tmp * ratio
        quantity_layers[:-1] += flux_diff * DT
        quantity_layers[1:] -= flux_diff * DT
        
        conc_layers = np.where(volume_layers[:,np.newaxis] > 0, quantity_layers / volume_layers[:,np.newaxis], [0]*len(MW_HNS))

        bench_transfer_layer += time.perf_counter() - start_bench
        start_bench = time.perf_counter()
        
        #dissolution
        if np.sum(active_part) > 0:
            r_bubble = radius_sphere(volume_bubble_part)
            area_bubble = area_sphere(volume_bubble_part)

            k_dis_part = []

            for l in range(len(NAME_HNS)):
                k_dis_part.append(k_dis_vecto(DC_HNS[l], r_bubble*2, ws_bubble_part, clean = CLEAN_W))
            k_dis_part = np.transpose(np.array(k_dis_part))

            dis_flux = np.zeros((np.sum(mask_active), len(NAME_HNS)))
            #layer on which the particle is located
            for j in np.where(np.sum(l_part_in_layer, axis=0)>0)[0]:
                # particles
                for index, k in enumerate(np.where(l_part_in_layer[:,j]>0)[0]):
                    for l in range(len(NAME_HNS)):
                        if water_conc_hns[l] > 0:
                            conc_layers[j,l] = water_conc_hns[l]
                        q_trans =  DT * number_bubble_part[k] * l_part_in_layer[k,j] * flux(k_dis_part[index,l],area_bubble[index],pressure_part[index] * H_HNS[l],conc_layers[j,l]) #TODO And the molar Fraction???? AND THE LENGTH SHOULD BE NORMALIZED!!!!!
                        if q_trans < 0:
                            q_trans = -min(quantity_layers[j,l], -q_trans)
                        else:
                            q_trans = min(number_bubble_part[k]*mol_bubble_part[k,l], q_trans) #TODO this will cause the flux to be too high, there are many layers!!!! qtrans must be smaller (/nbr bbl)
                        dis_flux[index,l] += q_trans
                        
                        if water_conc_hns[l] == 0:
                            quantity_layers[j,l] += q_trans

            mol_bubble_part[mask_active] -= dis_flux / number_bubble_part[mask_active, np.newaxis]
            active_part[np.where(np.sum(mol_bubble_part, axis=1) <1e-24)] = 0
            np.where(mol_bubble_part <0, 0, mol_bubble_part)
            
            #mol_dissolved[i+1] += np.sum(dis_flux[:,id_to_show])
            
        bench_diss += time.perf_counter() - start_bench
        start_bench = time.perf_counter()

        mask_active = (active_part ==1)


        mol_vola.append(mol_vola[-1])
        # volatilization
        for l in range(len(NAME_HNS)):
            if water_conc_hns[l] == 0:
                nd_h = nondimensional_henry(101325/H_HNS[l], T_W)
                l_ec = liquid_phase_exchange_coef(MW_HNS[l], WIND_SPEED)
                g_ec = gas_phase_exchange_coef(MW_HNS[l], WIND_SPEED, NORM_VEL_W)
                k_vola = mass_transfer_coefficient(nd_h, g_ec, l_ec)
                
                fl = mass_flux_lyman(k_vola, conc_layers[-1,0], 101325/H_HNS[l], MW_HNS[l])
                vola_mol = area_layers[-1] * fl * DT / MW_HNS[l]
                quantity_layers[-1,l] -= vola_mol
                vola_since_last_out[l] += vola_mol
                if l == id_to_show:
                    mol_vola[-1] += vola_mol
            # volatilizations_mol[i,l] = vola_mol
            
            # if l == id_to_show:
            #     mol_dissolved[i+1] -= vola_mol
            #     mol_vola[i+1] += vola_mol

        bench_vola += time.perf_counter() - start_bench
        start_bench = time.perf_counter()
        if len(release_events) >0:
            data_release ={
                "current_time": current_time,
                "releases":[]
            }
            for release in release_events:
                id_part = release[0]
                data_release["releases"].append({
                    "quantity_released":release[1].tolist(),
                    "x": x_part[id_part],
                    "y": y_part[id_part],
                    "radius_part": radius_part[id_part]
                })
            release_events = []

            with open(path_release_event, "a") as f:
                if first_release:
                    first_release = False
                else:
                    f.write(',\n')
                json.dump(data_release, f, indent=4)

        if current_time - last_output >= image_every:
            last_output = current_time
            
            data_part = {
                "current_time": current_time,
                "z_part": z_part.tolist(),
                "x_part": x_part.tolist(),
                "y_part": y_part.tolist(),
                "active_part": active_part.tolist(),
                "phase_part": phase_part.tolist(),
                "phase_change_t": phase_change_t.tolist(),
                "radius_part": radius_part.tolist(),
                "height_part": height_part.tolist(),
                "fully_inside_water_part": fully_inside_water_part.tolist(),
                "number_bubble_part": number_bubble_part.tolist(),
                "part_id": part_id.tolist()
            }
            for j in range(len(NAME_HNS)):
                data_part["mol_bubble_part_"+NAME_HNS[j]]= mol_bubble_part[:,j].tolist()


            with open(path_json_part, "a") as f:
                if first_output:
                    first_output = False
                else:
                    f.write(',\n')
                json.dump(data_part, f, indent=4)
                

            with open(path_csv_layers, 'a') as csvfile:
                csv_writer = csv.writer(csvfile)
                for j in range(N_LAYERS):
                    if activated_layers[j]:
                        row = [current_time, j]
                        for x in x_layers[j,:]:
                            row.append(x)
                        for y in y_layers[j,:]:
                            row.append(y)
                        for q in range(len(quantity_layers[j,:])):
                            if water_conc_hns[q] ==0: #do not add the quantity assumed constant
                                row.append(quantity_layers[j,q])
                        for k in range(len(vola_since_last_out)):
                            if water_conc_hns[k] ==0:
                                if j == N_LAYERS-1:
                                        row.append(vola_since_last_out[k])
                                        vola_since_last_out[k] = 0
                                else:
                                    row.append(0)
                        csv_writer.writerow(row)    
                    

        bench_json += time.perf_counter() - start_bench

        start_bench = time.perf_counter()
        if current_time - last_image >= image_every and plot_gif:
            
            nbr_img += 1
            # if nbr_img%10 !=0:
            #     continue
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            if animate:
                elev = 30
                azim = 45
                #ax.view_init(elev=elev*((i)/EVENT_DURATION)/3, azim=azim*((i)/EVENT_DURATION)/3)
                ax.view_init(elev=elev*0.25, azim=azim*0.25)
            ax.set_xlim([-5, 5])
            ax.set_ylim([-5, 5])
            ax.set_zlim([-TOTAL_DEPTH, 0])
            ax.set_title(f"{current_time:.2f} s")
            ax.xaxis.set_major_locator(plt.MaxNLocator(3))
            ax.set_xlabel('[m]')
            ax.set_ylabel('[m]')    
            ax.set_zlabel('Depth [m]')
            for j in np.where(active_part ==1)[0]:
                Xc,Yc,Zc = data_for_cylinder_along_z(x_part[j],y_part[j],z_part[j],radius_part[j],height_part[j])
                sum_bubble = np.sum(mol_bubble_part[j,:])

                ax.plot_surface(Xc, Yc, Zc, alpha=clamp(min(0.4, (sum_bubble*number_bubble_part[j])/mol_gas_init_hns)), color=(clamp(mol_bubble_part[j,1]/sum_bubble),clamp(mol_bubble_part[j,0]/sum_bubble),clamp(mol_bubble_part[j,2]/sum_bubble)))
            
            #sc = ax.scatter(part_x[mask_active], part_y[mask_active], z_part[mask_active], marker='o')
            for j in np.where(activated_layers ==1)[0]:
                j = int(j)
                X = [x_layers[j,0], x_layers[j,3], x_layers[j,2], x_layers[j,1], x_layers[j,0], x_layers[j,1], x_layers[j,2], x_layers[j,3]]
                Y = [y_layers[j,0], y_layers[j,3], y_layers[j,2], y_layers[j,1], y_layers[j,0], y_layers[j,1], y_layers[j,2], y_layers[j,3]]
                Z = [z_layers[j], z_layers[j], z_layers[j], z_layers[j], z_layers[j+1], z_layers[j+1], z_layers[j+1], z_layers[j+1]]
                vertices = [[0,1,2,3],[1,5,6,2],[3,2,6,7],[4,0,3,7],[5,4,7,6],[4,5,1,0]]

                tupleList = list(zip(X, Y, Z))

                poly3d = [[tupleList[vertices[ix][iy]] for iy in range(len(vertices[0]))] for ix in range(len(vertices))]
                ax.add_collection3d(Poly3DCollection(poly3d, facecolors='g', linewidths=1, alpha=min(0.1,10**conc_layers[j,0]/1000)))

            plt.savefig(f'tmp/{nbr_img:05d}.png')
            plt.close(fig)
            
            last_image = current_time
        bench_img += time.perf_counter() - start_bench
        if i % 100 == 0 and simulation_duration < 1e10:
            print(f"{int(current_time/simulation_duration*100+1)}%")

        mol_dissolved.append(0)
        for j in range(N_LAYERS):
            mol_dissolved[-1] += quantity_layers[j,id_to_show]

        mol_fresh.append(np.sum(mol_bubble_part[mask_active, id_to_show] * number_bubble_part[mask_active]))
        

        current_time += DT
        i += 1
        if GUI is not None:
            if GUI.in_post:
                simulation_duration = current_time-1
        last_tmstps = np.roll(last_tmstps,1)
        last_tmstps[-1] = time.perf_counter() - time_st_tsmt
        
    print(f"Select vars + arr: {bench_select_var}")
    print(f"Element lagr : {bench_elemnt_lagr}")
    print(f"Jet surface : {bench_jet_srfc}")
    print(f"Set layer : {bench_set_layer}")
    print(f"Transfer layer : {bench_transfer_layer}")
    print(f"Dissolved : {bench_diss}")
    print(f"Vola : {bench_vola}")
    print(f"Json : {bench_json}")
    print(f"Image : {bench_img}")

    with open(path_json_part, "a") as f:
        f.write(']}')
    with open(path_release_event, "a") as f:
        f.write(']}')

    # all = [mol_dissolved, mol_srfc, mol_vola, mol_fresh]
    # for i in range(len(all)):
    #     for j in range(len(all[i])):
    #         all[i][j] = all[i][j] * MW_HNS[id_to_show]

    # plt.stackplot(times, all, labels = ['Dissolved', 'Surface', 'Volatilized', 'Fresh'])
    # plt.plot((times[0], times[-1]), (0,MASSIC_FLOW_RATE * times[-1]))
    # plt.legend(loc='upper left')
    # plt.show()

    return (path_json_part, path_csv_layers, path_metadata, path_release_event)