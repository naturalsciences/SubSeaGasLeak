import os
from datetime import datetime
import logging
import logging.handlers
import subprocess
import subprocess
import traceback
import sys
import argparse
import json
import math
import netCDF4 as nc
import numpy as np
from bubble_model import run_model
from postproc import run_postproc
import random
import string

PATH_GEBCO = "data/GEBCO_2024.nc"


def random_char(y):
   """
   Return y random ascii_letters
   """
   return ''.join(random.choice(string.ascii_letters) for x in range(y))

def get_depth(lon_pt, lat_pt):
    """
    Read in the bathymetry the depth at the given location with a bilinear interpolation

    Parameters
    --------
    lon_pt : longitude of the point
    lat_pt : latitude of the point
    """
    if not os.path.exists(PATH_GEBCO):
        logging.error(f"File {PATH_GEBCO} does not exist. Please only use \"m above seabed\" as location_depth_units")

    with nc.Dataset(PATH_GEBCO, 'r') as ds:
        lon = ds.variables['lon'][:]
        lat = ds.variables['lat'][:]
        elevation = ds.variables['elevation'][:]
        id_lon = 1
        id_lat = 1

        while lon[id_lon] < lon_pt:
            id_lon += 1
        
        while lat[id_lat] < lat_pt:
            id_lat += 1

        coef_lon = (lon_pt - lon[id_lon-1])/(lon[id_lon]-lon[id_lon-1])
        coef_lat = (lat_pt - lat[id_lat-1])/(lat[id_lat]-lat[id_lat-1])

        lon_low = min(0,elevation[id_lat,id_lon-1]) * coef_lat + min(0,elevation[id_lat-1, id_lon-1] * (1-coef_lat))
        lon_high = min(0,elevation[id_lat,id_lon]) * coef_lat + min(0,elevation[id_lat-1, id_lon] * (1-coef_lat))

        return lon_high * coef_lon + lon_low * (1 - coef_lon)

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


def update_status(sim_ID, status, path_api, msg = ""):
    """
    Update the status of the simulation sim_ID with status and the message msg
    """
    proc = subprocess.Popen(f"python3 {path_api} \"{sim_ID}\" \"{status}\" \"{msg}\"", shell=True)
    proc.wait()

def load_and_check(sim_name, data, name, units, bypass_verif = 0, type = None):
    """
    Return the value of the name inside of type from data, check if units
    are matching. Change "^2" and "^3" to "²" and "³"

    params
    -----
    data: data from the json
    type: where is the name
    name: name of the value
    units: units of the value
    bypass_verif: if this value is at 1, the type of unit will not be verified
            if the value of interest is = -1.
    """
    if type == None:
        value = data[name]["value"]
        units_json = data[name]["units"]
    else:
        value = data[type][name]["value"]
        units_json = data[type][name]["units"]
    units_json = replace_exponents(units_json)
    units = replace_exponents(units)
    if units != units_json and (bypass_verif != 1 or bypass_verif == 1 and value != -1):
        logging.error(f'{type} : {name} has a units issue : {units_json} != {units}')
        update_status(sim_name, "ERROR",f'{type} : {name} has a units issue : {units_json} != {units}')
        raise Exception(f'{type} : {name} has a units issue : {units_json} != {units}')

    return value

def replace_exponents(text):
    return text.replace("^2", "²").replace("^3", "³")

def test_range(sim_name, to_test, min, max, name):
    """
    Test if the value "to_test" is smaller than max and bigger than min
    If not, an exception occurs, and take into account the name for the message
    """
    if to_test < min or to_test > max:
        logging.error(f"{name} Value out of range : {to_test} not inside [{min}, {max}]")
        update_status(sim_name, "ERROR",f"{name} Value out of range : {to_test} not inside [{min}, {max}]")
        raise Exception(f"{name} Value out of range : {to_test} not inside [{min}, {max}]")


def get_first_free_sim_id():
    sim_id = 0
    while os.path.exists(f'requests/{sim_id:05d}_raw_request.json'):
        sim_id += 1
    return sim_id

def create_raw_json(config, sim_id):
    os.makedirs('requests', exist_ok=True)
    with open(f'requests/{sim_id:05d}_raw_request.json', 'w') as f:
        json.dump(config, f, indent=2)

def get_results_folder(sim_id):
    return f'results/{sim_id:05d}'

def breach_flow_rate(rho_water, heigth, rho_hns, has_two_hole = False, g = 9.81):
    if has_two_hole:
        k = 0.703
    else:
        k = 0.246
    return k * np.sqrt(2*g*heigth * (rho_water - rho_hns)/rho_water)

def pipeline_flow_rate(depth, rho_water, T, mw_hns, pipeline_p, R = 8.314, g = 9.81, atm_P = 101325):
    p = depth*rho_water*g + atm_P
    return np.sqrt(2*R*T/(p*mw_hns)*(pipeline_p - p))

def run(gui, sim_name, path_api, already_preprocessed):
        if gui is not None:
            gui.is_running = True
        logging.info(f'Starting the simulation {sim_name}')
        
        if sim_name > 95000:
            logging.warning(f"Only {99999-sim_name} simulations remaining before running out of ids")

        #update_status(sim_name, "COMPUTING", path_api)


        try:
            if not already_preprocessed:
                # Preparing the request
                try:
                    f = open(f'requests/{sim_name:05d}_raw_request.json')
                except OSError:
                    logging.error(f"Unable to find requests/{sim_name:05d}_raw_request.json")
                    update_status(sim_name, "ERROR", f"Unable to find requests/{sim_name:05d}_raw_request.json")
                    sys.exit()
                data = json.load(f)
                f.close()

                #must read the bathymmetry
                if data["location_depth_units"].lower() == "m above seabed":
                    tot_depth = get_depth(data["location"]["coordinates"][0], data["location"]["coordinates"][1])
                    if tot_depth > 0:
                        logging.error(f"{sim_name:05d} location is not at sea")
                        update_status(sim_name, "ERROR", f"{sim_name:05d} location is not at sea")
                        sys.exit()
                    depth = -(tot_depth + data["location"]["coordinates"][2])
                    depth = max(0, depth)

                    data["location_depth_units"] = "m below sea surface"
                    data["location"]["coordinates"][2] = depth

                # switch water temperature to K if needed
                if data["environment"]["water_properties"]["temperature"]["units"] == "°C":
                    data["environment"]["water_properties"]["temperature"]["value"] += 273.15
                    data["environment"]["water_properties"]["temperature"]["units"] = "K"




                data["components"][0]["partial_pressure"] = {
                    "value": 0,
                    "units":"Pa"
                }
                data["components"][0]["initial_fraction_mol"] = {
                    "value":100,
                    "units":"%"
                }

                #compute the molar volume and diffusion coefficient
                molar_vol = load_and_check(sim_name, data["components"][0], "molar_volume", "m³/mol")
                molar_mass = load_and_check(sim_name, data["components"][0], "molar_mass", "g/mol")
                air_diffusivity = load_and_check(sim_name, data["components"][0], "air_diffusivity", "m²/s")
                if molar_vol == -1:
                    data["components"][0]["molar_volume"]["value"] = molar_volume_Le_Bas(molar_mass/1000, organic = False)
                
                if air_diffusivity == -1:
                    w_visc = load_and_check(sim_name, data["environment"], "viscosity", "Pa s", type ="water_properties")
                    data["components"][0]["air_diffusivity"]["value"] = diffusion_coefficient(data["components"][0]["molar_volume"]["value"], wat_viscosity = w_visc)
                

                #try reading the oxygen and nitrogen json files
                try:
                    f = open(f'gases/oxygen.json')
                except OSError:
                    logging.error(f"Unable to find gases/oxygen.json")
                    update_status(sim_name, "ERROR", f"Unable to find gases/oxygen.json")
                    sys.exit()
                oxygen = json.load(f)
                f.close()

                try:
                    f = open(f'gases/nitrogen.json')
                except OSError:
                    logging.error(f"Unable to find gases/nitrogen.json")
                    update_status(sim_name, "ERROR", f"Unable to find gases/nitrogen.json")
                    sys.exit()
                nitrogen  = json.load(f)
                f.close()


                data["components"].append(oxygen)
                data["components"].append(nitrogen)


                try:
                    with open(f"requests/{sim_name:05d}_request.json", "w") as f:
                        json.dump(data, f, indent = 4)
                except OSError:
                    logging.error(f"Unable to write to requests/{sim_name:05d}_request.json")
                    update_status(sim_name, "ERROR", f"Unable to write to requests/{sim_name:05d}_request.json")
                    sys.exit()

            # Reading the prepared data

            try:
                f = open(f'requests/{sim_name:05d}_request.json')
            except OSError:
                logging.error(f"Unable to find requests/{sim_name:05d}_request.json")
                update_status(sim_name, "ERROR", f"Unable to find requests/{sim_name:05d}_request.json")
                sys.exit()
            data = json.load(f)
            f.close()


            duration = load_and_check(sim_name, data, "duration", "s", type ="time_info")
            if duration == -1:
                duration = 1e20
            lat = data["location"]["coordinates"][1]
            test_range(sim_name, lat, -360, 360, "latitude")
            lon = data["location"]["coordinates"][0]
            test_range(sim_name, lat, -360, 360, "longitude")
            if data["location_depth_units"].lower() == "m below sea surface":
                depth = data["location"]["coordinates"][2]
            else:
                logging.error(f"Error of depth units")
                sys.exit(0)

            wind_speed = load_and_check(sim_name, data, "wind_speed", "m/s", type ="environment")
            w_t = load_and_check(sim_name, data["environment"], "temperature", "K", type ="water_properties")
            w_visc = load_and_check(sim_name, data["environment"], "viscosity", "Pa s", type ="water_properties")
            w_density = load_and_check(sim_name, data["environment"], "density", "kg/m³", type ="water_properties")
            current = load_and_check(sim_name, data["environment"], "current_speed", "m/s", type ="water_properties")
            K_z = load_and_check(sim_name, data["environment"], "vertical_turbulent_diffusivity", "m^2/s", type ="water_properties")
            K_xy = load_and_check(sim_name, data["environment"], "horizontal_turbulent_diffusivity", "m^2/s", type ="water_properties")
            is_clean = data["environment"]["water_properties"]["is_clean"]

            flow_rate = load_and_check(sim_name, data, "flow_rate", "kg/s", type ="release")
            breach_area = load_and_check(sim_name, data, "breach_area", "m²", type ="release")

            name = []
            partial_pressure = []
            initial_fraction = []
            molar_mass = []
            critical_temperature = []
            henry = []
            molar_volume = []
            air_diffusivity = []
            solubility = []

            for comp in data["components"]:
                name.append(comp["name"])
                partial_pressure.append(load_and_check(sim_name, comp, "partial_pressure", "Pa"))
                initial_fraction.append(load_and_check(sim_name, comp, "initial_fraction_mol", "%")/100)
                molar_mass.append(load_and_check(sim_name, comp, "molar_mass", "g/mol") / 1000) #to kg/mol
                if comp["critical_temperature"]["units"] == "°C":
                    logging.warning(f"Critical temperature of {comp['name']} is in °C, it will be treated as if it was K. Only K is supported")
                    comp["critical_temperature"]["units"] = "K"
                critical_temperature.append(load_and_check(sim_name, comp, "critical_temperature", "K"))
                henry.append(load_and_check(sim_name, comp, "henry", "mol/(m³ Pa)"))
                molar_volume.append(load_and_check(sim_name, comp, "molar_volume", "m³/mol"))
                air_diffusivity.append(load_and_check(sim_name, comp, "air_diffusivity", "m²/s"))

            chara_height = load_and_check(sim_name, data, "jet_characteristic_height", "m", type = "release")
            total_mass = load_and_check(sim_name, data, "total_mass", "kg", type = "release")
            avg_diam = load_and_check(sim_name, data, "average_bubble_diameter", "m", type = "release")
            mu_bubble = math.log(avg_diam)

            max_change = data["advanced"]["max_shape_change"]
            timestep_per_part = data["advanced"]["timestep_per_part"]
            number_layers = data["advanced"]["number_layers"]
            sigma_bubble = data["advanced"]["sigma_bubble"]
            x_res = data["advanced"]["output_grid_x"]
            y_res = data["advanced"]["output_grid_y"]
            alias_factor = data["advanced"]["aliasing_factor"]

            atm_pressure = sum(partial_pressure)
            air_conc = []
            for i in range(len(name)):
                air_conc.append(partial_pressure[i] * atm_pressure)

            time = data["time_info"]["simulation_start_time"]
            date_format = "%Y-%m-%dT%H:%M:%SZ"
            datetime_obj = datetime.strptime(time, date_format)

            # Convert the datetime object to seconds since epoch
            epoch_time = int(datetime_obj.timestamp())

            if not os.path.exists(f'results/') :
                os.mkdir(f'results/')

            if not os.path.exists(f'results/{sim_name:05d}/') :
                os.mkdir(f'results/{sim_name:05d}/')

            out_files = f"results/{sim_name:05d}/{sim_name:05d}"
            out_zip = f"results/SubSeaGasLeak_results_{sim_name:05d}_{random_char(3)}.zip"
            
            if not os.path.exists(f'results/{sim_name:05d}/timeseries/') :
                os.mkdir(f'results/{sim_name:05d}/timeseries/')
            out_timeseries = f'results/{sim_name:05d}/timeseries/{sim_name:05d}_'

            path_json_part, path_csv_layers, path_metadata, path_release = run_model(
                depth, wind_speed, w_density, w_t, w_visc, current, is_clean, atm_pressure, breach_area, flow_rate,
                lon, lat, epoch_time, duration, number_layers, timestep_per_part, max_change, chara_height,
                K_z, K_xy, name, air_conc, initial_fraction, molar_mass, critical_temperature, henry, molar_volume, air_diffusivity,
                mu_bubble, sigma_bubble, False, out_files, GUI = gui
            )
            if gui is not None:
                gui.in_post = True
                gui.model_status.setText(f"Simulation in postprocessing")
            id_gas_interest = 0 
            #print(path_json_part, path_csv_layers, path_metadata, path_release)
            run_postproc(path_metadata, 
                            path_json_part, 
                            path_csv_layers, 
                            path_release, 
                            id_gas_interest,
                            x_res, 
                            y_res, 
                            alias_factor,
                            molar_mass,
                            out_files,
                            out_zip,
                            out_timeseries,
                            lon, 
                            lat,
                            epoch_time
                            )

        

        except Exception as e:
            logging.error(f"Error in the json :{traceback.format_exc()}")
            update_status(sim_name, "ERROR", f"Error in the json :{e}")
            sys.exit()



        logging.info(f"Simulation {sim_name} is finished")
        if gui is not None:
            gui.is_running = False
            gui.in_post = False
            gui.model_status.setText(f"Simulation {sim_name} is finished")

if __name__ == "__main__":

    #initialize the logger
    log_name = datetime.today().strftime('%Y_%m')

    if not os.path.exists(f'log/') :
        os.mkdir(f'log/')

    handler = logging.handlers.WatchedFileHandler(os.environ.get("LOGFILE", f"log/{log_name}.log"))
    formatter = logging.Formatter(fmt='%(asctime)s %(pathname)s:%(lineno)s %(levelname)-8s %(message)s',datefmt='%d %H:%M:%S')
    handler.setFormatter(formatter)
    logging.basicConfig(format='%(asctime)s %(pathname)s:%(lineno)s %(levelname)-8s %(message)s',
                        level=logging.INFO, datefmt='%d %H:%M:%S')
    root = logging.getLogger()
    root.addHandler(handler)

    try:

        #grab the simulation name
        parser = argparse.ArgumentParser()
        parser.add_argument("--no_gui", action="store_true", help="If the simulation prepared request is already there")
        parser.add_argument("--sim_ID", type=int, help="ID of the simulation")
        parser.add_argument("--path_API", help="path of the API")
        parser.add_argument("--already_preprocessed", action="store_true", help="If the simulation prepared request is already there")
        args = parser.parse_args()

        sim_name = args.sim_ID
        path_api = args.path_API
        already_preprocessed = args.already_preprocessed
        
        gui_window = None

        if not args.no_gui:
            import gui
            gui_window = gui.create_window()
        else:
            run(gui_window, sim_name, path_api, already_preprocessed)

    except Exception as e:
        logging.error(f"{traceback.format_exc()}")
        update_status(sim_name, "ERROR", path_api, f"{traceback.format_exc()}")
        sys.exit()
    

