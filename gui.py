import sys
import json
import run_bubble
import os
from datetime import datetime
from PySide6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                               QHBoxLayout, QFormLayout, QLineEdit, QSpinBox, 
                               QDoubleSpinBox, QCheckBox, QPushButton, QGroupBox, 
                               QScrollArea, QMessageBox, QFileDialog, QDateTimeEdit,
                               QLabel, QTabWidget, QGridLayout, QDialog, QDialogButtonBox)
from PySide6.QtCore import Qt, QDateTime
from PySide6.QtGui import QFont, QIcon

from threading import Thread



import subprocess
import platform

def open_file_manager(path=None):
    system = platform.system()
    
    if system == "Windows":
        if path:
            subprocess.run(['explorer', path])
        else:
            subprocess.run(['explorer'])
    elif system == "Linux":
        if path:
            subprocess.run(['xdg-open', path])
        else:
            subprocess.run(['xdg-open', '.'])  # Opens current directory
    else:
        print(f"Unsupported operating system: {system}")



class ComponentWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setup_ui()

    
    def setup_ui(self):
        layout = QFormLayout(self)
        
        self.name_edit = QLineEdit("methane")
        self.molar_mass_edit = QDoubleSpinBox()
        self.molar_mass_edit.setRange(0.1, 1000)
        self.molar_mass_edit.setValue(16.0)
        self.molar_mass_edit.setSuffix(" g/mol")
        
        self.critical_temp_edit = QDoubleSpinBox()
        self.critical_temp_edit.setRange(50, 1000)
        self.critical_temp_edit.setValue(190.56)
        self.critical_temp_edit.setSuffix(" K")
        
        self.henry_edit = QDoubleSpinBox()
        self.henry_edit.setRange(0.001, 10)
        self.henry_edit.setValue(0.015)
        self.henry_edit.setDecimals(6)
        self.henry_edit.setSuffix(" mol/(m³ Pa)")
        
        self.molar_volume_edit = QDoubleSpinBox()
        self.molar_volume_edit.setRange(-1, 1)
        self.molar_volume_edit.setValue(-1.0)
        self.molar_volume_edit.setDecimals(3)
        self.molar_volume_edit.setSuffix(" m³/mol")
        
        self.air_diffusivity_edit = QDoubleSpinBox()
        self.air_diffusivity_edit.setRange(-1, 1)
        self.air_diffusivity_edit.setValue(-1.0)
        self.air_diffusivity_edit.setDecimals(3)
        self.air_diffusivity_edit.setSuffix(" m²/s")
        
        layout.addRow("Component Name:", self.name_edit)
        layout.addRow("Molar Mass:", self.molar_mass_edit)
        layout.addRow("Critical Temperature:", self.critical_temp_edit)
        layout.addRow("Henry Constant:", self.henry_edit)
        layout.addRow("Molar Volume (-1 means no value):", self.molar_volume_edit)
        layout.addRow("Air Diffusivity (-1 means no value):", self.air_diffusivity_edit)
    
    def get_data(self):
        return {
            "name": self.name_edit.text(),
            "molar_mass": {
                "value": self.molar_mass_edit.value(),
                "units": "g/mol"
            },
            "critical_temperature": {
                "value": self.critical_temp_edit.value(),
                "units": "K"
            },
            "henry": {
                "value": self.henry_edit.value(),
                "units": "mol/(m³ Pa)"
            },
            "molar_volume": {
                "value": self.molar_volume_edit.value(),
                "units": "m³/mol"
            },
            "air_diffusivity": {
                "value": self.air_diffusivity_edit.value(),
                "units": "m²/s"
            }
        }


class UnderwaterGasReleaseGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Underwater Gas Release Model Configuration")
        self.setGeometry(100, 100, 800, 700)
        self.setup_ui()

        self.current_id = -1
        self.is_running = False
        self.in_post = False
    
    def setup_ui(self):
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        # Main layout
        main_layout = QVBoxLayout(central_widget)
        
        # Create tab widget
        tab_widget = QTabWidget()
        main_layout.addWidget(tab_widget)
        
        # Time & Location Tab
        time_location_tab = QWidget()
        tab_widget.addTab(time_location_tab, "Time & Location")
        self.setup_time_location_tab(time_location_tab)
        
        # Environment Tab
        environment_tab = QWidget()
        tab_widget.addTab(environment_tab, "Environment")
        self.setup_environment_tab(environment_tab)
        
        # Components Tab
        components_tab = QWidget()
        tab_widget.addTab(components_tab, "Chemical properties")
        self.setup_components_tab(components_tab)

        # Release Tab
        release_tab = QWidget()
        tab_widget.addTab(release_tab, "Release Parameters")
        self.setup_release_tab(release_tab)
        
        # Advanced Tab
        advanced_tab = QWidget()
        tab_widget.addTab(advanced_tab, "Advanced Settings")
        self.setup_advanced_tab(advanced_tab)
        
        # Buttons
        bottom_layout = QVBoxLayout()

        status_layout = QHBoxLayout()
        self.model_status = QLabel("")
        status_layout.addWidget(self.model_status)
        

        button_layout = QHBoxLayout()
        
        self.load_button = QPushButton("Load JSON")
        self.save_button = QPushButton("Save JSON")
        self.run_button = QPushButton("Run")
        self.stop_button = QPushButton("Stop Simulation to do post-processing")
        self.open_folder = QPushButton("Open results folder")
        self.generate_button = QPushButton("Generate & Preview JSON")
        
        self.load_button.clicked.connect(self.load_json)
        self.save_button.clicked.connect(self.save_json)
        self.run_button.clicked.connect(self.run_model)
        self.stop_button.clicked.connect(self.stop_model)
        self.open_folder.clicked.connect(self.open_results_folder)
        self.generate_button.clicked.connect(self.preview_json)
        
        # button_layout.addWidget(self.load_button)
        # button_layout.addWidget(self.save_button)
        button_layout.addWidget(self.run_button)
        button_layout.addWidget(self.stop_button)
        button_layout.addWidget(self.open_folder)
        button_layout.addStretch()
        button_layout.addWidget(self.generate_button)
        
        bottom_layout.addLayout(status_layout)
        bottom_layout.addLayout(button_layout)

        main_layout.addLayout(bottom_layout)
    
    def setup_time_location_tab(self, tab):
        layout = QVBoxLayout(tab)
        
        # Time Information Group
        time_group = QGroupBox("Time Information")
        time_layout = QFormLayout(time_group)
        
        self.start_time_edit = QDateTimeEdit()
        self.start_time_edit.setDateTime(QDateTime.currentDateTime())
        self.start_time_edit.setDisplayFormat("yyyy-MM-ddThh:mm:ssZ")
        
        self.duration_edit = QDoubleSpinBox()
        self.duration_edit.setRange(-1, 999999)
        self.duration_edit.setValue(3600)
        self.duration_edit.setSuffix(" s")
        
        time_layout.addRow("Start Time:", self.start_time_edit)
        time_layout.addRow("Maximum duration:", self.duration_edit)
        
        # Location Group
        location_group = QGroupBox("Location")
        location_layout = QFormLayout(location_group)
        
        self.latitude_edit = QDoubleSpinBox()
        self.latitude_edit.setRange(-90, 90)
        self.latitude_edit.setValue(51.2)
        self.latitude_edit.setDecimals(6)
        self.latitude_edit.setSuffix(" °")
        
        self.longitude_edit = QDoubleSpinBox()
        self.longitude_edit.setRange(-180, 180)
        self.longitude_edit.setValue(2.5)
        self.longitude_edit.setDecimals(6)
        self.longitude_edit.setSuffix(" °")
        
        self.depth_edit = QDoubleSpinBox()
        self.depth_edit.setRange(0, 10000)
        self.depth_edit.setValue(10)
        self.depth_edit.setSuffix(" m")
        
        location_layout.addRow("Latitude:", self.latitude_edit)
        location_layout.addRow("Longitude:", self.longitude_edit)
        location_layout.addRow("Depth:", self.depth_edit)
        
        layout.addWidget(time_group)
        layout.addWidget(location_group)
        layout.addStretch()
    
    def setup_environment_tab(self, tab):
        layout = QVBoxLayout(tab)
        
        # Wind Group
        wind_group = QGroupBox("Wind")
        wind_layout = QFormLayout(wind_group)
        
        self.wind_speed_edit = QDoubleSpinBox()
        self.wind_speed_edit.setRange(0, 100)
        self.wind_speed_edit.setValue(10)
        self.wind_speed_edit.setSuffix(" m/s")
        
        wind_layout.addRow("Wind Speed:", self.wind_speed_edit)
        
        # Water Properties Group
        water_group = QGroupBox("Water Properties")
        water_layout = QFormLayout(water_group)
        
        self.temperature_edit = QDoubleSpinBox()
        self.temperature_edit.setRange(-10, 50)
        self.temperature_edit.setValue(15)
        self.temperature_edit.setSuffix(" °C")
        
        self.viscosity_edit = QDoubleSpinBox()
        self.viscosity_edit.setRange(0.0001, 0.01)
        self.viscosity_edit.setValue(0.001)
        self.viscosity_edit.setDecimals(6)
        self.viscosity_edit.setSuffix(" Pa s")
        
        self.density_edit = QDoubleSpinBox()
        self.density_edit.setRange(900, 1100)
        self.density_edit.setValue(1013)
        self.density_edit.setSuffix(" kg/m³")
        
        # Current speed (3 components)
        current_layout = QHBoxLayout()
        self.current_x_edit = QDoubleSpinBox()
        self.current_x_edit.setRange(-10, 10)
        self.current_x_edit.setValue(0.05)
        self.current_x_edit.setDecimals(3)
        
        self.current_y_edit = QDoubleSpinBox()
        self.current_y_edit.setRange(-10, 10)
        self.current_y_edit.setValue(0.05)
        self.current_y_edit.setDecimals(3)
        
        self.current_z_edit = QDoubleSpinBox()
        self.current_z_edit.setRange(-10, 10)
        self.current_z_edit.setValue(0)
        self.current_z_edit.setDecimals(3)
        
        current_layout.addWidget(QLabel("X:"))
        current_layout.addWidget(self.current_x_edit)
        current_layout.addWidget(QLabel("Y:"))
        current_layout.addWidget(self.current_y_edit)
        current_layout.addWidget(QLabel("Z:"))
        current_layout.addWidget(self.current_z_edit)
        current_layout.addWidget(QLabel("m/s"))
        
        self.vertical_diff_edit = QDoubleSpinBox()
        self.vertical_diff_edit.setRange(0.001, 10)
        self.vertical_diff_edit.setValue(0.1)
        self.vertical_diff_edit.setDecimals(3)
        self.vertical_diff_edit.setSuffix(" m²/s")
        
        self.horizontal_diff_edit = QDoubleSpinBox()
        self.horizontal_diff_edit.setRange(0.001, 10)
        self.horizontal_diff_edit.setValue(1.5)
        self.horizontal_diff_edit.setDecimals(3)
        self.horizontal_diff_edit.setSuffix(" m²/s")
        
        self.is_clean_check = QCheckBox()
        self.is_clean_check.setChecked(False)
        
        water_layout.addRow("Temperature:", self.temperature_edit)
        water_layout.addRow("Viscosity:", self.viscosity_edit)
        water_layout.addRow("Density:", self.density_edit)
        water_layout.addRow("Current Speed:", current_layout)
        water_layout.addRow("Vertical Turbulent Diffusivity:", self.vertical_diff_edit)
        water_layout.addRow("Horizontal Turbulent Diffusivity:", self.horizontal_diff_edit)
        water_layout.addRow("Is the water pure? (not salty?):", self.is_clean_check)
        
        layout.addWidget(wind_group)
        layout.addWidget(water_group)
        layout.addStretch()
    
    def setup_release_tab(self, tab):
        layout = QVBoxLayout(tab)
        
        release_group = QGroupBox("Release Parameters")
        release_layout = QFormLayout(release_group)
        
        self.flow_rate_edit = QDoubleSpinBox()
        self.flow_rate_edit.setRange(0.001, 10000)
        self.flow_rate_edit.setValue(0.01)
        self.flow_rate_edit.setDecimals(4)
        self.flow_rate_edit.setSuffix(" kg/s")
        
        self.breach_area_edit = QDoubleSpinBox()
        self.breach_area_edit.setRange(0.001, 100)
        self.breach_area_edit.setValue(0.01)
        self.breach_area_edit.setDecimals(4)
        self.breach_area_edit.setSuffix(" m²")
        
        self.jet_height_edit = QDoubleSpinBox()
        self.jet_height_edit.setRange(0.1, 100)
        self.jet_height_edit.setValue(3)
        self.jet_height_edit.setSuffix(" m")
        
        self.bubble_diameter_edit = QDoubleSpinBox()
        self.bubble_diameter_edit.setRange(0.001, 1)
        self.bubble_diameter_edit.setValue(0.005)
        self.bubble_diameter_edit.setDecimals(6)
        self.bubble_diameter_edit.setSuffix(" m")
        
        self.total_mass_edit = QDoubleSpinBox()
        self.total_mass_edit.setRange(1, 10000)
        self.total_mass_edit.setValue(150)
        self.total_mass_edit.setSuffix(" kg")
        
        self.compute_flow_rate = QPushButton("Estimate Flow Rate")
        self.compute_flow_rate.clicked.connect(self.open_flow_rate_dialog)
        flow_rate_widget = QWidget()
        flow_rate_layout = QHBoxLayout(flow_rate_widget)
        flow_rate_layout.addWidget(self.flow_rate_edit)
        flow_rate_layout.addWidget(self.compute_flow_rate)
        flow_rate_layout.setContentsMargins(0, 0, 0, 0)  # Optional: remove margins

        # Add the combined widget to the form layout
        release_layout.addRow("Breach Area:", self.breach_area_edit)
        release_layout.addRow("Jet Characteristic Height:", self.jet_height_edit)
        release_layout.addRow("Average Bubble Diameter:", self.bubble_diameter_edit)
        release_layout.addRow("Total Mass:", self.total_mass_edit)
        release_layout.addRow("Flow Rate:", flow_rate_widget)
        
        layout.addWidget(release_group)
        layout.addStretch()
    
    def setup_components_tab(self, tab):
        layout = QVBoxLayout(tab)
        
        components_group = QGroupBox("Chemical properties")
        components_layout = QVBoxLayout(components_group)
        
        # For now, just one component. Can be extended for multiple components
        self.component_widget = ComponentWidget()
        components_layout.addWidget(self.component_widget)
        
        layout.addWidget(components_group)
        layout.addStretch()
    
    def setup_advanced_tab(self, tab):
        layout = QVBoxLayout(tab)
        
        advanced_group = QGroupBox("Advanced Simulation Settings")
        advanced_layout = QFormLayout(advanced_group)
        
        self.max_shape_change_edit = QSpinBox()
        self.max_shape_change_edit.setRange(1, 1000)
        self.max_shape_change_edit.setValue(100)
        
        self.timestep_per_part_edit = QSpinBox()
        self.timestep_per_part_edit.setRange(1, 100)
        self.timestep_per_part_edit.setValue(2)
        
        self.number_layers_edit = QSpinBox()
        self.number_layers_edit.setRange(1, 100)
        self.number_layers_edit.setValue(20)
        
        self.sigma_bubble_edit = QDoubleSpinBox()
        self.sigma_bubble_edit.setRange(0.1, 2)
        self.sigma_bubble_edit.setValue(0.635)
        self.sigma_bubble_edit.setDecimals(3)
        
        self.output_grid_x_edit = QSpinBox()
        self.output_grid_x_edit.setRange(1, 100)
        self.output_grid_x_edit.setValue(10)
        
        self.output_grid_y_edit = QSpinBox()
        self.output_grid_y_edit.setRange(1, 100)
        self.output_grid_y_edit.setValue(10)
        
        self.aliasing_factor_edit = QSpinBox()
        self.aliasing_factor_edit.setRange(1, 10)
        self.aliasing_factor_edit.setValue(2)
        
        advanced_layout.addRow("Max Shape Change:", self.max_shape_change_edit)
        advanced_layout.addRow("Timestep Per Part:", self.timestep_per_part_edit)
        advanced_layout.addRow("Number of Layers:", self.number_layers_edit)
        advanced_layout.addRow("Sigma Bubble:", self.sigma_bubble_edit)
        advanced_layout.addRow("Output Grid X:", self.output_grid_x_edit)
        advanced_layout.addRow("Output Grid Y:", self.output_grid_y_edit)
        advanced_layout.addRow("Aliasing Factor:", self.aliasing_factor_edit)
        
        layout.addWidget(advanced_group)
        layout.addStretch()
    
    def generate_json(self):
        """Generate the JSON configuration from current GUI values."""
        # Format datetime to match the expected format
        start_time = self.start_time_edit.dateTime().toString("yyyy-MM-ddThh:mm:ss") + "Z"
        
        config = {
            "time_info": {
                "simulation_start_time": start_time,
                "duration": {
                    "value": self.duration_edit.value(),
                    "units": "s"
                }
            },
            "location": {
                "type": "Point",
                "coordinates": [
                    self.longitude_edit.value(),
                    self.latitude_edit.value(),
                    self.depth_edit.value()
                ]
            },
            "location_depth_units": "m below sea surface",
            "environment": {
                "wind_speed": {
                    "value": self.wind_speed_edit.value(),
                    "units": "m/s"
                },
                "water_properties": {
                    "temperature": {
                        "value": self.temperature_edit.value(),
                        "units": "°C"
                    },
                    "viscosity": {
                        "value": self.viscosity_edit.value(),
                        "units": "Pa s"
                    },
                    "density": {
                        "value": self.density_edit.value(),
                        "units": "kg/m³"
                    },
                    "current_speed": {
                        "value": [
                            self.current_x_edit.value(),
                            self.current_y_edit.value(),
                            self.current_z_edit.value()
                        ],
                        "units": "m/s"
                    },
                    "vertical_turbulent_diffusivity": {
                        "value": self.vertical_diff_edit.value(),
                        "units": "m^2/s"
                    },
                    "horizontal_turbulent_diffusivity": {
                        "value": self.horizontal_diff_edit.value(),
                        "units": "m^2/s"
                    },
                    "is_clean": self.is_clean_check.isChecked()
                }
            },
            "release": {
                "flow_rate": {
                    "value": self.flow_rate_edit.value(),
                    "units": "kg/s"
                },
                "breach_area": {
                    "value": self.breach_area_edit.value(),
                    "units": "m²"
                },
                "jet_characteristic_height": {
                    "value": self.jet_height_edit.value(),
                    "units": "m"
                },
                "average_bubble_diameter": {
                    "value": self.bubble_diameter_edit.value(),
                    "units": "m"
                },
                "total_mass": {
                    "value": self.total_mass_edit.value(),
                    "units": "kg"
                }
            },
            "components": [self.component_widget.get_data()],
            "advanced": {
                "max_shape_change": self.max_shape_change_edit.value(),
                "timestep_per_part": self.timestep_per_part_edit.value(),
                "number_layers": self.number_layers_edit.value(),
                "sigma_bubble": self.sigma_bubble_edit.value(),
                "output_grid_x": self.output_grid_x_edit.value(),
                "output_grid_y": self.output_grid_y_edit.value(),
                "aliasing_factor": self.aliasing_factor_edit.value()
            }
        }
        
        return config
    
    def preview_json(self):
        """Generate and preview the JSON in a message box."""
        try:
            config = self.generate_json()
            json_str = json.dumps(config, indent=2)
            
            msg = QMessageBox()
            msg.setWindowTitle("JSON Preview")
            msg.setText("Generated JSON Configuration:")
            msg.setDetailedText(json_str)
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to generate JSON: {str(e)}")
    
    def save_json(self):
        """Save the JSON configuration to a file."""
        try:
            config = self.generate_json()
            
            filename, _ = QFileDialog.getSaveFileName(
                self, "Save JSON Configuration", 
                "gas_release_config.json", 
                "JSON Files (*.json);;All Files (*)"
            )
            
            if filename:
                with open(filename, 'w') as f:
                    json.dump(config, f, indent=2)
                QMessageBox.information(self, "Success", f"Configuration saved to {filename}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save JSON: {str(e)}")
    
    def load_json(self):
        """Load a JSON configuration from a file."""
        try:
            filename, _ = QFileDialog.getOpenFileName(
                self, "Load JSON Configuration", 
                "", 
                "JSON Files (*.json);;All Files (*)"
            )
            
            if filename:
                with open(filename, 'r') as f:
                    config = json.load(f)
                
                self.load_config(config)
                QMessageBox.information(self, "Success", f"Configuration loaded from {filename}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to load JSON: {str(e)}")
    
    def load_config(self, config):
        """Load configuration data into the GUI."""
        try:
            # Time info
            start_time_str = config["time_info"]["start_time"].replace("Z", "")
            start_time = QDateTime.fromString(start_time_str, "yyyy/MM/dd;hh:mm:ss")
            self.start_time_edit.setDateTime(start_time)
            self.duration_edit.setValue(config["time_info"]["duration"]["value"])
            
            # Location
            self.latitude_edit.setValue(config["location"]["latitude"])
            self.longitude_edit.setValue(config["location"]["longitude"])
            self.depth_edit.setValue(config["location"]["depth"]["value"])
            
            # Environment
            self.wind_speed_edit.setValue(config["environment"]["wind_speed"]["value"])
            
            water_props = config["environment"]["water_properties"]
            self.temperature_edit.setValue(water_props["temperature"]["value"])
            self.viscosity_edit.setValue(water_props["viscosity"]["value"])
            self.density_edit.setValue(water_props["density"]["value"])
            
            current = water_props["current_speed"]["value"]
            self.current_x_edit.setValue(current[0])
            self.current_y_edit.setValue(current[1])
            self.current_z_edit.setValue(current[2])
            
            self.vertical_diff_edit.setValue(water_props["vertical_turbulent_diffusivity"]["value"])
            self.horizontal_diff_edit.setValue(water_props["horizontal_turbulent_diffusivity"]["value"])
            self.is_clean_check.setChecked(water_props["is_clean"])
            
            # Release
            release = config["release"]
            self.flow_rate_edit.setValue(release["flow_rate"]["value"])
            self.breach_area_edit.setValue(release["breach_area"]["value"])
            self.jet_height_edit.setValue(release["jet_characteristic_height"]["value"])
            self.bubble_diameter_edit.setValue(release["average_bubble_diameter"]["value"])
            self.total_mass_edit.setValue(release["total_mass"]["value"])
            
            # Components (first component only)
            if config["components"]:
                comp = config["components"][0]
                self.component_widget.name_edit.setText(comp["name"])
                self.component_widget.molar_mass_edit.setValue(comp["molar_mass"]["value"])
                self.component_widget.critical_temp_edit.setValue(comp["critical_temperature"]["value"])
                self.component_widget.henry_edit.setValue(comp["henry"]["value"])
                self.component_widget.molar_volume_edit.setValue(comp["molar_volume"]["value"])
                self.component_widget.air_diffusivity_edit.setValue(comp["air_diffusivity"]["value"])
            
            # Advanced
            advanced = config["advanced"]
            self.max_shape_change_edit.setValue(advanced["max_shape_change"])
            self.timestep_per_part_edit.setValue(advanced["timestep_per_part"])
            self.number_layers_edit.setValue(advanced["number_layers"])
            self.sigma_bubble_edit.setValue(advanced["sigma_bubble"])
            self.output_grid_x_edit.setValue(advanced["output_grid_x"])
            self.output_grid_y_edit.setValue(advanced["output_grid_y"])
            self.aliasing_factor_edit.setValue(advanced["aliasing_factor"])
            
        except KeyError as e:
            QMessageBox.warning(self, "Warning", f"Missing field in JSON: {str(e)}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to load configuration: {str(e)}")

    def run_model(self):
        self.model_status.setText("Model starting...")


        if self.is_running:
            QMessageBox.warning(self, "Warning", "A simulation is already running.")
            return

        self.current_id = run_bubble.get_first_free_sim_id()
        run_bubble.create_raw_json(self.generate_json(), self.current_id)

        thread = Thread(target = run_bubble.run, args = (self, self.current_id, None, False))
        thread.daemon = True
        thread.start()

    def stop_model(self):
        if self.is_running:
            if not self.in_post:
                self.in_post = True
            else:
                QMessageBox.warning(self, "Warning", "The post-processing is already running.")
        else:
            QMessageBox.warning(self, "Warning", "No simulation is running.")
    
    def open_results_folder(self):
        if self.current_id == -1:
            QMessageBox.warning(self, "Warning", "No simulation started yet.")
        else: #open the system file manager
            open_file_manager(run_bubble.get_results_folder(self.current_id))

    def open_flow_rate_dialog(self):
        dialog = FlowRateDialog(self)
        if dialog.exec() == QDialog.Accepted:
            value = dialog.get_value()
            self.flow_rate_edit.setValue(value)


class FlowRateDialog(QDialog):
    def __init__(self, main_window, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Enter Flow Rate Parameters")
        self.setMinimumWidth(400)

        # Main layout
        layout = QVBoxLayout(self)




        # Tabs
        self.tabs = QTabWidget()
        self.single_hole_tab = self.create_single_or_double_tab(main_window)
        self.double_hole_tab = self.create_single_or_double_tab(main_window)
        self.pipeline_tab = self.create_pipeline_tab(main_window)

        self.tabs.addTab(self.single_hole_tab, "Single Hole")
        self.tabs.addTab(self.double_hole_tab, "Double Hole")
        self.tabs.addTab(self.pipeline_tab, "Pipeline")

        layout.addWidget(self.tabs)

        # Dialog buttons
        self.buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        layout.addWidget(self.buttons)

    def create_single_or_double_tab(self, main_window):
        tab = QWidget()
        layout = QFormLayout(tab)

        density_edit = QDoubleSpinBox()
        density_edit.setRange(900, 1100)
        density_edit.setValue(main_window.density_edit.value())
        density_edit.setSuffix(" kg/m³")
        layout.addRow("Ambient Water Density:", density_edit)

        height_edit = QDoubleSpinBox()
        height_edit.setRange(0.01, 1000)
        height_edit.setValue(5)
        height_edit.setSuffix(" m")
        layout.addRow("Height of the hole above the tank bottom:", height_edit)

        breach_area_edit = QDoubleSpinBox()
        breach_area_edit.setRange(0.001, 100)
        breach_area_edit.setValue(main_window.breach_area_edit.value())
        breach_area_edit.setDecimals(4)
        breach_area_edit.setSuffix(" m²")
        layout.addRow("Breach Area:", breach_area_edit)

        density_hns_edit = QDoubleSpinBox()
        density_hns_edit.setRange(1, 2000)
        density_hns_edit.setValue(660)
        density_hns_edit.setDecimals(4)
        density_hns_edit.setSuffix(" kg/m³")
        layout.addRow("HNS Density:", density_hns_edit)

        return tab

    def create_pipeline_tab(self, main_window):
        tab = QWidget()
        layout = QFormLayout(tab)

        depth_edit = QDoubleSpinBox()
        depth_edit.setRange(0.01, 12000)
        depth_edit.setValue(main_window.depth_edit.value())
        depth_edit.setSuffix(" m")
        layout.addRow("Depth:", depth_edit)

        density_edit = QDoubleSpinBox()
        density_edit.setRange(900, 1100)
        density_edit.setValue(main_window.density_edit.value())
        density_edit.setSuffix(" kg/m³")
        layout.addRow("Ambient Water Density:", density_edit)

        temperature_edit = QDoubleSpinBox()
        temperature_edit.setRange(-10, 50)
        temperature_edit.setValue(main_window.temperature_edit.value())
        temperature_edit.setSuffix(" °C")
        layout.addRow("Temperature:", temperature_edit)

        pressure_edit = QDoubleSpinBox()
        pressure_edit.setRange(1, 1500)
        pressure_edit.setValue(50)
        pressure_edit.setDecimals(2)
        pressure_edit.setSuffix(" Bar")
        layout.addRow("Pipeline pressure:", pressure_edit)

        breach_area_edit = QDoubleSpinBox()
        breach_area_edit.setRange(0.001, 100)
        breach_area_edit.setValue(main_window.breach_area_edit.value())
        breach_area_edit.setDecimals(4)
        breach_area_edit.setSuffix(" m²")
        layout.addRow("Breach Area:", breach_area_edit)

        density_hns_edit = QDoubleSpinBox()
        density_hns_edit.setRange(1, 2000)
        density_hns_edit.setValue(660.0)
        density_hns_edit.setDecimals(4)
        density_hns_edit.setSuffix(" kg/m³")
        layout.addRow("HNS Density:", density_hns_edit)

        molar_weight_edit = QDoubleSpinBox()
        molar_weight_edit.setRange(1, 500)
        molar_weight_edit.setValue(main_window.component_widget.molar_mass_edit.value())
        molar_weight_edit.setDecimals(4)
        molar_weight_edit.setSuffix(" g/mol")
        layout.addRow("Molar Mass:", molar_weight_edit)

        return tab

    def get_value(self):
        current_index = self.tabs.currentIndex()

        if current_index == 0:
            layout = self.single_hole_tab.layout()
        elif current_index == 1:
            layout = self.double_hole_tab.layout()
        else:
            layout = self.pipeline_tab.layout()

        params = []
        for i in range(layout.rowCount()):
            params.append(layout.itemAt(i, QFormLayout.FieldRole).widget().value())

        flow_speed = 0
        area = 0
        density = 0
        if current_index == 0:
            flow_speed = run_bubble.breach_flow_rate(params[0], params[1], params[3], has_two_hole = False)
            area = params[2]
            density = params[3]
        elif current_index == 1:
            flow_speed = run_bubble.breach_flow_rate(params[0], params[1], params[3], has_two_hole = True)
            area = params[2]
            density = params[3]
        else:
            flow_speed = run_bubble.pipeline_flow_rate(params[0], params[1], params[2], params[6]*1000, params[3]*100000)
            area = params[4]
            density = params[5]

        return flow_speed * area * density
    


def create_window():
    app = QApplication()
    app.setStyle('Fusion')
    
    window = UnderwaterGasReleaseGUI()
    window.show()

    sys.exit(app.exec())
