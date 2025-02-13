"""
____________________________________________________________________
Function: Input Initialization (Input) 

  Inputs: 
  
        none

  Outputs: 
  
        Inputs : a dictionary containing datetime information for light calculation
                 and information about previous flux calculations
        start_date : a datetime object representing the starting date and time, will be
                        changed throughout simulation

    This function initializes the "Inputs" dictionary, which is used throughout dFBA simulation
    to pass parameters and saved solutions among the functions that perform optimization and ODE solving.
"""

def Input():

    import numpy as np
    import scipy, Light, Extras
    from datetime import datetime

    # Initialize a dictionary to hold all inputs to the class.
    Inputs = {}

    # User defined inputs here!!! Names will be more descriptive in this user-facing section.
    Inputs['Year'] = 2024
    Inputs['DayNumber'] = 0 # Daynumber 0 is a flag for not using DAYNUM/YY mode.
    Inputs['Month'] = 3 # Month 0 is a flag for not using MM/DD/YYYY mode.
    Inputs['Day'] = 1 # Day 0 is a flag for not using MM/DD/YYYY mode.
    Inputs['Hour'] = 12 # 24-hour time
    Inputs['Minute'] = 0
    Inputs['Second'] = 0
    Inputs['Interval'] = 0 # Defaults to zero. Do not change, interval will be handed in dFBA.
    Inputs['Latitude'] = 75 # Bottom of Arctic Circle
    Inputs['Longitude'] = 0 # Prime Meridian
    Inputs['Timezone'] = 0 # We are at the prime meridian. (This will be automatically calculated.)
    Inputs['Pressure'] = 1013 # DEFAULT: 1013 mb = 1 atmosphere
    Inputs['Temperature'] = 0 # DEFAULT: 10 deg C (arctic will be much colder)

    Inputs['Turbidity'] = 0.084
    Inputs['Water Vapor'] = 1.4164
    Inputs['Ozone'] = 0.3438 # GET CITATIONS FOR THESE FROM ORIGINAL SMARTS2
    Inputs['Albedo'] = 0.1 # Average albedo of ocean surface

    # Input wavelength-dependent data tables.
    Inputs['WVL-ETR'] = np.loadtxt('.\\Data\\WVL-ETR.csv',delimiter=',',dtype='float64', encoding='utf-8-sig')
    Inputs['WVL-ABS'] = np.loadtxt('.\\Data\\WVL-ABS.csv',delimiter=',',dtype='float64', encoding='utf-8-sig')

    # Create an object that will store data for wavelength-dependent light transmission through the atmosphere.
    ICalc = Light.TotSolEng(Inputs)
    # Calculate the total incident extraterrestrial radiation (ETR) based on date, time and location.
    ICalc.CalcTotSolarEng(Inputs)
    # Calculate wavelength-dependent light transmission based on ETR spectrum.
    ICalc.CalcSolarSpec(Inputs)

    # Calculate wavelength-dependent light absorption/scattering by the ocean at a certain depth.
    ICalc.waterabs(0)
    ICalc.DTOTphoto = np.asarray(ICalc.DTOTphoto).squeeze()

    Inputs['Itot'] = scipy.integrate.simpson(ICalc.DTOTphoto,ICalc.allWVLphoto)

    # Integrate only the photoactive region of the incident spectrum to get total photoactive light intensity.
    Inputs['Iphoto'] = Extras.par(ICalc.DTOTphoto,ICalc.allWVLphoto)

    start = str(Inputs['Month']).zfill(2) + "/" + str(Inputs['Day']).zfill(2) + "/" + str(Inputs['Year']) + " " + str(Inputs['Hour'])
    start_date = datetime.strptime(start, "%m/%d/%Y %H")

    return Inputs, start_date