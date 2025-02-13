"""
____________________________________________________________________
Function: Total Solar Energy (TotSolEng) calculation class 

  Inputs: 
  
        Inputs : a dictionary containing datetime information for light calculation
                 and information about previous flux calculations        

    This class is a combined port of two solar light intensity calculators:
        1) "Solar Position and Intensity" (SOLPOS 2.0), written in C by Martin Rymes from NREL:
                https://www.nrel.gov/grid/solar-resource/solpos.html
        2) "Simple Model of the Atmospheric Radiative Transfer of Sunshine" (SMARTS2),
            written in FORTRAN by Christian Gueymard from the Florida Solar Energy Center:
                http://publications.energyresearch.ucf.edu/wp-content/uploads/2018/06/FSEC-PF-270-95.pdf

    SOLPOS 2.0 is used to calculate the integrated total solar light intensity in W/m^2, while SMARTS2
    is used to calculate wavelength-dependent absorption based on the extraterrestrial (top of atmosphere)
    emission spectrum of the sun.

    Ultimately, wavelength-dependent light intensity is very important to calculate the light absorption
    by diatom and cyanobacteria chlorophyll. 

    We take an object-oriented approach as the spectrum and calculated values will need to be
    self-referenced throughout the entirety of the dFBA simulation.
"""

class TotSolEng:
    
    """
    ____________________________________________________________________
    Method: initialization method (init)

    Inputs: 

            Inputs : a dictionary containing datetime information for light calculation
                 and information about previous flux calculations 

        Here, the inputs dictionary is converted into units used in the calculator,
        and stored within the object. Constants and flags are also defined.
    """  
    def __init__(self,Inputs):
        
        import math
        import pandas as pd
        import numpy as np
        import matplotlib.pyplot as plt
        import scipy
        import time
        
        # All constants here.
        self.raddeg = 0.0174532925 # radians -> degrees
        self.degrad = 57.295779513 # degrees -> radians
        self.solcon = 1367.0 # From NASA data (cite?)
        self.tilt = 0 # Will always be hitting a horizontal surface
        self.aspect = 180 # Facing does not matter, surface is horizontal
        self.sbwid = 7.6 # shadow band width
        self.sbrad = 31.7 # shadow band radius
        self.sbsky = 0.04 # partly cloudy factor
        self.omega = 0.945 # single-scattering albedo constant
        self.omegap = 0.095
        self.GG = 0.65 # aerosol asymmetry factor
        self.RR = 22 # maximum ozone height
        # Cumulative number of days in the year on the 1st of each month. Index 0 is a placeholder, since month is 1<x<12.
        self.month_days = [0, 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]
    
        self.turb = Inputs['Turbidity']
        self.watervap = Inputs['Water Vapor']
        self.ozone = Inputs['Ozone']
        self.RHO = Inputs['Albedo']
    
        self.year = Inputs['Year']
        self.month = Inputs['Month']
        self.day = Inputs['Day']
        self.daynum = Inputs['DayNumber']
        self.leap = False
        
        # Define leap parameter for error checking. Include all Julian calendar exceptions.
        if self.year % 4 == 0 and (self.year % 100 != 0 or self.year % 400 == 0):
            self.leap == True
        
        # Search for flags:
        # If no month and day reported, use total number of days.
        if Inputs['Month'] == 0 and Inputs['Day'] == 0:
            self.daynum = Inputs['DayNumber']
        # If total number of days not reported, use MM/DD/YYYY.
        elif Inputs['DayNumber'] == 0:
            self.month = Inputs['Month']
            self.day = Inputs['Day']
            self.dom2doy(self.month,self.day)
            
        self.hour = Inputs['Hour']
        self.minute = Inputs['Minute']
        self.second = Inputs['Second']
        self.interval = Inputs['Interval']
        self.latitude = Inputs['Latitude']
        self.longitude = Inputs['Longitude']
        self.timezone = Inputs['Timezone']
        self.press = Inputs['Pressure']
        self.temp = Inputs['Temperature']
        self.allWVL = Inputs['WVL-ETR'][:,0]
        self.allETR = Inputs['WVL-ETR'][:,1]
        self.tblWVL = Inputs['WVL-ABS'][:,0]
        self.allAW = Inputs['WVL-ABS'][:,1]
        self.allAO = Inputs['WVL-ABS'][:,2]
        self.allAU = Inputs['WVL-ABS'][:,3]

    """
    ____________________________________________________________________
    Method: main light intensity calculation (CalcTotSolarEng)

    Inputs: 

            Inputs : a dictionary containing datetime information for light calculation
                 and information about previous flux calculations 

        This is the main function of the ported SOLPOS 2.0 code. This method calls all SOLPOS 2.0
        submethods, which calculate geometries, then calculates various parameters based on these 
        geometries. Finally, refraction, absorption, and absorption by airmass are calculated.

    """     
    def CalcTotSolarEng(self,Inputs):
        # Begin by checking for errors.
        self.ErrorCheck(Inputs)
        
        self.geometry()
        self.localtrig()
        self.zen_no_ref()
        self.sshas()
        self.sbcfs()
        self.tsts()
        self.srsss()
        self.sazms()
        self.refrac()
        self.amasss()
        self.primes()
        self.etr()

    """
    ____________________________________________________________________
    Method: checking inputs for errors (ErrorCheck)

    Inputs: 

            Inputs : a dictionary containing datetime information for light calculation
                 and information about previous flux calculations 

        This function makes sure that all input values are within reasonable ranges for SOLPOS 2.0
        and SMARTS2. For example, dates and times are forced to have realistic values, and be within 
        the timespan for which the algorithm was fit. This function will throw an error if you attempt
        to run CalcTotSolarEng with an incorrect input value.

    """    

    def ErrorCheck(self,Inputs):
        # Algorithm is bounded only between 1950 and 2050.
        if Inputs['Year'] < 1950 or Inputs['Year'] > 2050:
            raise ValueError('Year must be 1950<x<2050.') 
        # Must give a real month. Month = 0 allowed in case of flag.
        if Inputs['Month'] < 0 or Inputs['Month'] > 12:
            raise ValueError('Month must be 1<x<12.') 
        # Must give a real day. Day = 0 allowed in case of flag.
        # Additional checks must be performed for months with less days and leap years.
        if Inputs['Day'] < 0 or Inputs['Day'] > 31:
            raise ValueError('Day must be 1<x<31.')
        elif (Inputs['Month'] in [4,6,9,11]) and Inputs['Day'] > 30:
            raise ValueError('Day must be 1<x<30.')
        elif Inputs['Month'] == 2 and self.leap and Inputs['Day'] > 29:
            raise ValueError('Day must be 1<x<29.')
        elif Inputs['Month'] == 2 and not self.leap and Inputs['Day'] > 28:
            raise ValueError('Day must be 1<x<28.')
        # Must input a valid number of days. 
        # Additional check must be performed for leap years.
        if Inputs['DayNumber'] < 0 or (self.leap and Inputs['DayNumber'] > 366):
            raise ValueError('DayNumber must be 0<x<366 this year.')
        if Inputs['DayNumber'] < 0 or (not self.leap and Inputs['DayNumber'] > 365):
            raise ValueError('DayNumber must be 0<x<365 this year.')
        # Now, define the hour, minute, second bounds.
        if Inputs['Hour'] < 0 or Inputs['Hour'] > 23:
            raise ValueError('Hour must be 0<x<23.')
        if Inputs['Minute'] < 0 or Inputs['Minute'] > 59:
            raise ValueError('Minute must be 0<x<59.')
        if Inputs['Second'] < 0 or Inputs['Second'] > 59:
            raise ValueError('Second must be 0<x<59.')
        # Define limits on the time zone.
        if abs(Inputs['Timezone']) > 12:
            raise ValueError('Timezone must be -12<x<12.')
        # Define limits on the spatial coordinates.
        if abs(Inputs['Latitude']) > 90:
            raise ValueError('Latitude must be -90<x<90.')
        if abs(Inputs['Longitude']) > 180:
            raise ValueError('Longitude must be -180<x<180.')

    """
    ____________________________________________________________________
    Method: wavelength-dependent light intensity calculation (CalcSolarSpec)

    Inputs: 

            Inputs : a dictionary containing datetime information for light calculation
                 and information about previous flux calculations 

        This is the main function of the ported SMARTS2 code. This method first runs SOLPOS 2.0
        by called CalcTotSolarEng, then expands upon the geometries calculated by SOLPOS 2.0 to
        calculate wavelength-dependent transmission of light through the atmosphere.

        SMARTS2 considers scattering, absorption, and reflection in its model.
    """   

    def CalcSolarSpec(self,Inputs):
        # Must run previous program to acquire geometry and other relevant parameters.
        self.CalcTotSolarEng(Inputs)
        
        self.fortrangeom()
        self.forwardscat()
        self.ozonemass()
        self.waveloop()

    """
    ____________________________________________________________________
    Method: pure water absorbance calculation (waterabs)

    Inputs: 

            depth : distance below ocean surface (meters)

        This method builds upon the calculated incident spectrum at the ocean's surface by
        using wavelength-dependent water absorption coefficients to attenuate the light. The
        absorption coefficient data can be found here: https://doi.org/10.1364/AO.36.008710.
        Pure water absorbs more strongly at longer wavelengths, so the deeper one looks, the
        less red light is available.
    """  
    def waterabs(self,depth):
        import math
        import numpy as np
        # If a spectrum for light passing through the atmosphere does not exist, throw error.
        if not hasattr(self,'DTOT'):
            raise Exception('Global Irradiance Spectrum at sea level must be calculated first: .CalcSolarSpec(Inputs)')
    
        # (Pope; Fry, 1997) Data.
        WVLabscons = np.loadtxt('.\\Data\\PureWaterAbs.csv',delimiter=',',dtype='float64', encoding='utf-8-sig')
        
        absWVL = WVLabscons[:,0]
        abscons = WVLabscons[:,1]
        
        photoactmin = int(min(absWVL))
        photoactmax = int(max(absWVL))
        
        self.allWVLphoto = self.allWVL[200:568]
        self.DTOTphoto = self.DTOT[200:568]
        
        absorb = np.interp(self.allWVLphoto,absWVL,abscons)
        
        # Preallocate for underwater measurements.
        self.subDTOTphoto = np.zeros((len(self.DTOTphoto),1))
        
        for idx, Io in enumerate(self.DTOTphoto):
            self.subDTOTphoto[idx] = self.DTOTphoto[idx]*math.exp(-depth*absorb[idx])
            
        return self.subDTOTphoto


    """
    ___________________________________________________________________________________________
    Method: water absorbance calculation (SeaAP_Abs) including Attenuation by Particulates (AP)

    Inputs: 

            depth : distance below ocean surface (meters)
            X : current biomass concentration (mg dry weight/L)

        The author of the paper cited above also collected absorption coefficient data for impure
        sea water. They split this into two categories, particulates and dissolved carbon. Here, we
        introduce additional attenuation by particulates (organic biomass). We are able to make the
        that the particulates Pope et al. measure are mostly photosynthetic organisms because of the
        characteristic absorption peaks of chlorophyll present at 430 and 670 nm in Pope's dissertation:
        https://www.proquest.com/openview/62db766a95a8d4ec3276c2a4fd8b1184/1?pq-origsite=gscholar&cbl=18750&diss=y
    """     
       
    def SeaAP_Abs(self,depth,X):
        Xmax = 2.1 # Maximum biomass concentration at 29degN, where data is pulled from.
        
        import math
        import numpy as np
        
        # (Pope; Fry) Particulate (biomass) Data.
        WVLabscons = np.loadtxt('.\Data\SeaAP.csv',delimiter=',',dtype='float64', encoding='utf-8-sig')
        
        absWVL = WVLabscons[:,0]
        abscons = WVLabscons[:,1]
        
        photoactmin = int(min(absWVL))
        photoactmax = int(max(absWVL))
        
        absorb = np.interp(self.allWVLphoto,absWVL,abscons)
        
        for idx, Io in enumerate(self.DTOTphoto):
            self.subDTOTphoto[idx] = self.DTOTphoto[idx]*math.exp(-depth*absorb[idx]*(X/Xmax))
            
        return self.DTOTphoto
        
    def dom2doy(self,month,day):
        # Add the cumulative days at the start of the month to the current day of the month.
        self.daynum = day + self.month_days[month]
        
        if self.leap:
            self.daynum += 1
    
    def geometry(self):
        import math
        
        ### Day angle (fraction of 360 degrees traveled around the sun)
        self.dayang = 360.0 * ((self.daynum - 1) / 365.0)
        
        ### Earth radius vector (ERV). Typically, ERV*Solar constant = solar energy.
        # Adjusts for position in orbit, distance from sun.
        # The following variables are only needed in the geometry function, keep local.
        
        # Sine of the day angle. 
        sd = math.sin(self.raddeg*self.dayang)
        # Cosine of the day angle.
        cd = math.cos(self.raddeg*self.dayang)
        # Twice the day angle.
        d2 = 2.0*self.dayang
        c2 = math.cos(self.raddeg*d2)
        s2 = math.sin(self.raddeg*d2)
        
        # Calculate the Earth Radius Vector
        self.erv = 1.000110 + 0.034221*cd + 0.001280*sd
        # Extra correction factors.
        self.erv += 0.000719*c2 + 0.000077*s2
        
        ### Universal Coordinated Time (Greenwich Standard) in seconds
        self.utime = self.hour*3600 + self.minute*60 + self.second - self.interval/2
        # Convert back into hours and adjust for timezone.
        self.utime = (self.utime/3600) - self.timezone
        
        ### Julian Day beginning at year 1950.
        # Years since 1950:
        delta = self.year - 1949
        # Number of leap years since 1950. Subtracting 0.01 so that 1951 is not considered a leap year, for example.
        leap = round((delta/4)) 
        self.julday = 32916.5 + delta*365 + leap + self.daynum + self.utime/24
        
        ### Ecliptic coordinate "time" since year 2000.
        # Subtract 50 years worth of days.
        self.ectime = self.julday - 51545
        
        ### Mean longitude (ecliptic coordinate)
        self.mnlong = 280.460 + 0.9856474*self.ectime
        # Remove multiples of 360, until reaching a coterminal angle 0<x<360.
        self.mnlong -= 360*math.floor(self.mnlong/360)
        
        ### Mean anomaly (ecliptic coordinate)
        self.mnanom = 357.528 + 0.9856003*self.ectime
        # Remove multiples of 360, until reaching a coterminal angle 0<x<360.
        self.mnanom -= 360*math.floor(self.mnanom/360)
        
        ### Ecliptic longitude (ecliptic coordinate)
        self.eclong = self.mnlong + 1.915*math.sin(self.mnanom*self.raddeg) + 0.020*math.sin(2*self.mnanom*self.raddeg)
        # Remove multiples of 360, until reaching a coterminal angle 0<x<360.    
        self.eclong -= 360*math.floor(self.eclong/360)
        
        ### Obliquity of the Ecliptic (tilt of the Earth). This equation may have a SIGN ERROR!
        self.ecobli = 23.439 - 4e-7*self.ectime
        
        ### Declination
        self.declin = self.degrad * math.asin((math.sin(self.ecobli*self.raddeg)*math.sin(self.eclong*self.raddeg)))
        
        ### Right Ascension
        # This is a fraction, calculate it separately for the arctangent function.
        top = math.cos(self.raddeg*self.ecobli)*math.sin(self.raddeg*self.eclong)
        bottom = math.cos(self.raddeg*self.eclong)
        
        # Put it all together. "arctangent2" returns a value pi<x<-pi. "arctangent" returns -pi/2<x<pi/2.
        self.rascen = self.degrad*math.atan2(top,bottom)
        # Make the angle positive (0<x<360)
        if self.rascen < 0:
            self.rascen += 360
        
        ### Greenwich mean sidereal time
        self.gmst = 6.697375 + 0.0657098242*self.ectime + self.utime
        # Remove multiples of 24 hours to get within 0<x<24.
        self.gmst -= 24*math.floor(self.gmst/24)
        
        ### Local mean sidereal time
        self.lmst = 15*self.gmst + self.longitude
        # Remove multiples of 360 degrees for 0<x<360.
        self.lmst -= 360*math.floor(self.lmst/360)
        
        ### Hour angle
        self.hrang = self.lmst - self.rascen
        
        # Two checks to force the angle between -180<x<180
        if self.hrang < -180:
            self.hrang += 360
        elif self.hrang > 180:
            self.hrang -= 360
        
    ### Quick calculations of important trig values that will be used in successive functions.    
    def localtrig(self):
        import math
        # cosine of the declination
        self.cd = math.cos(self.raddeg*self.declin)
        # cosine of the hour angle
        self.ch = math.cos(self.raddeg*self.hrang)
        # cosine of latitutde
        self.cl = math.cos(self.raddeg*self.latitude)
        # sine of the declination
        self.sd = math.sin(self.raddeg*self.declin)
        # sin of latitude
        self.sl = math.sin(self.raddeg*self.latitude)
    
    ### Zenith calculation with no refraction.
    def zen_no_ref(self):
        import math
        self.cz = self.sd*self.sl + self.cd*self.cl*self.ch
        
        # Checks for rounding errors that would give unreasonable values.
        if abs(self.cz) > 1:
            if self.cz > 1:
                self.cz = 1
            else:
                self.cz = -1
        
        ### Extraterrestrial solar zenith angle:
        self.zenetr = math.acos(self.cz)*self.degrad
        
        # Limit degrees below horizon at night to avoid weirdness.
        if self.zenetr > 99:
            self.zenetr = 99
            
        ### Extraterrestrial solar elevation angle (how high is sun in sky?)
        self.elevetr = 90 - self.zenetr
        
    ### Sunset hour angle, where in the sky does the sun set?
    def sshas(self):
        import math
        cdcl = self.cd*self.cl
        
        if abs(cdcl) >= 0.001:
            cssha = (-self.sl*self.sd)/cdcl
            # Catch cases for if the sun is setting directly over the equator.
            if cssha < -1:
                self.ssha = 180
            elif cssha > 1:
                self.ssha = 0
            else:
                self.ssha = self.degrad*math.acos(cssha)
        # More catch cases. If cdcl is very small, sun sets very close to equator.    
        elif (self.declin >= 0 and self.latitude > 0) or (self.declin < 0 and self.latitude < 0):
            self.ssha = 180
        else:
            self.ssha = 0
        
    ### Shadowband correction factor    
    def sbcfs(self):
        import math
        p = (0.6366198*self.sbwid)/(self.sbrad*(self.cd**3))
        t1 = self.sl*self.sd*self.ssha*self.raddeg
        t2 = self.cl*self.cd*math.sin(self.ssha*self.raddeg)
        
        self.sbcf = self.sbsky + 1/(1-p*(t1+t2))
        
    ### True solar time, minutes from midnight
    def tsts(self):
        self.tst = (180+self.hrang)*4
        
        ### Separate out fixed true solar time to use in equation of time.
        self.tstfix = self.tst - self.hour*60 - self.minute - self.second/60 + self.interval/120
        
        # bind the fixed true solar time within a single day -720<x<720 minutes.
        while self.tstfix > 720:
            self.tstfix -= 1440
        while self.tstfix < -720:
            self.tstfix += 1440
            
        ### Equation of time.
        self.eqntim = self.tstfix + 60*self.timezone - 4*self.longitude
    
    ### Sunrise and Sunset times (minutes from midnight)
    def srsss(self):
        # Permanent night
        if self.ssha <= 1:
            self.sretr = 2999
            self.ssetr = -2999
        # Permanent daylight
        elif self.ssha >= 179:
            self.sretr = -2999
            self.ssetr = 2999
        else:
            self.sretr = 720 - 4*self.ssha - self.tstfix
            self.ssetr = 720 + 4*self.ssha - self.tstfix
            
    ### Solar azimuth angle
    def sazms(self):
        import math
        self.ce = math.cos(self.raddeg*self.elevetr)
        self.se = math.sin(self.raddeg*self.elevetr)
        
        self.azim = 180
        cecl = self.ce*self.cl
        
        if abs(cecl) >= 0.001:   
            ca = ((self.se*self.sl)-self.sd)/cecl
            if ca > 1.0:
                ca = 1.0
            elif ca < -1.0:
                ca = -1.0
                
            self.azim = 180.0 - (math.acos(ca)*self.degrad)
            
            if self.hrang > 0:
                self.azim = 360 - self.azim

    # Refraction of light based on previous angle calculations.        
    def refrac(self):
        import math
        # If the angle is sufficiently steep, assume no refraction.
        if self.elevetr > 85:
            self.refcor = 0    
        else:
            # Three different models depending on the steepness of the sun's elevation.
            tanelev = math.tan(self.raddeg*self.elevetr)
            if self.elevetr >= 5:
                self.refcor = (58.1/tanelev) - (0.07/(tanelev**3)) + (0.000086/(tanelev**5))
            elif self.elevetr >= -0.575:
                self.refcor = 1735 + self.elevetr*(-518.2 + self.elevetr*(103.4 + self.elevetr*(-12.79 + self.elevetr*0.711)))
            else:
                self.refcor = -20.774/tanelev
                
            
        self.prestemp = (self.press*283)/(1013*(273+self.temp))
        self.refcor *= self.prestemp/3600
        
        self.elevref = self.elevetr + self.refcor
        
        if self.elevref < 0:
            self.elevref = 0

        # Convert between zenith and elevation.    
        self.zenref = 90 - self.elevref
        self.coszen = math.cos(self.raddeg*self.zenref)

    # Calculation of airmass in the way of light intensity.    
    def amasss(self):
        import math
        if self.zenref > 93:
            self.amass = -1
            self.ampress = -1
            
        else:
            self.amass = 1/(math.cos(self.raddeg*self.zenref)+0.50572*((96.07995-self.zenref)**-1.6364))
            self.ampress = self.amass*(self.press/1013)

    # Normalize and unnormalize Kt values.    
    def primes(self):
        import math
        self.unprime = 1.031*math.exp(-1.4/(0.9+(9.4/self.amass))) + 0.1
        self.prime = 1/self.unprime

    # Extraterrestrial radiation given angle to sun.   
    def etr(self):
        if self.coszen > 0:
            self.etrn = self.solcon*self.erv
            self.etrang = self.etrn*self.coszen
            
        else:
            self.etrn = 0
            self.etrang = 0
    
    ### NEXT: Switching to ported SMARTS2 FORTRAN code. ###
    
    def fortrangeom(self):
        import math
        self.CI = math.cos(self.zenref*self.raddeg)
        self.zsin = math.sin(self.zenref*self.raddeg)


    def forwardscat(self):
        import math
        ALG = math.log(1-self.GG)
        AFS = ALG*(1.459+(ALG*(0.1595+(ALG*0.4129))))
        BFS = ALG*(0.0783+(ALG*(-0.3824-(ALG*0.5874))))
        self.FSP = 1 - (0.5*math.exp((AFS+(BFS/1.8))/1.8))
        self.FS = 1 - (0.5*math.exp((AFS+(BFS*self.coszen))*self.coszen))



    def ozonemass(self):
        self.AMO = (1+self.RR)/(((self.coszen**2)+(2*self.RR))**0.5)
        
    

    def waveloop(self):
        import math
        import numpy as np
        indexfinder = 0
        self.DTOT = np.zeros((len(self.allWVL),1))
        self.DIR = np.zeros((len(self.allWVL),1))
        self.DIF = np.zeros((len(self.allWVL),1))
    
        if self.zenref > 90:
            return
    
        for WVLOld in self.allWVL:
            # Change to nanometers.
            WVL = WVLOld/1000
            
            H0 = self.allETR[np.where(self.allWVL == WVLOld)[0]]
            
            AW = np.interp(WVLOld,self.tblWVL,self.allAW)
            AO = np.interp(WVLOld,self.tblWVL,self.allAO)
            AU = np.interp(WVLOld,self.tblWVL,self.allAU)
            
            H0 = self.erv*H0
            
            omegal = self.omega*math.exp(-self.omegap*(math.log((WVL/0.4))**2)) 
            
            # Wavelength dependent albedo, requires textbook parameters...
            
            ### Transmittance Section
            # 
            TR = math.exp(-self.ampress/((WVL**4)*(115.6406-(1.335/(WVL**2))))) 
            # Ozone        
            TO = math.exp(-AO*self.ozone*self.AMO)
            # Water vapor
            TW = math.exp((-0.2385*AW*self.watervap*self.amass)/((1+(20.07*AW*self.watervap*self.amass))**0.45))
            # Uniformly mixed gases
            TU = math.exp((-1.41*AU*self.ampress)/((1+118.93*AU*self.ampress)**0.45))
            
            # Airmass absorption
            DELA = self.turb*((WVL/0.5)**(-1.14))
            TAS = math.exp(-omegal*DELA*self.amass)
            TAA = math.exp(-(1-omegal)*DELA*self.amass)
            TA = math.exp(-DELA*self.amass)
            
            # Direct transmission of light, percentage
            DIR = H0*TR*TO*TW*TU*TA
            self.DIR[np.where(self.allWVL == WVLOld)[0]] = (DIR*self.coszen)
            
            # Correlations from SMARTS2 paper, see link above for link.
            DRAY = H0*self.coszen*TO*TW*TU*TAA*(1-(TR**0.95))*0.5
            DAER = H0*self.coszen*TO*TW*TU*TAA*(TR**1.5)*(1-TAS)*self.FS
            TRP = math.exp(-1.8/((WVLOld**4)*115.6406-(1.335/(WVLOld**2))))
            TWP = math.exp((-0.2385*AW*self.watervap*1.8)/(1+((20.07*AW*self.watervap*1.8)**0.45)))
            TUP = math.exp((-1.41*AU*1.8)/((1+((118.93*AU*1.8)**0.45))))
            TASP = math.exp(-omegal*DELA*1.8)
            TAAP = math.exp(-(1-omegal)*DELA*1.8)
            RHOA = TUP*TWP*TAAP*((0.5*(1-TRP))+((1-self.FSP)*TRP*(1-TASP)))
            DRGD = (((DIR*self.coszen)+(DRAY+DAER))*self.RHO*RHOA)/(1-(self.RHO*RHOA))
            DIF = DRAY+DAER+DRGD
            
            # Correction factor at low wavelengths.
            CRC = 1
            if WVL <= 0.45:
                CRC = (WVL+0.55)**1.8
            
            # Dividing by 341.8 to convert to mmol/(m^2*s).
            self.DIF[np.where(self.allWVL == WVLOld)[0]] = (DIF*CRC)           
            self.DTOT[np.where(self.allWVL == WVLOld)[0]] = (DIR*self.coszen+DIF)