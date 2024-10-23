import numpy as np
import matplotlib.pyplot as plt

def RH_Dew_point(t_db, t_wb, P, verbose = False):
    '''
    Function for calculating the relative humidity (%) and dew-point in (degrees Celsius)
    from the values of atmospheic pressure in (kPa) and dry-bulb and wet-bulb temperatures
    in degrees Celsius.
    
    # Example Inputs:
    t_wb = 20#10 # Celsius
    t_db = 30#15 # Celsius
    P = 101.325 # kPa

    Equations from:
        Environmental Engineering in South African Mines, The Mine Ventilation 
        Society of South Africa, 1989, pp 451-455. 
        (https://earthscience.stackexchange.com/questions/19151/how-is-relative-humidity-determined-from-a-wet-and-dry-bulb-readings)
    '''

    # The psychrometric equations 
    
    # Constants:
    M_w = 18.016
    M_a = 28.9664
    f = 1.0048
    
    Pp_ws = 0.6105*np.exp(17.27*t_wb/(237.3 + t_wb)) # kPa
    
    # 1. Calculate r_o, the moisture content of saturated air at the wet-bulb temperature.
    
    r_o = (M_w/M_a)*(f*Pp_ws/(P - f*Pp_ws)) # kg/kg
    
    # 2. Calculate the following enthalpy terms for dry air, liquid water and water vapour
    
    # Dry air
    H_ao = 1.005*t_wb # kJ/kg
    H_ai = 1.005*t_db # kJ/kg
    # Liquid water
    Hp_wl = 6.3e-6*t_wb**3 - 7.27e-4*t_wb**2 + 4.2058*t_wb + 0.03 # kJ/kg
    # Water vapour
    Hp_wo = -6.62e-6*t_wb**3 - 1.94e-4*t_wb**2 + 1.8375*t_wb + 2500.83 # kJ/kg
    Hp_wi = -6.62e-6*t_db**3 - 1.94e-4*t_db**2 + 1.8375*t_db + 2500.83 # kJ/kg
    
    # 3. Calculate r, the actual moisture content of the air:
    
    r = (r_o*(Hp_wo - Hp_wl) - (H_ai - H_ao))/(Hp_wi - Hp_wl) # kg/kg
    
    # 4. Calculate P_w, the vapour pressure:
    
    P_w = r*P/(f*((M_w/M_a) + r)) # kPa
    
    # 5. Calculate v, the specific volume
    R_Ma = (8.31436/M_a)*(1 - ((5.307e-6*P + 9.49e-6) - (8.115e-8*P + 2.794e-6)*t_db)) #kJ/kg.K 
    T = 273.15 + t_db # K
    v = R_Ma*T/(P - P_w) # m^3/kg    
    
    # 6. Calculate w the density
    
    w = (1 + r)/v # kg/m^3
    
    # 7. Calculate H the enthalpy
    
    H = H_ai + r*Hp_wi # kJ/kg
    
    # 8. Calculate S, the sigma heat
    
    S = H - r*Hp_wl # kJ/kg
    
    # 9. Calculate phi the relative humidity
    Ppp_ws = 0.6105*np.exp(17.27*t_db/(237.3 + t_db)) # kPa
    
    phi = 100*P_w/Ppp_ws
    
    # 10. Calculate t_dp the dew point temperature 
    x = np.log(P_w/0.6105)
    t_dp = 237.3*x/(17.27 - x) # Celsius

    if verbose:
        print('Absolute humidity %3.4f kg/kg'  % r)
        print('Relative humidity %3.2f percent' % phi)
        print('Dew point temperature %3.2f degrees Celsius' % t_dp)
    return phi, t_dp, r

def Wet_bulb_temp(t_db, RH, verbose = False):
    '''Estimate the wet-bulb temperature using Stull`s formula
    Accurate for 5% < RH < 99%
    Input:
        t_db = Dry-bulb temperature (degrees Celsius)
          RH = Relative humidity (%) 
    Source: https://www.omnicalculator.com/physics/wet-bulb

    Stull, R. (2011). Wet-bulb temperature from relative humidity and air 
    temperature. Journal of applied meteorology and climatology, 50(11), 
    2267-2269.
    '''
    t_wb = t_db*np.arctan(0.151977*np.sqrt(RH + 8.313659)) + \
            0.00391838*np.sqrt(RH**3)*np.arctan(0.023101*RH) \
            -np.arctan(RH - 1.676331) + np.arctan(t_db + RH) - 4.686035
    if verbose:
        print('Wet-bulb temperature %3.2f degrees Celsius' % t_wb)
    return t_wb

"""
Lawrence, Mark G., 2005: The relationship between relative humidity and the 
dewpoint temperature in moist air: A simple conversion and applications.
Bull. Amer. Meteor. Soc., 86, 225-233. 
doi: http;//dx.doi.org/10.1175/BAMS-86-2-225
"""

def magnus_formula(T):
    '''Magnus formla for saturation vapor pressure
    valid within the range -40~C < T <= 50~C for 
    dry-bulb temperatures in degrees Celsius.
    
    Output in hPa.
    '''
    A_1 = 17.625
    B_1 = 243.04 # C
    C_1 = 610.94 # Pa
    e_s = C_1*np.exp(A_1*T/(B_1 + T))
    return 0.01*e_s

def dew_point(T, RH):
    '''Dew-point temperature given the dry-bulb temperature T (Celsius) and
    RH in percentage'''
    A_1 = 17.625
    B_1 = 243.04 # C
    T_d = B_1*(np.log(RH/100) + A_1*T/(B_1 + T))/(A_1 - np.log(RH/100) - A_1*T/(B_1 + T))
    return T_d

############################################################

def dewpoint_linear(T,RH):
    '''Linear approximation valid for RH > 50% '''
    return T - ((100 - RH)/5)
