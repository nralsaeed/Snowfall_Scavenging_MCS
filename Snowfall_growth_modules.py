##Binning modules used for the MCS profile data:

import numpy as np
import pandas as pd
from scipy.stats import binned_statistic_2d as bin2d
from scipy.interpolate import pchip
import statsmodels.api as sm
from planetThermo import pco2
from planetThermo import tco2


def closest_nonzero(lst, start_index):
    nonzeros = [(i, x , start_index-i) for i, x in enumerate(lst) if x != 0]
    sorted_nonzeros = sorted(nonzeros, key=lambda x: abs(x[0] - start_index))
    return sorted_nonzeros[0][1]

#Altitude and Latitude binning (Zonal mean with altitude)
def Alt_Lat_binning(dfc,choice, latbinsize):
    '''
    bins the data by LS (per 1 LS) and then bins the temp, dust and H2Oice quantities
    in a 2D meshgrid of altitude (per km) and latitude(per 1 degrees).
    
    input: MCS dataframe that includes: TEMPERATURE, PRESSURE, DUST
    ,ALTITUDE, LATITUDE, LONGITUDE, H2OICE, LS, LOCALTIME
    
    output: 
    - altitude bins array
    - latitude bins array
    - a dictionary structure with each entry representing 1 LS containing a matrix of
    temperature (or dust or h2oice) quantities in each meshpoint
    '''

    #creating the base bins (won't change for the different LS's)
    
    #setting altitude bins
    #using overall max and min of all data read 
    #step size is set to 1km 
    A_bins = np.linspace(0,np.int(np.ceil(dfc.ALTITUDE.max())),np.int(np.ceil(dfc.ALTITUDE.max()))+1)
    
    #setting latitude bins
    #bin size determine by input
    Lmin = np.int(np.floor(dfc.LATITUDE.min()))
    Lmax = np.int(np.ceil(dfc.LATITUDE.max()))

    if latbinsize == 1: 
        L_bins = np.linspace(Lmin, Lmax,((Lmax - Lmin))+1)
    elif latbinsize == 2:
        L_bins = np.linspace(Lmin, Lmax,np.int(((Lmax - Lmin)/2))+1)
    else:
        print('error: choose bin size either 1 or 2')

    
    
    #Setting Ls bins
    #determining LS span (minimum and maximum)
    LSmin = np.int(round(dfc.LS.min()))
    LSmax = np.int(round(dfc.LS.max()))
    #How many Ls bins total: 
    #this will be used to create array of LSs
    LSsiz  = (LSmax - LSmin)
    #creating array of Ls with spacing = binsize
    LS_arr = np.linspace(LSmin, LSmax, LSsiz+1)
    LS_arr[LS_arr!= LSmax] #removing extra entry from entry
    
    #empty dictionary for storing purposes
    results = dict() 
    
    #looping through the different LS bins to and performing the 2d binning with pressure and lat for each LS
    for i in range(int(LSsiz)):
        #pick out only entries with that specific Ls
        dfc1 = dfc[(dfc.LS >= LS_arr[i]) & (dfc.LS < LS_arr[i] + 1)]
        
        #ignore if empty
        if  dfc1.size == 0:
            continue
        
        #binning temperature
        if choice == 'temp':
            temp = bin2d(dfc1.LATITUDE, dfc1.ALTITUDE, dfc1.TEMPERATURE, 
                         bins = [L_bins,A_bins], statistic = 'mean')[0]
            ans = temp.T
            
        #binning dust
        elif choice == 'dust':
            dust = bin2d(dfc1.LATITUDE, dfc1.ALTITUDE, dfc1.DUST,
                         bins = [L_bins,A_bins], statistic = 'mean')[0]
            ans = dust.T
            
        #binning H2O ice
        elif choice == 'h2oice':
            h2oice = bin2d(dfc1.LATITUDE, dfc1.ALTITUDE, dfc1.H2OICE,
                                        bins = [L_bins,A_bins], statistic = 'mean')[0]
            ans = h2oice.T

        #binning Pressure:
        elif choice == 'pressure':
            press = bin2d(dfc1.LATITUDE, dfc1.ALTITUDE, dfc1.PRESSURE,
                                        bins = [L_bins,A_bins], statistic = 'mean')[0]
            ans = press.T
            
        #error message if no choice was made
        else: 
            print('invalid variable input please choose from the following: temp, dust, h2oice ')
           
        #save current Ls binned data to dictionary    
        results[str(LS_arr[i])] = ans
    return L_bins, A_bins, results



    #Pressure and Latitude binning (Zonal mean with pressure)
def Pres_Lat_binning(dfc,choice):
    '''
    bins the data by LS (per 1 LS) and then bins the temp, dust and H2Oice quantities
    in a 2D meshgrid of pressure (logged) and latitude(per 2 degrees).
    
    input: MCS dataframe that includes: TEMPERATURE, PRESSURE, DUST
    ,ALTITUDE, LATITUDE, LONGITUDE, H2OICE, LS, LOCALTIME
    
    output: 
    - pressure bins array
    - latitude bins array
    - a dictionary structure with each entry representing 1 LS containing a matrix of
    temperature (or dust or h2oice) quantities in each meshpoint
    '''

    #creating the base bins (won't change for the different LS's)
    
    #setting pressure bins
    #using overall max and min of all data read  
    P_bins = np.logspace(np.log10(dfc.PRESSURE.min()), np.log10(dfc.PRESSURE.max()), base=10, num = 100)
    
    #setting latitude bins
    #looking at 2 degree bins for now 
    Lmin = np.floor(dfc.LATITUDE.min())
    Lmax = np.ceil(dfc.LATITUDE.max())
    L_bins = np.linspace(Lmin, Lmax,((Lmax - Lmin)/2)+1)
    
    
    #Setting Ls bins
    #determining LS span (minimum and maximum)
    LSmin = round(dfc.LS.min())
    LSmax = round(dfc.LS.max())
    #How many Ls bins total: 
    #this will be used to create array of LSs
    LSsiz  = (LSmax - LSmin)
    #creating array of Ls with spacing = binsize
    LS_arr = np.linspace(LSmin, LSmax, LSsiz+1)
    LS_arr[LS_arr!= LSmax] #removing extra entry from entry
    
    #empty dictionary for storing purposes
    results = dict() 
    
    #looping through the different LS bins to and performing the 2d binning with pressure and lat for each LS
    for i in range(int(LSsiz)):
        #pick out only entries with that specific Ls
        dfc1 = dfc[(dfc.LS >= LS_arr[i]) & (dfc.LS < LS_arr[i] + 1)]
        
        #ignore if empty
        if  dfc1.size == 0:
            continue
        
        #binning temperature
        if choice == 'temp':
            temp = bin2d(dfc1.LATITUDE, dfc1.PRESSURE, dfc1.TEMPERATURE, 
                         bins = [L_bins,P_bins], statistic = 'mean')[0]
            ans = temp.T
            
        #binning dust
        elif choice == 'dust':
            dust = bin2d(dfc1.LATITUDE, dfc1.PRESSURE, dfc1.DUST,
                         bins = [L_bins,P_bins], statistic = 'mean')[0]
            ans = dust.T
            
        #binning H2O ice
        elif choice == 'h2oice':
            h2oice = bin2d(dfc1.LATITUDE, dfc1.PRESSURE, dfc1.H2OICE,
                                        bins = [L_bins,P_bins], statistic = 'mean')[0]
            ans = h2oice.T
            
        #error message if no choice was made
        else: 
            print('invalid variable input please choose from the following: temp, dust, h2oice ')
            
        
        #save current Ls binned data to dictionary    
        results[str(LS_arr[i])] = ans
        
    return L_bins, P_bins, results

def plot_limits(temp):
    ''' define the limits for the colorbars for the lat alt plots'''
    
    #first create the array of names to input as the list key

    names = list(temp.keys())

    #now create the empty arrays to store the minimum temp/dust/ice recorded in each frame (each time unit)
    tmin = np.zeros(len(temp))
    tmax = np.zeros(len(temp))

    #actually find the minimum temp/dust/ice recorded in each frame and store it 
    for i in range(len(temp)):
        t = temp[names[i]]
        tmin[i] = np.nanmin(t)
        tmax[i] = np.nanmax(t)

    #find the minimums and maximums of all the temp/dust/ice values for the entire movie 
    tmins = np.nanmin(tmin)
    tmaxs = np.nanmax(tmax)

    return t,tmins, tmaxs

def Planck(T, waveno):
    '''
    Calculates the blackbody radiation of an object of temperature T at wavenumber waveno
    inputs: 
    T: temperature of blackbody
    waveno: wavenumber
    output:
    Planck's spectral radiance in W/m2·sr−1·cm−1 
    '''
    # plancks constant
    alpha1 = 1.191042e-8 #W/m2 sr cm-1
    alpha2 = 1.4387752 #K cm
    B = (alpha1*(waveno**3))/(np.exp((alpha2*waveno)/(T))-1)
    return B

def Median_Ls(Opc,Press, ptcl_size,dz): 
    '''Finds the median Pressure(z) and Number_density(z) in a specific Ls by 
    taking the values for all lats and calculating a median'''
    
    #define csts
    rhoc = 1600 #kg/m3 #solid density of CO2 ice (Find reference)
    Kb = 1.38064852e-23 #m2 kg/ s2 K #boltzmann constant
    dc = 3.3e-10 #m CO2 molecule diameter
    Tav = 140 #K #average temp at the poles
    mu = 8.5e-6 #(kg/ m s) OR (Pa s) # molecular viscosity of atmosphere at 150 K
    gmars = 3.711 #m/s² mars gravity
    
    fctr = np.int(1e3/dz)

 
    Exp_density = np.zeros(shape =(fctr*Opc.shape[0]+1,Opc.shape[1]))
    Exp_press = np.zeros(shape =(fctr*Press.shape[0]+1,Press.shape[1]))
    for j in range(Press.shape[1]):
        opc = Opc[:,j]
        press = Press[:,j]
        if len(np.where(press == 0)[0]) > 0:
            altmax = np.where(press == 0)[0][0] - 1
        else: 
            altmax = len(press) - 1
        if opc[:altmax].any() == True:
            opc_smooth, z_array = Z_interpolater(opc,dz,altmax)
            press_smooth, alts = Z_interpolater(press,dz,altmax)
        
            lat_rho,vels,num_den_exp = find_Nd_vels(opc_smooth,press_smooth,ptcl_size,dz)
        
            indexor = len(num_den_exp)
            Exp_density[:indexor,j] = num_den_exp
            indexor2 = len(press_smooth)
            Exp_press[:indexor2,j] = press_smooth
        else:
            continue
            
    #Now, take the median of all the values along the latitudes
    #mask all the rows where values = 0
    msk_press = np.ma.masked_equal(Exp_press, 0)
    msk_den = np.ma.masked_equal(Exp_density, 0)

    #calculate median
    M_den = np.ma.median(msk_den, axis=1).filled(0)
    M_press  = np.ma.median(msk_press, axis=1).filled(0)
    M_frost  = tco2(0.95*M_press)
    
    #index of highest alt with pressure reading:
    if len(np.where(M_press == 0)[0]) > 0:
        altmax2 = np.where(M_press == 0)[0][0] - 1
    else: 
        altmax2 = len(M_press) - 1
    
    M_den = M_den[:altmax2]
    M_press = M_press[:altmax2]
    M_frost = M_frost[:altmax2]
    
    #empty array to store the settling velocity w
    M_vels = np.zeros(len(M_den))
    for z in range(len(M_den)):
        mfp = (Kb*Tav)/(np.sqrt(2)*np.pi*(dc**2)*M_press[z]) #mean free path in air
        M_vels[z] = ((2*rhoc*gmars*(ptcl_size**2))/(9*mu))*(1+ (mfp/ptcl_size)) #m/s
    
        
    return M_den, M_press, M_frost, M_vels


def Z_interpolater(Array,dz,altmax):
    ''' function to smoothen the temp,press,variable(z) arrays in altitude. Essentially to decrease the altitude step between layers
    input: Array(z)
    dz = step size
    altmax = maximum altitude with pressure reading
    '''
    array = Array.copy()
    
    array = array[:altmax]
    
    #the factor by which you are multiplying the altitude seperation (assuming initial spacing is 1000 meters)
    fctr = np.int(1e3/dz)
    
    #altitude limit
    Z = int(len(array))*1e3 #m
    #original altitude separation
    Idz = 1000 #m

    Nz = int(round(Z/float(Idz)))
    z = np.linspace(0, Nz*Idz, Nz+1) # Mesh points in space

    #smoothen the plots (because as is, it is too jumpy and causes issues with the finite diff method)
    #plus this way I can have a smaller spacial step in altitude than 1km
    #set up the smoothing with Piecewise Cubic Hermite Interpolating Polynomial
    pch = pchip(z[:-1],array)

    #create more finely space array in z axis (has to be at least 2*Nz)
    xs = np.linspace(0, Z, (fctr*Nz)+1)
    xs2 = np.append(xs, Z+(Idz/dz))
    z_array = xs2


    #interpolate for both no density u and settling velocity w
    array_smooth = pch(xs2)

    #make sure there are no negative values accidentally reached from the interpolation:
    array_smooth[array_smooth<0]=0
    
    return array_smooth, z_array


def IR_cooling_radii_vels(NumDen,Tfrost,Press, radii, vels,dz,dt,T,ptcl_size):
    ''' Function to caclulate the radius distribution of CO2 ice
    particles in the atmosphere r(z) to be used as an input in the cloud model
    We do this by taking the case of a single particle, and following it to the surface
    calculating its growth and new size at each altitude step
    Inputs: 
    Outputs:
    r(z)
    w(z)
    T total time 
    '''
    
    Nt = int(round(T/float(dt)))
    #setting up some constants for later:
    ac = ptcl_size #m starting radius of CO2 ice particles (10 microns)
    WN = 667.3 #(wavenumber center frequency of 15 micron band)  
    mubar = 0.488 #see Pollack et al 1981 for this and the following
    do = 54.4 
    To = 300
    qd = 0.879
    fo = 0.153
    qf = -0.256
    aa = 0.323
    bb = 0.566
    nomw = 2000 # m-1 nominal width
    gmars = 3.711 #m/s² mars gravity
    Cp = 700 #J/kg K specific heat of CO2 atmosphere.. an assumption (ask Paul?)
    Latent = 5.9e5 #J/kg latent heat of CO2 sublimation (hayne 2014)
    Qc = 3.0249 #efficiency of extinction for CO2 
    rhoc = 1600 #kg/m3 #solid density of CO2 ice (Find reference)
    Kb = 1.38064852e-23 #m2 kg/ s2 K #boltzmann constant
    dc = 3.3e-10 #m CO2 molecule diameter
    Tav = 140 #K #average temp at the poles
    mu = 8.5e-6 #(kg/ m s) OR (Pa s) # molecular viscosity of atmosphere at 150 K
    zlength = len(Press)-1
 

    #First lets define the lookup matrix for the flux equivalent widths 
    #that we will use in the growth calculation
    E_matrix = np.zeros((zlength+1,zlength+1))
    for i in range(zlength):
        for j in range(zlength): 
            if i < j:
                Pbar = sum(Press[k] for k in range(i,j+1))/len(range(i,j+1))
                Tbar = sum(Tfrost[k]*Press[k] for k in range(i,j+1))/sum(Press[k] for k in range(i,j+1))
            elif i > j:
                Pbar = sum(Press[k] for k in range(j,i+1))/len(range(j,i+1))
                Tbar = sum(Tfrost[k]*Press[k] for k in range(j,i+1))/sum(Press[k] for k in range(j,i+1))
            elif i == j: 
                Pbar = Press[i]
                Tbar = Tfrost[i]
            Pbarml = 0.01*Pbar #convert av press to milibars for calculation
            apl = 1305.5*0.01*np.abs(Press[i] - Press[j]) #absorber path length in atmos cm
            E_matrix[i,j] = mubar*do*((Tbar/To)**qd)*np.log(1+(fo*((Tbar/To)**qf)*(apl**bb)*(Pbarml**aa)))
            
    #calculation of IR flux change for a single layer i of the atmosphere: 
    #note: the gas flux Fg, doesn't need to be recalculated every time step.. it will be the same throughout
    #the aerosol flux will change as the number density changes so we need to include those calculations inside the time loop
    
    ## IR flux from gas:
    
    #define surface temp to be used in gas flux calcs
    Ts = Tfrost[0]
    #Create arrays to store the IR gas flux from each layer
    Fg_down = np.zeros(zlength+1)
    Fg_up = np.zeros(zlength+1)
    
    #now we calculate the IR gas flux coming off from each layer, in both the upward and downward direction
    for i in range(1,zlength+1):
        #from gas cooling:
        Fg_down[i] = 2*np.pi*sum((E_matrix[i,j-1]-E_matrix[i,j])*Planck(Tfrost[j-1],WN) for j in range(1,i+1))
        Fg_up[i] = 2*np.pi*Planck(Ts,WN)*((0.5*nomw) - E_matrix[i,zlength]) \
                + 2*np.pi*sum((E_matrix[i,j+1]-E_matrix[i,j])*Planck(Tfrost[j+1],WN) for j in range(i+1,zlength))
    
    #boundary conditions:
    Fg_down[0] = 0 #nothing penetrating the surface
    Fg_up[zlength] = np.pi*Planck(Ts,WN)*nomw #total leaving the top of the atmosphere is whatever heat came from the surface to begin with
    
    #create empty matrix to store radius at each altitude 
    r_update = np.zeros(zlength+1)
    #create empty array to store settling velocity at each altitude
    vels_update = np.zeros(zlength+1)
    
    #set starting radius 
    rold = ac
    
    #create empty matrix to store new_radius and new vels at each altitude for every time_step
    r_evolution = np.zeros(shape =(zlength+1,Nt+1))
    v_evolution = np.zeros(shape =(zlength+1,Nt+1))
    
    r_evolution[:,0] = radii
    v_evolution[:,0] = vels
    #start timer:
    n=0
    #start counter
    j = 0
    
    #track new height reached by particle:
    new_h = dt*vels

    while n < T:
        
        
        for idx in range(zlength-1):
            
            #if dt <= dz/vels[idx]:
            if new_h[idx] <= dz: #if particle is still in that layer, compute new particle size
                #cooling from gas W/m2
                delF = ((np.abs(Fg_up[idx]) - np.abs(Fg_down[idx])) - (np.abs(Fg_up[idx+1]) - np.abs(Fg_down[idx+1])))
                massgrowth = np.abs(delF*dt/Latent) #where Latent is the latent heat of CO2 condesation
                if NumDen[idx] < 500:
                    particlegrowth = (massgrowth)/(800*dz) #mass growth per particle
                else:
                    particlegrowth = (massgrowth)/(NumDen[idx]*dz) #mass growth per particle

                rold = radii[idx]
                rnew = ((3*particlegrowth/(4*np.pi*rhoc))+(rold**3))**(1/3)
                r_update[idx] = rnew

                #calculate the settling velocity associated with this new particle size
                mfp = (Kb*Tav)/(np.sqrt(2)*np.pi*(dc**2)*Press[idx]) #mean free path in air
                vels_update[idx] = ((2*rhoc*gmars*(rnew**2))/(9*mu))*(1+ (mfp/rnew)) #m/s
                
            else: #if particle has moved to layer below, calculate the new radius reached but now 
                #using info of layer below + average it with the radius already reached in that layer
                #cooling from gas W/m2
                if idx > 0: #as long as its not the surface layer do this
                    nidx = idx-1
                    delF = ((np.abs(Fg_up[nidx]) - np.abs(Fg_down[nidx])) - (np.abs(Fg_up[nidx+1]) - np.abs(Fg_down[nidx+1])))
                    massgrowth = np.abs(delF*dt/Latent) #where Latent is the latent heat of CO2 condesation
                    if NumDen[nidx] < 100:
                        particlegrowth = (massgrowth)/(100*dz) #mass growth per particle
                    else:
                        particlegrowth = (massgrowth)/(NumDen[nidx]*dz) #mass growth per particle

                    rold = radii[idx]
                    rnew = ((((3*particlegrowth/(4*np.pi*rhoc))+(rold**3))**(1/3)) + r_update[nidx])/2
                    r_update[nidx] = rnew

                    #calculate the settling velocity associated with this new particle size
                    mfp = (Kb*Tav)/(np.sqrt(2)*np.pi*(dc**2)*Press[nidx]) #mean free path in air
                    vels_update[nidx] = ((2*rhoc*gmars*(rnew**2))/(9*mu))*(1+ (mfp/rnew)) #m/s
                    
                    #reset the start height of this particle now that it is in a new layer
                    new_h[idx] = 0
                    #print('error, timestep too large, particle has moved to layer below')
                else:
                    #print('error, surface reached')
                    continue
        r_update[0] =  r_update[1]
        vels_update[0] = vels_update[1]
        
        #for the lowest section of the atmosphere (<10 km), any small fluctuations can grow to be really large
        #so we replace any outliers with the median value of radius, and recalculate the velocity.
        
        no_outliers_r = remove_outliers_radius(r_update[:10])
        no_outliers_v = remove_outliers_radius(vels_update[:10])
        r_update[:10] = no_outliers_r
        vels_update[:10] = no_outliers_v
        
        #store the new radius and vels reached after dt time elapsed
        r_evolution[:,j+1] = r_update
        v_evolution[:,j+1] = vels_update
        
        
        #update radii and vels:
        radii = r_update
        vels = vels_update
        
        #update height of ptcl moving down 
        new_h =  new_h + dt*vels
        #update_counter
        j = j+1
        #update time
        n = n+dt
        
    
    
    return r_evolution, v_evolution



def find_Nd_vels(Opacity, Press,ptcl_size,dz):
    '''Function that calculates the number density + settling velocity of particles in the atmosphere for a single column 
    (i.e. in this case, for one latitude bin in one Ls) given a specific particle size distribution with altitude.
    inputs: Opacity vs alt array, Press vs alt array, particle radius vs alt.
    outputs: Number density vs alt array (m-3), settling velocities vs alt array (m/s)'''
    #define constants:
    h = dz #m # vertical height of 1 pixel
    Qc = 3.0249 #efficiency of extinction for CO2
    rhoc = 1600 #kg/m3 #solid density of CO2 ice (Find reference)
    R_m = 3.3895e6 #radius of Mars in m
    h = 1000 #m # vertical height of 1 pixel
    mu = 8.5e-6 #(kg/ m s) OR (Pa s) # molecular viscosity of atmosphere at 150 K
    G = 6.67408e-11 #m3/kg s  # gravitational constant
    Kb = 1.38064852e-23 #m2 kg/ s2 K #boltzmann constant
    dc = 3.3e-10 #m CO2 molecule diameter
    M_m = 6.39e23 #kg # mass of Mars in Kgrams
    gr = 3.711 #m/s²
    Tav = 140 #K #average temp at the poles
    fctr = np.int(1e3/dz)
    C = Opacity.copy() #m^-1 
    #step 1: convert the co2 ice opacities to optical depth
    C = C*h
    
    #step 2: calculate the no. density of CO2 ice particles in each pixel
    Rho =  C/(Qc*np.pi*(ptcl_size**2)*h)
    #Next: extrapolate to the surface

    #calculate the scale height of the atmosphere:
    Scale_H = (Kb*Tav)/(44*1.67e-27*gr)
    
    #make a copy of the array to create a new one with filled in gap
    NRho = Rho.copy()
    NRhoexp = NRho.copy()
    
    #altitude limit
    Z = int(len(NRho))*1e3 #m
    Idz = 1000 #m
    Nz = int(round(Z/float(Idz)))
    ##create array
    xs = np.linspace(0, Z, (fctr*Nz)+1)
    z_array = np.append(xs, Z+(Idz/dz))
    
    nonnan = np.argwhere(~np.isnan(NRho)&~(NRho==0))
    #don't do anything for arrays with all nan entries
    if nonnan.size == 0:
        pass
    #for arrays with non nan entries, extrapolate to the bottom using the bottom most non nan entry
    else:
        bottomindex = np.int(np.min(nonnan))
        #take the av of the bottom two points for the extrapolation
        if NRho[bottomindex+1] > 0: 
            AvNRho = (NRho[bottomindex] + NRho[bottomindex+1])/2
            #set the bottom most entry to be the av
            NRho[bottomindex] = AvNRho
        else:
            AvNRho = NRho[bottomindex]
            
        #using n = n_o exp(-z/H) find the no. density at each altitude below the bottom most entry
        for z in range(bottomindex):
            #using n = n_o exp(-z/H) find the no. density at each altitude below the bottom most entry
                
            N_d = (AvNRho/(np.exp(-z_array[bottomindex]/Scale_H)))*(np.exp(-z_array[z]/Scale_H)) 
            NRho[z] = N_d 
            NRhoexp[z] = N_d 
        #calculate the exponential shape of the no density to be used in the radius calculation
        for z in range(bottomindex, len(C)):
            N_d = (AvNRho/(np.exp(-z_array[bottomindex]/Scale_H)))*(np.exp(-z_array[z]/Scale_H)) 
            NRhoexp[z] = N_d
    Num_den = NRho
    Num_den_exp = NRhoexp
    #lastly calculate the settling velocity 
    
    #empty array to store the settling velocity w
    vels = np.zeros(len(C))
    for z in range(len(C)):
        mfp = (Kb*Tav)/(np.sqrt(2)*np.pi*(dc**2)*Press[z]) #mean free path in air
        vels[z] = ((2*rhoc*gr*(ptcl_size**2))/(9*mu))*(1+ (mfp/ptcl_size)) #m/s
    
    
    
    return Num_den, vels, Num_den_exp


    
def Cloud_Growth_CN(I, Pressure, R_evolution,V_evolution, T, dz, dt, K):
    """Solve dudt =k*d2u_dx2 + w*dudz + n*dwdz on (0,Z)x(0,T]
    using Crank Nicolson finite difference scheme
    where u is the vertical diffusion of a snow cloud
    And also incorporating a distribution of starting particle sizes
    with each time step
    I, Temp, Press, and w are all unique to each profile
    I = n(z) initial number density profile at t=0
    Temp = T(z) temperature profile
    Press = P(z) pressure profile
    w = w(z) settling velocity
    T = time boundary in seconds
    dz = altitude bin (stepsize) in meters
    dt = time bin (stepsize) in seconds
    K = diffusion coefficient
    """
    
    Nt = int(round(T/float(dt))) 
    t = np.linspace(0, Nt*dt, Nt+1) # Mesh points in time
    
    #the factor by which you are multiplying the altitude seperation
    fctr = np.int(1e3/dz)
    
    #index of highest alt with pressure reading:
    if len(np.where(Pressure == 0)[0]) > 0:
        altmax = (np.where(Pressure == 0)[0][0])-1 #(Multiply by 1000/dz because the pressure array isnt already interpolated to smaller dz)
    else: 
        altmax = len(Pressure)-1
    
    #You must now truncate the input array so that it ends at the maximum altitude:
    I_smooth = I[:altmax]
    old_radii = R_evolution[:altmax,0]
    ws = V_evolution[:altmax,0]

    zs = (np.linspace(0, len(I_smooth)*1e3, len(I_smooth)+1))/fctr

    #create empty arrays to store
    u = np.zeros((len(I_smooth))) # unknown u at new time level
    u_n = np.zeros((len(I_smooth))) # u at the previous time level
    
    #create empty matrix to store no. den flux at each altitude for every point in time
    u_evolution = np.zeros(shape =((len(I_smooth)),Nt+1))
    #create empty array to store CO2 being removed at each time step at the surface
    sediment_evolution = np.zeros(Nt+1)

    # Set initial condition u(z,0) = I(z)
    for f in range(0, len(I_smooth)-1):
        u_n[f] = I_smooth[f]
        
    u_evolution[:,0] = u_n
    
    
    #set up the tridiagonal matrix for Crank Nicholson
    #first extract the LHS coefficients:
    A = np.zeros(len(I_smooth)+1)
    B = np.zeros(len(I_smooth)+1)
    C = np.zeros(len(I_smooth)+1)
    E = np.zeros(len(I_smooth)+1)
    F = np.zeros(len(I_smooth)+1)
    G = np.zeros(len(I_smooth)+1)
    
    for i in range(0, len(I_smooth)-1):
        # extract coefficients of Crank Nicholson matrix
        WW = ws[i+1] - ws[i]
        A[i] = (-2*dt*K) - (ws[i]*dt*dz)
        B[i] = (4*(dz**2)) + (4*dt*K)
        C[i] = (-2*dt*K) + (ws[i]*dt*dz)
        E[i] = (2*dt*K) + (ws[i]*dt*dz) 
        F[i] = (4*(dz**2)) - (4*dt*K) + (2*dt*dz*WW)
        G[i] = (2*dt*K) - (ws[i]*dt*dz)
        
    ##set boundary conditions at surface, diffusion terms go away for coefficients 
    #of surface (z = 0) + coefficients of z=-1 are equal to 0 since it doesnt exist: 
    
    A[0] = 0
    B[0] = (4*(dz**2))
    C[0] = (-2*dt*K) + (ws[0]*dt*dz)
    E[0] =  0
    WW0 =  ws[1] - ws[0]
    F[0] = (4*(dz**2)) + (2*dt*dz*WW0)
    G[0] = (2*dt*K) - (ws[0]*dt*dz)
    
    
    # create tridiagonal matrix for Crank Nicholson coefficients
    #left hand side
    M = np.diagflat([C[i] for i in range(len(I_smooth)-1)], -1) +\
            np.diagflat([B[i] for i in range(len(I_smooth))]) +\
                np.diagflat([A[i] for i in range(len(I_smooth)-1)], 1)
    #right hand side
    L = np.diagflat([G[i] for i in range(len(I_smooth)-1)], -1) +\
            np.diagflat([F[i] for i in range(len(I_smooth))]) +\
                np.diagflat([E[i] for i in range(len(I_smooth)-1)], 1)
    
    #initialize the matrix:
    Uinit = np.zeros((len(I_smooth)))
    Uinit = u_n
    
    
    #start time
    n = 0
    #start counter
    j = 1
    
    #start the time loop to advance the number density:
    while n < T:

        # Compute u at inner mesh points
        # by solving the Crank Nicholson matrix:
        Usol = np.linalg.solve(M, L.dot(Uinit))
        
        #boundary condition at zmax, u = 0
        Usol[len(I_smooth)-1] = 0
        
        
        #also at the surface the flux/removal of snow = velocity*density, so we set that in the boundary conditions as well:
        #boundary condition at surface:
        #No diffusion at the surface (second derivative goes to 0)
        ww1 = (ws[0]*(Uinit[1] - Uinit[0]))/(dz)
        uu1 = (Uinit[0]*(ws[1] - ws[0]))/(dz)
        
        Usol[0] = Uinit[0] + (uu1 + ww1)*dt 
        #Usol[0] = Uinit[0] - (ws[0]*Uinit[0]*dt)/dz
        
        #update solution 
        Uinit = Usol
        
        #saving the column no density at each timestep
        u_evolution[:,j] = Uinit
        
        
        #store sedimentation rate:
        sediment = (ws[0]*Uinit[0])
        #saving the sedimentation at each time step (unit = m-2 s-1)
        sediment_evolution[j] = sediment
        
        new_radii = R_evolution[:altmax,j]
        new_ws = V_evolution[:altmax,j]
    
        #update radii and velocities:
        old_radii = new_radii
        ws = new_ws
        
        #recalculate the coefficients for the crank nicholson matrix:
        for i in range(0, len(I_smooth)-1):
            # extract coefficients of Crank Nicholson matrix
            WW = ws[i+1] - ws[i]
            A[i] = (-2*dt*K) - (ws[i]*dt*dz)
            B[i] = (4*(dz**2)) + (4*dt*K)
            C[i] = (-2*dt*K) + (ws[i]*dt*dz)
            E[i] = (2*dt*K) + (ws[i]*dt*dz) 
            F[i] = (4*(dz**2)) - (4*dt*K) + (2*dt*dz*WW)
            G[i] = (2*dt*K) - (ws[i]*dt*dz)

        ##update boundary conditions at surface, diffusion terms go away for coefficients 
        #of surface (z = 0) + coefficients of z=-1 are equal to 0 since it doesnt exist: 

        A[0] = 0
        B[0] = (4*(dz**2))
        C[0] = (-2*dt*K) + (ws[0]*dt*dz)
        E[0] =  0
        WW0 =  ws[1] - ws[0]
        F[0] = (4*(dz**2)) + (2*dt*dz*WW0)
        G[0] = (2*dt*K) - (ws[0]*dt*dz)


        # create tridiagonal matrix for Crank Nicholson coefficients
        #left hand side
        M = np.diagflat([C[i] for i in range(len(I_smooth)-1)], -1) +\
                np.diagflat([B[i] for i in range(len(I_smooth))]) +\
                    np.diagflat([A[i] for i in range(len(I_smooth)-1)], 1)
        #right hand side
        L = np.diagflat([G[i] for i in range(len(I_smooth)-1)], -1) +\
                np.diagflat([F[i] for i in range(len(I_smooth))]) +\
                    np.diagflat([E[i] for i in range(len(I_smooth)-1)], 1)

        #update counter
        j = j+1
        #update time
        n = n + dt
        
    surface_radius = new_radii[0]
        
        
    return Uinit, zs, t, u_evolution,sediment_evolution, surface_radius # Uinit holds latest u



def All_lat_snowfall(Opc,Press,Tfrost,T,dz,dt,K,LatV,ptcl_size):

    fctr = np.int(1e3/dz)

    #create matrix to store updated no density
    #the dimensions of this matrix depends on the chosen spacial step dz
    New_density = np.zeros(shape =(fctr*Opc.shape[0],Opc.shape[1])) 
    #create 3D matrix to store time evolution of number density (1 dimension for. space (altitude), 1d for no. density, and 1d for time)
    N_time_evolution = np.zeros(shape =(fctr*Opc.shape[0],np.int(T/dt)+1,Opc.shape[1])) 
    #first index is altitude, second is time, third is latitude

    #create matrix to store time evolution of the sedimentation at surface for each latitude:
    Sed_evol_lat = np.zeros(shape =(np.int(T/dt)+1,Opc.shape[1])) #first index is time, second index is latitude
    
    M_den, M_press, M_frost, M_vels = Median_Ls(Opc,Press,ptcl_size,dz)

    start_radii = ptcl_size*np.ones(len(M_den))

    Radius_evol, Vel_evol = IR_cooling_radii_vels(M_den,M_frost,M_press, start_radii, M_vels,dz,dt,T,ptcl_size)

    for j in range(len(LatV)-1): 
        opc = Opc[:,j]
        press = Press[:,j]
        tfrost = Tfrost[:,j]
        if len(np.where(press == 0)[0]) > 0:
            altmax = np.where(press == 0)[0][0]-1
        else:
            altmax = len(press) -1

        if opc[:altmax].any() == True:
            #also calculate the smoothened arrays for the specific latitude
            opc_smooth, z_array = Z_interpolater(opc,dz,altmax)
            press_smooth, z_array = Z_interpolater(press,dz,altmax)
        
            #calculate the number density and settling velocity for the calculated radius distribution with altitude
            lat_rho,vels,num_den_exp = find_Nd_vels(opc_smooth,press_smooth,ptcl_size,dz)
        
            #Run the number density in the cloud settling model
            N_new, zs, tt, N_evolution,flux,surface_radius= Cloud_Growth_CN(lat_rho,press_smooth,Radius_evol,Vel_evol,T, dz, dt, K)
            indexor = N_evolution.shape[0]
            N_time_evolution[:indexor,:,j] = N_evolution
            New_density[:indexor,j] = N_new
            Sed_evol_lat[:,j] = flux
        else:
            surface_radius = 0

    return New_density,N_time_evolution,Sed_evol_lat,surface_radius



def N_Water_buildup(Sed_evol_lat,LatV,dt): 
    #for NORTH POLE   
    #finding the amount of water being removed (measured in equivalent deposited thickness):
    #define constants:
    ah = 4e-6 #m molecule cross section
    rhoh = 1000 #individual particle mass density of water ice (Cotton et al 2013 JRMS) kg/m3
    R_m = 3.3895e6 #radius of Mars in m
    m_p = (4/3)*rhoh*np.pi*(ah**3) #average mass of H2O ice particle

    #first of all we have to calculate the amount of CO2 being removed from the atmosphere,, this is going to be our 
    #boundary condition that we have selected for the finite difference (at z=0, downward flux = w(0)*n(0))

    #create empty matrix to store water deposited per latitude, evolving through time:
    H2O_deposit = np.zeros(shape = (Sed_evol_lat.shape[0],Sed_evol_lat.shape[1])) #first index is time, second is lat
    H2O_mass = np.zeros(shape = (Sed_evol_lat.shape[0],Sed_evol_lat.shape[1]))
    for i in range(Sed_evol_lat.shape[0]):

        # define the sedimentation flux of CO2 at each time step
        F_co2 = Sed_evol_lat[i,:]  #particles m^-2 s^-1

        for j in range(Sed_evol_lat.shape[1]):

            #find the mass flux of H2O assuming nCO2 = nH2O and hence flux CO2 = flux H2O
            
            sin = np.abs(np.sin(np.deg2rad(LatV[j])))
            sin2 = np.abs(np.sin(np.deg2rad(LatV[j+1])))
        
            height = R_m*(sin2-sin)
            SA = 2*np.pi*R_m*height
            
            Ndep =F_co2[j]*dt*SA # no of CO2 particles
            
            water_leeched_mass = Ndep*m_p #kg of H2O deposited

            water_leeched_thickness =F_co2[j]*dt*((4/3)*np.pi*((ah)**3)) #m

            H2O_deposit[i,j] = water_leeched_thickness
            H2O_mass[i,j] = water_leeched_mass

    #cumulative sum of amount of water being deposited so that the plot shows the build up        
    H2O_cum_thick = np.cumsum(H2O_deposit, axis = 0)
    H2O_cum_mass = np.cumsum(H2O_mass, axis = 0)

    return H2O_cum_thick, H2O_cum_mass

def N_CO2_buildup(Sed_evol_lat,surface_radius,LatV,dt):  
    # for NORTH POLE
    #finding the amount of CO2 being removed (measured in equivalent deposited thickness):
    #define constants:
    rhoc = 1600 #kg/m3 #solid density of CO2 ice (Find reference)
    R_m = 3.3895e6 #radius of Mars in m
    
    #calculate average particale mass
    m_p = (4/3)*rhoc*np.pi*(surface_radius**3) #average mass of CO2 ice particle


    #first of all we have to calculate the amount of CO2 being removed from the atmosphere,, this is going to be our 
    #boundary condition that we have selected for the finite difference (at z=0, downward flux = w(0)*n(0))

    #create empty matrix to store water deposited per latitude, evolving through time:
    CO2_deposit = np.zeros(shape = (Sed_evol_lat.shape[0],Sed_evol_lat.shape[1])) #first index is time, second is lat
    CO2_mass = np.zeros(shape = (Sed_evol_lat.shape[0],Sed_evol_lat.shape[1]))
    for i in range(Sed_evol_lat.shape[0]):

        # define the sedimentation flux of CO2 at each time step
        F_co2 = Sed_evol_lat[i,:]  #particles m^-2 s^-1

        for j in range(Sed_evol_lat.shape[1]):
            #find the mass flux of CO2:
            sin = np.abs(np.sin(np.deg2rad(LatV[j])))
            sin2 = np.abs(np.sin(np.deg2rad(LatV[j+1])))
        
            height = R_m*(sin2-sin)
            SA = 2*np.pi*R_m*height
            
            Ndep =F_co2[j]*dt*SA # no of CO2 particles
            
            co2_leeched_mass = Ndep*m_p #kg of CO2 deposited

            co2_leeched_thickness =(F_co2[j]*dt*((4/3)*np.pi*((surface_radius)**3))) #m

            CO2_deposit[i,j] = co2_leeched_thickness
            CO2_mass[i,j] = co2_leeched_mass

    #cumulative sum of amount of CO2 being deposited so that the plot shows the build up        
    CO2_cum_thick = np.cumsum(CO2_deposit, axis = 0)
    CO2_cum_mass = np.cumsum(CO2_mass, axis = 0)

    
    return CO2_cum_thick, CO2_cum_mass


def S_Water_buildup(Sed_evol_lat,LatV,dt): 
    #for SOUTH POLE
    #finding the amount of water being removed (measured in equivalent deposited thickness):
    #define constants:
    ah = 4e-6 #m molecule cross section
    rhoh = 1000 #individual particle mass density of water ice (Cotton et al 2013 JRMS) kg/m3
    R_m = 3.3895e6 #radius of Mars in m
    m_p = (4/3)*rhoh*np.pi*(ah**3) #average mass of H2O ice particle

    #first of all we have to calculate the amount of CO2 being removed from the atmosphere,, this is going to be our 
    #boundary condition that we have selected for the finite difference (at z=0, downward flux = w(0)*n(0))

    #create empty matrix to store water deposited per latitude, evolving through time:
    H2O_deposit = np.zeros(shape = (Sed_evol_lat.shape[0],Sed_evol_lat.shape[1])) #first index is time, second is lat
    H2O_mass = np.zeros(shape = (Sed_evol_lat.shape[0],Sed_evol_lat.shape[1]))
    for i in range(Sed_evol_lat.shape[0]):

        # define the sedimentation flux of CO2 at each time step
        F_co2 = Sed_evol_lat[i,:]  #particles m^-2 s^-1

        for j in range(Sed_evol_lat.shape[1]):

            #find the mass flux of H2O assuming nCO2 = nH2O and hence flux CO2 = flux H2O
            
            sin = np.abs(np.sin(np.deg2rad(np.abs(LatV[j]))))
            sin2 = np.abs(np.sin(np.deg2rad(np.abs(LatV[j+1]))))
        
            height = R_m*np.abs(sin2-sin)
            SA = 2*np.pi*R_m*height
            
            Ndep =F_co2[j]*dt*SA # no of CO2 particles
            
            water_leeched_mass = Ndep*m_p #kg of H2O deposited

            water_leeched_thickness =F_co2[j]*dt*((4/3)*np.pi*((ah)**3)) #m

            H2O_deposit[i,j] = water_leeched_thickness
            H2O_mass[i,j] = water_leeched_mass

    #cumulative sum of amount of water being deposited so that the plot shows the build up        
    H2O_cum_thick = np.cumsum(H2O_deposit, axis = 0)
    H2O_cum_mass = np.cumsum(H2O_mass, axis = 0)

    return H2O_cum_thick, H2O_cum_mass

def S_CO2_buildup(Sed_evol_lat,surface_radius,LatV,dt): 
    #for SOUTH POLE   
    #finding the amount of CO2 being removed (measured in equivalent deposited thickness):
    #define constants:
    rhoc = 1600 #kg/m3 #solid density of CO2 ice (Find reference)
    R_m = 3.3895e6 #radius of Mars in m
    
    #calculate average particale mass
    m_p = (4/3)*rhoc*np.pi*(surface_radius**3) #average mass of CO2 ice particle

    #first of all we have to calculate the amount of CO2 being removed from the atmosphere,, this is going to be our 
    #boundary condition that we have selected for the finite difference (at z=0, downward flux = w(0)*n(0))

    #create empty matrix to store water deposited per latitude, evolving through time:
    CO2_deposit = np.zeros(shape = (Sed_evol_lat.shape[0],Sed_evol_lat.shape[1])) #first index is time, second is lat
    CO2_mass = np.zeros(shape = (Sed_evol_lat.shape[0],Sed_evol_lat.shape[1]))
    for i in range(Sed_evol_lat.shape[0]):

        # define the sedimentation flux of CO2 at each time step
        F_co2 = Sed_evol_lat[i,:]  #particles m^-2 s^-1

        for j in range(Sed_evol_lat.shape[1]):

            #find the mass flux of CO2:
            sin = np.abs(np.sin(np.deg2rad(np.abs(LatV[j]))))
            sin2 = np.abs(np.sin(np.deg2rad(np.abs(LatV[j+1]))))
        
            height = R_m*np.abs(sin2-sin)
            SA = 2*np.pi*R_m*height
            
            Ndep =F_co2[j]*dt*SA # no of CO2 particles
            
            co2_leeched_mass = Ndep*m_p #kg of CO2 deposited

            co2_leeched_thickness =(F_co2[j]*dt*((4/3)*np.pi*((surface_radius)**3))) #m

            CO2_deposit[i,j] = co2_leeched_thickness
            CO2_mass[i,j] = co2_leeched_mass

    #cumulative sum of amount of CO2 being deposited so that the plot shows the build up        
    CO2_cum_thick = np.cumsum(CO2_deposit, axis = 0)
    CO2_cum_mass = np.cumsum(CO2_mass, axis = 0)

    return CO2_cum_thick, CO2_cum_mass


def remove_outliers(input_array):
    '''function that sets all values greater that 3 std dev from the mean to the cutoff value
    also removes all zero values'''
    Array = input_array.copy()
    #remove negative values
    Array[Array<0]=0
    #calculate original mean
    MnVal = np.median(Array)
    #remove values too large
    cutoff = np.median(Array)+(3*np.std(Array))
    #and too small
    cutoff2 = np.median(Array)-(3*np.std(Array))
    
    Array[Array>=cutoff]=MnVal
    Array[Array<=cutoff2]=MnVal
    return Array

def remove_outliers_radius(input_array):
    '''function that sets all values greater that 3 std dev from the mean to the cutoff value
    also removes all zero values'''
    Array = input_array.copy()
    #remove negative values
    Array[Array<0]=0
    #calculate original mean
    MnVal = np.median(Array)
    #remove values too large
    cutoff = np.median(Array)+(1.5*np.std(Array))
    #and too small
    cutoff2 = np.median(Array)-(1.5*np.std(Array))
    
    Array[Array>=cutoff]=MnVal
    Array[Array<=cutoff2]=MnVal
    return Array

def Ls_filler(Matrix, missing_Ls, LatV, intnames,dups):
    '''Fills in the gaps in the cumulative water deposit
    Matrix: either Water_deposit, or Mass_deposit for either water or CO2
    missing_Ls = array of missing Ls values (as strongs matching the names)
    LatV: latitude array
    intnames = list of Ls's as integers
    dups = array that dictates how many rows to duplicate(if gap is 3 Ls long then duplicate three times)'''
    count = 0
    for i in range(len(missing_Ls)):
        idx = np.abs(intnames-missing_Ls[i]).argmin() + count
        
        if idx < Matrix.shape[0]-1:
            insert = (Matrix[idx]+ Matrix[idx+1])/2
            lowess = sm.nonparametric.lowess(insert, LatV[:-1],frac=0.3)
            Matrix = np.insert(Matrix, np.repeat([idx+1],dups[i]),lowess[:, 1],axis=0)
            count = dups[i] +count
        else:
            insert = Matrix[idx]
            lowess = sm.nonparametric.lowess(insert, LatV[:-1],frac=0.3)
            Matrix = np.insert(Matrix, np.repeat([idx+1],dups[i]),lowess[:, 1],axis=0)
            count = dups[i] +count
    return Matrix
