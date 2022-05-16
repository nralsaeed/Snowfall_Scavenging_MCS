##Binning modules used for the MCS profile data:

import numpy as np
import pandas as pd
from scipy.stats import binned_statistic_2d as bin2d
from scipy.interpolate import pchip

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

def radius_profile(NumDen, Temp, Tfrost,Press, Z, dz,Nz):
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
    
    #setting up some constants for later:
    ac = 1e-5 #m starting radius of CO2 ice particles (10 microns)
    WN = 667.3 #(wavenumber center frequency of 15 micron band)  
    sigma = 5.67e-8 #W/m2 K4 Stephen Boltzmann constant
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
    emm = 0.80686 #placeholder for now, using only Mie params for r = 10
    Cp = 700 #J/kg K specific heat of CO2 atmosphere.. an assumption (ask Paul?)
    Latent = 5.9e5 #J/kg latent heat of CO2 sublimation (hayne 2014)
    Qc = 3.0249 #efficiency of extinction for CO2 
    rhoc = 1600 #kg/m3 #solid density of CO2 ice (Find reference)
    Kb = 1.38064852e-23 #m2 kg/ s2 K #boltzmann constant
    dc = 3.3e-10 #m CO2 molecule diameter
    Tav = 140 #K #average temp at the poles
    mu = 8.5e-6 #(kg/ m s) OR (Pa s) # molecular viscosity of atmosphere at 150 K
    layerheight = dz
    fctr = np.int(1e3/dz)
    
    #radiis = np.linspace(20,80,26)*1e-6
    #ftimes = np.zeros(len(radiis))
 

    #First lets define the lookup matrix for the flux equivalent widths 
    #that we will use in the growth calculation
    E_matrix = np.zeros(((fctr*Nz)+1,(fctr*Nz)+1))
    for i in range((fctr*Nz)+1):
        for j in range((fctr*Nz)+1): 
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
    Fg_down = np.zeros((fctr*Nz)+1)
    Fg_up = np.zeros((fctr*Nz)+1)
    
    #now we calculate the IR gas flux coming off from each layer, in both the upward and downward direction
    for i in range(1,(fctr*Nz)+1):
        #from gas cooling:
        Fg_down[i] = 2*np.pi*sum((E_matrix[i,j-1]-E_matrix[i,j])*Planck(Tfrost[j-1],WN) for j in range(1,i+1))
        Fg_up[i] = 2*np.pi*Planck(Ts,WN)*((0.5*nomw) - E_matrix[i,(fctr*Nz)]) \
                + 2*np.pi*sum((E_matrix[i,j+1]-E_matrix[i,j])*Planck(Tfrost[j+1],WN) for j in range(i+1,(fctr*Nz)))
    
    #boundary conditions:
    Fg_down[0] = 0 #nothing penetrating the surface
    Fg_up[(fctr*Nz)] = np.pi*Planck(Ts,WN)*nomw #total leaving the top of the atmosphere is whatever heat came from the surface to begin wtih
    
    #create empty matrix to store radius at each altitude 
    r_update = np.zeros((fctr*Nz)+1)
    #create empty array to store settling velocity at each altitude
    ws = np.zeros((fctr*Nz)+1)
    #create empty array to store tracker of fall time
    fall_time_track = np.zeros((fctr*Nz)+1)
    #create empty array to store tracker of cumulative time
    cum_track = np.zeros((fctr*Nz)+1)
    
        
    
    #set starting radius and starting altitude, this will be updated in the loop
    rold = ac
    
    
    idx = fctr*Nz -1 
    
        
    
    r_update[idx] = rold
    
    #calculate the settling velocity associated with this starting particle size
    mfp = (Kb*Tav)/(np.sqrt(2)*np.pi*(dc**2)*Press[idx]) #mean free path in air
    ws[idx] = ((2*rhoc*gmars*(rold**2))/(9*mu))*(1+ (mfp/rold)) #m/s
    #start time counter
    Time_total = 0
    #start the loop to track the particle all the way to the surface:
    while idx > 0:

        #Now to calculate the new radius at the end of this layer
        #cooling from gas W/m2
        delFgas = ((np.abs(Fg_up[idx]) - np.abs(Fg_down[idx])) - (np.abs(Fg_up[idx+1]) - np.abs(Fg_down[idx+1])))
        
        #cooling from aerosol itself W/m2
        delFaer = 0.5*emm*sigma*(Tfrost[idx]**4)
                
        delF = delFgas 
        dt = dz/ws[idx]
        massgrowth = np.abs(delF*dt/Latent) #where Latent is the latent heat of CO2 condesation

        if NumDen[idx] < 1000:
            particlegrowth = (massgrowth)/(1000*dz) #mass growth per particle
        else:
            particlegrowth = (massgrowth)/(NumDen[idx]*dz) #mass growth per particle
        
        rnew = ((3*particlegrowth/(4*np.pi*rhoc))+(rold**3))**(1/3)
        r_update[idx-1] = rnew

        #calculate the settling velocity associated with this new particle size
        mfp = (Kb*Tav)/(np.sqrt(2)*np.pi*(dc**2)*Press[idx]) #mean free path in air
        ws[idx-1] = ((2*rhoc*gmars*(rnew**2))/(9*mu))*(1+ (mfp/rnew)) #m/s
        
        #save how long it took for this particle to fall dz height
        fall_time_track[idx-1] = dt
        
        #update time counter
        Time_total = Time_total + dt
        
        #save how long it took to get to this radius
        cum_track[idx-1] = Time_total
            

        #update index counter
        idx = idx-1
        
        #update radius
        rold = rnew
    
    #for rad in range(len(radiis)):
        #vels = ((2*rhoc*gmars*(radiis[rad]**2))/(9*mu))*(1+ (mfp/radiis[rad]))
        #ftimes[rad] = dz/vels
    return r_update, ws



def Ls_radius_dist(rho, Temperature, FrostTemp,Pressure, dz): 
    '''Finds the radius distribution of the CO2 ice in a specific Ls by 
    taking the median CO2 ice density distribution within the polar vortex area (all lats) 
    and inputting it into the radius profile finder function'''
    
    #first take the median of all the values along the latitudes
    msk_rho = np.ma.masked_equal(rho, 0)
    msk_temp = np.ma.masked_equal(Temperature, 0)
    msk_press = np.ma.masked_equal(Pressure, 0)
    msk_tfrost = np.ma.masked_equal(FrostTemp, 0)
    
    M_rho    = np.ma.median(msk_rho, axis=1).filled(0)
    M_temp   = np.ma.median(msk_temp, axis=1).filled(0)
    M_press  = np.ma.median(msk_press, axis=1).filled(0)
    M_tfrost = np.ma.median(msk_tfrost, axis=1).filled(0)

    #index of highest alt with pressure reading:
    if len(np.where(M_press == 0)[0]) > 0:
        altmax = np.where(M_press == 0)[0][0] -1
    else: 
        altmax = len(M_rho) -1

    #altitude limit
    Z = int(altmax)*1e3 #m
    #original altitude separation
    Idz = 1000 #m

    Nz = int(round(Z/float(Idz)))
    z = np.linspace(0, Nz*Idz, Nz+1) # Mesh points in space

    
    #the factor by which you are multiplying the altitude seperation
    fctr = np.int(1e3/dz)
    
    
    #You must now truncate all the input arrays so that they end at the maximum altitude:
    M_rho = M_rho[:altmax]
    M_temp = M_temp[:altmax]
    M_tfrost = M_tfrost[:altmax]
    M_press = M_press[:altmax]
    
    
    #smoothen the plots (because as is, it is too jumpy and causes issues with the finite diff method)
    #plus this way I can have a smaller spacial step in altitude than 1km
    #set up the smoothing with Piecewise Cubic Hermite Interpolating Polynomial
    pch = pchip(z[:-1],M_rho)
    pchT = pchip(z[:-1],M_temp)
    pchP = pchip(z[:-1],M_press)
    pchFT = pchip(z[:-1],M_tfrost)

    #create more finely space array in z axis (has to be at least 2*Nz)
    xs = np.linspace(0, Z, (fctr*Nz)+1)
    
    #interpolate for both no density u and settling velocity w
    Rho_smooth = pch(xs)
    Press = pchP(xs)
    Temp = pchT(xs)
    Tfrost = pchFT(xs)
    #make sure there are no negative values accidentally reached from the interpolation:
    Rho_smooth[Rho_smooth<0]=0
    Press[Press<0]=0
    Temp[Temp<0]=0
    Tfrost[Tfrost<0]=0
    
    #calculate the radius distribution in atm, and associated velocities
    radius_z, ws = radius_profile(Rho_smooth, Temp, Tfrost,Press, Z, dz,Nz)
    
    return radius_z, ws

    
def Cloud_Model_CN(I,radii, vels, Pressure, T, dz, dt, K):
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
    
    #index of highest alt with pressure reading:
    if len(np.where(Pressure == 0)[0]) > 0:
        altmax = np.where(Pressure == 0)[0][0] - 1
    else: 
        altmax = len(I) - 1

    #altitude limit
    Z = int(altmax)*1e3 #m
    #original altitude separation
    Idz = 1000 #m

    Nz = int(round(Z/float(Idz)))
    z = np.linspace(0, Nz*Idz, Nz+1) # Mesh points in space
    
    #the factor by which you are multiplying the altitude seperation
    fctr = np.int(1e3/dz)
    
    
    #You must now truncate the input array so that it ends at the maximum altitude:
    I = I[:altmax]
    
    #smoothen the plots (because as is, it is too jumpy and causes issues with the finite diff method)
    #plus this way I can have a smaller spacial step in altitude than 1km
    #set up the smoothing with Piecewise Cubic Hermite Interpolating Polynomial
    pch = pchip(z[:-1],I)
    
    #create more finely space array in z axis (has to be at least 2*Nz)
    xs = np.linspace(0, Z, (fctr*Nz)+1)
    
    #interpolate 
    I_smooth = pch(xs)
    
    #make sure there are no negative values accidentally reached from the interpolation:
    I_smooth[I_smooth<0]=0
    
    #Make sure to truncate the radius profile and velocities to the maximum altitude as well
    altmax2 = len(I_smooth)
    if altmax2 > len(radii):
        print('it dont match', len(radii),altmax2, len(I))
        I_smooth = I_smooth[:len(radii)]
        radius_z = radii
        ws = vels
    else:
        radius_z = radii[:altmax2]
        ws = vels[:altmax2]


    
    #create empty arrays to store
    u = np.zeros((fctr*Nz)) # unknown u at new time level
    u_n = np.zeros((fctr*Nz)) # u at the previous time level
    
    #create empty matrix to store mass flux  and radius growth at each altitude for every point in time
    u_evolution = np.zeros(shape =((fctr*Nz),Nt+1))
    #create empty array to store CO2 being removed at each time step at the surface
    sediment_evolution = np.zeros(Nt+1)

    # Set initial condition u(z,0) = I(z)
    for f in range(0, (fctr*Nz)):
        u_n[f] = I_smooth[f]
        
    u_evolution[:,0] = u_n
    
    
    #set up the tridiagonal matrix for Crank Nicholson
    #first extract the LHS coefficients:
    A = np.zeros(fctr*Nz+1)
    B = np.zeros(fctr*Nz+1)
    C = np.zeros(fctr*Nz+1)
    E = np.zeros(fctr*Nz+1)
    F = np.zeros(fctr*Nz+1)
    G = np.zeros(fctr*Nz+1)
    for i in range(0, (fctr*Nz)):
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
    M = np.diagflat([C[i] for i in range(fctr*Nz-1)], -1) +\
            np.diagflat([B[i] for i in range(fctr*Nz)]) +\
                np.diagflat([A[i] for i in range(fctr*Nz-1)], 1)
    #right hand side
    L = np.diagflat([G[i] for i in range(fctr*Nz-1)], -1) +\
            np.diagflat([F[i] for i in range(fctr*Nz)]) +\
                np.diagflat([E[i] for i in range(fctr*Nz-1)], 1)
    
    #initialize the matrix:
    Uinit = np.zeros((fctr*Nz))
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
        Usol[fctr*Nz-1] = 0
        
        
        #also at the surface the flux/removal of snow = velocity*density, so we set that in the boundary conditions as well:
        #boundary condition at surface:
        #No diffusion at the surface (second derivative goes to 0)
        ww1 = (ws[0]*(Uinit[1] - Uinit[0]))/(dz)
        uu1 = (Uinit[0]*(ws[1] - ws[0]))/(dz)
        
        Usol[0] = Uinit[0] + (uu1 + ww1)*dt 
        
        
        #update solution 
        Uinit = Usol
        
        #saving the column no density at each timestep
        u_evolution[:,j] = Uinit
        
        
        #store sedimentation rate:
        sediment = (ws[0]*Uinit[0])
        
        #saving the sedimentation at each time step (unit = m-2 s-1)
        sediment_evolution[j] = sediment 
        
        #update counter
        j = j+1
        #update time
        n = n + dt
        
        
    return Uinit, z, t, u_evolution, xs, sediment_evolution, radius_z # Uinit holds latest u


def All_lat_snowfall(N_c,Temp,Tfrost,Press,T,dz,dt,K,LatV):

    fctr = np.int(1e3/dz)

    #create matrix to store updated no density
    #the dimensions of this matrix depends on the chosen spacial step dz
    New_density = np.zeros(shape =(fctr*N_c.shape[0],N_c.shape[1])) 
    #create 3D matrix to store time evolution of number density (1 dimension for. space (altitude), 1d for no. density, and 1d for time)
    N_time_evolution = np.zeros(shape =(fctr*N_c.shape[0],np.int(T/dt)+1,N_c.shape[1])) 
    #first index is altitude, second is time, third is latitude
    #create matrix to store the radius distribution for each lat in here
    Radii_lats = np.zeros(shape =(fctr*N_c.shape[0],N_c.shape[1])) 
    #first index is altitude, second is latitude

    #create matrix to store time evolution of the sedimentation at surface for each latitude:
    Sed_evol_lat = np.zeros(shape =(np.int(T/dt)+1,N_c.shape[1])) #first index is time, second index is latitude
    radiuses, vels = Ls_radius_dist(N_c,Temp,Tfrost,Press, dz)
    for j in range(len(LatV)-1):        
        rho = N_c[:,j]
        press = Press[:,j]
        counting = sum(i > 0 for i in press)

        
        if rho.any() == True:
        #if counting > 3:
            N_new, zz, tt, N_evolution, zs,flux, radiiuses= Cloud_Model_CN(rho,radiuses,vels,press,T,dz,dt,K)
        else:
            continue
        indexor = N_evolution.shape[0]
        N_time_evolution[:indexor,:,j] = N_evolution
        New_density[:indexor,j] = N_new
        Sed_evol_lat[:,j] = flux
        Radii_lats[:indexor,j] = radiiuses[:-1]

    return New_density,N_time_evolution,Sed_evol_lat,Radii_lats


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

def N_CO2_buildup(Sed_evol_lat,Radii,LatV,dt):  
    # for NORTH POLE
    #finding the amount of CO2 being removed (measured in equivalent deposited thickness):
    #define constants:
    rhoc = 1600 #kg/m3 #solid density of CO2 ice (Find reference)
    R_m = 3.3895e6 #radius of Mars in m
    
    #define the final particle sizes at the surface level
    #its the median of the particle sizes at the surface across all latitudes for that Ls
    ac = np.nanmean(Radii[0,:]) #m ice particle cross section on surface #index 0 is at surface
    #calculate average particale mass
    m_p = (4/3)*rhoc*np.pi*(ac**3) #average mass of CO2 ice particle


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

            co2_leeched_thickness =(F_co2[j]*dt*((4/3)*np.pi*((ac)**3))) #m

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

def S_CO2_buildup(Sed_evol_lat,Radii,LatV,dt): 
    #for SOUTH POLE   
    #finding the amount of water being removed (measured in equivalent deposited thickness):
    #define constants:
    rhoc = 1600 #kg/m3 #solid density of CO2 ice (Find reference)
    R_m = 3.3895e6 #radius of Mars in m

    #define the final particle sizes at the surface level
    #its the median of the particle sizes at the surface across all latitudes for that Ls
    ac = np.nanmean(Radii[0,:]) #m ice particle cross section on surface #index 0 is at surface
    #calculate average particale mass
    m_p = (4/3)*rhoc*np.pi*(ac**3) #average mass of CO2 ice particle

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

            co2_leeched_thickness =(F_co2[j]*dt*((4/3)*np.pi*((ac)**3))) #m

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
    MnVal = np.mean(Array)
    #remove values too large
    cutoff = np.mean(Array)+(3*np.std(Array))
    #and too small
    cutoff2 = np.mean(Array)-(3*np.std(Array))
    
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
