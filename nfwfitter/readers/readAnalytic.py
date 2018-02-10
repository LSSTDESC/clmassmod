#######################
# Creates 2D shear maps using Diemer density profile
#######################

import numpy as np
import scipy.interpolate
import yaml

import nfwutils
import nfwmodeltools
import catalog
import colossusMassCon as cmc

import colossus.cosmology.cosmology as cCosmo
import colossus.halo.concentration as chc
import colossus.halo.profile_dk14 as dk14prof
import colossus.defaults as cDefaults


#######################

    

#################################

class AnalyticSimReader(object):

    def __init__(self, *args, **keywords):

        self._sim = None
        self._filebase = None
        self.profile = None

    #########

    def configure(self, config):
        self.profile = config['analytic_profile']
        

    #########

    def getCosmology(self):

        return nfwutils.Cosmology(omega_m = 0.25, omega_l = 0.75, h=0.73)

    #########

    def load(self, filebase):

        nfwutils.global_cosmology.set_cosmology(self.getCosmology())
        cmc.matchCosmo()

        if self._sim is None or self._filebase != filebase:
            self._sim = self.profile(filebase)
            self._filebase = filebase
        
        return self._sim.copy()
    
###############################
###############################

def calcLensingTerms(density_profile, rmax, max_r_integrate = cDefaults.HALO_PROFILE_SURFACE_DENSITY_MAX_R_INTEGRATE):

    assert(rmax < cDefaults.HALO_PROFILE_SURFACE_DENSITY_MAX_R_INTERPOLATE)

    log_min_r = np.log(cDefaults.HALO_PROFILE_DELTA_SIGMA_MIN_R_INTERPOLATE)
    log_max_r = np.log(np.max(rmax) * 1.01)
    table_log_r = np.arange(log_min_r, log_max_r + 0.01, 0.01)
    print table_log_r, len(table_log_r)
    table_r = np.exp(table_log_r)
    table_log_Sigma = np.log(density_profile.surfaceDensity(table_r, interpolate=False, max_r_integrate = max_r_integrate))

    print table_log_Sigma

    log_surface_density_interp = scipy.interpolate.InterpolatedUnivariateSpline(table_log_r, table_log_Sigma)

    kappa_enc_integrand = np.exp(table_log_Sigma)*(table_r**2)

    kappa_enc_traprule_element = (table_log_r[1:] - table_log_r[:-1])*(kappa_enc_integrand[1:] + kappa_enc_integrand[:-1])

    kappa_enc_traprule = 0.5*np.cumsum(kappa_enc_traprule_element)

    log_kappa_enc_interp = scipy.interpolate.InterpolatedUnivariateSpline(table_log_r[1:], np.log(kappa_enc_traprule))

    surface_density = lambda r: np.exp(log_surface_density_interp(np.log(r)))

    deltaSigma = lambda r: np.exp(log_kappa_enc_interp(np.log(r)))*2./r**2 - surface_density(r)

    return surface_density, deltaSigma

###

def calcLensing(density_profile, r_kpch, zcluster, max_r_integrate = cDefaults.HALO_PROFILE_SURFACE_DENSITY_MAX_R_INTEGRATE):

    surfacedensity_func, deltaSigma_func = calcLensingTerms(density_profile, np.max(r_kpch), max_r_integrate = max_r_integrate)
    surfacedensity = surfacedensity_func(r_kpch)
    deltaSigma = deltaSigma_func(r_kpch)


        

    curcosmo = nfwutils.global_cosmology
    Dl = curcosmo.angulardist(zcluster)
    beta_inf = curcosmo.beta([1e6], zcluster)
    sigma_crit = (curcosmo.v_c**2/(4*np.pi*curcosmo.G))/(Dl*beta_inf)  #units are M_dot / Mpc^2
    convert_units = 1./(curcosmo.h*1e6)
    converted_sigma_crit = sigma_crit * convert_units
       
    gamma_t_inf = deltaSigma/converted_sigma_crit
    kappa_inf = surfacedensity/converted_sigma_crit

    return gamma_t_inf, kappa_inf


##########

class DiemerAnalyticSim(catalog.Catalog):


    #########

    def __init__(self, filebase):

        print 'Loading %s' % filebase

        #the file in this case is a configuration file with information about mass, concentration, redshift, etc

        super(AnalyticSim, self).__init__()

        with open(filebase) as input:
            config = yaml.load(input)

        max_dist = config['max_dist']
        gridlength = config['gridlength']
        zcluster = config['zcluster']
        m200 = config['m200']
        c200 = config['c200']
    

        xs = np.linspace(-max_dist, max_dist, gridlength)
        ys = np.linspace(-max_dist, max_dist, gridlength)

        x_mpc, y_mpc = [grid.flatten() for grid in np.meshgrid(xs, ys, indexing='ij')]
        r_mpc = np.sqrt(x_mpc**2 + y_mpc**2)

        
        Dl = nfwutils.global_cosmology.angulardist(zcluster)
        x_arcmin = (x_mpc/Dl)*(180./np.pi)*60.
        y_arcmin = (y_mpc/Dl)*(180./np.pi)*60.

        #Diemer radii are in units of kpc/h\n",
        r_kpch = (r_mpc*1000*nfwutils.global_cosmology.h)

        density_profile = dk14prof.getDK14ProfileWithOuterTerms(M = m200, c = c200, z = zcluster, mdef = '200c')


        gamma_t_inf, kappa_inf = calcLensing(density_profile, r_kpch, zcluster)

        posangle = np.arctan2(y_mpc, x_mpc)

        cos2phi = np.cos(2*posangle)
        sin2phi = np.sin(2*posangle)
    
        gamma1_inf = -gamma_t_inf*cos2phi
        gamma2_inf = -gamma_t_inf*sin2phi
            
        self.zcluster = zcluster

        self.x_mpc = x_mpc
        self.y_mpc = y_mpc
        self.x_arcmin = x_arcmin
        self.y_arcmin = y_arcmin

        self.gamma1_inf = gamma1_inf
        self.gamma2_inf = gamma2_inf
        self.kappa_inf = kappa_inf

######################
######################
 
class NFWAnalyticSim(catalog.Catalog):


    #########

    def __init__(self, filebase):

        print 'Loading %s' % filebase

        #the file in this case is a configuration file with information about mass, concentration, redshift, etc

        super(AnalyticSim, self).__init__()

        with open(filebase) as input:
            config = yaml.load(input)

        max_dist = config['max_dist']
        gridlength = config['gridlength']
        zcluster = config['zcluster']
        m200 = config['m200']
        c200 = config['c200']
        rscale = nfwutils.rscaleConstM(m200, c200, zcluster, 200.)
    

        xs = np.linspace(-max_dist, max_dist, gridlength)
        ys = np.linspace(-max_dist, max_dist, gridlength)

        x_mpc, y_mpc = [grid.flatten() for grid in np.meshgrid(xs, ys, indexing='ij')]
        r_mpc = np.sqrt(x_mpc**2 + y_mpc**2)

        
        Dl = nfwutils.global_cosmology.angulardist(zcluster)
        x_arcmin = (x_mpc/Dl)*(180./np.pi)*60.
        y_arcmin = (y_mpc/Dl)*(180./np.pi)*60.

        rho_c_over_sigma_c = 1.5 * nfwutils.global_cosmology.angulardist(zcluster) * nfwutils.global_cosmology.beta([1e6], zcluster)[0] * nfwutils.global_cosmology.hubble2(zcluster) / nfwutils.global_cosmology.v_c**2


        gamma_t_inf = nfwmodeltools.NFWShear(r_mpc,
                                             c200,
                                             rscale,
                                             rho_c_over_sigma_c)
        kappa_inf = nfwmodeltools.NFWKappa(r_mpc,
                                           c200,
                                           rscale,
                                           rho_c_over_sigma_c)


        posangle = np.arctan2(y_mpc, x_mpc)

        cos2phi = np.cos(2*posangle)
        sin2phi = np.sin(2*posangle)
    
        gamma1_inf = -gamma_t_inf*cos2phi
        gamma2_inf = -gamma_t_inf*sin2phi
            
        self.zcluster = zcluster

        self.x_mpc = x_mpc
        self.y_mpc = y_mpc
        self.x_arcmin = x_arcmin
        self.y_arcmin = y_arcmin

        self.gamma1_inf = gamma1_inf
        self.gamma2_inf = gamma2_inf
        self.kappa_inf = kappa_inf

######################


        
        
    


    
    
    

    
    
    
