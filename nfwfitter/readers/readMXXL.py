#######################
# Reads in MXXL files for processing
#######################

import numpy as np
import astropy.io.ascii as asciireader
import binaryutils

import nfwutils
import catalog


#######################


class MXXLBinary(object):

    def __init__(self, filename):

        self.filename = filename
        self.parseBinary()

    ####

    def parseBinary(self):

        #The files are binary (little endian). As pseudo-C struct, the file looks like:
        with open(self.filename, 'rb') as input:
            input = open(self.filename, 'rb')

            #double lower_bound[2];                          // lower bound of area represented by plane (in this case in comoving Mpc/h)
            #NOTE: For "mass_map", units of Mpc/h
            #NOTE: For convergence_map, units of arcsec

            self.lower_bound = binaryutils.readArray(input, 'd', (2,))

            #double upper_bound[2];                          // upper bound of area represented by plane (in this case in comoving Mpc/h)

            #NOTE: For "mass_map", units of Mpc/h
            #NOTE: For convergence_map, units of arcsec
            self.upper_bound = binaryutils.readArray(input, 'd', (2,))

            #double plane_angle;                             // angle of coordinate axes of plane (should be 0.5 * M_PI)
            self.plane_angle = binaryutils.readVal(input, 'd')

            #double redshift;                                // redfhift of the plane (should be 0.9887)
            self.redshift = binaryutils.readVal(input, 'd')

            #double particle_mass;                           // mass of simulation particle
            self.particle_mass = binaryutils.readVal(input, 'f')

            #int    N_pixels[2];                             // number of pixels in each dimension (should both be 1024 in this case)
            self.npixels = binaryutils.readArray(input, 'i', (2,))

            #float  mass_density[N_pixels[0] * N_pixels[1]]; // surface mass density (in simulation units, i.e. 10^10 M_solar/h per comoving (Mpc/h)^2)
            self.data = binaryutils.readArray(input, 'f', self.npixels)

    def grid(self):

        gridDelta = (self.upper_bound - self.lower_bound)/(self.npixels)  # grid delta in arcsec

        X2,X1 = np.meshgrid(np.arange(self.npixels[1]), np.arange(self.npixels[0]))
        
        deltaX1_arcmin = (self.lower_bound[0] + (X1 + 0.5)*gridDelta[0])/60.
        deltaX2_arcmin = (self.lower_bound[1] + (X2 + 0.5)*gridDelta[1])/60.

        dL = nfwutils.global_cosmology.angulardist(self.redshift)

        deltaX1_mpc = (deltaX1_arcmin * dL * np.pi)/(180.*60)
        deltaX2_mpc = (deltaX2_arcmin * dL * np.pi)/(180.*60)
        return (deltaX1_mpc, deltaX2_mpc), (deltaX1_arcmin, deltaX2_arcmin)
        
        

        


    

#################################

class MXXLSimReader(object):

    def __init__(self, *args, **keywords):

        self._sim = None
        self._filebase = None

    #########

    def getCosmology(self):

        return nfwutils.Cosmology(omega_m = 0.25, omega_l = 0.75, h=0.73)

    #########

    def load(self, filebase):

        if self._sim is None or self._filebase != filebase:
            self._sim = MXXLSim(filebase)
            self._filebase = filebase
        
        return self._sim.copy()
    
###############################


class MXXLSim(catalog.Catalog):


    #########

    def __init__(self, filebase):

        print 'Loading %s' % filebase

        super(MXXLSim, self).__init__()

        # returns set of angle and mpc distances from the true center (flattened)
        # for each position, provide reduced shear g1 and g2, as well as z and beta for each source

        kappafile = '{0}.convergence_map'.format(filebase)
        gamma1file = '{0}.shear_1_map'.format(filebase)
        gamma2file = '{0}.shear_2_map'.format(filebase)
        answerfile = '{0}.answer'.format(filebase)

        kappa = MXXLBinary(kappafile)
        gamma1 = MXXLBinary(gamma1file)
        gamma2 = MXXLBinary(gamma2file)


        self.zcluster = kappa.redshift


        delta_mpc, delta_arcmin = kappa.grid()
        delta_mpc = [x.flatten() for x in delta_mpc]
        delta_arcmin = [x.flatten() for x in delta_arcmin]


        self.x_mpc = delta_mpc[0]
        self.y_mpc = delta_mpc[1]
        self.x_arcmin = delta_arcmin[0]
        self.y_arcmin = delta_arcmin[1]



        beta_inf = nfwutils.global_cosmology.beta([1e6], self.zcluster)

        # Three components that I can plot
        self.gamma1_inf = beta_inf*gamma1.data.flatten() # Shear wrt x axis
        self.gamma2_inf = beta_inf*gamma2.data.flatten() # Shear wrt 45 deg
        self.kappa_inf = beta_inf*kappa.data.flatten()

    def reshape_components(self) :
        import numpy as np
        for component in [ self.gamma1_inf, self.gamma2_inf, self.kappa_inf ] :
            component_shape = int(np.sqrt(component.shape))

            component = component.reshape((component_shape, component_shape))


        
        
    


    
    
    

    
    
    
