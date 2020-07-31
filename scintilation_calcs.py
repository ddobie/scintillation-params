import pickle
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np

def extract_ne2001_params(filename='ne2001_vals_array.pickle'):
    vals_array = pickle.load(open(filename,'rb'))
    
    ls = vals_array[0,:,:]
    bs = vals_array[1,:,:]
    nu0 = vals_array[2,:,:]
    DM = vals_array[3,:,:]
    SM = vals_array[4,:,:]
    theta_F0 = vals_array[5,:,:]
    D = vals_array[6,:,:]

    ls *= u.deg
    bs *= u.deg
    nu0 *= u.GHz
    theta_F0 *= 1e-3
    theta_F0 *= u.mas
    D *= u.kpc
    
    return ls, bs, nu0, theta_F0, D
	
def find_nearest(v, arr):
    idx = (np.abs(arr - v)).argmin()
    return idx
	
def diffractive_params(nu, nu0, theta_F0, theta_s=1e-5*u.mas):
    theta_d = theta_F0 * (nu/nu0)**(6/5) #eq. 16
    point_source_correction =  np.minimum(np.ones(theta_d.shape), theta_d/theta_s) 
    
    m = np.ones(nu0.shape) * point_source_correction #eq. 14
    delta_nu = nu * (nu/nu0)**(17/5) #eq. 15
    
    t_d = 0.72 * 2*u.hour *(nu/nu0)**(6/5) / point_source_correction #eq. 17
    
    
    return m, delta_nu, t_d, theta_d

def refractive_params(nu, nu0, theta_F0, theta_s=1e-5*u.mas):
    theta_r = theta_F0 * (nu0/nu)**(11/5) #eq. 12
    point_source_correction = np.minimum(np.ones(theta_r.shape), theta_r/theta_s)
    m = (nu/nu0)**(17/30) * (point_source_correction)**(7/6) #eq. 11
    t_r = 0.72 * 2*u.hour * (nu0/nu)**(11/5) / point_source_correction #eq. 13
    
    
    weak = nu > nu0
    m[weak] = 0
    t_r[weak] = np.nan
    theta_r[weak] = np.nan

    return m, t_r, theta_r
	
def get_scint_info(skycoord, nu):
    nearest_l = find_nearest(skycoord.l.wrap_at(180*u.deg), ls[:,0])
    nearest_b = find_nearest(skycoord.b, bs[0,:])

    el = (nearest_l, nearest_b)
    
    nu0_val = nu0[el]
    theta_F0_val = theta_F0[el]*1e3
    
    print("nu0: {}".format(nu0_val))
    print()
    
    diff_params = diffractive_params(nu, nu0_val, theta_F0_val)
    ref_params = refractive_params(nu, nu0_val, theta_F0_val)
    
    print("Diffractive Scintillation")
    print("Modulation Index: {:.3f}".format(diff_params[0]))
    print("Scintillation Bandwidth: {:.1f}".format(diff_params[1].to(u.MHz)))
    print("Diffractive Timescale: {:.1f}".format(diff_params[2].to(u.min)))
    
    print()
    
    print("Refractive Scintillation")
    print("Modulation Index: {:.3f}".format(ref_params[0]))
    print("Refractive Timescale: {:.1f}".format(ref_params[1].to(u.day)))
	
ls, bs, nu0, theta_F0, D = extract_ne2001_params()

teresa_source = SkyCoord(71.0822, -64.5155, unit=(u.degree, u.degree), frame='fk5').transform_to('galactic')

get_scint_info(teresa_source, 0.9*u.GHz)