dir_path = '/Users/mholman/Dropbox/support/'

import os
import bisect
import pickle

from collections import defaultdict
from collections import Counter

import scipy
import numpy as np
import healpy as hp
import spiceypy as spice

import MPC_library
import kepcart as kc
import tracklets as tr

# Load a few spice kernels
dir_path = '/Users/mholman/Dropbox/support/'
spice.furnsh(dir_path+'/kernels/MetaK_jwst.txt')

au2m = 149597870700
au_km = au2m/1000.

def convert_Mars(line, dt):
    ''' Converting the format in Mars's file
        good_preds_loc.txt to simulated tracklets
        based on the angular rates
        
        dt is a +/- time offset in seconds
    '''
    objName, MJD, RA, Dec, rate_ra, rate_dec, flux = line.rstrip().split()
    RA, Dec = float(RA), float(Dec)
    rate_ra, rate_dec = float(rate_ra), float(rate_dec)

    # Change RA rate and Dec rate units. Incorporate cos(Dec) term
    # in RA rate.
    # Units are being converted from arcsec per hour to
    # degrees per second.
    dra = rate_ra/np.cos(Dec*np.pi/180)*(dt/3600)*(1/3600)
    ddec = rate_dec*(dt/3600)*(1/3600)
    
    ra = RA+dra
    dec = Dec+ddec
    
    et = spice.str2et('JD '+ str(float(MJD)+2400000.5+dt/(24*60*60)))
    jd_tdb = spice.j2000() + et/(24*60*60)
    flux = float(flux)
    return objName, jd_tdb, ra, dec, flux

def convert_Mars_v2(line, dt):
    ''' Converting the format in Mars's file
        good_preds_loc.txt to simulated tracklets
        based on the angular rates
        
        dt is a +/- time offset in seconds
    '''
    objName, MJD, RA, Dec, x, y, rate_ra, rate_dec, flux = line.rstrip().split()
    RA, Dec = float(RA), float(Dec)
    x, y = float(x), float(y)
    rate_ra, rate_dec = float(rate_ra), float(rate_dec)

    # Change RA rate and Dec rate units. Incorporate cos(Dec) term
    # in RA rate.
    # Units are being converted from arcsec per hour to
    # degrees per second.
    dra = rate_ra/np.cos(Dec*np.pi/180)*(dt/3600)*(1/3600)
    ddec = rate_dec*(dt/3600)*(1/3600)
    
    ra = RA+dra
    dec = Dec+ddec
    
    et = spice.str2et('JD '+ str(float(MJD)+2400000.5+dt/(24*60*60)))
    jd_tdb = spice.j2000() + et/(24*60*60)
    flux = float(flux)
    return objName, jd_tdb, ra, dec, flux

def read_Mars_detection(line):
    ''' Reading the data in Mars's file good_preds_loc.txt and
    converting MJD to jd_tdb.
    '''
    objName, MJD, RA, Dec, rate_ra, rate_dec, flux = line.rstrip().split()
    ra, dec = float(RA), float(Dec)
    rate_ra, rate_dec = float(rate_ra), float(rate_dec)

    et = spice.str2et('JD '+ str(float(MJD)+2400000.5))
    jd_tdb = spice.j2000() + et/(24*60*60)
    flux = float(flux)
    return objName, jd_tdb, ra, dec, rate_ra, rate_dec, flux



def simplify_names(infilename, outfilename):
    names_dict={}
    with open(infilename) as infile, open(outfilename, 'w') as outfile:
        line=infile.readline().replace('filename', '')
        outfile.write(line)
        for i, line in enumerate(infile):
            objName, filename, MJD, RA, Dec, rate_ra, rate_dec, flux = line.rstrip().split()
            trk_name = '%06d' % (i)
            names_dict[trk_name]=objName
            outstring = '%s %s %s %s %s %s %s\n' % (trk_name, MJD, RA, Dec, rate_ra, rate_dec, flux)
            outfile.write(outstring)
    return names_dict

def simplify_names_v2(infilename, outfilename):
    names_dict={}
    with open(infilename) as infile, open(outfilename, 'w') as outfile:
        line=infile.readline().replace('filename', '')
        outfile.write(line)
        for i, line in enumerate(infile):
            print(len(line.rstrip().split()))
            print(line.rstrip())
            objName, filename, MJD, RA, Dec, x, y, rate_ra, rate_dec, likelihood, flux = line.rstrip().split()
            trk_name = '%06d' % (i)
            names_dict[trk_name]=objName
            outstring = '%s %s %s %s %s %s %s %s %s %s\n' % (trk_name, MJD, RA, Dec, x, y, rate_ra, rate_dec, likelihood, flux)
            outfile.write(outstring)
    return names_dict

def make_tracklets_from_rate_detections(infilename, outfilename, dt=600):
    ''' 
    Convert RA/Dec/rate detections from Mars's file KBMOD
    results to simulated tracklets with three RA/Dec detections.
    The observations are offset by a time dt in seconds.
    '''
    outfile = open(outfilename, 'w')
    with open(infilename) as file:
        file.readline()
        for i, line in enumerate(file):

            objName, jd_tdb, ra, dec, flux = convert_Mars(line, -dt)
            outstring = '%s %.11lf %.13lf %.13lf %.5lf\n' % (objName, jd_tdb, ra, dec, flux)
            outfile.write(outstring)
        
            objName, jd_tdb, ra, dec, flux = convert_Mars(line, 0)
            outstring = '%s %.11lf %.13lf %.13lf %.5lf\n' % (objName, jd_tdb, ra, dec, flux)
            outfile.write(outstring)
        
            objName, jd_tdb, ra, dec, flux = convert_Mars(line, +dt)
            outstring = '%s %.11lf %.13lf %.13lf %.5lf\n' % (objName, jd_tdb, ra, dec, flux)
            outfile.write(outstring)
    outfile.close()

def make_tracklets_from_rate_detections_v2(infilename, outfilename, dt=600):
    ''' 
    Convert RA/Dec/rate detections from Mars's file KBMOD
    results to simulated tracklets with three RA/Dec detections.
    The observations are offset by a time dt in seconds.
    '''
    outfile = open(outfilename, 'w')
    with open(infilename) as file:
        file.readline()
        for i, line in enumerate(file):

            objName, jd_tdb, ra, dec, flux = convert_Mars_v2(line, -dt)
            outstring = '%s %.11lf %.13lf %.13lf %.5lf\n' % (objName, jd_tdb, ra, dec, flux)
            outfile.write(outstring)
        
            objName, jd_tdb, ra, dec, flux = convert_Mars_v2(line, 0)
            outstring = '%s %.11lf %.13lf %.13lf %.5lf\n' % (objName, jd_tdb, ra, dec, flux)
            outfile.write(outstring)
        
            objName, jd_tdb, ra, dec, flux = convert_Mars_v2(line, +dt)
            outstring = '%s %.11lf %.13lf %.13lf %.5lf\n' % (objName, jd_tdb, ra, dec, flux)
            outfile.write(outstring)
    outfile.close()


def read_astrometry(line):
    ''' Reading the reformatted files after
        generating simulated tracklets.
    '''
    objName, jd_tdb, ra, dec, flux = line.rstrip().split()
    jd_tdb = float(jd_tdb)
    ra = float(ra)
    dec = float(dec)
    flux = float(flux)
    return objName, jd_tdb, ra, dec, flux

def transform_astrometry_Mars(in_filename, out_filename, vec, t_ref, readfunc=read_astrometry, obsCode='274', ecliptic=False):

    outfile = open(out_filename, 'w')
    outstring = '# Input file: %s\n' % (in_filename)
    outfile.write(outstring)

    outstring = '# Tangent vector (equatorial):\n# vec = %.16lf %.16lf %.16lf\n' % (vec[0], vec[1], vec[2])
    outfile.write(outstring)

    if ecliptic:
        outstring = '# ecliptic projection frame\n'
        outfile.write(outstring)
    else:
        outstring = '# equatorial projection frame\n'
        outfile.write(outstring)
    
    vec /= np.linalg.norm(vec)
    if ecliptic:
        vec = tr.equatorial_to_ecliptic(vec)

    mat = tr.xyz_to_proj_matrix(vec)
    
    et = (t_ref-spice.j2000())*24*60*60
    pos, _ = spice.spkpos('JWST', et, 'J2000', 'NONE', 'SSB')

    ref_pos = spice.convrt(pos, 'KM', 'AU')
    if ecliptic:
        ref_pos = tr.equatorial_to_ecliptic(ref_pos)
        
    ref_pos = np.dot(mat, ref_pos)
    ref_pos *= 0.0

    outstring = '# reference time:\n# t_ref = %.16lf\n' % (t_ref)
    outfile.write(outstring)
    
    outstring = '# Barycenter position in output frame:\n# bary_pos = %.16lf %.16lf %.16lf\n#\n' % (-ref_pos[0], -ref_pos[1], -ref_pos[2])
    outfile.write(outstring)
        
    outstring = "#trackletID  obsCode flux     jd_tdb             x               y               z               xe              ye              ze    \n"
    outfile.write(outstring)
    with open(in_filename, 'r') as f:
        for i, line in enumerate(f):

            if line.startswith('#'):
                continue

            objName, jd_tdb, raDeg, decDeg, flux = readfunc(line)

            xt = np.cos(decDeg*np.pi/180.)*np.cos(raDeg*np.pi/180.)
            yt = np.cos(decDeg*np.pi/180.)*np.sin(raDeg*np.pi/180.)  
            zt = np.sin(decDeg*np.pi/180.)
                
            r_target = np.array((xt, yt, zt))
                
            if ecliptic:
                r_target = tr.equatorial_to_ecliptic(r_target)
                
            r_target = np.dot(mat, r_target)
                
            xt, yt, zt = r_target
                
            theta_x = xt/zt
            theta_y = yt/zt

            et = (jd_tdb-spice.j2000())*24*60*60
                
            pos, _ = spice.spkpos('JWST', et, 'J2000', 'NONE', 'SSB')
            bary_obs = spice.convrt(pos, 'KM', 'AU')

            if ecliptic:
                bary_obs = tr.equatorial_to_ecliptic(bary_obs)
                
            bary_obs = np.dot(mat, bary_obs)
            bary_obs -= ref_pos
                
            xo, yo, zo = bary_obs
                

            outstring = "%11s %4s %9.5lf %15.12lf %15.12lf %15.12lf %15.12lf %15.12lf %15.12lf %15.12lf\n"% \
                (objName, obsCode, flux, jd_tdb, xt, yt, zt, xo, yo, zo)
                
            outfile.write(outstring)

    outfile.close()

    return -ref_pos

def transform_astrometry_Mars_v2(in_filename, out_filename, vec, t_ref, readfunc=read_astrometry, obsCode='274', ecliptic=False):

    outfile = open(out_filename, 'w')
    outstring = '# Input file: %s\n' % (in_filename)
    outfile.write(outstring)

    outstring = '# Tangent vector (equatorial):\n# vec = %.16lf %.16lf %.16lf\n' % (vec[0], vec[1], vec[2])
    outfile.write(outstring)

    if ecliptic:
        outstring = '# ecliptic projection frame\n'
        outfile.write(outstring)
    else:
        outstring = '# equatorial projection frame\n'
        outfile.write(outstring)
    
    vec /= np.linalg.norm(vec)
    if ecliptic:
        vec = tr.equatorial_to_ecliptic(vec)

    mat = tr.xyz_to_proj_matrix(vec)
    
    et = (t_ref-spice.j2000())*24*60*60
    pos, _ = spice.spkpos('JWST', et, 'J2000', 'NONE', 'SSB')

    ref_pos = spice.convrt(pos, 'KM', 'AU')
    if ecliptic:
        ref_pos = tr.equatorial_to_ecliptic(ref_pos)
        
    ref_pos = np.dot(mat, ref_pos)
    ref_pos *= 0.0

    outstring = '# reference time:\n# t_ref = %.16lf\n' % (t_ref)
    outfile.write(outstring)
    
    outstring = '# Barycenter position in output frame:\n# bary_pos = %.16lf %.16lf %.16lf\n#\n' % (-ref_pos[0], -ref_pos[1], -ref_pos[2])
    outfile.write(outstring)
        
    outstring = "#trackletID  obsCode flux     jd_tdb             x               y               z               xe              ye              ze    \n"
    outfile.write(outstring)
    with open(in_filename, 'r') as f:
        for i, line in enumerate(f):

            if line.startswith('#'):
                continue

            objName, jd_tdb, raDeg, decDeg, flux = readfunc(line)

            xt = np.cos(decDeg*np.pi/180.)*np.cos(raDeg*np.pi/180.)
            yt = np.cos(decDeg*np.pi/180.)*np.sin(raDeg*np.pi/180.)  
            zt = np.sin(decDeg*np.pi/180.)
                
            r_target = np.array((xt, yt, zt))
                
            if ecliptic:
                r_target = tr.equatorial_to_ecliptic(r_target)
                
            r_target = np.dot(mat, r_target)
                
            xt, yt, zt = r_target
                
            theta_x = xt/zt
            theta_y = yt/zt

            et = (jd_tdb-spice.j2000())*24*60*60
                
            pos, _ = spice.spkpos('JWST', et, 'J2000', 'NONE', 'SSB')
            bary_obs = spice.convrt(pos, 'KM', 'AU')

            if ecliptic:
                bary_obs = tr.equatorial_to_ecliptic(bary_obs)
                
            bary_obs = np.dot(mat, bary_obs)
            bary_obs -= ref_pos
                
            xo, yo, zo = bary_obs
                

            outstring = "%11s %4s %9.5lf %15.12lf %15.12lf %15.12lf %15.12lf %15.12lf %15.12lf %15.12lf\n"% \
                (objName, obsCode, flux, jd_tdb, xt, yt, zt, xo, yo, zo)
                
            outfile.write(outstring)

    outfile.close()

    return -ref_pos

def transform_detections_Mars(in_filename, out_filename, vec, t_ref, readfunc=read_Mars_detection, obsCode='274', ecliptic=False):

    outfile = open(out_filename, 'w')
    outstring = '# Input file: %s\n' % (in_filename)
    outfile.write(outstring)

    outstring = '# Tangent vector (equatorial):\n# vec = %.16lf %.16lf %.16lf\n' % (vec[0], vec[1], vec[2])
    outfile.write(outstring)

    if ecliptic:
        outstring = '# ecliptic projection frame\n'
        outfile.write(outstring)
    else:
        outstring = '# equatorial projection frame\n'
        outfile.write(outstring)
    
    vec /= np.linalg.norm(vec)
    if ecliptic:
        vec = tr.equatorial_to_ecliptic(vec)

    mat = tr.xyz_to_proj_matrix(vec)
    
    et = (t_ref-spice.j2000())*24*60*60
    pos, _ = spice.spkpos('JWST', et, 'J2000', 'NONE', 'SSB')

    ref_pos = spice.convrt(pos, 'KM', 'AU')
    if ecliptic:
        ref_pos = tr.equatorial_to_ecliptic(ref_pos)
        
    ref_pos = np.dot(mat, ref_pos)
    ref_pos *= 0.0

    ref_vel = ref_pos*0.0

    outstring = '# reference time:\n# t_ref = %.16lf\n' % (t_ref)
    outfile.write(outstring)
    
    outstring = '# Barycenter position in output frame:\n# bary_pos = %.16lf %.16lf %.16lf\n' % (-ref_pos[0], -ref_pos[1], -ref_pos[2])
    outfile.write(outstring)
        
    outstring = '# Barycenter velocity in output frame:\n# bary_vel = %.16lf %.16lf %.16lf\n#\n' % (-ref_vel[0], -ref_vel[1], -ref_vel[2])
    outfile.write(outstring)
        
    outstring = "#trackletID  obsCode flux     jd_tdb             x               y               z               dx              dy              dz               xe              ye              ze              vxe             vye             vze    \n"
    outfile.write(outstring)
    with open(in_filename, 'r') as f:
        f.readline()
        for i, line in enumerate(f):

            objName, jd_tdb, raDeg, decDeg, rate_ra, rate_dec, flux = readfunc(line)

            cosd = np.cos(decDeg*np.pi/180.)
            sind = np.sin(decDeg*np.pi/180.)            
            cosa = np.cos(raDeg*np.pi/180.)
            sina = np.sin(raDeg*np.pi/180.)            

            xt = cosd * cosa
            yt = cosd * sina 
            zt = sind

            dx = -zt*xt/cosd*rate_dec - yt/cosd*rate_ra
            dy = -zt*yt/cosd*rate_dec + xt/cosd*rate_ra
            dz =  cosd*rate_dec

            dx = -sind * cosa * rate_dec - cosd * sina * rate_ra
            dy = -sind * sina * rate_dec + cosd * cosa * rate_ra
            dz = cosd * rate_dec
                
            r_target = np.array((xt, yt, zt))
            rate_target = np.array((dx, dy, dz))            
                
            if ecliptic:
                r_target = tr.equatorial_to_ecliptic(r_target)
                rate_target = tr.equatorial_to_ecliptic(rate_target)                
                
            r_target = np.dot(mat, r_target)
            xt, yt, zt = r_target
                
            rate_target = np.dot(mat, rate_target)
            dx, dy, dz = rate_target            

            et = (jd_tdb-spice.j2000())*24*60*60

            state,_ = spice.spkezr('JWST', et, 'J2000', 'NONE', 'SSB')            
            #pos, _ = spice.spkpos('JWST', et, 'J2000', 'NONE', 'SSB')
            #bary_obs = spice.convrt(pos, 'KM', 'AU')

            bary_obs_pos = spice.convrt(state[0:3], 'KM', 'AU')
            bary_obs_vel = spice.convrt(state[3:6], 'KM', 'AU')*24*60*60            

            if ecliptic:
                bary_obs_pos = tr.equatorial_to_ecliptic(bary_obs_pos)
                bary_obs_vel = tr.equatorial_to_ecliptic(bary_obs_vel)                
                
            bary_obs_pos = np.dot(mat, bary_obs_pos)
            bary_obs_pos -= ref_pos
                
            xo, yo, zo = bary_obs_pos

            bary_obs_vel = np.dot(mat, bary_obs_vel)
            bary_obs_vel -= ref_vel

            vxo, vyo, vzo = bary_obs_vel

            outstring = "%11s %4s %9.5lf %15.12lf %15.12lf %15.12lf %15.12lf %15.12lf %15.12lf %15.12lf %15.12lf %15.12lf %15.12lf %15.12lf %15.12lf %15.12lf\n"% \
                (objName, obsCode, flux, jd_tdb, xt, yt, zt, dx, dy, dz, xo, yo, zo, vxo, vyo, vzo)
                
            outfile.write(outstring)

    outfile.close()

    return -ref_pos

def format_astrometry_Mars(filename, h_filename, readfunc=read_astrometry, obsCode='274'):
    with open(h_filename, 'w') as outfile:

        outstring = "#trackletID obsCode flux filter  jd_tdb    x_target      y_target      z_target       x_obs           y_obs           z_obs     \n"
        outfile.write(outstring)
        with open(filename, 'r') as f:
            for i, line in enumerate(f):
                
                objName, jd_tdb, raDeg, decDeg, flux = readfunc(line)

                xt = np.cos(decDeg*np.pi/180.)*np.cos(raDeg*np.pi/180.)
                yt = np.cos(decDeg*np.pi/180.)*np.sin(raDeg*np.pi/180.)  
                zt = np.sin(decDeg*np.pi/180.)

                RA_sig = 0.001
                Dec_sig = 0.001

                et = (jd_tdb-spice.j2000())*24*60*60
                
                pos, _ = spice.spkpos('JWST', et, 'J2000', 'NONE', 'SSB')
                bary_obs = spice.convrt(pos, 'KM', 'AU')

                xo, yo, zo = bary_obs
                filt = '_'
                
                outstring = "%11s %4s %9.5lf %s %15.12lf %15.12lf %15.12lf %15.12lf %15.12lf %15.12lf %15.12lf\n"% \
                            (objName, obsCode, flux, filt, jd_tdb, xt, yt, zt, xo, yo, zo)
                
                outfile.write(outstring)


def format_astrometry_Mars_orbfit(in_filename, out_filename, tracklet_names, readfunc=read_astrometry, obsCode='274', RA_sig=0.03, Dec_sig=0.03):
    with open(out_filename, 'w') as outfile:
        with open(in_filename, 'r') as f:
            outstring = "#trackletID jd_tdb                ra               ra_unc         dec               dec_unc           x_obs           y_obs           z_obs        obsCode flux \n"
            outfile.write(outstring)
            for i, line in enumerate(f):
                
                objName, jd_tdb, raDeg, decDeg, flux = readfunc(line)

                xt = np.cos(decDeg*np.pi/180.)*np.cos(raDeg*np.pi/180.)
                yt = np.cos(decDeg*np.pi/180.)*np.sin(raDeg*np.pi/180.)  
                zt = np.sin(decDeg*np.pi/180.)
                
                et = (jd_tdb-spice.j2000())*24*60*60
                
                pos, _ = spice.spkpos('JWST', et, 'J2000', 'NONE', 'SSB')
                bary_obs = spice.convrt(pos, 'KM', 'AU')

                xe, ye, ze = bary_obs

                #sscanf(inbuff, "%s %lf %lf %lf %lf %lf %lf %lf %lf %s", desig, &jd_tdb, &ra, &ra_unc, &dec, &dec_unc, &xe, &ye, &ze, obsCode);
                
                outstring = "%11s %15.12lf %15.12lf %10.2le %15.12lf %10.2le  %15.12lf %15.12lf %15.12lf %4s %9.5lf\n"% \
                    (tracklet_names[objName], jd_tdb, raDeg, RA_sig, decDeg, Dec_sig, xe, ye, ze, obsCode, flux)
                
                outfile.write(outstring)



def phi_xy(det, t_ref, g, gdot, acc, speed_of_light=MPC_library.Constants.speed_of_light):
    theta_x, theta_y, xe, ye, ze = det.theta_x, det.theta_y, det.xe, det.ye, det.ze
    tp = det.t - t_ref - ze/speed_of_light
    gx, gy, gz = 0.5*acc*tp*tp
    fac = 1 + gdot*tp + g*gz - g*ze
    phi_x = theta_x*fac + g*xe - g*gx
    phi_y = theta_y*fac + g*ye - g*gy
    
    return phi_x, phi_y, tp, fac, ze, gz
