#!/usr/bin/env python

"""Process and link the JWST rate detections.
"""

# Import standard packages
import os
import sys
import random
from collections import defaultdict
from itertools import combinations

# Import third-party packages
import numpy as np
import healpy as hp
import sqlite3
import pickle
import pandas as pd
from scipy import spatial
from astropy.io import fits
import spiceypy as spice
import scipy as sp

import rebound
import assist 

# Import homegrown libraries
# Import homegrown libraries
# These are all things in the support directory
import sys
sys.path.append(".")
cwd = os.getcwd()
os.chdir('/Users/mholman/Dropbox/support')
import MPC_library # for a small number of routines
import tracklets as tr
import JWST as jw
os.chdir('kepcart_dir')
import kepcart as kc
os.chdir(cwd)

first=lambda x: x[0]
second=lambda x: x[1]
third=lambda x: x[2]
fourth=lambda x: x[3]

from collections import namedtuple
Detection = namedtuple("Detection", ("t", "theta_x", "theta_y", "xe", "ye", "ze", "flux"))
RateDetection = namedtuple("RateDetection", ("t", "theta_x", "theta_y", "dtheta_x", "dtheta_y", "xe", "ye", "ze", "vxe", "vye", "vze", "flux"))

def convert_Mars_v2(line, dt):
    ''' Converting the format in Mars's file
        good_preds_loc.txt to simulated tracklets
        based on the angular rates
        
        dt is a +/- time offset in seconds
    '''
    objName, MJD, RA, Dec, x, y, rate_ra, rate_dec, likelihood, flux, _ = line.rstrip().split()
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

def convert_implant(line, dt):
    ''' Converting the format in Mars's file
        good_preds_loc.txt to simulated tracklets
        based on the angular rates
        
        dt is a +/- time offset in seconds
    '''
    items = line.rstrip().split()
    objName = items[2]
    mag, RA, Dec, rate_ra, rate_dec, MJD, likelihood = items[10:17]
    flux = np.power(10, -0.4*(mag-27))

    #_, _, objName, MJD, RA, Dec, x, y, rate_ra, rate_dec, likelihood, flux = line.rstrip().split()

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

def make_tracklets_from_rate_detections_v3(infilename, outfilename, dt=600):
    ''' 
    Convert RA/Dec/rate detections from Mars's file KBMOD
    results to simulated tracklets with three RA/Dec detections.
    The observations are offset by a time dt in seconds.
    '''
    outfile = open(outfilename, 'w')
    with open(infilename) as file:
        file.readline()
        for i, line in enumerate(file):

            objName, jd_tdb, ra, dec, flux = convert_implant(line, -dt)
            outstring = '%s %.11lf %.13lf %.13lf %.5lf\n' % (objName, jd_tdb, ra, dec, flux)
            outfile.write(outstring)
        
            objName, jd_tdb, ra, dec, flux = convert_implant(line, 0)
            outstring = '%s %.11lf %.13lf %.13lf %.5lf\n' % (objName, jd_tdb, ra, dec, flux)
            outfile.write(outstring)
        
            objName, jd_tdb, ra, dec, flux = convert_implant(line, +dt)
            outstring = '%s %.11lf %.13lf %.13lf %.5lf\n' % (objName, jd_tdb, ra, dec, flux)
            outfile.write(outstring)
    outfile.close()

def read_Mars_detection_v3(line):
    ''' Reading the data in Mars's file good_preds_Model0_0.9_names.txt and
    converting MJD to jd_tdb.
    '''
    objName, MJD, RA, Dec, x, y, rate_ra, rate_dec, likelihood, flux = line.rstrip().split()
    ra, dec = float(RA), float(Dec)
    x, y = float(x), float(y)
    rate_ra, rate_dec = float(rate_ra), float(rate_dec)

    et = spice.str2et('JD '+ str(float(MJD)+2400000.5))
    jd_tdb = spice.j2000() + et/(24*60*60)
    
    likelihood = float(likelihood)
    flux = float(flux)
    return objName, jd_tdb, ra, dec, x, y, rate_ra, rate_dec, likelihood, flux

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

def transform_astrometry_Mars_positions(in_filename, out_filename, vec, t_ref, readfunc=read_astrometry, obsCode='274', ecliptic=False):

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

def transform_astrometry_Mars_states(in_filename, out_filename, vec, t_ref, readfunc=read_Mars_detection_v3, obsCode='274', ecliptic=False):

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
    
    ref_vel = 0.0*ref_pos

    outstring = '# reference time:\n# t_ref = %.16lf\n' % (t_ref)
    outfile.write(outstring)
    
    outstring = '# Barycenter position in output frame:\n# bary_pos = %.16lf %.16lf %.16lf\n' % (-ref_pos[0], -ref_pos[1], -ref_pos[2])
    outfile.write(outstring)
    
    outstring = '# Barycenter velocity in output frame:\n# bary_vel = %.16lf %.16lf %.16lf\n#\n' % (-ref_vel[0], -ref_vel[1], -ref_vel[2])
    outfile.write(outstring)
        
    #outstring = "#trackletID  obsCode flux     jd_tdb             x               y               z               xe              ye              ze    \n"
    
    outstring = "#trackletID  obsCode flux     jd_tdb             x               y               z               dx                  dy                  dz                   xe              ye              ze              vxe                 vye                 vze    \n"
    outfile.write(outstring)
    with open(in_filename, 'r') as f:
        for i, line in enumerate(f):

            if line.startswith('#'):
                continue

            objName, jd_tdb, raDeg, decDeg, x, y, rate_ra, rate_dec, likelihood, flux = readfunc(line)
            
            # Units are being converted from arcsec per hour to
            # radians per day.
            dra_dt  = rate_ra*(1/3600)*(np.pi/180)*24.0
            ddec_dt = rate_dec*(1/3600)*(np.pi/180)*24.0
            
            cd = np.cos(decDeg*np.pi/180.)
            sd = np.sin(decDeg*np.pi/180.)
            ca = np.cos(raDeg*np.pi/180.)
            sa = np.sin(raDeg*np.pi/180.)
            
            xt = cd*ca
            yt = cd*sa
            zt = sd
            
            vx = -sd*ca*ddec_dt - cd*sa*dra_dt
            vy = -sd*sa*ddec_dt + cd*ca*dra_dt
            vz = cd*ddec_dt
            
            r_target = np.array((xt, yt, zt))
            v_target = np.array((vx, vy, vz))
            
            if ecliptic:
                r_target = tr.equatorial_to_ecliptic(r_target)
                v_target = tr.equatorial_to_ecliptic(v_target)
                
            r_target = np.dot(mat, r_target)
            v_target = np.dot(mat, v_target)
                
            xt, yt, zt = r_target
            vx, vy, vz = v_target
                
            theta_x = xt/zt
            theta_y = yt/zt
            
            theta_x_dot = vx/zt - xt/(zt*zt)*vz
            theta_y_dot = vy/zt - yt/(zt*zt)*vz

            et = (jd_tdb-spice.j2000())*24*60*60
            
            state, _ = spice.spkezr('JWST', et, 'J2000', 'NONE', 'SSB')
            pos, vel = state[0:3], state[3:6]
            bary_obs_pos = spice.convrt(pos, 'KM', 'AU')
            bary_obs_vel = spice.convrt(vel, 'KM', 'AU')
            bary_obs_vel *= 24*60*60
            
            if ecliptic:
                bary_obs_pos = tr.equatorial_to_ecliptic(bary_obs_pos)
                bary_obs_vel = tr.equatorial_to_ecliptic(bary_obs_vel)
                
            bary_obs_pos = np.dot(mat, bary_obs_pos)
            bary_obs_pos -= ref_pos
            
            bary_obs_vel = np.dot(mat, bary_obs_vel)
            bary_obs_vel -= ref_vel
                
            xo, yo, zo = bary_obs_pos
            vxo, vyo, vzo = bary_obs_vel
                

            outstring = "%11s %4s %9.5lf %15.12lf %15.12lf %15.12lf %15.12lf %19.12le %19.12le %19.12le %15.12lf %15.12lf %15.12lf %19.12le %19.12le %19.12le\n"% \
                (objName, obsCode, flux, jd_tdb, xt, yt, zt, vx, vy, vz, xo, yo, zo, vxo, vyo, vzo)
                
            outfile.write(outstring)

    outfile.close()

    return -ref_pos

def parse_detection(line):
    trackletID, obsCode, flux, jd_tdb, x, y, z, xe, ye, ze = line.rstrip().split()
    jd_tdb = float(jd_tdb)
    x, y, z = float(x), float(y), float(z)
    flux = float(flux)
    xe, ye, ze = float(xe), float(ye), float(ze)
    return trackletID, obsCode, (jd_tdb, x, y, z, xe, ye, ze, flux)

def parse_rate_detection(line):
    trackletID, obsCode, flux, jd_tdb, x, y, z, dx, dy, dz, xe, ye, ze, vxe, vye, vze = line.rstrip().split()
    jd_tdb = float(jd_tdb)
    x, y, z = float(x), float(y), float(z)
    dx, dy, dz = float(dx), float(dy), float(dz)
    flux = float(flux)
    xe, ye, ze = float(xe), float(ye), float(ze)
    vxe, vye, vze = float(vxe), float(vye), float(vze)
    return trackletID, obsCode, (jd_tdb, x, y, z, dx, dy, dz, xe, ye, ze, vxe, vye, vze, flux)

def solve_rate_detection(GMtotal, t_ref, g, gdot, obs, speed_of_light=MPC_library.Constants.speed_of_light):

    t_emit = (obs.t-t_ref - obs.ze/speed_of_light)

    acc_z = -GMtotal*g*g
    
    f = 1.0 + gdot*t_emit + 0.5*g*acc_z*t_emit*t_emit - g*obs.ze
    fdot = gdot + g*acc_z*t_emit - g*obs.vze
    
    adot = obs.dtheta_x*f + obs.theta_x*fdot + g*obs.vxe
    bdot = obs.dtheta_y*f + obs.theta_y*fdot + g*obs.vye   
    
    alpha = -adot*t_emit + obs.theta_x*f + g*obs.xe
    beta  = -bdot*t_emit + obs.theta_y*f + g*obs.ye
    
    return alpha, adot, beta, bdot, t_emit, obs.flux

def fit_tracklet_rms(t_ref, g, gdot, v, GM=MPC_library.Constants.GMsun, speed_of_light=MPC_library.Constants.speed_of_light):
    # Here's a version that incorporates radial gravitational
    # acceleration
    
    t_emit = [(obs.t-t_ref - obs.ze/speed_of_light) for obs in v]
    flux = [obs.flux for obs in v]        
    acc_z = -GM*g*g
    fac =[(1.0 + gdot*t + 0.5*g*acc_z*t*t - g*obs.ze) for obs, t in zip(v, t_emit)]
                
    A = np.vstack([t_emit, np.ones(len(t_emit))]).T 
    
    x = [obs.theta_x*f + g*obs.xe for obs, f in zip(v, fac)]                 
    (mx, cx), res_x = np.linalg.lstsq(A, x, rcond=None)[0:2]
    res_x = np.sqrt(res_x[0]/len(v))
            
    y = [obs.theta_y*f + g*obs.ye for obs, f in zip(v, fac)]                 
    (my, cy), res_y = np.linalg.lstsq(A, y, rcond=None)[0:2]
    res_y = np.sqrt(res_y[0]/len(v))
    
    return (cx, mx, cy, my, res_x, res_y, t_emit[0], flux[0])

def make_graph(visit_trees, visit_labels, rad):
    graph=defaultdict(list)
    for (i, j) in combinations(visit_trees.keys(), 2):
        matches = visit_trees[i].query_ball_tree(visit_trees[j], rad)
        matches_dict={k:match for (k, match) in enumerate(matches) if len(match)>0}
        for k, matches in matches_dict.items():
            for m in matches:
                graph[visit_labels[i][k]].append(visit_labels[j][m])
    return graph

def get_paths(node, graph):
    if node not in graph:
        return [[node]]
    else:
        paths=[]
        for neighbor in graph[node]:
            for path in get_paths(neighbor, graph):
                paths.append(path)
            paths.append([])
        return [[node]+path for path in paths]

def get_all_paths(graph):
    paths=[]
    for node in graph:
        for path in get_paths(node, graph):
            paths.append(path)
    return paths

p = [0.0, 0.0, 0.0, 0.0, 0.0]
def fit_cluster(cluster_key, tracklets_dict, t_ref, gdot=0.0):
    v=[]
    for k in cluster_key:
        for t in tracklets_dict[k]:
            v.append(t)
    soln = sp.optimize.least_squares(tr.residuals, p, jac=tr.partials, args=(v, t_ref, gdot), method='lm')
    rms = np.sqrt(soln.cost*2/(2*len(v)-len(p)))*206265
    return rms, soln

def visit_number(t):
    if t<-2:
        return 0
    if t<1:
        return 1
    else: 
        return 2

def format_cluster_orbfit(cluster_key, tracklet_lines, out_filename, readfunc = read_astrometry, obsCode='274', RA_sig=0.03, Dec_sig=0.03):

    outfile = open(out_filename, 'w')
    for k in cluster_key:
        for line in tracklet_lines[k]:
            objName, jd_tdb, raDeg, decDeg, flux = readfunc(line)
            et = (jd_tdb-spice.j2000())*24*60*60
                
            pos, _ = spice.spkpos('JWST', et, 'J2000', 'NONE', 'SSB')
            bary_obs = spice.convrt(pos, 'KM', 'AU')

            xe, ye, ze = bary_obs

            outstring = "%11s %15.12lf %15.12lf %10.2le %15.12lf %10.2le  %15.12lf %15.12lf %15.12lf %4s %9.5lf\n"% \
                (objName, jd_tdb, raDeg, RA_sig, decDeg, Dec_sig, xe, ye, ze, obsCode, flux)
            outfile.write(outstring)
    outfile.close()
                
def main(in_filename):

    # Create some output file names
    out_filename =  in_filename.replace('names.txt', 'rate_detections.tng')
    tracklets_filename =  in_filename.replace('names.txt', 'tracklets.txt')
    tracklets_tangent_filename =  in_filename.replace('names.txt', 'tracklets.tng')

    orbit_fits_filename =  in_filename.replace('names.txt', 'orbit_fits.txt')    

    # Load a few spice kernels
    dir_path = '/Users/mholman/Dropbox/support/'
    spice.furnsh(dir_path+'/kernels/MetaK_jwst.txt')

    # Load ephemeris files for ASSIST.
    ephem = assist.Ephem("/Users/mholman/assist/data/linux_p1550p2650.440", "/Users/mholman/assist/data/sb441-n16.bsp")
    #GMsun = ephem.get_particle('Sun', 0).m

    # Calculate GMtotal
    GMtotal = 0
    GMs = [ephem.get_particle(i, 0).m for i in range(27)]
    GMs_sorted = sorted(GMs)
    for GM in GMs_sorted:
        GMtotal += GM

    # Unit vector to JWST survey center, in equatorial coordinates
    vec = np.array((-0.8556287151668741, -0.4820832782601389, -0.18839314183602954))
    vec /= np.linalg.norm(vec)

    # reference time
    t_ref = 2459974.5

    # Synthetize 3-detection tracklets from Mars's rate detections
    make_tracklets_from_rate_detections_v2(in_filename, tracklets_filename, dt=600)

    tracklet_lines = defaultdict(list)
    with open(tracklets_filename) as infile:
        for line in infile:
            objName, jd_tdb, raDeg, decDeg, flux = read_astrometry(line)
            tracklet_lines[objName].append(line)

    # Transform the 3-detection tracklets to tangent plane coordinates 
    transform_astrometry_Mars_positions(tracklets_filename, tracklets_tangent_filename, vec, t_ref, readfunc=read_astrometry, obsCode='274', ecliptic=True)

    # Read 3-detection tracklets into a dictionary
    known_tracklets=set()
    tracklets = defaultdict(list)
    with open(tracklets_tangent_filename) as file:
        for line in file:
            if line.strip().startswith('#'):
                continue
            trackletID, obsCode, data = parse_detection(line)
            if trackletID not in known_tracklets:
                jd_tdb, x, y, z, xe, ye, ze, flux = data
                det = Detection(jd_tdb, x/z, y/z, xe, ye, ze, flux)
                tracklets[trackletID].append(det)

    # Transform the rate detections, with velocities, to tangent plane coordinates
    #transform_astrometry_Mars_states(in_filename, out_filename, vec, t_ref, readfunc=read_Mars_detection_v3, obsCode='274', ecliptic=True)

    gs = np.round(np.linspace(0., 0.06, 6*20+1), 7)
    #gdots = np.linspace(-1, 1, 9)*1e-4
    gdots = np.round(np.linspace(-1, 1, 9)*1e-4, 7)    
    #gdots = [0.0]
    g_gdots = [(x,y) for x in gs for y in gdots]

    master_results = {}
    for g_gdot in g_gdots:
        print(g_gdot)
        g, gdot = g_gdot
        results = []
        for i, (k_i, tracklet_i) in enumerate(tracklets.items()):
            v = tracklet_i
            a, adot, b, bdot, res_x, res_y, t_emit, flux = fit_tracklet_rms(t_ref, g, gdot, v)
            result = (k_i, a, adot, b, bdot, t_emit, flux)
            results.append(result)
        master_results[g_gdot] = results

    with open('master_implant.pkl', 'wb') as outfile:
        pickle.dump(master_results, outfile)

    '''
    # Read the processed rate detections
    rate_detections = {}
    with open(out_filename) as file:
        for line in file:
            if line.strip().startswith('#'):
                continue
            objID, obsCode, data = parse_rate_detection(line)
            jd_tdb, x, y, z, dx, dy, dz, xe, ye, ze, vxe, vye, vze, flux = data
        
            theta_x = x/z
            theta_y = y/z

            theta_x_dot = dx/z - x/(z*z)*dz
            theta_y_dot = dy/z - y/(z*z)*dz
        
            rate_det = RateDetection(jd_tdb, theta_x, theta_y, theta_x_dot, theta_y_dot, xe, ye, ze, vxe, vye, vze, flux)
        
            rate_detections[objID] = rate_det
    
    master_results = {}
    for g_gdot in g_gdots:
        g, gdot = g_gdot
        results = []
        for objID in sorted(rate_detections):
            obs = rate_detections[objID]
            alpha, adot, beta, bdot, t_emit, flux = solve_rate_detection(GMtotal, t_ref, g, gdot, obs)
            result = (objID, alpha, adot, beta, bdot, t_emit, flux)
            results.append(result)
        master_results[g_gdot] = results

    # At this point, start clustering arrows

    with open('master_rates.pkl', 'wb') as outfile:
        pickle.dump(master_results, outfile)

    '''

    dt = 5
    rad = 3.0e-4
    graphs={}
    for (g, gdot), arrows in master_results.items():

        #Separate arrows by visit number
        visit_arrows=defaultdict(list)
        for arrow in arrows:
            vn = visit_number(arrow[5])
            visit_arrows[vn].append(arrow)

        # Make a separate tree for each visit,
        # along with dictionaries that retain
        # the names/labels by index.
        visit_trees={}
        visit_labels={}
        for vn, arrows in visit_arrows.items():
            i = 0
            label_dict={}

            points=[]
            for k, cx, mx, cy, my, t, flux in arrows:
                label_dict[i] = k
                points.append([cx, mx*dt, cy, my*dt])
                i +=1
            
            # Skip if there are no arrows
            if len(points)==0:
                continue
            
            points=np.array(points)
            tree = sp.spatial.cKDTree(points)
            visit_trees[vn] = tree
            visit_labels[vn] = label_dict
        
        # Make the corresponding graph for this
        # g, gdot pair
        graphs[g, gdot]=make_graph(visit_trees, visit_labels, rad)

    all_triples = {}
    all_pairs = {}
    for k, grph in graphs.items():
        paths = get_all_paths(grph)
        triples = [m for m in sorted(list(set([tuple(path) for path in paths if len(path)>1]))) if len(m)==3]
        pairs = [m for m in sorted(list(set([tuple(path) for path in paths if len(path)>1]))) if len(m)==2]
        all_triples[k]=triples
        all_pairs[k]=pairs

    triples_set = set()
    for k, triples in all_triples.items():
        for triple in triples:
            triples_set.add((triple))

    pairs_set = set()
    for k, pair in all_pairs.items():
        for pair in pairs:
            pairs_set.add((pair))

    print(len(triples_set), len(pairs_set))


    rms_values = []
    rms_thresh = 0.2
    mag_range_thresh = 3.0
    tracklet_clusters = defaultdict(list)
    i=0
    all_triples = []
    orbit_fits_file = open(orbit_fits_filename, 'w')

    outstring  = '# id0    id1    id2   rms   sig_mag      dist    flux0  flux1  flux2\n'
    orbit_fits_file.write(outstring)
    for triple in triples_set:
        rms, soln = fit_cluster(triple, tracklets, t_ref)
        rms_values.append(rms)
        flux = np.array([tracklets[d][0].flux for d in triple])
        inst_mag = -2.5*np.log10(flux)
        mag_range = inst_mag.max()-inst_mag.min()
        cluster_filename = 'astrometry/'+'_'.join(triple)+'.ast'
        if rms<rms_thresh and mag_range<mag_range_thresh and 1/soln.x[4]>0:
            all_triples.append(triple)
            ids = tuple([d for d in triple])
            fluxes = tuple([f for f in flux])
            outstring  = '%s %s %s' % (ids[0], ids[1], ids[2])
            outstring += '%6.3f %6.3f %12.2f   ' % (rms, inst_mag.std(), 1/soln.x[4])
            outstring += '%6.2f %6.2f %6.2f' % (fluxes[0], fluxes[1], fluxes[2])
            outstring += '\n'
            orbit_fits_file.write(outstring)

            format_cluster_orbfit(triple, tracklet_lines, cluster_filename)
            i += 1
            for d in triple:
                tracklet_clusters[d].append((triple, rms, inst_mag.std()))

    orbit_fits_file.close()
    #rms_values = np.array(rms_values)
        
    print('finished')
    return

if __name__ == "__main__":

    if len(sys.argv)==2:    
        in_filename  = sys.argv[1]
    else:
        exit(-1)

    main(in_filename)

