# -*- Copyright (c) 2018, Bethany Ann Ludwig, All rights reserved. -*-
"""
NAME:
    Prune
PURPOSE:
    To remove everything outside of a 22" radius from E0102
"""
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np

def Prune(data,hdr):
	# Constants
	c_dec = -72.03125; c_ra = 16.00875
	d = 58 #kpc
	arcs = 22 #"

	w = WCS(hdr)
	# For each pixel we need to know it's ra/dec coordinates
	a,b = np.shape(data)
	row = np.arange(0,a); col = np.arange(0,b)
	row,col=np.meshgrid(row,col)
	row=row.flatten(); col=col.flatten()

	all_ra, all_dec = w.all_pix2world(col,row,1)
	# Numbers here are from Karin's paper.
	c1 = SkyCoord(c_ra*u.deg, c_dec*u.deg, distance=d*u.kpc, frame='icrs')
	c2 = SkyCoord(all_ra*u.deg, all_dec*u.deg, distance=d*u.kpc, frame='icrs')

	sep = c1.separation_3d(c2)
	radius = d*u.kpc*arcs/206265

	look =np.where(sep > radius)
	data[row[look],col[look]] = np.nan

	return data

def Annulus(data,hdr,thickness):
	# Constants
	c_dec = -72.03125; c_ra = 16.00875
	d = 58 #kpc
	arcs = 22 #"
	
	w = WCS(hdr)
	# For each pixel we need to know it's ra/dec coordinates
	a,b = np.shape(data)
	row = np.arange(0,a); col = np.arange(0,b)
	row,col=np.meshgrid(row,col)
	row=row.flatten(); col=col.flatten()

	all_ra, all_dec = w.all_pix2world(col,row,1)
	# Numbers here are from Karin's paper.
	c1 = SkyCoord(c_ra*u.deg, c_dec*u.deg, distance=d*u.kpc, frame='icrs')
	c2 = SkyCoord(all_ra*u.deg, all_dec*u.deg, distance=d*u.kpc, frame='icrs')

	sep = c1.separation_3d(c2)
	radius = d*u.kpc*arcs/206265
	ring = d*u.kpc*thickness/206265

	# Clear out SNR 
	look =np.where(sep < radius)

	data[row[look],col[look]] = np.nan

	# Clear out everything outside annulus 
	look =np.where(sep > ring)
	data[row[look],col[look]] = np.nan

	return data

def SelectRegionalData(data,hdr,origin_ra,origin_dec,radius):
	# Constants
	w = WCS(hdr)
	c_ra = origin_ra; c_dec = origin_dec
	d = 58 #kpc
	arcs = radius
	data = np.copy(data)

	# For each pixel we need to know it's ra/dec coordinates
	a,b = np.shape(data)
	row = np.arange(0,a); col = np.arange(0,b)
	row,col=np.meshgrid(row,col)
	row=row.flatten(); col=col.flatten()

	all_ra, all_dec = w.all_pix2world(col,row,1)
	
	c1 = SkyCoord(c_ra*u.deg, c_dec*u.deg, distance=d*u.kpc, frame='icrs')
	c2 = SkyCoord(all_ra*u.deg, all_dec*u.deg, distance=d*u.kpc, frame='icrs')

	sep = c1.separation_3d(c2)
	radius = d*u.kpc*arcs/206265

	look =np.where(sep > radius)
	data[row[look],col[look]] = np.nan

	return data