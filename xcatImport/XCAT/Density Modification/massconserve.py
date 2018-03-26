# -*- coding: utf-8 -*-
"""
massconserve.py
Copyright 2013, Christopher L. Williams

Please cite Williams et al., Med. Phys. 40 (7), July 2013 if you use
this software in your work.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Requires numpy and scipy as well as the (included) mvf_lite
packages

Usage:
Run "python massconserve.py --help" to see a usage message.
 

Update history:
6/1/2013 - Initial version CW
7/8/2013 - Fixed tumor attenuation CW
"""
import numpy as np
import mvf_lite as mvf
import scipy.ndimage
import logging
import argparse

def loadxcatbin(fname,shape):
   '''
   Helper function to load an XCAT binary file.
   '''
   nx,ny,nslices=shape
   return ((np.fromfile(fname,dtype=np.float32)).reshape([nx,ny,nslices],order='F'))

def writexcatbin(arr,fname):
    '''
    Writes out an XCAT binary file
    '''
    np.array(arr.transpose(),dtype=np.float32).tofile(fname) # The XCAT is F-ordered...
    
def getphantomsize(logfile):
    '''
    Reads the phantom geometry out of an XCAT log file
    '''
    nx,ny=None,None
    startslice=None
    endslice=None
    f=open(logfile,'r')
    try:
        for line in f:
            ls=line.strip()
            if ls.startswith("array_size ="):
                    l,v=line.split('=')
                    array_size=int(v.split()[0].strip())
                    nx,ny=(array_size,array_size)
            elif ls.startswith("starting slice number"):
                l,v=line.split('=')
                startslice=int(v.split()[0].strip())
            elif ls.startswith("ending slice number"):
                l,v=line.split('=')
                endslice=int(v.split()[0].strip())
            if nx is not None and startslice is not None and endslice is not None:
                break
        for var in [startslice,endslice,nx,ny]:
            if var is None:
                raise KeyError("Missing XCAT geometry information in logfile!")      
        nslices=endslice-startslice+1
        return nx,ny,nslices
    finally:
        f.close()
        
def getphantomscale(parfile):
    '''
    Reads the phantom scale out of an XCAT parameter file
    '''
    xscale,yscale=None,None
    zscale=None
    f=open(parfile,'r')
    try:
        for line in f:
            ls=line.strip()
            if ls.startswith("pixel_width ="):
                    l,v=line.split('=')
                    xscale=float(v.split()[0].strip())
                    yscale=xscale
            elif ls.startswith("slice_width ="):
                l,v=line.split('=')
                zscale=float(v.split()[0].strip())
            if xscale is not None and zscale is not None:
                break
        for var in [xscale,zscale]:
            if var is None:
                raise KeyError("Missing XCAT geometry information in logfile!")      
        return xscale,yscale,zscale
    finally:
        f.close()
        
def gettumorgeom(parfile):
    '''
    Reads the tumor geometry out of an XCAT parameter file
    '''
    x,y,z,dia=None,None,None,None
    f=open(parfile,'r')
    try:
        for line in f:
            ls=line.strip()
            if ls.startswith("x_location ="):
                    l,v=line.split('=')
                    x=float(v.split()[0].strip())
            elif ls.startswith("y_location ="):
                    l,v=line.split('=')
                    y=float(v.split()[0].strip())
            elif ls.startswith("z_location ="):
                    l,v=line.split('=')
                    z=float(v.split()[0].strip())
            elif ls.startswith("lesn_diameter ="):
                l,v=line.split('=')
                dia=float(v.split()[0].strip())
        for var in [z,y,z,dia]:
            if var is None:
                raise KeyError("Missing Tumor geometry information in parfile!")      
        return x,y,z,dia/10.0
    finally:
        f.close()

def getlungatten(logfile):
    '''
    Reads the lung attenuation (in 1/pixel units) out of an XCAT log file
    ''' 
    nextLung=False
    f=open(logfile,'r')
    try:
        for line in f:       
            if nextLung is True:
                if line.strip().startswith("Lung"):
                    l,v=line.split('=')
                    density=float(v.split()[0].strip())
                    return density
            elif line.startswith("Linear Attenuation Coefficients (1/pixel):"):
                nextLung=True    
        raise KeyError("No lung density found!")
    finally:
        f.close()

def gettumoratten(logfile):
    '''
    Reads the tumor attenuation (in 1/cm units) out of an XCAT log file.  Note this is in 1/cm as opposed
	to 1/pixel.
    ''' 
    nextLung=False
    f=open(logfile,'r')
    try:
        for line in f:       
            if nextLung is True:
                if line.strip().startswith("Body (water)"):
                    l,v=line.split('=')
                    density=float(v.split()[0].strip())
                    return density
            elif line.startswith("Linear Attenuation Coefficients (1/cm):"):
                nextLung=True    
        raise KeyError("No tumor (Body) density found!")
    finally:
        f.close()
		
def mvftumor(p0,pp,mvec,lungatt,tpos,pscale,tdia,tatt=0.1837,dodensity=True,filt=None):
    '''
    Generates a tumor based on an MVF.  The tumor deforms with the lung
    '''
    #lungatt=xcat.getLungAtten(self.logfile)
    w=np.where(pp>0)
    ww=np.array(w).T
    vinv=mvec.vec_cinverse_multi(ww)
    opos=ww+vinv
    dist=np.sqrt( (pscale[0]*(opos[:,0]-tpos[0]))**2 + (pscale[1]*(opos[:,1]-tpos[1]))**2 + (pscale[2]*(opos[:,2]-tpos[2]))**2 )
    tw=np.where(dist <= (tdia/2.0))
    tt=np.zeros_like(pp)
    
    if dodensity is True:
        wnl=np.where(abs(p0-lungatt)>0.0001)
        dvolume=mvec.displaced_volume_map_floor()
        dvolume[wnl]=1.0  
        if filt is not None:            
            dvolume[np.isnan(dvolume)]=1.0 # get rid of edge issues
            dvolume_filt=scipy.ndimage.median_filter(dvolume,size=filt)
        else:
            dvolume_filt=dvolume
        
        ofloor=np.array(np.floor(opos[tw[0],:]),dtype=np.int)            
        tt[(w[0][tw[0]],w[1][tw[0]],w[2][tw[0]])]=(tatt*pscale[0])/dvolume_filt[(ofloor[:,0],ofloor[:,1],ofloor[:,2])]            
    else:
        tt[(w[0][tw],w[1][tw],w[2][tw])]=(tatt*pscale[0])
   
    return tt


        
def massconserve(p0,pp,mvec,lungatt,dofilter=True):
        ''' 
        Take an xcat phantom and make it mass conserving
        '''
        # This is how the volume of the reference voxels have changed
        logging.info("Calculating volume deformations")
        dvolume=mvec.displaced_volume_map_floor() 

        wnl=np.where(abs(p0-lungatt)>0.0001)
        w=np.where(abs(pp-lungatt)<=0.0001)   
        
        dvolume[wnl]=1.0  
        dvolume[np.isnan(dvolume)]=1.0 # get rid of edge issues
        if dofilter:
            logging.info("Filtering the volume deformations")
            dvolume_filt=scipy.ndimage.median_filter(dvolume,size=3)                             
        else:
            dvolume_filt=dvolume
        wa=np.array(w)
        logging.info("Inverting motion vector field")        
        opix=np.floor(wa+mvec.vec_cinverse_multi(wa.T).T).astype(np.int)

        wz=np.where(  (opix[0,:]<0) | (opix[0,:]>=pp.shape[0])
                    | (opix[1,:]<0) | (opix[1,:]>=pp.shape[1])
                    | (opix[2,:]<0) | (opix[2,:]>=pp.shape[2]) )
        if wz[0].size > 0:
            logging.warning("Lung-like regions detected at the phantom boundary!")
            logging.warning("Check to make sure that the lung is entirely in the phantom")
            logging.warning("and set si_air_flag=0 and li_air_flag=0 when you make your XCAT")
            logging.warning("Mass MAY NOT be conserved correctly in these boundary areas")
            opix[0,(opix[0,:]>=pp.shape[0])]=pp.shape[0]-1
            opix[1,(opix[1,:]>=pp.shape[1])]=pp.shape[1]-1
            opix[2,(opix[2,:]>=pp.shape[2])]=pp.shape[2]-1
            opix[0,(opix[0,:]<0)]=0
            opix[1,(opix[1,:]<0)]=0
            opix[2,(opix[2,:]<0)]=0
            
        out=pp.copy()
        logging.info("Applying mass correction to lung")
        out[w]=pp[w]/dvolume_filt[(opix[0,:],opix[1,:],opix[2,:])]
        return out
        
        
if __name__ == "__main__":
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Modify an XCAT phantom to '+
            'conserve mass.')
    parser.add_argument('referencephantomfile', metavar='reference_phantom',
                        type=str,help='Path to the reference (Frame=1) phantom attenuation file.  This frame must be the same one as the reference frame for the motion vector file.')
    parser.add_argument('phantomfile', metavar='phantom_file',type=str,
                        help='Path to the target frame phantom attenuation file')
    parser.add_argument('vectorfile', metavar='vector_file',type=str,
                        help='Path to the XCAT motion vector file between the reference frame and the target frame.  If the string "NONE" is given, it will not assume any motion (useful if you just want to insert a tumor into the phantom)')
    parser.add_argument('parfile', metavar='parfile',type=str,
                        help='Path of the XCAT parameter (".par") file used to generate the phantom')
    parser.add_argument('logfile', metavar='log_file',type=str,
                        help='file name of the XCAT log file (usually <phantomname>_log)')
    parser.add_argument('outphantom', metavar='out_phantom',type=str,
                        help='File name of the mass-corrected phantom (this file will be created)')
    parser.add_argument('-q','--quiet', action='store_true',
                        help='Suppress nonessential messages')
    parser.add_argument('-n','--nofilter', action='store_true',
                        help='Do not perform a median filter on the volume displacement (the median filter removes some edge effects due to the voxelization)')
    parser.add_argument('-b','--binary', action='store_true',
                        help='Use a binary motion vector file (this is NOT standard)')
    parser.add_argument('-c','--constant', action='store_true',
                        help='Do not conserve mass (keep the lung density constant)')
    parser.add_argument('-t','--tumor', action='store_true',
                        help='Add the tumor into the phantom at the location specified in the parameter file')
    
    args=parser.parse_args()
    
    # Should we apply a filter to the volume deformation field?
    if args.nofilter:
        filt=None
    else:
        filt=3
    
    # Display logging information?
    if not args.quiet:
        logging.basicConfig(level=logging.INFO)
    
    # Load some metadata from the XCAT logfile
    logging.info("Getting metadata from parameter file: %s"%args.parfile)    
    pscale=getphantomscale(args.parfile)
    tgeom=gettumorgeom(args.parfile)
    logging.info("Getting metadata from logfile: %s"%args.logfile)    
    s=getphantomsize(args.logfile)
    lungatt=getlungatten(args.logfile)
    tumoratt=gettumoratten(args.logfile)
    
    tpos=tgeom[0:3]
    tdia=tgeom[3]
    
    # Load the motion vector file
    logging.info("Loading MVF: %s"%args.vectorfile)

    # If no vector file is given, assume the identity
    # This is a convenient way to just put a tumor in
    if args.vectorfile.upper() == "NONE":
        mvec=mvf.mvf(None,dim=s,identity=True)
    else:
        mvec=mvf.mvf(args.vectorfile,dim=s,binary=args.binary)
    
        
    logging.info("MVF Loaded")
    logging.info("Loading phantom: %s"%args.phantomfile)
    pp=loadxcatbin(args.phantomfile,s)
    logging.info("Loading reference frame phantom: %s"%args.referencephantomfile)
    p0=loadxcatbin(args.referencephantomfile,s)
    
    if args.constant:
        out=pp
    else:
        out=massconserve(p0,pp,mvec,lungatt,dofilter=filt)
    
    if args.tumor:
        logging.info("Adding Tumor")
        tt=mvftumor(p0,pp,mvec,lungatt,tpos,pscale,tdia,tatt=tumoratt,dodensity=True,filt=filt)    
        w=np.where(tt > 0)
        out[w]=tt[w]
    
    
    logging.info("Writing output file: %s"%args.outphantom)
    writexcatbin(out,args.outphantom)
    logging.info("Done!")
                        

    