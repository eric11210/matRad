# -*- coding: utf-8 -*-
"""
mvf_lite.py
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


Requires numpy and scipy.  It uses the scipy weave library to
compile C code, and you additionally need to specify an
installed compiler to use (default is GCC)

This contains a helper library to deal with the XCAT motion vector files

Update history:
6/1/2013 - Initial version CW
"""

import numpy as np
import scipy.weave.ext_tools

#WCOMPILER='msvc'
WCOMPILER='gcc'

def build_cfunc():
    '''
    This is a helper function to force the building of the weave c functions
    '''
    mod=scipy.weave.ext_tools.ext_module('cmvec',compiler=WCOMPILER)
    
    ax,ay,az,bx,by,bz,cx,cy,cz=1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0
    cmd=\
    """
    return_val=fabs( (ax*by*cz)-(ax*bz*cy)+(ay*bz*cx)-(ay*bx*cz)+(az*bx*cy)-(az*by*cx))/6;
    """
    fn = scipy.weave.ext_tools.ext_function('pyramidvol',cmd,['ax','ay','az','bx','by','bz','cx','cy','cz'])
    mod.add_function(fn)
    
    nx,ny,nz=1,1,1
    marr=np.array((1,1,1),dtype=np.float64)
    out=np.array((1,1,1),dtype=np.float64)
    cmd=\
    """
    int x,y,z;
    int dx,dy,dz;
    int xp,yp,zp;
    double total;
    
    double ax,ay,az,bx,by,bz,cx,cy,cz;
       
    for (x=0;x<int(nx)-1;x++) 
        for (y=0;y<int(ny)-1;y++) 
            for (z=0;z<int(nz)-1;z++){
                total=0;
                xp=x+1;
                yp=y+1;
                zp=z+1;
                
                ax=*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+z*marr_array->strides[2]+0*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+0*marr_array->strides[3]); //a=P0-P3
                ay=*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+z*marr_array->strides[2]+1*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+1*marr_array->strides[3]);
                az=*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+z*marr_array->strides[2]+2*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+2*marr_array->strides[3]);
                
                //bx=marr(xp,y,z,0)-marr(x,y,zp,0); //b=P1-P3
                //by=marr(xp,y,z,1)-marr(x,y,zp,1);
                //bz=marr(xp,y,z,2)-marr(x,y,zp,2);
                bx=*(double*)(marr_array->data+xp*marr_array->strides[0]+y*marr_array->strides[1]+z*marr_array->strides[2]+0*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+0*marr_array->strides[3]); 
                by=*(double*)(marr_array->data+xp*marr_array->strides[0]+y*marr_array->strides[1]+z*marr_array->strides[2]+1*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+1*marr_array->strides[3]); 
                bz=*(double*)(marr_array->data+xp*marr_array->strides[0]+y*marr_array->strides[1]+z*marr_array->strides[2]+2*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+2*marr_array->strides[3]); 
                
                
                //cx=marr(x,yp,z,0)-marr(x,y,zp,0); //c=P2-P3
                //cy=marr(x,yp,z,1)-marr(x,y,zp,1);
                //cz=marr(x,yp,z,2)-marr(x,y,zp,2);   
                cx=*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+z*marr_array->strides[2]+0*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+0*marr_array->strides[3]); 
                cy=*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+z*marr_array->strides[2]+1*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+1*marr_array->strides[3]); 
                cz=*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+z*marr_array->strides[2]+2*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+2*marr_array->strides[3]); 
                
                
                total+=fabs( (ax*by*cz)-(ax*bz*cy)+(ay*bz*cx)-(ay*bx*cz)+(az*bx*cy)-(az*by*cx))/6;
                
                //ax=marr(xp,yp,z,0)-marr(x,y,zp,0); //a=P5-P3
                //ay=marr(xp,yp,z,1)-marr(x,y,zp,1);
                //az=marr(xp,yp,z,2)-marr(x,y,zp,2);
                ax=*(double*)(marr_array->data+xp*marr_array->strides[0]+yp*marr_array->strides[1]+z*marr_array->strides[2]+0*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+0*marr_array->strides[3]); 
                ay=*(double*)(marr_array->data+xp*marr_array->strides[0]+yp*marr_array->strides[1]+z*marr_array->strides[2]+1*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+1*marr_array->strides[3]); 
                az=*(double*)(marr_array->data+xp*marr_array->strides[0]+yp*marr_array->strides[1]+z*marr_array->strides[2]+2*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+2*marr_array->strides[3]); 
                
                //b=P1-P3 already set
                //c=P2-P3 already set
                total+=fabs( (ax*by*cz)-(ax*bz*cy)+(ay*bz*cx)-(ay*bx*cz)+(az*bx*cy)-(az*by*cx))/6;
                
                //a=P5-P3
                //bx=marr(x,yp,zp,0)-marr(x,y,zp,0); //b=P4-P3
                //by=marr(x,yp,zp,1)-marr(x,y,zp,1);
                //bz=marr(x,yp,zp,2)-marr(x,y,zp,2);
                bx=*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+0*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+0*marr_array->strides[3]); 
                by=*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+1*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+1*marr_array->strides[3]); 
                bz=*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+2*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+2*marr_array->strides[3]); 
                //c=P2-P3
                total+=fabs( (ax*by*cz)-(ax*bz*cy)+(ay*bz*cx)-(ay*bx*cz)+(az*bx*cy)-(az*by*cx))/6;
                
                //ax=marr(xp,yp,z,0)-marr(x,yp,zp,0); //a=P5-P4
                //ay=marr(xp,yp,z,1)-marr(x,yp,zp,1);
                //az=marr(xp,yp,z,2)-marr(x,yp,zp,2);
                ax=*(double*)(marr_array->data+xp*marr_array->strides[0]+yp*marr_array->strides[1]+z*marr_array->strides[2]+0*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+0*marr_array->strides[3]); 
                ay=*(double*)(marr_array->data+xp*marr_array->strides[0]+yp*marr_array->strides[1]+z*marr_array->strides[2]+1*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+1*marr_array->strides[3]); 
                az=*(double*)(marr_array->data+xp*marr_array->strides[0]+yp*marr_array->strides[1]+z*marr_array->strides[2]+2*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+2*marr_array->strides[3]);                                             
                
                
                //bx=marr(xp,y,zp,0)-marr(x,yp,zp,0); //b=P6-P4
                //by=marr(xp,y,zp,1)-marr(x,yp,zp,1);
                //bz=marr(xp,y,zp,2)-marr(x,yp,zp,2);
                bx=*(double*)(marr_array->data+xp*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+0*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+0*marr_array->strides[3]); 
                by=*(double*)(marr_array->data+xp*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+1*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+1*marr_array->strides[3]); 
                bz=*(double*)(marr_array->data+xp*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+2*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+2*marr_array->strides[3]);                                             
                                                              
                //cx=marr(xp,yp,zp,0)-marr(x,yp,zp,0); //c=P7-P4
                //cy=marr(xp,yp,zp,1)-marr(x,yp,zp,1);
                //cz=marr(xp,yp,zp,2)-marr(x,yp,zp,2);                                                              
                cx=*(double*)(marr_array->data+xp*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+0*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+0*marr_array->strides[3]); 
                cy=*(double*)(marr_array->data+xp*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+1*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+1*marr_array->strides[3]); 
                cz=*(double*)(marr_array->data+xp*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+2*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+2*marr_array->strides[3]);                                             
                
                total+=fabs( (ax*by*cz)-(ax*bz*cy)+(ay*bz*cx)-(ay*bx*cz)+(az*bx*cy)-(az*by*cx))/6;
                
                //a=P5-P4
                //b=P6-P4
                //cx=marr(xp,y,z,0)-marr(x,yp,zp,0); //c=P1-P4
                //cy=marr(xp,y,z,1)-marr(x,yp,zp,1);
                //cz=marr(xp,y,z,2)-marr(x,yp,zp,2);  
                cx=*(double*)(marr_array->data+xp*marr_array->strides[0]+y*marr_array->strides[1]+z*marr_array->strides[2]+0*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+0*marr_array->strides[3]); 
                cy=*(double*)(marr_array->data+xp*marr_array->strides[0]+y*marr_array->strides[1]+z*marr_array->strides[2]+1*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+1*marr_array->strides[3]); 
                cz=*(double*)(marr_array->data+xp*marr_array->strides[0]+y*marr_array->strides[1]+z*marr_array->strides[2]+2*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+2*marr_array->strides[3]);                                             
                
                total+=fabs( (ax*by*cz)-(ax*bz*cy)+(ay*bz*cx)-(ay*bx*cz)+(az*bx*cy)-(az*by*cx))/6;
                
                //ax=marr(x,y,zp,0)-marr(x,yp,zp,0); //a=P3-P4
                //ay=marr(x,y,zp,1)-marr(x,yp,zp,1);
                //az=marr(x,y,zp,2)-marr(x,yp,zp,2);
                ax=*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+0*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+0*marr_array->strides[3]); 
                ay=*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+1*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+1*marr_array->strides[3]); 
                az=*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+2*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+2*marr_array->strides[3]);                                             
                
                //c=P1-P4
                //b=P6-P4
            
                
                total+=fabs( (ax*by*cz)-(ax*bz*cy)+(ay*bz*cx)-(ay*bx*cz)+(az*bx*cy)-(az*by*cx))/6;

                //out(x,y,z)=total;
                *(double *)(out_array->data+x*out_array->strides[0]+y*out_array->strides[1]+z*out_array->strides[2])=total;

            }
          
    """
    fn = scipy.weave.ext_tools.ext_function('displaced_volume_map_floor',cmd,['marr','out','nx','ny','nz'])
    mod.add_function(fn)
    
    npos=1
    maxiter=1
    nx,ny,nz=1,1,1
    parr=np.array((1,1,1,1),dtype=np.float64)
    cmd=\
    """
    int i,p;        
    double x,y,z,f_x,f_y,f_z,pos_x,pos_y,pos_z;
    int xf,yf,zf,xp,yp,zp;
    int warned;    
    double weight;    
      
    for (p=0; p < int(npos);p++){
        //x=guesses(p,0);
        //y=guesses(p,1);
        //z=guesses(p,2);
        x=(*(double*)(parr_array->data+p*parr_array->strides[0]+0*parr_array->strides[1]));
        y=(*(double*)(parr_array->data+p*parr_array->strides[0]+1*parr_array->strides[1]));
        z=(*(double*)(parr_array->data+p*parr_array->strides[0]+2*parr_array->strides[1]));   
        pos_x=x;
        pos_y=y;
        pos_z=z;
        
        for (i=0;i<maxiter;i++){
        
            /* Calculate forward vector using linear interpolation */
    
            xf=(int)floor(x);
            yf=(int)floor(y);
            zf=(int)floor(z);      
            f_x=0;
            f_y=0;
            f_z=0;
            for (xp=0;xp<=1;xp++)
                for (yp=0;yp<=1;yp++)
                    for (zp=0;zp<=1;zp++) {
                        weight=1.0;
                        if (xp==0)
                            weight*=(1-(x-xf));
                        else
                            weight*=(x-xf);
                        if (yp==0)
                            weight*=(1-(y-yf));
                        else
                            weight*=(y-yf);
                        if (zp==0)
                            weight*=(1-(z-zf));
                        else
                            weight*=(z-zf);
                           
                        if (! ( (xf+xp >= marr_array->dimensions[0]) || (xf+xp <0) || (yf+yp >= marr_array->dimensions[1]) || (yf+yp <0) || (zf+zp >= marr_array->dimensions[2]) || (zf+zp <0) )) {
                            f_x+=(*(double *)(marr_array->data+(xf+xp)*marr_array->strides[0]+(yf+yp)*marr_array->strides[1]+(zf+zp)*marr_array->strides[2]+0*marr_array->strides[3])-(xf+xp))*weight;
                            f_y+=(*(double *)(marr_array->data+(xf+xp)*marr_array->strides[0]+(yf+yp)*marr_array->strides[1]+(zf+zp)*marr_array->strides[2]+1*marr_array->strides[3])-(yf+yp))*weight;
                            f_z+=(*(double *)(marr_array->data+(xf+xp)*marr_array->strides[0]+(yf+yp)*marr_array->strides[1]+(zf+zp)*marr_array->strides[2]+2*marr_array->strides[3])-(zf+zp))*weight;
                        }
                    }

            /* Okay, we've got the forward vector*/
        
            x=pos_x-f_x;
            y=pos_y-f_y;
            z=pos_z-f_z;
            
        }
        *(double*)(out_array->data+p*out_array->strides[0]+0*out_array->strides[1])=-1.0*f_x;
        *(double*)(out_array->data+p*out_array->strides[0]+1*out_array->strides[1])=-1.0*f_y;
        *(double*)(out_array->data+p*out_array->strides[0]+2*out_array->strides[1])=-1.0*f_z;
    }
    """
    fn = scipy.weave.ext_tools.ext_function('vec_cinverse_multi',cmd,['marr','out','maxiter','parr','npos'])
    mod.add_function(fn)

    mod.compile()

try:
    import cmvec
except:
    print "Building..."
    build_cfunc()
    import cmvec

def _pyramidvol(pt1,pt2,pt3,pt4):
    '''
    Helper function to calculate the volume of a pyramid (tetrahedron).  The four verticies are
    specified by p1, p2, p3 and p4, which are 3-element numpy arrays.
    '''
    
    ax=float(pt1[0]-pt4[0])
    ay=float(pt1[1]-pt4[1])
    az=float(pt1[2]-pt4[2])
    
    bx=float(pt2[0]-pt4[0])
    by=float(pt2[1]-pt4[1])
    bz=float(pt2[2]-pt4[2])
    
    cx=float(pt3[0]-pt4[0])
    cy=float(pt3[1]-pt4[1])
    cz=float(pt3[2]-pt4[2])

    return cmvec.pyramidvol(ax,ay,az,bx,by,bz,cx,cy,cz)

class mvf():
    '''
    The MVF class is used to manage XCAT motion vector files.
    '''
    
    mv=None    
    
    def __init__(self,mvfname,dim=(256,256,120),arr=None,binary=False,identity=False):
        '''
        Initializes an MVF.  All unspecified pixels in the XCAT MVF file will 
        be assumed to have no motion.
        '''
        
        if arr is not None:
            self.mv=arr
        elif identity is True:
            self.mv=np.recarray((dim[0],dim[1],dim[2]),dtype=[('x',np.float),('y',np.float),('z',np.float)])
            for i in range(dim[0]):
                self.mv[i,:,:].x=i
            for i in range(dim[1]):
                self.mv[:,i,:].y=i
            for i in range(dim[2]):
                self.mv[:,:,i].z=i
        elif binary is True:
            
            self.mv=np.recarray((dim[0],dim[1],dim[2]),dtype=[('x',np.float),('y',np.float),('z',np.float)])
            for i in range(dim[0]):
                self.mv[i,:,:].x=i
            for i in range(dim[1]):
                self.mv[:,i,:].y=i
            for i in range(dim[2]):
                self.mv[:,:,i].z=i

            dt=np.dtype([('x0',np.int16),('y0',np.int16),('z0',np.int16),('x1',np.float64),('y1',np.float64),('z1',np.float64)])
            t=np.fromfile(mvfname,dtype=dt)
            self.mv.x[(t['x0'],t['y0'],t['z0'])]=t['x1']
            self.mv.y[(t['x0'],t['y0'],t['z0'])]=t['y1']
            self.mv.z[(t['x0'],t['y0'],t['z0'])]=t['z1']
        else:
            f=open(mvfname)
            line=f.readline() # Skip first line
            line=f.readline() # Read header line
            #frame=int(line.split()[7]) #read the frame that the MVF is for
        
            self.mv=np.recarray((dim[0],dim[1],dim[2]),dtype=[('x',np.float),('y',np.float),('z',np.float)])
            for i in range(dim[0]):
                self.mv[i,:,:].x=i
            for i in range(dim[1]):
                self.mv[:,i,:].y=i
            for i in range(dim[2]):
                self.mv[:,:,i].z=i
                
        
            for line in f:
                s=line.split()
                x0=int(s[2])
                y0=int(s[3])
                z0=int(s[4])
        
                self.mv[x0,y0,z0].x=float(s[6])        
                self.mv[x0,y0,z0].y=float(s[7])
                self.mv[x0,y0,z0].z=float(s[8])  
 
    
    def displaced_volume_map_floor(self):
        '''
        Quickly (using weave) computes the displaced volume average for all
        voxels in an MVF.  Returns a numpy array of the relative volume changes
        of the voxels.  The voxel volume is 
        '''

        marr=self.mv.view((np.float,3))
        
        out=np.empty(self.mv.shape)
        nx=marr.shape[0]
        ny=marr.shape[1]
        nz=marr.shape[2]
        cmvec.displaced_volume_map_floor(marr,out,nx,ny,nz)
        out[nx-1,:,:]=np.NaN
        out[:,ny-1,:]=np.NaN
        out[:,:,nz-1]=np.NaN
        return out

    def vec_cinverse_multi(self,positions,maxiter=15):
        '''
        Inverts the motion vector field for each position in the position list (list or array of 3-element tuples)
        '''
      
        parr=np.array(positions,dtype=np.float)
        npos=parr.shape[0]                 
        marr=self.mv.view((np.float,3))
        out=np.zeros_like(parr,dtype=np.float)
        cmvec.vec_cinverse_multi(marr,out,maxiter,parr,npos)
        return out