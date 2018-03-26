README for massconserve.exe and massconserve.py/mvf_lite.py

Copyright 2013, Christopher L. Williams

Please cite Williams et al., Med. Phys. 40 (7), July 2013 if you use this software in your work.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. 

#####################################################################

This software is designed to modify an existing XCAT phantom in order to locally conserve mass in the lung.  The software works by using a vector field generated with the XCAT software (mode=4) to calculate the local volume changes during respiration.  It relates each region in a "target" frame back to a "reference" frame, and corrects the volume of the target frame so that mass is conserved.

The software requires 5 inputs:

* A reference XCAT phantom. This is commonly the frame #1 of the XCAT phantom, generated with "mode=0" in the XCAT software (e.g. test_atn_1.bin).

* A target XCAT phantom frame.  This is some other frame of the XCAT phantom which you want to correct (e.g. test_atn_3.bin), also generated with "mode=0" in the XCAT.

* A vector file, generated with mode=4 that is valid for the previously given phantom frames (e.g. test_vec_frame1_to_frame3.txt)

* The XCAT parameter (".par") file.  This is used to extract metadata about the shape of the phantom, and lesion location.

* The XCAT LOG file.  This should be the log file that was produced when XCAT was run with mode=0.  This is used to determine attenuation values for various tissues (e.g. test_log).

The software can optionally insert a tumor into a mode=0 phantom, using the -t flag, which will place a tumor such that it moves in the same manner as the surrounding lung tissues.

Further details on usage can be viewed by invoking the "--help" flag.