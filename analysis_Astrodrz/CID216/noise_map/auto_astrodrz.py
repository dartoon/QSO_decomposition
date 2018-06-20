 #source activate iraf27
 #!/usr/bin/env python
import drizzlepac
from drizzlepac import astrodrizzle
unlearn astrodrizzle
astrodrizzle.AstroDrizzle('*flt.fits',output='final', build=yes, static=yes,
                          skysub=yes,driz_separate=yes, median=yes, blot=yes,
                          driz_cr=yes,driz_combine=yes,final_wcs=yes,
                          final_bits=576,final_scale=0.0642,final_pixfrac=0.8,
                          final_rot=0.0)

#final_kernel='gaussian',


