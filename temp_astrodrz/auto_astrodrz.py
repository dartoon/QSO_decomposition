 #!/usr/bin/env python
import drizzlepac
from drizzlepac import astrodrizzle
unlearn astrodrizzle
astrodrizzle.AstroDrizzle('*flt.fits',output='f160w', build=yes, static=no,
                          skysub=yes,driz_separate=no, median=no, blot=no,
                          driz_cr=no,driz_combine=yes,final_wcs=yes,
                          final_bits=576,final_scale=0.0642,final_pixfrac=0.8,
                          final_rot=0.0)

#final_kernel='gaussian',
