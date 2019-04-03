The catalogs contain the following columns:

# COL 1,  NUMBER: ID, identical to the official CANDELS SExtractor-based segmentation map ID, and identical in the catalogs for the various filters.
# COL 2,  RA: from SExtractor, measured in F160W
# COL 3,  DEC: idem.
# COL 4,  f: FLAG value (0: good fit; 1: suspicious fit; 2: bad fit; 3: no fit -- see data paper)
# COL 5,  mag: total AB magnitude from best-fitting Sersic model (GALFIT)
# COL 6,  dmag: 1-sigma uncertainty in mag
# COL 7,  re: semi-major axis in arcsec of the ellipse that contains half of the total light in the best fitting Sersic model
# COL 8,  1-sigma uncertainty on re
# COL 9,  n: Sersic index of the best-fitting Sersic model
# COL 10, dn: 1-sigma uncertainty on n
# COL 11, q: axis ratio of the best fitting Sersic model
# COL 12, dq: 1-sigma uncertainty on q
# COL 13, pa: position angle in degrees (0: North; 90: East)
# COL 14, dpa: 1-sigma uncertainty on pa
# COL 15, sn: S/N as measured in the F???W filter, integrated over the F160W segmentation region

Note that the uncertainties presented here are not computed by galfit, but
derived as described in the data paper.

The cutout tarballs are each several GB large and contain the galfit-produced
images (extension 1), Sersic model images (extension 2), and residual images
(extension 3). The names of the fits files correspond to the NUMBER (COL 1) in
the catalogs.

