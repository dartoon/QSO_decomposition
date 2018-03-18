# QSO_decomposition

    "''
    With this task, we can do the QSO decomposition. These are inclduing
    1. Giving the HST image, and doing the SWarp.
    2. Quickly the cut out and label the PSF.
    3. Do the SB compare. 
    4. Average the PSF,
    5. Doing MCMC
    '''
1. Do the SWarp. Go to folder 'fits_image' and run "python auto_indivi_4swarp.py". After inputing the 'target' name (such as 'CID70'), all the data images will be at "target/swarp", with a 'auto_swarp.sh' in that folder. One can do '. /auto_swarp.sh' and combining these images in name "coadd.fits".

2. Load the star.reg and see which is star and which is QSO. Selecting the QSO and PSF candidates. Run the cutout.py template in folder "target/analysis" to cut out the PSF and QSO. This will get the .png image with PSF id and QSO position labeled.

# Analysis the PSF with QSO. Select PSFs and average. MCMC fitting. Get results.
3. Use