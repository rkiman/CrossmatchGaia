from CrossmatchGaia import jdutil
import numpy as np
from astropy.table import Table

def calc_sep(ra1,ra2,dec1,dec2):
    '''
    input:
    ra1 (array): right ascention first group of stars (deg) 
    dec1 (array): declination first group of stars (deg) 
    ra2 (array): right ascention second group of stars (deg) 
    dec2 (array): declination second group of stars (deg)
    
    returns:
    sep (array): angular separation (arcsec)
    '''
    positions_gaiamatches = SkyCoord(ra1 * u.deg, dec1 * u.deg, frame='icrs')
    positions_sdss = SkyCoord(ra2 * u.deg, dec2 * u.deg, frame='icrs')

    sep_arcsec = positions_gaiamatches.separation(positions_sdss).arcsecond
    return sep_arcsec

def move_epoch_back(ra, dec, pmra, pmde, mjd, goodpm, yr_from):
    '''
    warning: this function needs python 3
    needs:
    ra (array/list): right ascention (deg)
    de (array/list): declination (deg)
    pmra (array/list): proper motion ra (mas/yr)
    pmde (array/list): proper motion dec (mas/yr)
    mjd (array/list):  Modified Julian Date (epoch to move stars to)
    goodpm (array/list): 1 if pm can be used, 0 if not
    yr_from (float): epoch to move the stars from (yr)

    returns:
    ra_new: right ascention in the epoch of mjd
    dec_new: declination in the epoch of mjd
    '''
    year_to = np.array([jdutil.jd_to_date(jdutil.mjd_to_jd(x))[0] 
        for x in mjd])
    ra_new = np.zeros(len(ra))
    dec_new = np.zeros(len(ra))

    mask_good_pm = goodpm == 1
    dt = year_to - yr_from
    pmra_deg_yr = pmra*(1/3600000)/np.cos(np.deg2rad(dec))
    pmde_deg_yr = pmde*(1/3600000)

    ra_new[mask_good_pm] = ra[mask_good_pm] + dt*pmra_deg_yr[mask_good_pm]
    dec_new[mask_good_pm] = dec[mask_good_pm] + dt*pmde_deg_yr[mask_good_pm]

    ra_new[~mask_good_pm] = ra[~mask_good_pm] 
    dec_new[~mask_good_pm] = dec[~mask_good_pm] 

    return ra_new,dec_new

def move_epoch_back_yr(ra, de, pmra, pmde, yr2, goodpm, yr):
    '''
    warning: this function needs python 3
    needs:
    ra (array/list): right ascention (deg)
    de (array/list): declination (deg)
    pmra (array/list): proper motion ra (mas/yr)
    pmde (array/list): proper motion dec (mas/yr)
    mjd (array/list):  Modified Julian Date (epoch to move stars to)
    goodpm (array/list): good proper motion flag from sdss
    yr (float): epoch to move the stars from (yr)

    returns:
    ra_new: right ascention in the epoch of mjd
    dec_new: declination in the epoch of mjd
    '''
    year = yr2
    ra_new = np.zeros(len(ra))
    de_new = np.zeros(len(ra))
    for i in range(len(ra)):
        if(goodpm[i] == 1):
            ra_new[i] = ra[i] + (year-yr)*pmra[i]*(1/3600000)/np.cos(np.deg2rad(de[i]))
            de_new[i] = de[i] + (year-yr)*pmde[i]*(1/3600000)
        else:
            ra_new[i] = ra[i] 
            de_new[i] = de[i] 
    return ra_new,de_new

def move_epoch(ra, de, pmra, pmde, mjd, goodpm, yr):
    '''
    warning: this function needs python 3
    needs:
    ra (array/list): right ascention (deg)
    de (array/list): declination (deg)
    pmra (array/list): proper motion ra (mas/yr)
    pmde (array/list): proper motion dec (mas/yr)
    mjd (array/list):  Modified Julian Date
    goodpm (array/list): good proper motion flag from sdss
    yr (float): epoch to move the stars to (yr)

    returns:
    ra_new: right ascention in 2015.5
    dec_new: declination in 2015.5
    '''
    year = np.array([jdutil.jd_to_date(jdutil.mjd_to_jd(x))[0] for x in mjd])
    ra_new = np.zeros(len(ra))
    de_new = np.zeros(len(ra))
    for i in range(len(ra)):
        if(goodpm[i] == 1):
            ra_new[i] = ra[i] + (yr-year[i])*pmra[i]*(1/3600000)/np.cos(np.deg2rad(de[i]))
            de_new[i] = de[i] + (yr-year[i])*pmde[i]*(1/3600000)
        else:
            ra_new[i] = ra[i] 
            de_new[i] = de[i] 
    return ra_new,de_new

def move_single_epoch(ra,de,pmra,pmde,mjd):
    '''
    warning: this function needs python 3
    needs:
    ra (float): right ascention (deg)
    de (float declination (deg)
    pmra (float): proper motion ra (mas/yr)
    pmde (float): proper motion dec (mas/yr)
    mjd (float):  Modified Julian Date
    
    returns:
    ra_new: right ascention in 2015
    dec_new: declination in 2015 
    '''
    year = jdutil.jd_to_date(jdutil.mjd_to_jd(mjd))[0]  

    ra_new = ra + (2015.0-year)*pmra*(1/3600000)/np.cos(de)
    de_new = de + (2015.0-year)*pmde*(1/3600000)

    return ra_new,de_new

def check_duplicates(table,extra_id):
    '''
    This function looks for duplicates in the sample. In case there are duplicates, use remove_duplicates.

    needs:
    table (astropy table) 
    extra_id (string) name of the column where the ids are repeated if it is a duplicated
    '''
    for i in range(len(table[extra_id])-1):
        num = table[extra_id][i+1] - table[extra_id][i]
        if(num == 0):
            print(table[extra_id][i])


def remove_duplicates(table,extra_id,separation):
    '''
    This function removes the duplicates after doing a match, keeping only the closest neighbours.

    needs:
    table (astropy table) 
    extra_id (string) name of the column where the ids are repeated if it is a duplicated
    separation (string) name of the column with the separation between the matches

    Note: if after running this functions ones there are still duplicates, run it again.
    '''
    remove = []
    for i in range(len(table[extra_id])-1):
        num = table[extra_id][i+1] - table[extra_id][i]
        if(num == 0):
            diff = table[separation][i+1]-table[separation][i]
            if(diff > 0):
                remove.append(i+1)
            else:
                remove.append(i)

    remove.sort()
    n=0
    for x in remove:
        table.remove_row(x-n) #Each time I remove a row, the number of all the following rows changes.
        n+=1
    return
