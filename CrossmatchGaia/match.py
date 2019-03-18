from CrossmatchGaia import jdutil
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np

def calc_sep(ra1,ra2,dec1,dec2):
    '''
    input:
    ra1 (float/array): right ascention first group of stars (deg) 
    dec1 (float/array): declination first group of stars (deg) 
    ra2 (float/array): right ascention second group of stars (deg) 
    dec2 (float/array): declination second group of stars (deg)
    
    returns:
    sep_arcsec (float/array): angular separation (arcsec)
    '''
    positions1 = SkyCoord(ra1 * u.deg, dec1 * u.deg, frame='icrs')
    positions2 = SkyCoord(ra2 * u.deg, dec2 * u.deg, frame='icrs')

    sep_arcsec = positions1.separation(positions2).arcsecond
    return sep_arcsec

def move_epoch_back(ra, dec, pmra, pmde, goodpm, mjd_from=None, 
                    yr_from=None, mjd_to=None, yr_to=None):
    '''
    warning: this function needs python 3
    needs:
    ra (array/list): right ascention (deg)
    de (array/list): declination (deg)
    pmra (array/list): proper motion ra (mas/yr)
    pmde (array/list): proper motion dec (mas/yr)
    mjd_to (array/list):  Modified Julian Date (epoch to move stars to)
    yr_to (array/list):  epoch to move the stars to (yr)
    goodpm (array/list): 1 if pm can be used, 0 if not
    mjd_from (array/list):  Modified Julian Date (epoch to move stars from)
    yr_from (float): epoch to move the stars from (yr)

    returns:
    ra_new (array/float): right ascention in the epoch of mjd
    dec_new (array/float): declination in the epoch of mjd
    '''
    if mjd_to is None:
        if yr_to is None:
            raise ValueError("One of mjd_to or yr_to has to be defined.")
        else:
            year_to = yr_to
    else:
        if isinstance(mjd_to, np.ndarray):
            year_to = np.array([jdutil.jd_to_date(jdutil.mjd_to_jd(x))[0] 
                for x in mjd_to])
        else:
            year_to = jdutil.jd_to_date(jdutil.mjd_to_jd(mjd_to))[0] 

    if mjd_from is None:
        if yr_from is None:
            raise ValueError("One of mjd_from or yr_from has to be defined.")
        else:
            year_from = yr_from
    else:
        if isinstance(mjd_from, np.ndarray):
            year_from = np.array([jdutil.jd_to_date(jdutil.mjd_to_jd(x))[0] 
                for x in mjd_from])
        else:
            year_from = jdutil.jd_to_date(jdutil.mjd_to_jd(mjd_from))[0]

    if isinstance(ra, np.ndarray):
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
    else:
        if goodpm == 1:
            dt = year_to - yr_from
            pmra_deg_yr = pmra*(1/3600000)/np.cos(np.deg2rad(dec))
            pmde_deg_yr = pmde*(1/3600000)

            ra_new = ra + dt*pmra_deg_yr
            dec_new = dec + dt*pmde_deg_yr
        else:
            raise ValueError("If you don't trust the pm don't move the epoch.")

    return ra_new,dec_new


def check_duplicates(table,extra_id):
    '''
    This function looks for duplicates in the sample. 
    In case there are duplicates, use remove_duplicates.

    needs:
    table (astropy table) 
    extra_id (string) name of the column where the ids are repeated 
    if it is a duplicated
    '''
    for i in range(len(table[extra_id])-1):
        num = table[extra_id][i+1] - table[extra_id][i]
        if(num == 0):
            print(table[extra_id][i])


def remove_duplicates(table,extra_id,separation):
    '''
    This function removes the duplicates after doing a match, 
    keeping only the closest neighbours.

    needs:
    table (astropy table) 
    extra_id (string) name of the column where the ids are 
                      repeated if it is a duplicated
    separation (string) name of the column with the separation 
                        between the matches

    Note: if after running this functions ones there are still 
    duplicates, run it again.
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
        table.remove_row(x-n) #Each time I remove a row, 
        n+=1                  #the number of all the following rows changes.
    return
