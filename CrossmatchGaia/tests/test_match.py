import numpy as np
from CrossmatchGaia import pre_query_gaia 

def test_move_epoch_back():
	ra = np.array([1, np.nan, 1.1])
	dec = np.array([3, np.nan, 2.3])
	pmra = np.array([3, np.nan, 8.5])
	pmdec = np.array([2, np.nan, 1.7])
	mjd = np.array([55568, 55614, 55530])
	goodpm = np.array([1,1,1])
	yr = 2015.5
	ra_new,dec_new = pre_query_gaia.move_epoch_back(ra, dec, pmra, pmdec, 
										   goodpm, mjd_to=mjd, yr_from=yr)
	mask1 = np.array([True if str(x)!= 'nan' else False for x in ra])
	mask2 = np.array([True if str(x)!= 'nan' else False for x in ra_new])

	assert len(ra_new[mask2]) == len(ra[mask1])
	assert len(ra_new[mask2]) == 2
	assert np.isclose(ra_new[0],0.99999624)
	assert np.isclose(dec_new[0],2.9999975)
	assert np.isclose(ra_new[2],1.099987)
	assert np.isclose(dec_new[2],2.2999974)
	assert np.isclose(pre_query_gaia.move_epoch_back(1, 3, 3, 2, 1, mjd_to=55568, 
		              yr_from=2015.5)[0],0.99999624)

def test_calc_sep():

	ra1 = 185.9891288578898
	dec1 = 24.891416254023895
	ra2 = 69.56196935050878
	dec2 = 26.194469466376514

	sep = pre_query_gaia.calc_sep(ra1,ra2,dec1,dec2)

	assert np.isclose(sep, 360589.5322959722)
