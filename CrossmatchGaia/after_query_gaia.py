
from astropy.table import Table, join

def remove_duplicates(table):
	remove = []
	for i in range(len(table['extra_id'])-1):
	    num = table['extra_id'][i+1] - table['extra_id'][i]
	    if(num == 0):
	        diff = table['dist'][i+1]-table['dist'][i]
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

def match_MLSDSS(table):
	MLSDSS = Table.read('../Catalogs/MLSDSS.fits')
	gaia_MLSDSS = join(MLSDSS,table, keys='extra_id', join_type='outer', uniq_col_name='flag')
	print(len(gaia_MLSDSS))
	gaia_MLSDSS.write('MLSDSS_gaia_5.fits', format='fits')
	return 

gaia_MLSDSS_pre = Table.read('async_20180425062521.vot',format='votable')
remove_duplicates(gaia_MLSDSS_pre)
print(len(gaia_MLSDSS_pre))

match_MLSDSS(gaia_MLSDSS_pre)


