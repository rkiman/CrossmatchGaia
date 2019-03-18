'''
This code is based on the documentation: 
https://astroquery.readthedocs.io/en/latest/gaia/gaia.html
'''
from astroquery.gaia import Gaia

def query_gaia_match(user_name,table_2_match,radius_arc):
	'''
	This function queries Gaia DR2 to find the match with table
	table_2_match. This table needs to be uploaded to the Gaia
	Archive in advance. The steps to upload the table can be found
	in http://gea.esac.esa.int/archive-help/index.html.

	Input:
	user_name(str): User name in the Gaia Archive.
	table_2_match(str): Name of the table uploaded to the Gaia 
	Archive.
	radius_arc(float): Radius of search in arcsec.

	Output:
	file.vot: Output file from the Archive with all the columns
	in table_2_match plus all Gaia DR2 columns. 
	'''

	Gaia.login_gui() 
	#Ask for userName and userPassword to authenticated access mode
	#This could be done with: 
	#Gaia.login(user='user_name', password='user_pass') 

	#Cross-match user table and gaia source
	job = Gaia.launch_job_async("""\
		SELECT crossmatch_positional(\
		'user_{}','{}',\
		'gaiadr2','gaia_source',\
		{},\
		'xmatch')\
		FROM dual;\
		""".format(user_name,table_2_match,radius_arc))

	#For the matches saved in test_xmatch, get the information from Gaia. 
	#The last line saves the information into a vot table
	job2 = Gaia.launch_job_async("""\
		SELECT c."dist", a.*, b.* \
		FROM user_{}.{} AS a, \
		gaiadr2.gaia_source AS b, \
		user_{}.xmatch AS c \
		WHERE (c.{}_{}_oid = a.{}_oid AND \
		c.gaia_source_source_id = b.source_id)\
		""".format(user_name,table_2_match,user_name,
			table_2_match,table_2_match,table_2_match), dump_to_file=True)

	Gaia.logout()

