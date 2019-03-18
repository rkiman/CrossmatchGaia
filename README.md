# CrossmatchGaia

This Package is useful to cross-match a sample with Gaia DR2. It will query the Gaia Archive so you need to have an account: http://gea.esac.esa.int/archive/

The table with ra and dec from your sample to match with Gaia DR2 needs to be uploaded before running this code. For the steps on how to upload the table go through the steps in this tutorial: http://gea.esac.esa.int/archive-help/index.html

Right now the main example makes the query to crossmatch. Need to add the next steps. I also want to add a function at the beginning that creates a csv file to upload to the Gaia Archive including the extra_id number that will be useful later. 