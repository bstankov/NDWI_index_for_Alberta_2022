#-------------------------------------------------------------------------------
# Name:        module1
# Identify files rfrom Near Real time i age depository and write them
#onto a list and then use wget to download those files
# Author:      bob.stankovic
#
# Created:     12/02/2021
# Copyright:   (c) bob.stankovic 2021
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import requests
import re, os, sys
from datetime import date
import numpy as np
#from bitBucke_func import download_hdf

def mkdir(dirname, remove=True, chdir=False):
    import shutil
    """create a directory dirnme.  if it iexists     , it is removed by shutil.rmtree
    """
    if os.path.isdir(dirname):
        if remove:
            shutil.rmtree(dirname)
        else:
            return False  # did not make new directory
    os.mkdir(dirname)
    return
def re_read_urlbase(url_base):
    '''  read the content from the web and get ready for parsing '''

    mod10_req = requests.get(url_base)  #read folders and sort them
    print(mod10_req.status_code)
    print(mod10_req.ok)
    mod10_sadrzaj = mod10_req.text
    #print(dir(req))
    mod10_stripped = re.sub('<[^<]+?>', '',mod10_sadrzaj)
    mod10_size = len(mod10_stripped)
    print("leng:", mod10_size)   #there are 73697
    return mod10_stripped
def write_file_out(targets, outF ):
    with open(outF, 'w') as out_file:
        for i in targets:
            out_file.write(str(i.strip()))

            out_file.write('\n')
    return

def get_MOD09GA1N_hdfs():
    '''read the NSIDC web content, parse it for daily folders
    identify the most recent folder and go and grab 4 HDF files covering AB  '''
    #downloading folder for MOD10A1F.061 data products

    #MODIS/Terra Surface Reflectance 8-Day L3 Global 500 m SIN Grid
    base = "https://nrt4.modaps.eosdis.nasa.gov/api/v2/content/archives/allData/61/MOD09A1N/2021/"

    folders_req = re_read_urlbase(base)

    #identify folders that start with 2021
    #pattern = re.compile(r'2021[.]\d{2}[.]\d{2}')
    pattern = re.compile(r'\d{3}')
    ##    matches = pattern.finditer(stripped)
    matches = pattern.findall(folders_req )
    jd = []
    for match in matches:
        match_ = match[0:3]
        if int(match_) < 365:
            jd.append(match_)
            print(match_)
    cur_fold = jd[-1]
    print("the most recent folder:",  cur_fold )
    folder_base = "%s%s" % (base,  cur_fold )

    print("folder_basee:", folder_base)

    # #create directory to download data into
    root = r"U:/RS_Task_Workspaces/NDWI/2021/data"
    dirname = root+'/'+ cur_fold +'/'
    print("dirname:",dirname)
    mkdir( dirname, remove=True, chdir=False)
    # #****************************************************

    # #read the most recent folder's content...
    mod09A1N_stripped = re_read_urlbase(folder_base)
    print("leng:",len(mod09A1N_stripped ))
    content_leng = len(mod09A1N_stripped )

    #identify necessary tiles covering Alberta
    hdf_dir=['h10v03','h10v04','h11v03','h12v03']
    key = ['.hdf','.xml']
    key_1 = 'hdf.xml'

    targets = []   #it contains hdf files to be downloaded
    if content_leng > 5000:
        #IN CASE THE MOST RECENT FOLDER MISSING NEEDED FILES WE download data from

        #print("sadrz:", mod09A1N_stripped)
        for line in mod09A1N_stripped.split('\n'):
           # if (key[0] in line and not key[1] in line):
            if ('hdf' in line and not 'hdf.met' in line ):
                line_1 = line.split('.hdf')
                final_name = line_1[0]+'.hdf'
               # print("fin_name:",final_name)
                for j in hdf_dir:
                    if j in line:
                       print("fname:",final_name)
                       targets.append(folder_base +'/'+ str(final_name.strip()))

        print(" target:",targets)
        # a textfile to list all files to be downloaded
       # outF = root + '/' + "mod09GA1N_list.txt"
        outF = dirname +  "mod09A1N_lista.txt"
        #write data for download into a text file ************************
        write_file_out(targets, outF)

    else:
        print("the most current folder is not populated yet")


    return outF


#generat a textfile to be downloaded
outfile = get_MOD09GA1N_hdfs()
print("outfile:",outfile)

#download files from the textfile using 'bitBucket.py" script
saveDir = os.path.dirname(outfile)
files = os.path.basename(outfile)
print("saveDir:",saveDir)
#Scripted downloads will need to use MODAPS app keys in order to be properly authorized.
#MODAPS app keys are string tokens that identify who you are. App keys get passed in the Authorization header of each HTTP GET request.
token = "Authorization: Bearer YnN0YW5rb3Y6WW05aUxuTjBZVzVyYjNacFkwQm5iM1l1WVdJdVkyRT06MTYxNTk5Mzg3ODowMTRkYmE3MjliOWY4MTQ1ODJmNWRiZjJkYzJlNzY5MTk3MzM1MzE4"
with open(outfile, 'r') as reader:
    # Read and print the entire file line by line
    for line in reader:
        # print(line, end='')
        com = 'wget -e robots=off -m -np -R .html,.tmp -nH --cut-dirs=4 "%s" --header "%s" -P %s' % (line.strip(), token, saveDir)
        print(com)
        os.system(com)




