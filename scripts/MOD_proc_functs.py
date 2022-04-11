import gdal
import sys
import shutil
import subprocess
import os
import glob
import datetime
import numpy as np
import osr
import ogr
import gdalconst
import datetime
gdal.UseExceptions()
#import gdalnumeric
#import builtins
# from pathlib import Path
print("hi all libraries are loaded....")

class proces_HDF(object):

    ''' the module process MODIS hdf file in .../AB_NDWI/hdf folder . First of all, it
    extracts the 2  and 5 neccessary image bands , reproject them to geographic and then to 10 tm projection.
    After that it creates NDWI image, which is then copied to the archive. there are two archives dating since January
    2017 onwards, one for path H11V03 and another for H12V03.  finally , all the images in the archive are used
    to create a time series mosaic in 'vrt.
    '''
    def save_vrt2tif(self, output_vrt, tifname):
        gdal_translate = r'C:/Program Files/GDAL/gdal_translate.exe'
        # save virt as tif
        cmd = '-of GTiff'
        tran_cmd = ' '.join([gdal_translate, cmd, output_vrt, tifname])
        print("transcom:", tran_cmd)
        subprocess.Popen(tran_cmd)
        # return tifname

    def gdalmerge(self ,output_vrt, direct_list):
        ''' create vrt band to update main file'''
        gdalmerge = r'C:/Program Files/GDAL/gdal_merge.py'
        cmd = "-separate -o " + output_vrt + ' -of GTiff ' + direct_list
        fullCmd = ' '.join([gdalmerge, cmd])
        print("gdalmerge:", fullCmd)
        subprocess.Popen(fullCmd)
        print("output file is done")
        print("\n")

    def img2array(self,raster_in):
        raster = gdal.Open(raster_in)
        band = raster.GetRasterBand(1)
        array = band.ReadAsArray()
        return array

    def array2raster(self,out_raster_name, raster_5, raster_2):
        ''' the function is fed by two images MOD_b2 and MOD_b5
        and it scales the image and spits out NDWI index out'''
        np.seterr(divide='ignore', invalid='ignore')
        ds2 = gdal.Open(raster_2)
        b2_array = ds2.GetRasterBand(1).ReadAsArray()
        ds5 = gdal.Open(raster_5)
        b5_array = ds5.GetRasterBand(1).ReadAsArray()
        b2_array = np.array(b2_array, dtype=float)
        b5_array = np.array(b5_array, dtype=float)
        # *************************
        ds_proj = ds2.GetProjection()
        geot = ds2.GetGeoTransform()
        width = ds2.RasterXSize
        height = ds2.RasterYSize
        # ***********************************
        ndwi = ((b2_array * 0.0001) - (b5_array * 0.0001)) / ((b2_array * 0.0001) + (b5_array * 0.0001))

        driver = gdal.GetDriverByName('GTiff')
        out_ndwi = driver.Create(out_raster_name, width, height, 1, gdal.GDT_Float32)
        out_ndwi.SetGeoTransform((geot))
        out_ndwi.SetProjection((ds_proj))
        out_band = out_ndwi.GetRasterBand(1).WriteArray(ndwi)
        ds5 = None
        ds2 = None
        return out_ndwi

    def mkdir(self ,dirname, remove=True, chdir=False):
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
    def clean_dir(self,apath):
        os.chdir(apath)
        for root, dirs, files in os.walk(".", topdown=False):
            for file in files:
                print(os.path.join(root, file))
                os.remove(os.path.join(root, file))
        return

    def append_date(adatum, date_file):
        ''' when we do updates with a new file
        it opens a file with dates *.dates and append
        a new date to .dates file'''
        with open(date_file, 'a') as f:

            f.write(f'\n{adatum}') #for python 3

    def save_raster(self ,output_name, raster_data, dataset, driver="GTiff"):
        """
        A function to save a 1-band raster using GDAL to the file indicated
        by ``output_name``.
        Parameters:
            output_name: str ........        The output filename, with full path and extension if required
        raster_data: array ........        The array that we want to save
        dataset: str.............        Filename of a GDAL-friendly dataset that we want to use to
            read geotransform & projection information
        driver: str .......        A GDAL driver string, like GTiff or HFA.
        """
        # Open the reference dataset
        g = gdal.Open(dataset)
        # Get the Geotransform vector
        geo_transform = g.GetGeoTransform()
        x_size = g.RasterXSize  # Raster xsize
        y_size = g.RasterYSize  # Raster ysize
        srs = g.GetProjectionRef()  # Projection
        # Need a driver object. By default, we use GeoTIFF
        driver = gdal.GetDriverByName(driver)
        dataset_out = driver.Create(output_name, x_size, y_size, 1, gdal.GDT_Float32)
        dataset_out.SetGeoTransform(geo_transform)
        dataset_out.SetProjection(srs)
        dataset_out.GetRasterBand(1).WriteArray(raster_data)
        dataset_out = None

    def export_bands(self ,src, dst):
        # cmd = "gdal_translate.exe -b 2"
        cmd = "gdal_translate.exe -of MEM -b 2"
        fullCmd = ' '.join([cmd, src, dst])
        print("com:", fullCmd)
        os.system(fullCmd)
        return dst

    def extport2sinus(self ,hdf_layer, dst_singrd):
        # gdalwarp - of GTiff HDF4_EOS: EOS_GRID:"MOD09A1.A2020001.h11v03.006.2020010223355.hdf": MOD_Grid_500m_Surface_Reflectance:sur_refl_b02 b2.tif'
        # cmd = 'gdalwarp.exe -of GTiff -tps -t_srs "EPSG:4326" -ts 2400 2400'
        cmd = 'gdalwarp.exe -of GTiff '
        # cmd = "gdal_translate.exe -of MEM -b 2"
        fullCmd = ' '.join([cmd, hdf_layer, dst_singrd])
        print("com:", fullCmd)
        os.system(fullCmd)
        return dst_singrd

    def modis2normsphere(self ,hdf_layer, dst_singrd):
        # gdalwarp - of GTiff HDF4_EOS: EOS_GRID:"MOD09A1.A2020001.h11v03.006.2020010223355.hdf": MOD_Grid_500m_Surface_Reflectance:sur_refl_b02 b2.tif'
        # cmd = 'gdalwarp.exe -of GTiff -tps -t_srs "EPSG:4326" -ts 2400 2400'
        cmd = 'gdalwarp.exe -of GTiff -t_srs "+proj=latlong +ellps=sphere"'
        # cmd = "gdal_translate.exe -of MEM -b 2"
        fullCmd = ' '.join([cmd, hdf_layer, dst_singrd])
        print("com:", fullCmd)
        os.system(fullCmd)
        return dst_singrd

    def sin_wgs84(self ,input_sin, out84):
        # gdalwarp -of GTiff -t_srs "EPSG:4326" -ts 2400 2400 "b2.tif" "b2_wgs84.tif"
        cmd = 'gdalwarp.exe -of GTiff -t_srs "EPSG:4326" -ts 2400 2400'
        #crs_str = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
        fullCmd = ' '.join([cmd, input_sin, out84])
        print("com:", fullCmd)
        os.system(fullCmd)
        return out84

    def wgs84_epsg3400(self ,out84, out_10tm):
        '''reproject image from Geographic to 10TM'''
        ##gdalwarp -t_srs EPSG:3400 -tr 500 500 -te 152000 5853000 860600 6660000 b2_wgs84.tif b2_wgs84_10tm_sub.tif
        cmd = 'gdalwarp.exe -t_srs EPSG:3400 -tr 500 500 -te 88259.837 5342900.976 944759.837 6744400.97'
        fullCmd = ' '.join([cmd, out84, out_10tm])
        print("com:", fullCmd)
        os.system(fullCmd)
        return out_10tm

    def run_gridClip(self ,Xmin, Ymin, Xmax, Ymax, src, dst):
        #cmd = "gdalwarp.exe -t_srs EPSG:3400 -tr 500 500 -te"
        Xmin = -120.6
        Ymin= 48.7
        Xmax = -109.5
        Ymax= 60.1
        cmd = "gdalwarp.exe -te"
        # gdal_Warp = 'C:/Program Files/GDAL/gdalwarp.exe'
        # fullCmd = ' '.join([gdal_Warp, cmd, str(Xmin), str(Ymin), str(Xmax), str( Ymax), "-dstnodata -9999.0 ", src, dst])
        fullCmd = ' '.join([cmd, str(Xmin), str(Ymin), str(Xmax), str(Ymax), "-dstnodata -9999.0 ", src, dst])
        print("com:", fullCmd)
        os.system(fullCmd)
        return dst
    # def get_subset(self,):
    # gdal_translate - of GTiff - a_ullr < top_left_lon > < top_left_lat
    # < bottom_right_lon > < bottom_right_lat > -a_srs EPSG: 4269 input.png output.tif

    def gridClip_shp(self , src, dst, shape):
        cmd = "gdalwarp.exe -dstnodata -9999 -cutline"
        gdal_Warp = 'C:/Program Files/GDAL/gdalwarp.exe'
        cropToCutline = True,
        # fullCmd = ' '.join([gdal_Warp, cmd, str(Xmin), str(Ymin), str(Xmax), str( Ymax), "-dstnodata -9999.0 ", src, dst])
        fullCmd = ' '.join([cmd, shape, src, dst])
        print("com:", fullCmd)
        os.system(fullCmd)
        return dst



    # def getdate_from_hdf_jd(self ,file):
    #     ''' convert ordinal day to a date and create name for a file'''
    #     import datetime
    #     file1 = file.split('.')[1][1:]
    #     jd = file1[2:]
    #     a = datetime.datetime.strptime(str(jd), '%y%j').date()
    #     adate = "A" + file1 + '_' + a.strftime('%Y%m%d')
    #     return adate

    def youCanQuoteMe(self ,item):
        return "\"" + item + "\""

    def img2array(raster_in):
        raster = gdal.Open(raster_in)
        band = raster.GetRasterBand(1)
        array = band.ReadAsArray()
        return array

    def array2img(self,dst_filename, img_in, array):

        driver = gdal.GetDriverByName('GTiff')
        dst_ds = driver.CreateCopy(dst_filename, img_in, strict=0)
        out_band = dst_ds.GetRasterBand(1).WriteArray(array)
        out_band = None
        dst_ds = None
        # return dst_ds
    #
    def getDate_from_hdfJD(self, file, index):
        ''' convert ordinal day to a date and create name for a file'''
        '''the following products have these indexes:
        MOD35_L2.A2021016.1530.061.2021017012106 ---> 2
        MOD09 -->   1
        MOD10 ----> 1
        '''
        # import datetime
        file1 = file.split('_')[index][1:]
        jd = file1[2:]
        a = datetime.datetime.strptime(str(jd), '%y%j').date()
        adate = "A" + file1 + '_' + a.strftime('%Y%m%d')
        return adate

    # def array_2_raster(self,out_raster_name, raster_in, array):
    #     import gdal
    #     raster = gdal.Open(raster_in)
    #     geotransform = raster.GetGeoTransform()
    #     origin_x = geotransform[0]
    #     origin_y = geotransform[3]
    #     pixel_width = geotransform[1]
    #     pixel_height = geotransform[5]
    #     cols = array.shape[1]
    #     rows = array.shape[0]
    #     driver = gdal.GetDriverByName('GTiff')
    #     out_img = driver.Create(out_raster_name, cols, rows, 1, gdal.GDT_Float32)
    #     out_img.SetGeoTransform((origin_x, pixel_width, 0, origin_y, 0, pixel_height))
    #     out_band = out_img.GetRasterBand(1)
    #     out_band.WriteArray(array)
    #     out_raster_srs = osr.SpatialReference()
    #     # copy projection information into output osr object
    #     out_raster_srs.ImportFromWkt(raster.GetProjectionRef())
    #     out_img.SetProjection(out_raster_srs.ExportToWkt())
    #
    #     out_band.FlushCache()
    #     return out_img
    # def array2raster(self,out_raster_name, raster_in, array):
    #     #raster = gdal.Open(raster_in)
    #     geotransform = raster_in.GetGeoTransform()
    #     origin_x = geotransform[0]
    #     origin_y = geotransform[3]
    #     pixel_width = geotransform[1]
    #     pixel_height = geotransform[5]
    #     cols = array.shape[1]
    #     rows = array.shape[0]
    #     driver = gdal.GetDriverByName('GTiff')
    #     out_img = driver.Create(out_raster_name, cols, rows, 1, gdal.GDT_Float32)
    #     out_img.SetGeoTransform((origin_x, pixel_width, 0, origin_y, 0, pixel_height))
    #     out_band = out_img.GetRasterBand(1)
    #     out_band.WriteArray(array)
    #     out_raster_srs = osr.SpatialReference()
    #     # copy projection information into output osr object
    #     out_raster_srs.ImportFromWkt(raster_in.GetProjectionRef())
    #     out_img.SetProjection(out_raster_srs.ExportToWkt())
    #     out_band.FlushCache()
    #     return out_img
    def interpolate_infs(self,X):
        """Overwrite INFs with column value interpolations."""
        for j in range(X.shape[1]):
            # mask_j = np.isnan(X[:,j])   #for nans
            mask_j = np.isinf(X[:, j])  # for infns
            X[mask_j, j] = np.interp(np.flatnonzero(mask_j), np.flatnonzero(~mask_j), X[~mask_j, j])
        return X

    def remove_nan(self,raster_in):
        """Overwrite NaNs with column value interpolations."""
        # overwrite 'inf' values with column mean
        raster_no_INF = self.interpolate_infs(raster_in)
        # get column means
        raster_no_INF_CM = np.nanmean(raster_no_INF, axis=0)
        # find indexes where we need replace 'nan' values with col means
        inds = np.where(np.isnan(raster_no_INF))
        #  place column means in the indices. align the arrays using take
        raster_no_INF[inds] = np.take(raster_no_INF_CM, inds[1])
        return raster_no_INF

    def copy2servers(self,  afile):
        '''
        :param listOfServers:
        :param afile:
        :return:
        it copies afile input onto different servers which paths are indicated in listOfServers
        '''

        servers = [
            r'\\GOA\Shared\AF\Wildfire_Geospatial\Data_Repository\MODIS\Snow_and_Ice_Composite\MODIS_TerraAqua_Tiffs\original',
            r'\\CAL-GOA-APM-419\NDSI_MODIS', r'\\C-GOA-APM-11648\NDSI_MODIS', r'\\CAL-GOA-APM-420\NDSI_MODIS',
            r'\\C-GOA-APM-11650\NDSI_MODIS', r'\\CAL-GOA-APM-422\NDSI_MODIS']
        serv_names = ['Wildfire', 'FEWS DEV server', 'FEWS UAT server Edmonton', 'FEWS UAT server Calgary',
                      'FEWS PROD server Edmonton (FSS)', 'FEWS PROD server Calgary (FSS)']

        try:
            for server in servers:
                print("s:", server)
                cgf_files = shutil.copy2(afile, os.path.join(server, os.path.basename(afile)))
                print("File copied successfully.")

        # If source and destination are same
        except shutil.SameFileError:
            print("Source and destination represents the same file.")

        # If there is any permission issue
        except PermissionError:
            print("Permission denied.")

        # For other errors
        except:
            print("Error occurred while copying file.")


        return


    def proc_MOD10A2_HDF(self, path):
        ''' processing _MOD10A2 hdf files'''
        os.chdir(path)
        # hdf_folder = path + '/hdf'
        #
        # cur_folder = path + '/current'
        # outdir = path + '/ARCHIVE'
        # create folder for rasters in Geographic grid
        path_geo = path + '/geo'
        self.mkdir(path_geo)
        # print("createing path 10tm: ", path_10tm)
        mosaic = path + '/mosaics'  # to hold mosaics of different layers
        self.mkdir(mosaic)
        # MOD10A2 layers......
        sc_8day_fold = path_geo + '/sc_8day'  # folder for NDSI
        ext_SC_fold = path_geo  + '/ext_sc'  # folder for snow_cover
        self.mkdir(sc_8day_fold)
        self.mkdir(ext_SC_fold)

        hdf_List = glob.glob("*.hdf")
        print("hdflst:", hdf_List)
        gdal.UseExceptions()
        x = 0
        for afile in hdf_List:

            print("i: ", afile)
            filename = afile
            g = gdal.Open(filename)
            # g should now be a GDAL dataset, but if the file isn't found
            if g is None:
                print("Problem opening file %s!" % filename)
            else:
                print("File %s opened fine" % filename)
            # get subdataset names...
            subdatasets = g.GetSubDatasets()
            for fname, name in subdatasets:
                print("\t", fname)
            #
            # Let's create a list with the selected layer names
            # selected_layers = ["NDSI_Snow_Cover", "NDSI_Snow_Cover_Basic_QA", "NDSI_Snow_Cover_Algorithm_Flags_QA","NDSI"] #MOD10A1
            selected_layers = ["Maximum_Snow_Extent", "Eight_Day_Snow_Cover"]  ##MOD10A1 layers
            # We will store the data in a dictionary
            # Initialise an empty dictionary
            data = {}
            file_template = 'HDF4_EOS:EOS_GRID:"%s":MOD_Grid_Snow_500m:%s'
            # HDF4_EOS:EOS_GRID:"MOD09A1N.A2020134.h10v03.061.NRT.hdf":MOD_Grid_500m_Surface_Reflectance:sur_refl_b01
            # HDF4_EOS:EOS_GRID:"MOD10A1.A2021015.h10v03.006.2021017023910.hdf":MOD_Grid_Snow_500m:NDSI_Snow_Cover
            # HDF4_EOS:EOS_GRID:"MOD10A2.A2021009.h10v03.006.2021018033643.hdf":MOD_Grid_Snow_500m:Maximum_Snow_Extent
            for i, layer in enumerate(selected_layers):
                this_file = file_template % (filename, layer)
                print("thisfile: ", this_file)
                print("Opening Layer %d: %s" % (i + 1, this_file))
                g = gdal.Open(this_file)
                if g is None:
                    raise IOError
                init_geot = g.GetGeoTransform()
                w = g.RasterXSize
                h = g.RasterYSize
                print("initgeo:", init_geot)
                data[layer] = g.ReadAsArray()
                print("lay_shape:", w, h)
                print("\t>>> Read %s!" % layer)
                print("proecessing: layer: i: ", i)
                parts = this_file.split(':')
                ime = parts[2][1:-5] + '_' + parts[4]
                norm_sph = ime.replace('.', '_') + '.tif'
                print("input 4 extport2sinu", this_file, norm_sph)
                # run function
                self.modis2normsphere(this_file, norm_sph)
                #
                # # run sinus to geogrpahic
                out_84 = norm_sph[:-4] + '_wgs84.tif'
               # print("Out_84 ime:", out_84)
                print("input 4 sin2wgs84", norm_sph,  out_84)
                # print("sin_inp_geo:", geo_tin)
                self.sin_wgs84(norm_sph, os.path.join(path_geo,out_84) ) # ISKLJUCI/UKLJUCI
                # extract date from the file
                adate = self.getdDate_from_hdfJD(out_84,1)  # Check this

                print("adate:", adate)

            try:
                gdal.UseExceptions()
                list_84 = glob.glob(path_geo + '/*tif')
                # copy the same products of 4 tiles into folders for mosaicking
                sc_8day_files = [shutil.copy2(x, os.path.join(sc_8day_fold, os.path.basename(x))) for x in list_84 if
                                 'Eight_Day_Snow_Cover_wgs84' in x]
                max_ext_files = [shutil.copy2(x, os.path.join(ext_SC_fold, os.path.basename(x))) for x in list_84 if
                                 'Maximum_Snow_Extent_wgs84' in x]

            except Exception as e:
                print("Exception: ", e)
                exit(1)
        #     #***************************************
        # # #GET DATE FROM CURRENT IMAGE AND UPDATE *.DATES FILE
        # # #***************************
        #
        print("update 'dates' file..")
        ymd = adate.split('_')[1]
        print("ymd:", ymd)
        date_file = 'C:/_LOCALdata/prj_2021/SNOW_modis/ARCHIVE.dates'
        ##ISLJUJCI/  UKLJUCI
        print("updating dates file....")
        with open(date_file, 'a') as f:
            # f.write(ymd + '\n')
            f.write(f'\n{ymd}')  # for python 3

        # mosaic all 4 tiles covering AB
        mos_8day = os.path.join(mosaic, ymd + '_8day.tif')  #
        mos_max_ext = os.path.join(mosaic, ymd + "_sc.tif")

        # gdal.BuildVRT(mos_8day, glob.glob("C:/_LOCALdata/prj_2021/SNOW_modis/MOD/MOD10A2/hdf/geo/sc_8day/*.tif"))
        # gdal.BuildVRT(mos_max_ext, glob.glob("C:/_LOCALdata/prj_2021/SNOW_modis/MOD/MOD10A2/hdf/geo/ext_sc/*.tif"))
        gdal.BuildVRT(mos_8day, glob.glob(sc_8day_fold +"/*.tif"))
        gdal.BuildVRT(mos_max_ext, glob.glob(ext_SC_fold +"/*.tif"))

        # clear the folder
        print("curwd:", os.getcwd())
        files = glob.glob(os.getcwd() + '/*.tif')
        for f in files:
            os.remove(f)
        files10 = glob.glob(path_geo+"/*.tif")
        for f10 in files10:
            os.remove(f10)

        return

    def proc_MOD10A1_HDF(self, path_MOD10A1):
        ''' processing _MOD10A1 hdf files'''
        os.chdir(path_MOD10A1)

        path_geo = path_MOD10A1 + '/geo'
        self.mkdir(path_geo)
        # print("createing path 10tm: ", path_10tm)
        mosaic = path_MOD10A1 + '/mosaics'  # to hold mosaics of different layers
        self.mkdir(mosaic)
        # MOD10A1 layers......
        ndsi_fold = path_geo + '/NDSI'  # folder for NDSI
        ndsi_SC_fold = path_geo + '/NDSI_SC'  # folder for snow_cover
        flags_fold = path_geo + '/flags'  # folder for NDSI
        qa_fold = path_geo + '/b_qa'  # folder for snow_cover
        self.mkdir(ndsi_fold)
        self.mkdir(ndsi_SC_fold)
        self.mkdir(flags_fold)
        self.mkdir(qa_fold)

        hdf_List = glob.glob(path_MOD10A1 +"/*.hdf")
        print("hdflst:", hdf_List)
        gdal.UseExceptions()
        x = 0
        for afile in hdf_List:

            print("i: ", afile)
            filename = afile
            g = gdal.Open(filename)
            # g should now be a GDAL dataset, but if the file isn't found
            if g is None:
                print("Problem opening file %s!" % filename)
            else:
                print("File %s opened fine" % filename)
            # get subdataset names...
            subdatasets = g.GetSubDatasets()
            for fname, name in subdatasets:
                print("\t", fname)
            #
            # Let's create a list with the selected layer names
            selected_layers = ["NDSI_Snow_Cover", "NDSI_Snow_Cover_Basic_QA", "NDSI_Snow_Cover_Algorithm_Flags_QA","NDSI"] #MOD10A1

            # We will store the data in a dictionary
            # Initialise an empty dictionary
            data = {}
            file_template = 'HDF4_EOS:EOS_GRID:"%s":MOD_Grid_Snow_500m:%s'
            # HDF4_EOS:EOS_GRID:"MOD09A1N.A2020134.h10v03.061.NRT.hdf":MOD_Grid_500m_Surface_Reflectance:sur_refl_b01
            # HDF4_EOS:EOS_GRID:"MOD10A1.A2021015.h10v03.006.2021017023910.hdf":MOD_Grid_Snow_500m:NDSI_Snow_Cover

            for i, layer in enumerate(selected_layers):
                this_file = file_template % (filename, layer)
                print("thisfile: ", this_file)
                print("Opening Layer %d: %s" % (i + 1, this_file))
                g = gdal.Open(this_file)
                if g is None:
                    raise IOError
                init_geot = g.GetGeoTransform()
                w = g.RasterXSize
                h = g.RasterYSize
                print("initgeo:", init_geot)
                data[layer] = g.ReadAsArray()
                print("lay_shape:", w, h)
                print("\t>>> Read %s!" % layer)
                print("proecessing: layer: i: ", i)
                parts = this_file.split(':')
                print("parts:",parts)
                ime = parts[2][1:-5] + '_' + parts[4]
                print("ime:", ime)
                norm_sph = ime.replace('.', '_') + '.tif'
                print("input 4 extport2sinu", this_file, norm_sph)
                # run function
               #  self.modis2normsphere(this_file, norm_sph)
               #  #
               #  # # run sinus to geogrpahic
               #  out_84 = norm_sph[:-4] + '_wgs84.tif'
               # # print("Out_84 ime:", out_84)
               #  print("input 4 sin2wgs84", norm_sph,  out_84)
               #  # print("sin_inp_geo:", geo_tin)
               #  self.sin_wgs84(norm_sph, os.path.join(path_geo,out_84) ) # ISKLJUCI/UKLJUCI
               #  # extract date from the file
               #
               # # adate = self.getDate_from_hdfJD(out_84,1)  # Check this
               #  adate = 'A2021015_20210115'
               #  print("adate:", adate)

            # try:
            #     gdal.UseExceptions()
            #     list_84 = glob.glob(path_geo + '/*tif')
            #     # copy the same products of 4 tiles into folders for mosaicking
            #
            #     l_10tm = glob.glob(path_10tm + '/*tif')
            #     # copy the same products of 4 tiles into folders for mosaicking
            #     ndsi_files = [shutil.copy2(x, os.path.join(ndsi_fold, os.path.basename(x))) for x in list_84 if
            #                   'NDSI_10tm' in x]
            #     sc_files = [shutil.copy2(x, os.path.join(ndsi_SC_fold, os.path.basename(x))) for x in list_84 if
            #                 'Cover_10tm' in x]
            #     flags_files = [shutil.copy2(x, os.path.join(flags_fold, os.path.basename(x))) for x in list_84 if
            #                    'Algorithm_Flags_QA' in x]
            #     qa_files = [shutil.copy2(x, os.path.join(qa_fold, os.path.basename(x))) for x in list_84 if 'Basic' in x]

        #     except Exception as e:
        #         print("Exception: ", e)
        #         exit(1)
        # #     #***************************************
        # # # #GET DATE FROM CURRENT IMAGE AND UPDATE *.DATES FILE
        # # # #***************************
        # #
        # print("update 'dates' file..")
        # ymd = adate.split('_')[1]
        # print("ymd:", ymd)
        # date_file = 'C:/_LOCALdata/prj_2021/SNOW_modis/ARCHIVE.dates'
        # ##ISLJUJCI/  UKLJUCI
        # print("updating dates file....")
        # with open(date_file, 'a') as f:
        #     # f.write(ymd + '\n')
        #     f.write(f'\n{ymd}')  # for python 3
        #
        # # mosaic all 4 tiles covering AB
        # mos_8day = os.path.join(mosaic, ymd + '_8day.tif')  #
        # mos_max_ext = os.path.join(mosaic, ymd + "_sc.tif")
        #
        # # gdal.BuildVRT(mos_8day, glob.glob("C:/_LOCALdata/prj_2021/SNOW_modis/MOD/MOD10A2/hdf/geo/sc_8day/*.tif"))
        # # gdal.BuildVRT(mos_max_ext, glob.glob("C:/_LOCALdata/prj_2021/SNOW_modis/MOD/MOD10A2/hdf/geo/ext_sc/*.tif"))
        # gdal.BuildVRT(mos_8day, glob.glob(sc_8day_fold +"/*.tif"))
        # gdal.BuildVRT(mos_max_ext, glob.glob(ext_SC_fold +"/*.tif"))
        #
        # # clear the folder
        # print("curwd:", os.getcwd())
        # files = glob.glob(os.getcwd() + '/*.tif')
        # for f in files:
        #     os.remove(f)
        # files10 = glob.glob(path_geo+"/*.tif")
        # for f10 in files10:
        #     os.remove(f10)

        return



    def proc_MOD09_HDF(self, path):
        ''' processing _MOD09 hdf files'''
        os.chdir(path)
        # hdf_folder = path + '/hdf'
        #
        # cur_folder = path + '/current'
        # outdir = path + '/ARCHIVE'
        # create folder for rasters in Geographic grid
        path_geo = path + '/geo'
        self.mkdir(path_geo)
        # print("createing path 10tm: ", path_10tm)
        mosaic = path + '/mosaics'  # to hold mosaics of different layers
        self.mkdir(mosaic)
        # MOD9 layers......

        b03 = path_geo + '/b03'  # folder for b4
        b02 = path_geo + '/b02'  # folder for b05
        b01= path_geo + '/b01'  # folder for b7

        b04 = path_geo + '/b04'  # folder for b4
        b05 = path_geo + '/b05'  # folder for b05
        b07 = path_geo + '/b07'  # folder for b7
        b06 = path_geo  + '/b06'  # folder for b06
        self.mkdir(b06)
        self.mkdir(b07)
        self.mkdir(b05)
        self.mkdir(b04)
        self.mkdir(b03)
        self.mkdir(b02)
        self.mkdir(b01)

        hdf_List = glob.glob("*.hdf")
        print("hdflst:", hdf_List)
        gdal.UseExceptions()
        x = 0
        for afile in hdf_List:

            print("i: ", afile)
            filename = afile
            g = gdal.Open(filename)
            # g should now be a GDAL dataset, but if the file isn't found
            if g is None:
                print("Problem opening file %s!" % filename)
            else:
                print("File %s opened fine" % filename)
            # get subdataset names...
            subdatasets = g.GetSubDatasets()
            for fname, name in subdatasets:
                print("\t", fname)
            #
            # Let's create a list with the selected layer names


            selected_layers = ["sur_refl_b01","sur_refl_b02","sur_refl_b03","sur_refl_b04","sur_refl_b05","sur_refl_b06","sur_refl_b07"]

             ##MOD10A1 layers
            # We will store the data in a dictionary
            # Initialise an empty dictionary
            data = {}
            file_template = 'HDF4_EOS:EOS_GRID:"%s":MOD_Grid_500m_Surface_Reflectance:%s'
            # HDF4_EOS:EOS_GRID:"MOD09A1N.A2020134.h10v03.061.NRT.hdf":MOD_Grid_500m_Surface_Reflectance:sur_refl_b01

            for i, layer in enumerate(selected_layers):
                this_file = file_template % (filename, layer)
                print("thisfile: ", this_file)
                print("Opening Layer %d: %s" % (i + 1, this_file))
        #         g = gdal.Open(this_file)
                if g is None:
                    raise IOError
                init_geot = g.GetGeoTransform()
                w = g.RasterXSize
                h = g.RasterYSize
                print("initgeo:", init_geot)
                data[layer] = g.ReadAsArray()
                print("lay_shape:", w, h)
                print("\t>>> Read %s!" % layer)
                print("proecessing: layer: i: ", i)
                parts = this_file.split(':')
                ime = parts[2][1:-5] + '_' + parts[4]
                norm_sph = ime.replace('.', '_') + '.tif'
                print("input 4 extport2sinu", this_file, norm_sph)
                # run function
                self.modis2normsphere(this_file, norm_sph)
                #
                # # run sinus to geogrpahic
                out_84 = norm_sph[:-4] + '_wgs84.tif'
               # print("Out_84 ime:", out_84)
                print("input 4 sin2wgs84", norm_sph,  out_84)
                # print("sin_inp_geo:", geo_tin)
                self.sin_wgs84(norm_sph, os.path.join(path_geo,out_84) ) # ISKLJUCI/UKLJUCI
                # extract date from the file
                #adate = self.getdate_from_hdf_jd(out_84)  # Check this
                adate = self.getDate_from_hdfJD(out_84,1)  # Check this
                ymd = adate.split('_')[1]
                print("ymd:", ymd)
                print("adate:", adate)
        #
            try:
                gdal.UseExceptions()
                list_84 = glob.glob(path_geo + '/*tif')
                # copy the same products of 4 tiles into folders for mosaicking
                b7_files = [shutil.copy2(x, os.path.join( b07 , os.path.basename(x))) for x in list_84 if
                                 'sur_refl_b07_wgs84' in x]
                b5_files = [shutil.copy2(x, os.path.join(b05, os.path.basename(x))) for x in list_84 if
                            'sur_refl_b05_wgs84' in x]

                b2_files = [shutil.copy2(x, os.path.join( b02 , os.path.basename(x))) for x in list_84 if
                                 'sur_refl_b02_wgs84' in x]

            except Exception as e:
                print("Exception: ", e)
                exit(1)
        # #     #***************************************
        # # # #GET DATE FROM CURRENT IMAGE AND UPDATE *.DATES FILE
        # # # #***************************
        # print("b5files:",b5_files )
        # print("b7files:", b7_files)
        print("update 'dates' file..")

        date_file = r'U:/RS_Task_Workspaces/NDWI/ARCHIVE.dates'
        ##ISLJUJCI/  UKLJUCI
        print("updating dates file....")
        with open(date_file, 'a') as f:
            # f.write(ymd + '\n')
            f.write(f'\n{ymd}')  # for python 3

        # mosaic all 4 tiles covering AB
        mos_b7 = os.path.join(mosaic, ymd + '_b7.tif')
        mos_b5 = os.path.join(mosaic, ymd + '_b5.tif')#
        mos_b2 = os.path.join(mosaic, ymd + "_b2.tif")

        # gdal.BuildVRT(mos_8day, glob.glob("C:/_LOCALdata/prj_2021/SNOW_modis/MOD/MOD10A2/hdf/geo/sc_8day/*.tif"))
        # gdal.BuildVRT(mos_max_ext, glob.glob("C:/_LOCALdata/prj_2021/SNOW_modis/MOD/MOD10A2/hdf/geo/ext_sc/*.tif"))
        gdal.BuildVRT(mos_b2, glob.glob(b02 +"/*.tif"))
        gdal.BuildVRT(mos_b5, glob.glob(b05 + "/*.tif"))
        gdal.BuildVRT(mos_b7, glob.glob(b07 +"/*.tif"))
        # *********************************************
        # now it is time to apply corrections to data file including Scalling factor  and additional offset
        # lst_day & LST_night (.02, NA)  --> 14054*0.02 = 281 K
        # day_view angle ( 1.0, -65)   --> (130*1) -65 = 65
        # emis_1, emis_2 (.002, 0.49) -->  240 * 0.002 + 0.49= 0.988
        #band reflectance --0.0001 scale for 500 m bands
        # Scale values of MODIS imagery stack
        #modis_bands_pre_scaled = modis_bands * 0.0001
        # ********************************
        # # reprojcet to 10 tm
        mos_b2_10tm = mos_b2.replace("b2","b2_10TM")
        mos_b5_10tm = mos_b5.replace("b5", "b5_10TM")
        self.wgs84_epsg3400(mos_b2, mos_b2_10tm)  # ISKLJUCI/UKLJUCI
        self.wgs84_epsg3400(mos_b5, mos_b5_10tm)  # ISKLJUCI/UKLJUCI
        #*****************************************************
        #create NDWI
        ndwi_name = os.path.join(os.path.dirname(mos_b2_10tm), "ndwi_" + os.path.basename(mos_b2_10tm).replace('_b2_', '_'))

        self.array2raster(ndwi_name, mos_b5_10tm, mos_b2_10tm)

        #**********************************
        # distribution of files to servers of Environment / Parks and Wildfire (KEvin Keats)
        # self.copy2servers(outclip)
        # # apply colors to a CFG file
        # self.aaply_colors(dst)
        # ******************************************
        # clear the folder
        # print("curwd:", os.getcwd())
        # files = glob.glob(os.getcwd() + '/*.tif')
        # for f in files:
        #     os.remove(f)
        # files10 = glob.glob(path_geo+"/*.tif")
        # for f10 in files10:
        #     os.remove(f10)

        return ndwi_name



    def proc_MOD10A1F_HDF(self, path):
        ''' processing _MOD10A1F hdf files'''
        os.chdir(path)
        # hdf_folder = path + '/hdf'

        path_geo = path + '/geo'
        self.mkdir(path_geo)
        # print("createing path 10tm: ", path_10tm)
        mosaic = path + '/mosaics'  # to hold mosaics of different layers
        self.mkdir(mosaic)
        # MOD9 layers......

        cgf = path_geo + '/CGF_NDSI_SC'  # folder for cloud gap filled
        cloud_pers = path_geo + '/cloud_persis'  # folder for cloud persistence
        qa = path_geo + '/qa'  # folder for basic qa
        flags = path_geo + '/flags'  # folder for flags
        sc = path_geo + '/sc'  # folder for daily snow cov
        self.mkdir(cgf)
        self.mkdir(cloud_pers)
        self.mkdir(qa)
        self.mkdir(flags)
        self.mkdir(sc)

        hdf_List = glob.glob("*.hdf")
        print("hdflst:", hdf_List)
        gdal.UseExceptions()
        x = 0
        for afile in hdf_List:
            print("i: ", afile)
            filename = afile
            g = gdal.Open(filename)
            # g should now be a GDAL dataset, but if the file isn't found
            if g is None:
                print("Problem opening file %s!" % filename)
            else:
                print("File %s opened fine" % filename)
            # get subdataset names...
            subdatasets = g.GetSubDatasets()
            for fname, name in subdatasets:
                print("\t", fname)
            #
            # Let's create a list with the selected layer names
            # RADI     selected_layers = ["CGF_NDSI_Snow_Cover","Cloud_Persistence","Basic_QA","Algorithm_Flags_QA","MYD10A1_NDSI_Snow_Cover"]
            selected_layers = ['CGF_NDSI_Snow_Cover', 'Cloud_Persistence', 'Basic_QA', 'Algorithm_Flags_QA']
            #

            ##MOD10A1F layers
            # We will store the data in a dictionary
            # Initialise an empty dictionary
            data = {}
            file_template = 'HDF4_EOS:EOS_GRID:"%s":MOD_Grid_Snow_500m:%s'
             #HDF4_EOS:EOS_GRID:"MOD10A1F.A2021015.h10v03.061.2021021032954.hdf":MOD_Grid_Snow_500m:CGF_NDSI_Snow_Cover
            for i, layer in enumerate(selected_layers):
                this_file = file_template % (filename, layer)
                print("thisfile: ", this_file)
                print("Opening Layer %d: %s" % (i + 1, this_file))
        #         g = gdal.Open(this_file)
                if g is None:
                    raise IOError
                init_geot = g.GetGeoTransform()
                w = g.RasterXSize
                h = g.RasterYSize
                print("initgeo:", init_geot)
                data[layer] = g.ReadAsArray()
                print("lay_shape:", w, h)
                print("\t>>> Read %s!" % layer)
                print("proecessing: layer: i: ", i)
                parts = this_file.split(':')
                ime = parts[2][1:-5] + '_' + parts[4]
                norm_sph = ime.replace('.', '_') + '.tif'
                print("input 4 extport2sinu", this_file, norm_sph)
                # run function
                self.modis2normsphere(this_file, norm_sph)
                #
                # # run sinus to geogrpahic
                out_84 = norm_sph[:-4] + '_wgs84.tif'
               # print("Out_84 ime:", out_84)
                print("input 4 sin2wgs84", norm_sph,  out_84)
                # print("sin_inp_geo:", geo_tin)
                self.sin_wgs84(norm_sph, os.path.join(path_geo,out_84) ) # ISKLJUCI/UKLJUCI
                # extract date from the file
                #adate = self.getdate_from_hdf_jd(out_84)  # Check this
                adate = self.getDate_from_hdfJD(out_84,1)  # Check this

                print("adate:", adate)

            try:
                gdal.UseExceptions()
                print("copy images into coresponding folders.....")
                list_84 = glob.glob(path_geo + '/*tif')
                # copy the same products of 4 tiles into folders for mosaicking

                cgf_files = [shutil.copy2(x, os.path.join( cgf , os.path.basename(x))) for x in list_84 if
                                 'CGF' in x]
                cloud_pers_files = [shutil.copy2(x, os.path.join(cloud_pers, os.path.basename(x))) for x in list_84 if
                          'Persistence_wgs84' in x]
                qa_files = [shutil.copy2(x, os.path.join(qa, os.path.basename(x))) for x in list_84 if
                          'Basic_QA_wgs84' in x]
                flags_files = [shutil.copy2(x, os.path.join(flags, os.path.basename(x))) for x in list_84 if
                           'Flags_QA_wgs84' in x]
                # sc_files = [shutil.copy2(x, os.path.join(sc, os.path.basename(x))) for x in list_84 if
                #             'MYD10A1_NDSI_Snow_Cover_wgs84' in x]
                # sc_files = [shutil.copy2(x, os.path.join(sc, os.path.basename(x))) for x in list_84 if
                #            '"%s"_NDSI_Snow_Cover_wgs84' in x]

            except Exception as e:
                print("Exception: ", e)
                exit(1)
        # #     #***************************************
        # # # #GET DATE FROM CURRENT IMAGE AND UPDATE *.DATES FILE
        # # # #***************************
        print("update 'dates' file..")
        ymd = adate.split('_')[1]
        print("ymd:", ymd)
        date_file = 'C:/_LOCALdata/prj_2021/SNOW_modis/ARCHIVE.dates'
        # ##ISLJUJCI/  UKLJUCI
        print("updating dates file....")
        with open(date_file, 'a') as f:
            # f.write(ymd + '\n')
            f.write(f'\n{ymd}')  # for python 3

        # mosaic all 4 tiles covering AB
        #name the output files
        mos_cgf = os.path.join(mosaic, "TerraAqua_CFG_" +ymd +'23.tif')  #
        mos_cloud_pers = os.path.join(mosaic, ymd + "_cloud_pers.tif")
        mos_qa = os.path.join(mosaic, ymd + '_qa.tif')  #
        mos_flags = os.path.join(mosaic, ymd + "_flags.tif")
        #mos_sc = os.path.join(mosaic, ymd + '_mod10_sc.tif')  #
                #build mosaics....
        gdal.BuildVRT(mos_cgf, glob.glob(cgf + "/*.tif"))
        gdal.BuildVRT(mos_cloud_pers, glob.glob(cloud_pers + "/*.tif"))
        gdal.BuildVRT(mos_qa, glob.glob(qa +"/*.tif"))
        gdal.BuildVRT(mos_flags, glob.glob(flags +"/*.tif"))
       # gdal.BuildVRT(mos_sc, glob.glob(sc + "/*.tif"))

        #now it is time to apply corrections to data file including Scalling factor  and additional offset
        #lst_day & LST_night (.02, NA)  --> 14054*0.02 = 281 K
        #day_view angle ( 1.0, -65)   --> (130*1) -65 = 65
        #emis_1, emis_2 (.002, 0.49) -->  240 * 0.002 + 0.49= 0.988
        #********************************
        #reprojcet to 10 tm
        # mos_cgf_10tm = mos_cgf.replace("Class","Class10t")
        # self.wgs84_epsg3400(mos_cgf, dst)  # ISKLJUCI/UKLJUCI

        # clip the product to AB boundaries
        dst =  mos_cgf.replace("CFG", "NDSI_Class" )
        # Xmin, Ymin, Xmax, Ymax = (48.7, -109.6, 60.08, -120.8)
        # self.run_gridClip(Xmin, Ymin, Xmax, Ymax, mos_cgf, dst)
        outclip = self.wgs84_epsg3400(mos_cgf, dst)

        #distribution of files to servers of Environment / Parks and Wildfire (KEvin Keats)
        self.copy2servers( outclip )
        # apply colors to a CFG file
        self.aaply_colors(dst)

        # clear the folder
        print("curwd:", os.getcwd())
        files = glob.glob(os.getcwd() + '/*.tif')
        for f in files:
            os.remove(f)
        files10 = glob.glob(path_geo+"/*.tif")
        for f10 in files10:
            os.remove(f10)

        return

    def proc_MOYD10A1F_HDF(self, path):
        ''' processing _MOD10A1F nad MYD10A1F hdf files'''
        os.chdir(path)
        # hdf_folder = path + '/hdf'

        path_geo = path + '/geo'
        self.mkdir(path_geo)
        # print("createing path 10tm: ", path_10tm)
        mosaic = path + '/mosaics'  # to hold mosaics of different layers
        self.mkdir(mosaic)
        # MOD9 layers......

        cgf = path_geo + '/CGF_NDSI_SC'  # folder for cloud gap filled
        cloud_pers = path_geo + '/cloud_persis'  # folder for cloud persistence
        qa = path_geo + '/qa'  # folder for basic qa
        flags = path_geo + '/flags'  # folder for flags
        sc = path_geo + '/sc'  # folder for daily snow cov
        self.mkdir(cgf)
        self.mkdir(cloud_pers)
        self.mkdir(qa)
        self.mkdir(flags)
        self.mkdir(sc)

        hdf_List = glob.glob("*.hdf")
        print("hdflst:", hdf_List)
        if "MOD" in hdf_List[0]:
            key = "MOD10A1"
        else:
            key = "MYD10A1"
        print("key:",key)
        gdal.UseExceptions()
        x = 0
        for afile in hdf_List:
            print("i: ", afile)
            filename = afile
            g = gdal.Open(filename)
            # g should now be a GDAL dataset, but if the file isn't found
            if g is None:
                print("Problem opening file %s!" % filename)
            else:
                print("File %s opened fine" % filename)
            # get subdataset names...
            subdatasets = g.GetSubDatasets()
            for fname, name in subdatasets:
                print("\t", fname)
            #
            # Let's create a list with the selected layer names
            # RADI     selected_layers = ["CGF_NDSI_Snow_Cover","Cloud_Persistence","Basic_QA","Algorithm_Flags_QA","MYD10A1_NDSI_Snow_Cover"]
            selected_layers = ['CGF_NDSI_Snow_Cover', 'Cloud_Persistence', 'Basic_QA', 'Algorithm_Flags_QA',key+"_NDSI_Snow_Cover"]
            #

            ##MOD10A1F layers
            # We will store the data in a dictionary
            # Initialise an empty dictionary
            data = {}
            file_template = 'HDF4_EOS:EOS_GRID:"%s":MOD_Grid_Snow_500m:%s'
             #HDF4_EOS:EOS_GRID:"MOD10A1F.A2021015.h10v03.061.2021021032954.hdf":MOD_Grid_Snow_500m:CGF_NDSI_Snow_Cover
            for i, layer in enumerate(selected_layers):
                this_file = file_template % (filename, layer)
                print("thisfile: ", this_file)
                print("Opening Layer %d: %s" % (i + 1, this_file))
        #         g = gdal.Open(this_file)
                if g is None:
                    raise IOError
                init_geot = g.GetGeoTransform()
                w = g.RasterXSize
                h = g.RasterYSize
                print("initgeo:", init_geot)
                data[layer] = g.ReadAsArray()
                print("lay_shape:", w, h)
                print("\t>>> Read %s!" % layer)
                print("proecessing: layer: i: ", i)
                parts = this_file.split(':')
                ime = parts[2][1:-5] + '_' + parts[4]
                norm_sph = ime.replace('.', '_') + '.tif'
                print("input 4 extport2sinu", this_file, norm_sph)
                # run function
                self.modis2normsphere(this_file, norm_sph)
                #
                # # run sinus to geogrpahic
                out_84 = norm_sph[:-4] + '_wgs84.tif'
               # print("Out_84 ime:", out_84)
                print("input 4 sin2wgs84", norm_sph,  out_84)
                # print("sin_inp_geo:", geo_tin)
                self.sin_wgs84(norm_sph, os.path.join(path_geo,out_84) ) # ISKLJUCI/UKLJUCI
                # extract date from the file
                #adate = self.getdate_from_hdf_jd(out_84)  # Check this
                adate = self.getDate_from_hdfJD(out_84,1)  # Check this

                print("adate:", adate)
       #
            try:
                gdal.UseExceptions()
                print("copy images into coresponding folders.....")
                list_84 = glob.glob(path_geo + '/*tif')
                # copy the same products of 4 tiles into folders for mosaicking

                cgf_files = [shutil.copy2(x, os.path.join( cgf , os.path.basename(x))) for x in list_84 if
                                 'CGF' in x]
                cloud_pers_files = [shutil.copy2(x, os.path.join(cloud_pers, os.path.basename(x))) for x in list_84 if
                          'Persistence_wgs84' in x]
                qa_files = [shutil.copy2(x, os.path.join(qa, os.path.basename(x))) for x in list_84 if
                          'Basic_QA_wgs84' in x]
                flags_files = [shutil.copy2(x, os.path.join(flags, os.path.basename(x))) for x in list_84 if
                           'Flags_QA_wgs84' in x]
                # sc_files = [shutil.copy2(x, os.path.join(sc, os.path.basename(x))) for x in list_84 if
                #             'MYD10A1_NDSI_Snow_Cover_wgs84' in x]
                sc_files = [shutil.copy2(x, os.path.join(sc, os.path.basename(x))) for x in list_84 if
                           key+ '_NDSI_Snow_Cover_wgs84' in x]

            except Exception as e:
                print("Exception: ", e)
                exit(1)
        # #     #***************************************
        # # # #GET DATE FROM CURRENT IMAGE AND UPDATE *.DATES FILE
        # # # #***************************
        print("update 'dates' file..")
        ymd = adate.split('_')[1]
        print("ymd:", ymd)
        date_file = 'C:/_LOCALdata/prj_2021/SNOW_modis/ARCHIVE.dates'
        # ##ISLJUJCI/  UKLJUCI
        print("updating dates file....")
        with open(date_file, 'a') as f:
            # f.write(ymd + '\n')
            f.write(f'\n{ymd}')  # for python 3

        # mosaic all 4 tiles covering AB
        #name the output files

        mos_cgf = os.path.join(mosaic, "TerraAqua_NDSI_CFG_" + ymd + '1423.tif')  #
        mos_cloud_pers = os.path.join(mosaic, key+ ymd + "_cloud_pers.tif")
        mos_qa = os.path.join(mosaic, key + ymd + '_qa.tif')  #
        mos_flags = os.path.join(mosaic, key + ymd + "_flags.tif")
        mos_sc = os.path.join(mosaic, key+ ymd +  '_sc.tif')  #
                #build mosaics....
        gdal.BuildVRT(mos_cgf, glob.glob(cgf + "/*.tif"))
        gdal.BuildVRT(mos_cloud_pers, glob.glob(cloud_pers + "/*.tif"))
        gdal.BuildVRT(mos_qa, glob.glob(qa +"/*.tif"))
        gdal.BuildVRT(mos_flags, glob.glob(flags +"/*.tif"))
        gdal.BuildVRT(mos_sc, glob.glob(sc + "/*.tif"))

        #now it is time to apply corrections to data file including Scalling factor  and additional offset
        #lst_day & LST_night (.02, NA)  --> 14054*0.02 = 281 K
        #day_view angle ( 1.0, -65)   --> (130*1) -65 = 65
        #emis_1, emis_2 (.002, 0.49) -->  240 * 0.002 + 0.49= 0.988

        # clip the product to AB boundaries
        dst = mos_cgf.replace("CFG", "CFG_Class")
        shape = "C:/_LOCALdata/prj_2021/SNOW_modis/FINAL/out/AB_prov.shp"

       # self.gridClip_shp(mos_cgf, dst,  shape)
        #outclip = self.run_gridClip(-120.6, 48.7, -109.5, 60.1, mos_cgf, dst)
        #Reporject ot 10 TM and clip into desired extent
        outclip = self.wgs84_epsg3400( mos_cgf, dst)
        #clip using provincial boundary

        #Distribute files to Environment/Parks and KEvin Keath
        #\\GOA\Shared\AF\Wildfire_Geospatial\Data_Repository\MODIS\Snow_and_Ice_Composite\MODIS_TerraAqua_Tiffs\original
        # # apply colors to a CFG file
        self.aaply_colors(outclip)

        # clear the folder
        print("curwd:", os.getcwd())
        files = glob.glob(os.getcwd() + '/*.tif')
        for f in files:
            os.remove(f)
        files10 = glob.glob(path_geo+"/*.tif")
        for f10 in files10:
            os.remove(f10)

        return

    # def run_gridClip(self, Xmin, Ymin, Xmax, Ymax, src, dst):
    #
    #     print("xmin, ymax,xmax,ymin", Xmin, Ymax, Xmax, Ymin)
    #     cmd = "-te xmin ymin xmax ymax"  # just for example!
    #     # gdalwarp -te xmin ymin xmax ymax -dstnodata "-9999.0" -dstnodata "-9999.0"
    #     gdal_Warp = r'C:\Program Files\GDAL\gdalwarp.exe'
    #     cmd = "-te "
    #     fullCmd = ' '.join(
    #         [gdal_Warp, cmd, str(Xmin), str(Ymin), str(Xmax), str(Ymax), " -dstnodata -9999.0 ", src, dst])
    #     print("full comanda", fullCmd)
    #     subprocess.Popen(fullCmd)
    #     return dst

    def aaply_colors(self,file_in):
        ''' apply colors to a raster file'''

        ds = gdal.Open(file_in, 1)
        band = ds.GetRasterBand(1)
        # create color table
        color = gdal.ColorTable()
        ## set color for each value
        # NO snow
        color.SetColorEntry(0, (0, 10, 0))
        #snow
        color.SetColorEntry(0, (0,255,0))
        color.SetColorEntry(1, (250,169,169))
        color.SetColorEntry(2, (250,169,168) )
        color.SetColorEntry(3, (250,169,167))
        color.SetColorEntry(4,(250,169,166))
        color.SetColorEntry(5,(250,169,165) )
        color.SetColorEntry(6,(250,169,164) )
        color.SetColorEntry(7, (250,169,163))
        color.SetColorEntry(8,(250,169,162))
        color.SetColorEntry(9, (250,169,161))
        color.SetColorEntry(10, (250,169,160))
        color.SetColorEntry(11, (250,214,160))
        color.SetColorEntry(12, (250,214,161))
        color.SetColorEntry(13, (250,214,162))
        color.SetColorEntry(14, (250,214,163))
        color.SetColorEntry(15, (250,214,164))
        color.SetColorEntry(16, (250,214,165))
        color.SetColorEntry(17, (250,214,166))
        color.SetColorEntry(18, (250,214,167))
        color.SetColorEntry(19, (250,214,168))
        color.SetColorEntry(20, (250,214,169))
        color.SetColorEntry(21, (250,214,171))
        color.SetColorEntry(22, (250,214,172))
        color.SetColorEntry(23, (250,214,173))
        color.SetColorEntry(24, (250,214,174))
        color.SetColorEntry(25, (250,214,175))
        color.SetColorEntry(26, (250,214,175))
        color.SetColorEntry(27, (250,214,174))
        color.SetColorEntry(28, (250,214,173))
        color.SetColorEntry(29, (250,214,172))
        color.SetColorEntry(30, (250,214,171))
        color.SetColorEntry(31, (250,214,187))
        color.SetColorEntry(32, (250,250,208))
        color.SetColorEntry(33, ( 250,250,209))
        color.SetColorEntry(34, (250,250,210))
        color.SetColorEntry(35, (250,250,211))
        color.SetColorEntry(36, (250,250,212))
        color.SetColorEntry(37, (250,250,213))
        color.SetColorEntry(38, (250,250,214))
        color.SetColorEntry(39, (250,250,215))
        color.SetColorEntry(40, (250,250,216))
        color.SetColorEntry(41, (250,250,217))
        color.SetColorEntry(42, (250,250,218))
        color.SetColorEntry(43, (250,250,219))
        color.SetColorEntry(44, (250,250,220))
        color.SetColorEntry(45, (250,250,221))
        color.SetColorEntry(46, (250,250,222))
        color.SetColorEntry(47, (250,250,223))
        color.SetColorEntry(48, (250,250,215))
        color.SetColorEntry(49, (241,241,144))
        color.SetColorEntry(50, (241,241,145))
        color.SetColorEntry(51, (241,241,146))
        color.SetColorEntry(52, (241,241,147))
        color.SetColorEntry(53, (241,241,148))
        color.SetColorEntry(54, (241,241,149))
        color.SetColorEntry(55, (241,241,150))
        color.SetColorEntry(56, (241,241,151))
        color.SetColorEntry(57, (241,241,152))
        color.SetColorEntry(58, (241,241,153))
        color.SetColorEntry(59, (241,241,154))
        color.SetColorEntry(60, (241,241,156))
        color.SetColorEntry(61, (241,241,155))
        color.SetColorEntry(62, (241,241,157))
        color.SetColorEntry(63, (241,241,158))
        color.SetColorEntry(64, (241,241,159))
        color.SetColorEntry(65, (255,255,0))
        color.SetColorEntry(66, (255,255,10))
        color.SetColorEntry(67, (255,255,11))
        color.SetColorEntry(68, (255,255,13))
        color.SetColorEntry(69, (255,255,15))
        color.SetColorEntry(70, (255,255,17))
        color.SetColorEntry(71, (90,217,250))
        color.SetColorEntry(72, (90,217,251))
        color.SetColorEntry(73, (90,217,252))
        color.SetColorEntry(74, (90,217,253))
        color.SetColorEntry(75, (90,217,254))
        color.SetColorEntry(76, (90,217,240))
        color.SetColorEntry(77, (90,217,241))
        color.SetColorEntry(78, (90,217,242))
        color.SetColorEntry(79, (90,217,243))
        color.SetColorEntry(80, (90,217,244))
        color.SetColorEntry(81, (90,217,245))
        color.SetColorEntry(82, (90,217,246))
        color.SetColorEntry(83, (90,217,247))
        color.SetColorEntry(84, (13,255,255))
        color.SetColorEntry(85, (13,255,253))
        color.SetColorEntry(86, (13,255,251))
        color.SetColorEntry(87, (13,255,249))
        color.SetColorEntry(88, (13,255,246))
        color.SetColorEntry(89, (13,255,244))
        color.SetColorEntry(90, (13,255,240))
        color.SetColorEntry(91, (13,254,238))
        color.SetColorEntry(92, (13,254,230))
        color.SetColorEntry(93, (13,254,228))
        color.SetColorEntry(94, (51,254,210))
        color.SetColorEntry(95, (51,253,210))
        color.SetColorEntry(96, (51,251,210))
        color.SetColorEntry(97, (51,249,210))
        color.SetColorEntry(98, (51,248,210))
        color.SetColorEntry(98, (51,247,210))
        color.SetColorEntry(99, (51,246,210))
        color.SetColorEntry(100, (51,245,210))
        color.SetColorEntry(101, (176, 173, 169))
        color.SetColorEntry(102, (176, 173, 169))
        color.SetColorEntry(103, (176, 173, 169))
        color.SetColorEntry(104, (176, 173, 169))
        color.SetColorEntry(105, (176, 173, 169))
        color.SetColorEntry(106, (176, 173, 169))
        color.SetColorEntry(107, (176, 173, 169))
        color.SetColorEntry(108, (176, 173, 169))
        color.SetColorEntry(109, (176, 173, 169))
        color.SetColorEntry(110, (176, 173, 169))
        color.SetColorEntry(111, (176, 173, 169))
        color.SetColorEntry(112, (176, 173, 169))
        color.SetColorEntry(113, (176, 173, 169))
        color.SetColorEntry(114, (176, 173, 169))
        color.SetColorEntry(115, (176, 173, 169))
        color.SetColorEntry(116, (176, 173, 169))
        color.SetColorEntry(117, (176, 173, 169))
        color.SetColorEntry(118, (176, 173, 169))
        color.SetColorEntry(119, (176, 173, 169))
        color.SetColorEntry(120, (176, 173, 169))
        color.SetColorEntry(121, (176, 173, 169))
        color.SetColorEntry(122, (176, 173, 169))
        color.SetColorEntry(123, (176, 173, 169))
        color.SetColorEntry(124, (176, 173, 169))
        color.SetColorEntry(125, (176, 173, 169))
        color.SetColorEntry(126, (176, 173, 169))
        color.SetColorEntry(127, (176, 173, 169))
        color.SetColorEntry(128, (176, 173, 169))
        color.SetColorEntry(129, (176, 173, 169))
        color.SetColorEntry(130, (176, 173, 169))
        color.SetColorEntry(131, (176, 173, 169))
        color.SetColorEntry(132, (176, 173, 169))
        color.SetColorEntry(133, (176, 173, 169))
        color.SetColorEntry(134, (176, 173, 169))
        color.SetColorEntry(135, (176, 173, 169))
        color.SetColorEntry(136, (176, 173, 169))
        color.SetColorEntry(137, (176, 173, 169))
        color.SetColorEntry(138, (176, 173, 169))
        color.SetColorEntry(139, (176, 173, 169))
        color.SetColorEntry(140, (176, 173, 169))
        color.SetColorEntry(141, (176, 173, 169))
        color.SetColorEntry(142, (176, 173, 169))
        color.SetColorEntry(143, (176, 173, 169))
        color.SetColorEntry(144, (176, 173, 169))
        color.SetColorEntry(145, (176, 173, 169))
        color.SetColorEntry(146, (176, 173, 169))
        color.SetColorEntry(147, (176, 173, 169))
        color.SetColorEntry(148, (176, 173, 169))
        color.SetColorEntry(149, (176, 173, 169))
        color.SetColorEntry(150, (176, 173, 169))
        color.SetColorEntry(151, (176, 173, 169))
        color.SetColorEntry(152, (176, 173, 169))
        color.SetColorEntry(153, (176, 173, 169))
        color.SetColorEntry(154, (176, 173, 169))
        color.SetColorEntry(155, (176, 173, 169))
        color.SetColorEntry(156, (176, 173, 169))
        color.SetColorEntry(157, (176, 173, 169))
        color.SetColorEntry(158, (176, 173, 169))
        color.SetColorEntry(159, (176, 173, 169))
        color.SetColorEntry(160, (176, 173, 169))
        color.SetColorEntry(161, (176, 173, 169))
        color.SetColorEntry(162, (176, 173, 169))
        color.SetColorEntry(163, (176, 173, 169))
        color.SetColorEntry(164, (176, 173, 169))
        color.SetColorEntry(165, (176, 173, 169))
        color.SetColorEntry(166, (176, 173, 169))
        color.SetColorEntry(167, (176, 173, 169))
        color.SetColorEntry(168, (176, 173, 169))
        color.SetColorEntry(169, (176, 173, 169))
        color.SetColorEntry(170, (176, 173, 169))
        color.SetColorEntry(171, (176, 173, 169))
        color.SetColorEntry(172, (176, 173, 169))
        color.SetColorEntry(173, (176, 173, 169))
        color.SetColorEntry(174, (176, 173, 169))
        color.SetColorEntry(175, (176, 173, 169))
        color.SetColorEntry(176, (176, 173, 169))
        color.SetColorEntry(177, (176, 173, 169))
        color.SetColorEntry(178, (176, 173, 169))
        color.SetColorEntry(179, (176, 173, 169))
        color.SetColorEntry(180, (176, 173, 169))
        color.SetColorEntry(181, (176, 173, 169))
        color.SetColorEntry(182, (176, 173, 169))
        color.SetColorEntry(183, (176, 173, 169))
        color.SetColorEntry(184, (176, 173, 169))
        color.SetColorEntry(185, (176, 173, 169))
        color.SetColorEntry(186, (176, 173, 169))
        color.SetColorEntry(187, (176, 173, 169))
        color.SetColorEntry(188, (176, 173, 169))
        color.SetColorEntry(189, (176, 173, 169))
        color.SetColorEntry(190, (176, 173, 169))
        color.SetColorEntry(191, (176, 173, 169))
        color.SetColorEntry(192, (176, 173, 169))
        color.SetColorEntry(193, (176, 173, 169))
        color.SetColorEntry(194, (176, 173, 169))
        color.SetColorEntry(195, (176, 173, 169))
        color.SetColorEntry(196, (176, 173, 169))
        color.SetColorEntry(197, (176, 173, 169))
        color.SetColorEntry(198, (176, 173, 169))
        color.SetColorEntry(199, (176, 173, 169))
        # missing data
        color.SetColorEntry(200, (255,255,255))
        # no decision
        color.SetColorEntry(201, (255,255,255))
        color.SetColorEntry(202, (237, 223, 232))
        color.SetColorEntry(203, (237, 223, 232))
        color.SetColorEntry(204, (237, 223, 232))
        color.SetColorEntry(205, (237, 223, 232))
        color.SetColorEntry(206, (237, 223, 232))
        color.SetColorEntry(207, (237, 223, 232))
        color.SetColorEntry(208, (237, 223, 232))
        color.SetColorEntry(209, (237, 223, 232))
        color.SetColorEntry(210, (237, 223, 232))
        # ningth
        color.SetColorEntry(211, (75, 75, 75))
        color.SetColorEntry(212, (237, 223, 232))
        color.SetColorEntry(213, (237, 223, 232))
        color.SetColorEntry(214, (237, 223, 232))
        color.SetColorEntry(215, (237, 223, 232))
        color.SetColorEntry(216, (237, 223, 232))
        color.SetColorEntry(217, (237, 223, 232))
        color.SetColorEntry(218, (237, 223, 232))
        color.SetColorEntry(219, (237, 223, 232))
        color.SetColorEntry(220, (237, 223, 232))
        color.SetColorEntry(221, (237, 223, 232))
        color.SetColorEntry(222, (237, 223, 232))
        color.SetColorEntry(223, (237, 223, 232))
        color.SetColorEntry(224, (237, 223, 232))
        color.SetColorEntry(225, (237, 223, 232))
        color.SetColorEntry(226, (237, 223, 232))
        color.SetColorEntry(227, (237, 223, 232))
        color.SetColorEntry(228, (237, 223, 232))
        color.SetColorEntry(229, (237, 223, 232))
        color.SetColorEntry(230, (237, 223, 232))
        color.SetColorEntry(231, (237, 223, 232))
        color.SetColorEntry(232, (237, 223, 232))
        color.SetColorEntry(233, (237, 223, 232))
        color.SetColorEntry(234, (237, 223, 232))
        color.SetColorEntry(235, (237, 223, 232))
        color.SetColorEntry(236, (237, 223, 232))
        # inland water
        color.SetColorEntry(237, (0,0,255))
        color.SetColorEntry(238, (121, 222, 19))
        # ocean
        color.SetColorEntry(239, (23, 34, 241))
        color.SetColorEntry(240, (237, 223, 232))
        color.SetColorEntry(241, (237, 223, 232))
        color.SetColorEntry(242, (237, 223, 232))
        color.SetColorEntry(243, (237, 223, 232))
        color.SetColorEntry(244, (237, 223, 232))
        color.SetColorEntry(245, (237, 223, 232))
        color.SetColorEntry(246, (237, 223, 232))
        color.SetColorEntry(247, (237, 223, 232))
        color.SetColorEntry(248, (237, 223, 232))
        color.SetColorEntry(249, (237, 223, 232))
        # cloud
        color.SetColorEntry(250, (153,51,255))
        color.SetColorEntry(251, (237, 223, 232))
        color.SetColorEntry(252, (237, 223, 232))
        color.SetColorEntry(253, (237, 223, 232))
        # detector saturated
        color.SetColorEntry(254, (255,0,255))
        # fill
        color.SetColorEntry(255, (255,0,255))



        band.SetRasterColorTable(color)
        band.SetRasterColorInterpretation(gdal.GCI_PaletteIndex)
        del band, ds

        print("it is done from function")
        return



def main():

    dir = r'U:/RS_Task_Workspaces/NDWI/2021/data'
    fold = [x[0] for x in os.walk(dir)]
    print("New data are in foler: ", fold[-1])
    nrt_refl = proces_HDF()
    ndwi = nrt_refl.proc_MOD09_HDF(fold[-1])
    #print("ADATE: ",ADATE)
    print("mosiaci is located in:",ndwi)

    print("doing composites...")
    nrt_refl = None
    shutil.copy2(ndwi, os.path.join(r'U:/RS_Task_Workspaces/NDWI/2021/ARCHIVE/ndwi', os.path.basename(ndwi)))
    os.chdir(dir)
    print("deleting the folowing folder:",fold[-1])
    shutil.rmtree(fold[-1])

if __name__ == "__main__":
    main()
    #send an email to interested parties........
    print("sending an email to various parties....")
    os.chdir(r'U:/RS_Task_Workspaces/NDWI/scripts')
    # com = "python final_goa_1.py"
    # os.system(com)


