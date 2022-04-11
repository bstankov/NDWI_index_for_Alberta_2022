# NDWI_index_for_Alberta_2022

## Daily NDWI index

This projects creates daily MODIS NDWI index raster files and weekly NDWI anomalies.
The NDWI is a remote sensing based indicator sensitive to the change in the water content of leaves (Gao, 1996).
NDWI is computed using the near infrared (NIR – MODIS band 2) and the short wave infrared (SWIR – MODIS band 5) reflectance.
From **https://nrt4.modaps.eosdis.nasa.gov/api/v2/content/archives/allData/61/MOD09A1N**

Using **"ident_download_MOD09A1N.py"**, we download following tiles of 8-day Modis composites: 
h10v03.061.NRT.hdf

h10v04.061.NRT.hdf

h11v03.061.NRT.hdf

h12v03.061.NRT.hdf

These hdf files are processed by **"MOD_proc_functs.py"**. The output of this proces is a current NDWI file which indicates vegetation crown water content.  
The NDWI product is dimensionless and varies between -1 to +1, depending on the leaf water content but also on the vegetation type and cover (Figure 1). 
High values of NDWI (in green) correspond to high vegetation water content and to high vegetation fraction cover. Low NDWI values (in red) correspond to low
vegetation water content and low vegetation fraction cover
![ndwi_20220407_10TM](https://user-images.githubusercontent.com/8118080/162656296-f129768c-c192-4e39-94dc-e9ed8b03f458.png)
Figure 1. Daily NDWI

### NDWI anomaly

To  produce NDWI index anomaly(Figure 2), we use historical series of Modis bands from January 2010 till end of 2021.
For each pixel we have Mean and Standard deviation for the period.

![image](https://user-images.githubusercontent.com/8118080/162657287-f9f34d3e-9189-429c-ad50-c5e3e0b5599f.png)

Figure 2. Calculation for NDWI anomaly

![image](https://user-images.githubusercontent.com/8118080/162658498-fa7ca0e6-b709-440e-aa9d-f68c2d4638af.png)

Figure 3. NDWI anomaly based on 10 year historical period.

The NDWI anomaly product is given in standard deviation units. It is commonly spread between -4 to +4 and from red to green.
Red shades portray negative anomalies.

REFERENCES:
Gao, B.-C. 1996. NDWI - A normalized difference water index for remote sensing of vegetation liquid water from space. Remote Sensing of Environment 58: 257-266.
