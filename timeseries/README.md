To use this package
=============
* cube file: time series displacement maps in a binary format of size N\*Ncol\*Nlines, where N is the number of acquisitions dates, Ncol is the number of columns per date and Nlines, the number of lines per dates. To create a cube, use create\_cube.py
* list\_images file: text file of list of images made of 5 columns containing for each images 1) number 2) date in YYYYMMDD format 3) what you want (not use) 4) numerical date 5) perpendicular baseline.
To convert date in YYYYMMDD format to numerical format, use date2dec.py.
* lect.in: text file containing the number of columns and lines (eg. 1420 4000). (Depricated: it has to be replace by cube.hdr) 


invers\_disp2coef.py
============
Spatial and temporal inversions of the time series delay maps based on an iteration procedure. Requiered at least a cube of images (BIP format) and a list of images with temporal and perpendicular baselines. 
At each iteration, (1) estimation of spatial ramps, (2) linear decomposition in time based on a library of temporal functions (linear, heaviside, logarithm, seasonal...) or a discretised vector file.  
(3) estimation of RMS of each maps that will be then used as weight for the next iteration. Possibility to also to correct for a term proportional to the topography.

```
invers_disp2coef.py -h | --help
```

correct\_ts\_from\_gacos.py
============
Correct InSAR Time Series data from Gacos atmospheric models (data to be download and cited on: ceg-research.ncl.ac.uk/v2/gacos/). 1) Convert .ztd files to .tif format, 2) crop, re-project and re-resample atmospheric models to data geometry 3) correct time series data.

```
correct_ts_from_gacos.py -h | --help
```


invers\_disp\_pixel.py
============
Temporal decomposition of the time series delay maps of selected pixels. 

```
Usage: invers_disp_pixel.py --cols=<values> --ligns=<values> [--cube=<path>] [--windowsize=<value>]  
[--lectfile=<path>] [--aps=<path>][--rmspixel=<path>] [--interseismic=<value>] [--threshold_rmsd=<value>] 
[--coseismic=<value>] [--postseismic=<value>] [--seasonal=<yes/no>] [--vector=<path>]
[--semianual=<yes/no>]  [--dem=<yes/no>] [--imref=<value>] [--cond=<value>] [--slowslip=<value>] [--ineq=<value>] 
[--name=<value>] [--rad2mm=<value>] [--plot=<yes/no>] [<iref>] [<jref>] [--bounds=<value>] 

invers_disp_pixel.py -h | --help
Options:

-h --help               Show this screen
--ncols VALUE           Pixel column numbers (eg. 200,400,450) 
--nligns VALUE          Pixel lign numbers  (eg. 1200,1200,3000) 
--cube PATH             Path to displacement file [default: depl_cumul_flat]
--windowsize VALUE      Number of pixels around the pixel defining the window [default: 0]
--lectfile PATH         Path to the lect.in file (output of invers_pixel) [default: lect.in]
--aps PATH              Path to the APS file giving the error associated to each dates [default: No weigthing]
--rmspixel PATH         Path to the RMS map giving the error associated to each pixel (e.g RMSpixel, output of invers_pixel) [default: None]        
--interseismic PATH     Add a linear function to the inversion
--threshold_rmsd VALUE  If interseismic = yes: first try inversion with ref/interseismic/dem only, if RMDS inversion > threshold_rmsd then add other 
basis functions [default: 1.] 
--coseismic PATH        Add heaviside functions to the inversion .Indicate coseismic time (e.g 2004.,2006.)
--postseismic PATH      Add logarithmic transients to each coseismic step. Indicate characteristic time of the log function, must be a serie of values 
of the same lenght than coseismic (e.g 1.,1.). To not associate postseismic function to a given coseismic step, put None (e.g None,1.) 
--slowslip   VALUE      Add slow-slip function in the inversion (as defined by Larson et al., 2004). Indicate median and characteristic time of the events
(e.g. 2004.,1,2006,0.5), default: None 
--vector PATH           Path to the vector text file containing a value for each dates [default: None]
--seasonal PATH         If yes, add seasonal terms in the inversion
--semianual PATH        If yes, add semianual  terms in the inversion
--dem PATH              If yes, add term proportional to the perpendicular baseline in the inversion
--imref VALUE           Reference image number [default: 1]
--cond VALUE            Condition value for optimization: Singular value smaller than cond x largest_singular_value are considered zero [default: 1.0e-10]
--ineq VALUE            If yes, add inequqlity constrained in the inversion: use least square result to iterate the inversion. Force postseismic to be the 
same sign than coseismic [default: no].       
--name Value            Name output figures [default: None] 
--rad2mm                Scaling value between input data (rad) and desired output [default: -4.4563]
--plot                  Display results [default: yes]            
--iref                  colum numbers of the reference pixel [default: None] 
--jref                  lign number of the reference pixel [default: None]
--bounds                Min,Max time series plots 
```

lect\_disp\_pixel.py
=============

Plot time series results obtained with invers_disp2coef.py for given pixels.

```
Usage: lect_cube_pixel.py --cols=<values> --ligns=<values> [--cube=<path>] [--lectfile=<path>] 
[--ref=<path>] [--slope=<path>] [--coseismic=<paths>] [--postseismic=<paths>] [--slowslip=<value>] 
 [--cos=<path>] [--sin=<path>] [--dem=<path>] [--aps=<path>] 
 [--rad2mm=<value>] [--plot=<yes/no>] [--name=<value>] [--imref=<value>] 
 [<iref>] [<jref>] [--bounds=<value>] 


Options:
-h --help           Show this screen
--ncols VALUE       Pixels column numbers (eg. 200,400,450) 
--nligns VALUE      Pixels lign numbers pixel (eg. 1200,1200,3000) 
--cube PATH         Path to displacement file [default: depl_cumul_flat]
--lectfile PATH     Path of the lect.in file [default: lect.in]
--ref PATH          Path to the reference map in r4 or tif format (e.g ref_coeff.r4) [default: ref_coeff.r4]
--slope PATH        Path to velocity map in r4 or tif format (e.g lin_coeff.r4) [default: None]
--cos PATH          Path to cosinus map in r4 format (e.g coswt_coeff.r4) [default: None]
--sin PATH          Path to sinus map in r4 or tif format (e.g sinwt_coeff.r4)  [default: None]
--coseismic         Time of the events (e.g 2013.2,2014.5) [default: None]
--postseismic       Characteristic time of the transient functions (e.g 1.,0.5) [default: None]
--slowslip   VALUE  Read slow-slip maps. Indicate median and characteristic time of the events
--ref PATH          Path to reference file (ref_coeff.r4) [default: None]
--dem PATH          Path to dem file error map in r4 or tif format (demerr_coeff.r4) [default: None]
--aps PATH          Path to aps file: 1 column file given aps for each dates (eg. aps.txt) [default: None]
--rad2mm                Scaling value between input data (rad) and desired output [default: -4.4563]
--name Value            Name output figures [default: None] 
--plot                  Display results [default: yes]            
--iref              colum numbers of the reference pixel [default: None] 
--jref              lign number of the reference pixel [default: None]
--imref VALUE           Reference image number [default: 1]
--bounds                Min,Max time series plots 
```

invers\_disp\_gps.py
============

Temporal decomposition of gps time series and projection into a LOS vector. 

```
Usage: invers_disp_gps.py --network=<path> --reduction=<path> [--dim=<value>] [--wdir=<path>] [--extension=<value>] [--coseismic=<value>][--postseismic=<value>] [--seasonal=<yes/no>] [--cond=<value>] [--ineq=<value>] [--proj=<value>] [--scale=<value>]

Options:
-h --help           Show this screen 
--network PATH      Text file of name lat, lon for each stations
--reduction PATH    Directory name containing the TS (each TS in a text file untitled by the GPS name: East, North, Down, sigma_East, sigma_North, sigma Down format)      
--dim VALUE         dimension of the gps time series [default:2]
--wdir PATH         Path to the network file and reduction directory [default: ./]
--extension Value   Extention of the the GPS TS files [default: .dat]
--coseismic PATH    Add heaviside functions to the inversion, indicate coseismic time (e.g 2004.,2006.)
--postseismic PATH  Add logarithmic transients to each coseismic step, indicate characteristic time of the log function, must be a serie of values of the same lenght than coseismic (e.g 1.,1.). To not associate postseismic function to a give coseismic step, put None (e.g None,1.) 
--seasonal PATH     If yes, add seasonal terms in the inversion
--cond VALUE        Condition value for optimization: Singular value smaller than cond x largest_singular_value are considered zero [default: 1.0e-10]
--ineq VALUE        If yes, add inequqlity constrained in the inversion: use lself.east square result to iterate the inversion. Force postseismic to be the same sign than coseismic [default: no].  
--proj VALUE        LOS projection [default for Envisat satellite: 0.38690,-0.09645,0.91706]
--scale VALUE       Scale factor for the three components [default: 1,1,1]
``` 
