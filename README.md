# TSdecomp
Tool package for iterative spatial and temporal decompositions of InSAR or GPS time series written in Python programming language.

 Requirements
=============
This project needs the following external components:
 * Python (2.7)
 * NumPy
 * SciPy
 * getopt
 * matplotlib
 * datetime


invers\_disp2coef.py
============
Spatial and temporal inversions of the time series delay maps (used depl_cumule (BIP format) and images_retenues, output of invers_pixel), based on an iteration procedure. 
At each iteration, (1) estimation of spatial ramps, (2) linear decomposition in time based on a library of temporal functions (linear, heaviside, logarithm, seasonal), 
(3) estimation of RMS that will be then used as weight for the next iteration. Possibility to also to correct for a term proportional to the topography.

```
Usage: invers_disp2coef.py [--cube=<path>] [--lectfile=<path>] [--list_images=<path>] [--aps=<path>] [--interseismic=<yes/no>] [--threshold_rmsd=<value>] \
[--coseismic=<values>] [--postseismic=<values>]  [--seasonal=<yes/no>] [--slowslip=<values>] [--semianual=<yes/no>]  [--dem=<yes/no>] [--vector=<path>] \
[--flat=<0/1/2/3/4/5/6/7/8/9>] [--nfit=<0/1>] [--ivar=<0/1>] [--niter=<value>]  [--spatialiter=<yes/no>]  [--sampling=<value>] [--imref=<value>] [--mask=<path>] \
[--rampmask=<yes/no>] [--threshold_mask=<value>] [--scale_mask=<value>] [--topofile=<path>] [--aspect=<path>] [--perc=<value>] \
[--tempmask=<yes/no>] [--cond=<value>] [--ineq=<value>] [--rmspixel=<path>] [--threshold_rms=<path>] \
[--crop=<values>] [--fulloutput=<yes/no>] [--geotiff=<path>] [--plot=<yes/no>] \
[<ibeg>] [<iend>] [<jbeg>] [<jend>] 

invers_disp2coef.py -h | --help

Options:
-h --help               Show this screen
--cube PATH             Path to displacement file [default: depl_cumul]
--lectfile PATH         Path to the lect.in file (output of invers_pixel) [default: lect.in]
--list_images PATH      Path to list images file made of 4 columns containing for each images 1) number 2) date in YYYYMMDD format 3) numerical date 4) perpendicular baseline [default: images_retenues] 
--aps PATH              Path to the APS file giving an input error to each dates [default: No weigthing if no spatial estimation or misfit spatial estimation used as input uncertianties]
--rmspixel PATH         Path to the RMS map that gives an error for each pixel (e.g RMSpixel, output of invers_pixel) [default: None]
--threshold_rms VALUE   Threshold on rmsmap for spatial estimations [default: 1.]
--interseismic YES/NO   Add a linear function in the inversion 
--threshold_rmsd VALUE  If interseismic = yes: first try inversion without coseismic and postseismic, 
: if RMDS inversion > threshold_rmsd then add other basis functions [default: 1.] 
--coseismic PATH        Add heaviside functions to the inversion, indicate coseismic time (e.g 2004.,2006.)
--postseismic PATH      Add logarithmic transients to each coseismic step, indicate characteristic time of the log function, must be a serie of values of the same lenght than coseismic (e.g 1.,1.). To not associate postseismic function to a give coseismic step, put None (e.g None,1.) 
--slowslip   VALUE      Add slow-slip function in the inversion (as defined by Larson et al., 2004). Indicate median and characteristic time of the events (e.g. 2004.,1,2006,0.5), default: None 
--vector PATH           Path to the vector text file containing a value for each dates [default: None]
--seasonal YES/NO       If yes, add seasonal terms in the inversion
--semianual YES/NO      If yes, add semianual terms in the inversion
--dem Yes/No            If yes, add term proportional to the perpendicular baseline in the inversion
--ivar VALUE            Define phase/elevation relationship: ivar=0 function of elevation, ivar=1 crossed function of azimuth and elevation
--nfit VALUE            Fit degree in azimuth or in elevation [0:linear (default), 1: quadratic]
--flat PATH             Remove a spatial ramp at each iteration. 
0: ref frame [default], 1: range ramp ax+b , 2: azimutal ramp ay+b, 3: ax+by+c, 
4: ax+by+cxy+d 5: ax**2+bx+cy+d, 6: ay**2+by+cx+d, 7: ay**2+by+cx**2+dx+e, 
8: ay**2+by+cx**3+dx**2+ex+f, 9: ax+by+cxy**2+dxy+e
--niter VALUE           Number of iterations. At the first iteration, image uncertainties is given by aps file or misfit spatial iteration, 
while for the next itarations, uncertainties are equals to the global RMS of the previous iteration for each map [default: 1]
--spatialiter  YES/NO   If yes iterate the spatial estimations at each iterations (defined by niter) on the maps minus the temporal terms (ie. interseismic, coseismic...) [default: no]
--sampling VALUE        Downsampling factor [default: 1]
--imref VALUE           Reference image number [default: 1]
--mask PATH             Path to mask file in r4 or tif format for the spatial estimations (Keep only values > threshold_mask for ramp estimation).
--rampmask YES/NO       Remove a quadratic ramp in range a linear in azimuth on the mask before computing threshold [default: no]
--threshold_mask VALUE  Threshold on mask: take only > values (use scale factor for convenience) [default: 0]
--scale_mask  VALUE     Scale factor to apply on mask
--tempmask YES/NO       If yes, also use the spatial mask for the temporal inversion [default: no]
--topofile PATH         Path to topographic file in r4 or tif format. If not None, add a phase-elevation relationship in the saptial estimation.
--aspect PATH           Path to aspect file in r4 or tif format: take into account the slope orientation in the phase/topo relationship [default: None].
--perc VALUE            Percentile of hidden LOS pixel for the spatial estimations to clean outliers [default:90.]
--crop VALUE            Define a region of interest for the temporal decomposition [default: 0,nlign,0,ncol]
--cond VALUE            Condition value for optimization: Singular value smaller than cond*largest_singular_value are considered zero [default: 1.0e-10]
--ineq VALUE            If yes, add ineguality constraints in the inversion: use least square result without post-seismic functions
as a first guess to iterate the inversion. Force postseismic to be the same sign and inferior than coseismic steps of the first guess [default: no].
--fulloutput YES/NO     If yes produce maps of models, residuals, ramps, as well as flatten cube without seasonal and linear term [default: no]  
--geotiff PATH          Path to Geotiff to save outputs in tif format. If None save output are saved as .r4 files [default: .r4] 
--plot YES/NO           Display plots [default: yes] 
--ibeg VALUE            Line numbers bounding the ramp estimation zone [default: 0]
--iend VALUE            Line numbers bounding the ramp estimation zone [default: nlign]
--jbeg VALUE            Column numbers bounding the ramp estimation zone [default: 0]
--jend VALUE            Column numbers bounding the ramp estimation zone [default: ncol]
```


invers\_disp\_pixel.py
============
Temporal inversions of the time series delays of selected pixels (used depl_cumule (BIP format) and images_retenues, output of invers_pixel). 

```
Usage: invers_disp_pixel.py --cols=<values> --ligns=<values> [--cube=<path>] [--windowsize=<value>]  [--lectfile=<path>] [--aps=<path>] \
[--rmspixel=<path>] [--interseismic=<value>] [--threshold_rmsd=<value>] [--coseismic=<value>] [--postseismic=<value>] [--seasonal=<yes/no>] [--vector=<path>]\
[--semianual=<yes/no>]  [--dem=<yes/no>] [--imref=<value>] [--cond=<value>] [--slowslip=<value>] [--ineq=<value>] \
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

Plot time series results obtained with invers_disp2coef.py for given pixels (reqiered image_retenues file)

```
Usage: lect_cube_pixel.py --cols=<values> --ligns=<values> [--cube=<path>] [--lectfile=<path>] \
[--ref=<path>] [--slope=<path>] [--coseismic=<paths>] [--postseismic=<paths>] [--slowslip=<value>] \
 [--cos=<path>] [--sin=<path>] [--dem=<path>] [--aps=<path>] \
 [--rad2mm=<value>] [--plot=<yes/no>] [--name=<value>] [--imref=<value>] \
 [<iref>] [<jref>] [--bounds=<value>] 


Options:
-h --help           Show this screen
--ncols VALUE       Pixels column numbers (eg. 200,400,450) 
--nligns VALUE      Pixels lign numbers pixel (eg. 1200,1200,3000) 
--cube PATH         Path to displacement file [default: depl_cumul_flat]
--lectfile PATH     Path of the lect.in file [default: lect.in]
--ref PATH          Path to the reference map in r4 format (e.g ref_coeff.r4) [default: ref_coeff.r4]
--slope PATH        Path to velocity map in r4 format (e.g lin_coeff.r4) [default: None]
--cos PATH          Path to cosinus map in r4 format (e.g coswt_coeff.r4) [default: None]
--sin PATH          Path to sinus map in r4 format (e.g sinwt_coeff.r4)  [default: None]
--coseismic         Time of the events (e.g 2013.2,2014.5) [default: None]
--postseismic       Characteristic time of the transient functions (e.g 1.,0.5) [default: None]
--slowslip   VALUE  Read slow-slip maps. Indicate median and characteristic time of the events
--ref PATH          Path to reference file (ref_coeff.r4) [default: None]
--dem PATH          Path to dem file error map in r4 format (demerr_coeff.r4) [default: None]
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

Temporal inversions of the gps time series displacements

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
