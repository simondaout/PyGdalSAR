Scripts specific to the NSBAS processing chain and FLATSIM products (https://flatsim.cnes.fr/doi/?lang=fr)  
=============

nsb\_unw\_ERAcorr\_stats.py / nsb\_unw\_GACOScorr\_stats.py
============

Tools to correct interferograms from ERA5 or GACOS atmospheric models

Reference:
* [Dodds, N., Daout, S., Walker, R. T., Begenjev, G., Bezmenov, Y., Mirzin, R., & Parsons, B. (2022). Interseismic deformation and strain-partitioning along the Main Köpetdag Fault, Turkmenistan, with Sentinel-1 InSAR time-series. Geophysical Journal International, 230(3), 1612-1629.](https://academic.oup.com/gji/article/230/3/1612/6568902)


nsb\_filtflatunw.py
============

Reference:
* [Daout, S., Sudhaus, H., Kausch, T., Steinberg, A., & Dini, B. (2019). Interseismic and postseismic shallow creep of the North Qaidam Thrust faults detected with a multitemporal InSAR analysis. Journal of Geophysical Research: Solid Earth, 124(7), 7259-7279.](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019JB017692)

Contrary to other processing chains, NSBAS provides various empirical corrections on the wrapped phase in order to reduce the variability of phase and thus enable unwrapping over large areas (Doin et al. (2015); Daout et al. (2017, 2018)). These corrections include (1) ramp corrections in range or azimuth to remove residual orbital errors, ionospheric delays, and effects of the satellite clock drift (Fig. processing.png), (2) empirical estimation of the phase delays due to differences in the stratified troposphere with a possible mask on deforming areas, (3) empirical removal of a fringe pattern (e.g co-seismic fringe pattern, seasonal ground deformation making the unwrapping difficult) (Fig. processing.png). All those corrections are estimated by searching on small moving sub-windows the proportionality between (1) phase and range, (2) phase and elevation, or (3) phase and the template that maximises the local phase coherence within the subwindows (Doin et al., 2011; Doin et al., 2015), respectively.

* The available jobs are check\_look, ecmwf, replace\_amp, filterSW, filterROI, flatr, flat\_az, flat\_atmo, flat\_model, colin, look\_int, unwrapping, add\_flatr\_back, add\_flata\_back, add\_model\_back, add\_atmo\_back. Bellow, some details about each of this processing step.

* The first thing you may want to do is multilook all your wrapped interferograms in order to reduce the size of your raster files. For this set Rlooks\_int in your procfile to a factor of 2, and run:

```
nsb\_filtflatunw.py --prefix= --suffix=\_sd --job=check\_look nsbas.proc
```

where, job corresponds to the list of post-processing jobs to be applied, and where prefix and suffix correspond to the –prefix and suffix of the interferogram names. Rlooks\_int can be 2/4/8/16/32. This averaging step reduces the size of the interferograms, smoothes the phase and reduce the noise. However, it might be dangerous to reduce the resolution of your interferograms if the deformation signal you are looking at is at high-frequency or if there is a high-frequency signal in your data that you do not want to smooth to avoid unwrapping errors. On the contrary, for long-wavelength deformation purposes like tectonic signal you may want to work at 4 or 8 looks on your wrapped interferograms.

* To replace the amplitude of the interferograms by the phase coherence, run:

```
nsb\_filtflatunw.py --prefix= --suffix=\_sd --job=replace\_amp nsbas.proc
```

This step will replace the amplitude of the complex interferogram by the coherence previously computed on 3x3 sliding windows by the program nsb\_gen\_phase\_par.pl. It allows to give more confidence to pixel with more phase coherence for the next steps instead of using the backscatter of the signal. Output interferograms will have the –prefix coh.

* To filter the phase of the interferograms, two options are available:

```
nsb\_filtflatunw.py --prefix=coh\_ --suffix=\_sd --job=filterROI nsbas.proc
nsb\_filtflatunw.py --prefix=coh\_ --suffix=\_sd --job=filterSW nsbas.proc
```

filterROI corresponds to the ROI\_PAC Goldstein adaptive filter. It requires the proc file parameter: filterStrength, which is the exponent applied to the power spectrum of the data and ranges between 0 and 4, and Filter\_method, which can be psfilt or adapt\_filt. Alternatively, filterSW is a coherence dependent filter that is a weighted average of the complex phase in sliding windows Doin et al. (2011), with the windows size defined by the parameterSWwindowsize in the proc file. A threshold on the amplitude (which is the coherence if ran replace\_amp job previously) is also defined by the parameter SWamplim. Several filter styles are available (SW: sliding window, SWc : cascade or SWcc: SWfilter\_casc\_grad, SWg: SWfilter\_grad, SWfilter\_gauss: gaussian windows). The cascade filter allows adaptative weighting of the filtered phase obtained with various window sizes using the coherence com- puted in each window Doin et al. (2011). This averaging filter, although crude, allows the retrieval of useful phase information from incoherent areas (mountain ranges, dunes...), but reduces phase gradients or jumps on the rapidly deforming areas, and therefore need to be applied carefully. As for the multilooking the choice of the filtering style and strength is a compromise between increasing the smoothing that will improve the unwrapping capability and removing the high-frequency content that will increase the unwrapping errors.

* To flatten the wrapped phase in range with an empirical function (See Doin et al. (2015)), use:

```
nsb\_filtflatunw.py --prefix=coh\_ --suffix=\_sd --job=flatr nsbas.proc
```

The estimation will be automatically done and the filterSW interferograms. The job filterSW will therefore be au- tomatically applied if the filterSW interferograms is not found. The required proc file parameters are: nfit\_range: median(-1),mean(0),linear=1,quadratic=2; thresh\_amp\_range: threshold on amplitude (coherence if replace\_amp has been applied, 0.2 might be a good choice). The output interferograms will have the additional suffix flatr.

*  Identically, to flatten the wrapped phase in azimuth with an empirical function (See Doin et al. (2015)), use:

```
nsb\_filtflatunw.py --prefix=coh\_ --suffix=\_sd\_flatr --job=flata nsbas.proc
```

The required proc file parameters are: nfit\_az: median(-1),mean(0),linear=1,quadratic=2; thresh\_amp\_az: thresh- old on amplitude (coherence if replace\_amp has been applied, 0.2 might be a good choice). The output interfero- grams will have the additional suffix flata.

*  To now flatr and flata remove an empirical phase-elevation relationship in order to reduce the stratified delays prior to unwrapping (See Doin et al. (2015)), run:

```
nsb\_filtflatunw.py --prefix=coh\_ --suffix=\_sd\_flatr\_flata --job=flat\_atmo nsbas.proc
```

The required proc file parameters are: nfit\_topo: order of the fit with the derivative phase with the topogra- phy: median(-1),mean(0),linear=1,quadratic=2; ivar: 1 function of elevation, 2:crossed function of azimuth and elevation, z\_ref: elevation reference for phase-elevation estimation (in meter), thresh\_amp\_topo: threshold on amplitude for phase-elevation estimation (if coherence, 0.2 might be a good choice). The estimation is also done on filterSW interferogram. Plot phase/topo relationship in date1-date2\_phase-topo.png file. The output interfero- grams will have the additional suffix flatz. Additionally, you may want to mask a given area for this step, if for example a lot of phase dispersion is observed on a dune area where the coherence if very low adding noise to the estimation. To do so, run:

```
nsb\_filtflatunw.py --prefix=coh\_ --suffix=\_sd\_flatr\_flata --job=flat\_atmo --ibeg\_mask=0 -- ibeg\_mask=1500 --jbeg\_mask=1000 --jend\_mask=3000 nsbas.proc
```

* An other flattening alternative to reduce the fringe rates and help for unwrapping, corresponds to remove a spatial model to each interferograms (See Daout et al. (2017)):

```
nsb\_filtflatunw.py --prefix=coh\_ --suffix=\_sd\_flatr\_flata\_flatz --job=flat\_model nsbas.proc
```

The model to be removed from the wrapped interfograms can be any patterns that helps you unwrapping (eg. coseismic fringes, permafrost, hydrological loading, subsidence...) and that you may have extracted from a PCA, ICA or a stack of few interferograms. Before applying the correction, it is important to crop this pattern around the deformation, remove long-wavelength ramps, and smooth strong gradients at the borders of the crop. To do that, you may use this command line:


*   After this flattening steps, you may want to replace the amplitude of the interferograms by the colinearity (See Pinel- Puyssegur et al. (2012)), running:

```
nsb\_filtflatunw.py --prefix=coh\_ --suffix=\_sd\_flatr\_flata\_flatz\_nomodel --job=colin nsbas.proc
```

this allows to re-estimate the variance of the interferogram of the reduced phase from ramp, phase/elevation or model empirical corrections. Also, contrary the to the coherence, the collinearity does not take into account the backscatter amplitude, which has been proven useful in natural environment, where backscatter variations are completely uncorrelated with the variance of the phase (i.e dark, smooth surfaces adjacent to bright, rough surfaces with identical stability). This step will change the –prefix of the interferogram to col.

* To multilook the interferograms to Rlooks\_unw, use the command line:

```
nsb\_filtflatunw.py --prefix=col\_ --suffix=\_sd\_flatz --job=look nsbas.proc
```

It required parameters Rlooks\_unw, which is the desired look number for unwrapping and can be 2,4,8,16,32.

*  Finally, the unwrapping is performed with the command line:

```
nsb\_filtflatunw.py --prefix=col\_ --suffix=\_sd\_flatz --job=unwrapping nsbas.proc
```

Several unwrapping methods are implemented. The unw\_method= mpd unwrapping method is performed using a specific scheme(López-Quiroz et al., 2009; Grandin et al., 2012). First, interferograms are multilook (Rlooks\_unw) and filtered with a (1) weighted average of the interferometric phase computed on sliding-windows (SW filter) and (2) an adaptative filter. Unwrapping is performed in adjoining sub-regions above a coherence threshold (threshold\_unw). Each newly unwrapped area is added to already unwrapped areas. The coherence threshold progressively decreases to propagate unwrapping further away. If necessary, high priority bridges defined in the text file bridges.in are used for the unwrapping. Afterwards, the phase difference between the smoothed unwrapped phase (SW) and the adaptative filtered interferograms is computed. The residual high-frequency signal is assumed to be between -π and π and added to the smooth unwrapped phase. Required parameters are: threshold\_unw, threshold\_unfilt, seedx, seedy. If unw\_method: roi, the program runs the ROIPAC algorthim. Required pa- rameters are: threshold\_unw, seedy, seedx. If unw\_method: icu: ICU algorithm.

- To launch automatically all the post-processing steps, run:

```
nsb\_filtflatunw.py --prefix= --suffix=\_sd --job=replace\_amp/flatr/flat\_atmo/flat\_model/colin/ look\_int/unwrapping --model=stack/stack\_clean.r4 nsbas.proc
```

- After unwrapping, you may want to visualise all unwrapped interferograms, running:

```
preview\_unw.py --int\_list=interf\_pair.rsc --radar={master\_date}/radar.hgt --outputdir=unw\_jpeg -- prefix=filt\_col --suffix=\_sd\_flatz --rlook=8
```

where, the radar argument points to the radar simulation of the Master image, outputdir is the output direc- tory where .jpeg are saved, int\_list points to the interferograms list file, and the interferograms are defined by prefixdate1-date2suffix\_rlookrlks.unw. The program also computed the list of successfull interferograms (interf\_pair\_sucess.txt) and the list of interferograms for which there are problems (interf\_pair\_problems.txt)

* Most likely, you will need to re-run unwrapping of several interferograms changing the filter strength/windows, the unwrapping threshold and imposing some unwrapping paths. To create bridges of a single interferogram, open the interferograms with MDX and right clic on the two points that need to be connected. Then copy the terminal output in a text file in the interferogram directory (eg. bridges.asc) and run:

```
bridge.py bridge.txt
```

The program will create the bridges.in file, which will be used for the unwrapping.

* To re-run the unwrapping on a sub-list of interferograms (e.g interf\_pair\_briges.txt), run:

```
nsb\_filtflatunw.py --prefix=col\_ --suffix=\_sd\_flatz --job=unwrapping --list\_int=interf\_pair\_bridges. txt nsbas.proc -f
```

The -f option forces the run of the jobs even if the output files are already available in the directories.

*  After unwrapping, nsb\_filtflatunw.py offers the possibility to re-introduce the previous corrections (orbital, atmo- spheric, model) previously removed to each interferogram to reconstruct the full unwrapped phase signal. Indeed, the the purpose of the corrections applied on the wrapped phase was not to separate the source of the signals but to flatten the interferogram to help for unwrapping. In addition, those corrections estimated on wrapped phase are possibly less accurate than those done on unwrapped phase and may contain some tectonic signals, that is why they can be reintroduced (Fig. 23d), running

```
nsb\_filtflatunw.py --prefix=col\_ --suffix=\_sd\_flatz\_flata\_flatz\_nomodel --job=add\_flatr\_back/ add\_flata\_back/add\_atmo\_back/add\_model\_back nsbas.proc
```

add\_flatr\_back: Add range ramp estimated on wrapped IFG back on unwrapped IFG. Required .flatr param- eter file containing polynomial fit produced by flatr. !!! Assumed flattening estimated on Rlooks\_int file. add\_flata\_back: Add azimutal ramp estimated on wrapped IFG back on unwrapped IFG. Required .flata pa- rameter file containing polynomial fit produced by flat\_az. !!! Assumed flattening estimated on Rlooks\_int file. add\_model\_back: Function adding model on unwrapped IFG previously removed on wrapped IFG. Required model file and parameter file .stack produced by flat\_model. add\_atmo\_back: Add back stratified model computed by flatten\_topo.


make\_fit\_az.sh / correct\_rgaz\_unw.py / correct\_ramp\_unw.py
============
Old NSBAS Tools: will be removed 

