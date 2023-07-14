Scripts for geocoded files

convertGNSS\_to\_LOS.py / plot\_los\_gps\_ts.py
============

Convert GNSS to LOS and plot GPS time series into LOS 

Example:
```
convertGNSS_to_LOS.py proj_to_T117.py
plot_los_gps.py --network=T117_gps_network.txt --reduction=time_series --extension=.Ty.neu --outputdir=ts_plots_t117
```

Reference:
* [Daout, S., D'Agostino, N., Pathier, E., Socquet, A., Lavé, J., Doin, M. P., ... & Benedetti, L. Along-Strike Variation of the Strain Partitioning within the Apennines as Seen from Large-Scale Multi-Temporal Insar Analysis. Available at SSRN 4429391.] (https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4429391)


tie\_from\_decomposition.py / tie\_to\_gps.py
============

Tie several LOS frame together either using the horizontal/vertical decompostion or external GNSS data

Reference:
* [Ou, Q., Daout, S., Weiss, J. R., Shen, L., Lazecky, M., Wright, T. J., & Parsons, B. E. (2022). Large-scale interseismic strain mapping of the NE Tibetan plateau from Sentinel-1 interferometry. Journal of Geophysical Research: Solid Earth, 127(6).](https://eprints.whiterose.ac.uk/187293/)


 References
============

* [Daout, Simon, et al. "Large‐scale InSAR monitoring of permafrost freeze‐thaw cycles on the Tibetan Plateau." Geophysical Research Letters 44.2 (2017): 901-909](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1002/2016GL070781)

* [Daout, S., Sudhaus, H., Kausch, T., Steinberg, A., & Dini, B. (2019). Interseismic and postseismic shallow creep of the North Qaidam Thrust faults detected with a multitemporal InSAR analysis. Journal of Geophysical Research: Solid Earth, 124(7), 7259-7279.](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019JB017692)

geo2r4.py
============
Convert geotiff or any raster made of 1 band into real4 format.

