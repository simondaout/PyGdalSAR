#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
invers_temp.py — Temporal decomposition of NSBAS cumulative time series.
Simplified version of invers_disp2coef.py (Simon Daout), modified from FORTRAN code model_temp_4var_seas.f90 (Marie-Pierre Doin)

Performs ONLY the temporal inversion (linear + seasonal ± semi-annual ±
bi-annual ± steps).  All spatial iterations, masks, ramps, topo corrections,
and empirical estimations are removed. Median susbstraction in the refernce zone is iterated as in FORTRAN code model_temp_4var_seas.f90.

Usage
-----
    # From the flatsim_validate directory, passing the track directory:
    python invers_temp.py data/Tienshan/D107_NORD

    # Explicit TS directory:
    python invers_temp.py data/Tienshan/D107_NORD/TS [options]

Arguments
---------
    track_dir              Track dir (with TS/ subdir) or TS directory directly
    --cube=<path>          GeoTIFF or ENVI cube  [default: auto-detect]
    --list_images=<path>   list_images.txt        [default: auto-detect]
    --rms=<path>           inrms.txt  [default: auto-detect; "none" to disable]
    --aps=<path>           inaps.txt  [default: auto-detect; "none" to disable]
    --niter=<n>            Number of iterations   [default: 2]
    --linear=<yes/no>      Linear velocity term   [default: yes]
    --seasonal=<yes/no>    Annual cos+sin terms   [default: yes]
    --semianual=<yes/no>   Semi-annual terms      [default: no]
    --bianual=<yes/no>     Bi-annual terms        [default: no]
    --steps=<t1,t2,...>    Heaviside step times   [default: None]
    --cond=<value>         SVD condition number   [default: 1e-5]
    --imref=<n>            Reference image (1-based) [default: 1]
    --dateslim=<dmin,dmax> Date limits e.g. 20141013,20220528 [default: all]
    --nproc=<n>            Number of CPU cores    [default: 4]
    --plot=<yes/no>        Show plots             [default: no]

Outputs (written in the TS directory, all as GeoTIFF)
------------------------------------------------------
    lin_coeff.tif / lin_sigcoeff.tif     linear velocity + uncertainty
    ref_coeff.tif / ref_sigcoeff.tif     constant term + uncertainty
    cos_coeff.tif / sin_coeff.tif        seasonal cos/sin amplitudes
    ampwt_coeff.tif / phiwt_coeff.tif   seasonal amplitude and phase
    *_sigcoeff.tif                       uncertainties for each map
    sigma_N.txt                          per-image RMS at each iteration
    inversion.eps                        summary map of all coefficients

Dependencies
------------
    numpy, scipy, gdal (osgeo), matplotlib
"""

import os, sys, glob, math, time, logging, warnings, argparse
import multiprocessing, gc

import numpy as np
from numpy.lib.stride_tricks import as_strided
import scipy.optimize as opt
from osgeo import gdal
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from datetime import datetime as dt

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

logging.basicConfig(level=logging.INFO,
                    format='line %(lineno)s -- %(levelname)s -- %(message)s')
logger = logging.getLogger('invers_temp.log')
start_time = time.time()



# ─────────────────────────────────────────────────────────────────────────────
#  Basis function classes
# ─────────────────────────────────────────────────────────────────────────────

class pattern:
    def __init__(self, name, reduction, date):
        self.name = name; self.reduction = reduction; self.date = date
    def info(self): print(f'  {self.name}  (t0={self.date})')

def Heaviside(t):
    h = np.zeros(len(t)); h[t >= 0] = 1.; return h

class reference(pattern):
    def g(self, t): return np.ones(t.size)

class linear(pattern):
    def __init__(self, name, reduction, date):
        super().__init__(name, reduction, date); self.to = date
    def g(self, t): return t - self.to

class cosvar(pattern):
    def __init__(self, name, reduction, date):
        super().__init__(name, reduction, date); self.to = date
    def g(self, t): return np.array([math.cos(2*math.pi*(ti-self.to)) for ti in t])

class sinvar(pattern):
    def __init__(self, name, reduction, date):
        super().__init__(name, reduction, date); self.to = date
    def g(self, t): return np.array([math.sin(2*math.pi*(ti-self.to)) for ti in t])

class cos2var(pattern):
    def __init__(self, name, reduction, date):
        super().__init__(name, reduction, date); self.to = date
    def g(self, t): return np.array([math.cos(4*math.pi*(ti-self.to)) for ti in t])

class sin2var(pattern):
    def __init__(self, name, reduction, date):
        super().__init__(name, reduction, date); self.to = date
    def g(self, t): return np.array([math.sin(4*math.pi*(ti-self.to)) for ti in t])

class cos5var(pattern):
    def __init__(self, name, reduction, date):
        super().__init__(name, reduction, date); self.to = date
    def g(self, t): return np.array([math.cos(math.pi*(ti-self.to)) for ti in t])

class sin5var(pattern):
    def __init__(self, name, reduction, date):
        super().__init__(name, reduction, date); self.to = date
    def g(self, t): return np.array([math.sin(math.pi*(ti-self.to)) for ti in t])

class steps(pattern):
    def __init__(self, name, reduction, date):
        super().__init__(name, reduction, date); self.to = date
    def g(self, t): return Heaviside(t - self.to)


# ─────────────────────────────────────────────────────────────────────────────
#  Utilities
# ─────────────────────────────────────────────────────────────────────────────

def date2dec(dates):
    times = []
    for d in np.atleast_1d(dates):
        x = dt.strptime(str(d), '%Y%m%d')
        times.append(float(x.strftime('%Y')) + float(x.strftime('%j')) / 365.1)
    return times


def write_envi_hdr(filename, shape, dtype='float32', interleave='bip'):
    dtype_map = {'float32': 4, 'int16': 2, 'float64': 5}
    lines, samples = shape[0], shape[1]
    bands = shape[2] if len(shape) == 3 else 1
    with open(filename + '.hdr', 'w') as f:
        f.write(f"ENVI\nsamples = {samples}\nlines   = {lines}\n"
                f"bands   = {bands}\ndata type = {dtype_map.get(dtype,4)}\n"
                f"interleave = {interleave}\nbyte order = 0\n")


def save_tif(arr, name, ts_dir, driver, gt, proj):
    """Write a 2-D float32 array as GeoTIFF."""
    outname = os.path.join(ts_dir, f'{name}.tif')
    ncols, nlines = arr.shape[1], arr.shape[0]
    ds = driver.Create(outname, ncols, nlines, 1, gdal.GDT_Float32)
    ds.GetRasterBand(1).WriteArray(arr.astype(np.float32))
    ds.SetGeoTransform(gt)
    ds.SetProjection(proj)
    ds.GetRasterBand(1).FlushCache()
    del ds
    logger.info(f'Saved: {outname}')


def _wls(A, b, w, cond=1e-5):
    """
    Weighted least-squares:  min ||W(Ax - b)||²
    w : 1-D weight array (not variances — weights applied directly).
    Returns (solution, sigmam).
    """
    Aw = A * w[:, np.newaxis]
    bw = b * w
    fsoln = np.linalg.lstsq(Aw, bw, rcond=cond)[0]
    try:
        varx   = np.linalg.pinv(A.T @ A)
        res2   = np.sum((b - A @ fsoln) ** 2)
        scale  = 1. / max(A.shape[0] - A.shape[1], 1)
        sigmam = np.sqrt(scale * res2 * np.diag(varx))
    except np.linalg.LinAlgError:
        sigmam = np.full(A.shape[1], np.nan)
    return fsoln, sigmam


# ─────────────────────────────────────────────────────────────────────────────
#  Multiprocessing workers
# ─────────────────────────────────────────────────────────────────────────────

_G = {}   # global state shared across workers in the same process

def _init_workers(N_, M_, dates_, basis_, Mbasis_, cond_, cte_coh_, avg_):
    global _G
    _G = dict(N=N_, M=M_, dates=dates_, basis=basis_, Mbasis=Mbasis_,
              cond=cond_, cte_coh=cte_coh_, avg=avg_)


def _temporal_decomp(disp, uncertainty):
    """
    Invert one pixel time series with an IRLS inner loop (2 iterations).
    uncertainty : per-image in_sigma (large = bad date).
    Converted to weight = 1/uncertainty before use.
    Inner-loop pixel weight per date k:
        pix_w_k = 1 / (cte_coh + |residual_k| / rms_pixel)
    """
    N, M, dates   = _G['N'], _G['M'], _G['dates']
    basis, Mbasis = _G['basis'], _G['Mbasis']
    cond, cte_coh = _G['cond'], _G['cte_coh']
    avg           = _G['avg']

    k  = np.flatnonzero(~np.isnan(disp))
    kk = len(k)
    m      = np.full(M, np.nan, dtype=np.float32)
    sigmam = np.full(M, np.nan, dtype=np.float32)
    mdisp  = np.full(N, np.nan, dtype=np.float32)

    if kk <= N / 6:
        return m, sigmam, mdisp

    tabx   = dates[k]
    # subtract per-date reference mean (Fortran: tab_deplac(k) = deplac(k) - avg(k))
    taby   = (disp[k] - avg[k]).astype(np.float64)
    # convert per-image uncertainty → weight (large uncertainty = small weight)
    weight_k = 1.0 / uncertainty[k].astype(np.float64)

    G = np.zeros((kk, M), dtype=np.float64)
    for l in range(Mbasis):
        G[:, l] = basis[l].g(tabx)

    # Start with outer weights only (iteration 0)
    pix_w = np.ones(kk)

    for inner_iter in range(3):   # iter 0 = init, iter 1-2 = IRLS (as in Fortran)
        #  W = 1/[(σ_APS+ε) * max(σm,ε) * (|r|+ε)]
        w_total = weight_k * pix_w
        bb, sigmam_tmp = _wls(G, taby, w_total, cond=cond)

        if inner_iter < 2:
            # compute weighted residuals → update pixel weights
            residuals = taby - G @ bb
            w2      = w_total ** 2
            rms_pix = np.sqrt(np.sum(residuals**2 * w2) / np.sum(w2))
            if rms_pix > 0:
                pix_w = 1.0 / (cte_coh + np.abs(residuals) / rms_pix)

    m      = bb.astype(np.float32)
    sigmam = sigmam_tmp.astype(np.float32)
    mdisp[k] = (G @ bb).astype(np.float32)

    return m, sigmam, mdisp


def _chunk_worker(args):
    chunk, uncertainty = args    # chunk: (N, P), uncertainty: per-image (N,)
    N, P = chunk.shape
    M = _G['M']
    m_all      = np.empty((P, M), dtype=np.float32)
    sigmam_all = np.empty((P, M), dtype=np.float32)
    models_all = np.empty((P, N), dtype=np.float32)
    for i in range(P):
        m, sm, md = _temporal_decomp(chunk[:, i], uncertainty)
        m_all[i], sigmam_all[i], models_all[i] = m, sm, md
    return m_all, sigmam_all, models_all


# ─────────────────────────────────────────────────────────────────────────────
#  Path helpers
# ─────────────────────────────────────────────────────────────────────────────

def _resolve_ts_dir(track_or_ts):
    p = os.path.abspath(track_or_ts)
    if os.path.basename(p) == 'TS' or not os.path.isdir(os.path.join(p, 'TS')):
        return p
    return os.path.join(p, 'TS')


def _find(ts_dir, *patterns):
    """Return first existing match among glob patterns, else None."""
    for pat in patterns:
        m = sorted(glob.glob(os.path.join(ts_dir, pat)))
        if m:
            return m[0]
    return None


# ─────────────────────────────────────────────────────────────────────────────
#  MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Temporal-only inversion of FLATSIM TS (no spatial iterations)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("track_dir",
                        help="Track directory (with TS/ subdir) or TS directory")
    parser.add_argument("--cube",        default=None,  help="Cube path (auto-detected)")
    parser.add_argument("--list_images", default=None,  help="list_images.txt path")
    parser.add_argument("--rms", default=None, help="inrms.txt path; pass none to use unit weights")
    parser.add_argument("--aps", default=None, help="inaps.txt path; pass none to use unit weights")
    parser.add_argument("--niter",       type=int,   default=2,    help="Iterations [2]")
    parser.add_argument("--linear",      default='yes', help="Linear term [yes]")
    parser.add_argument("--seasonal",    default='yes', help="Annual seasonal [yes]")
    parser.add_argument("--semianual",   default='no',  help="Semi-annual [no]")
    parser.add_argument("--bianual",     default='no',  help="Bi-annual [no]")
    parser.add_argument("--steps",       default=None,  help="Step times e.g. 2010.5,2015.2")
    parser.add_argument("--cond",        type=float, default=1e-5, help="SVD cond [1e-5]")
    parser.add_argument("--cte_coh",     type=float, default=0.4,
                        help="IRLS damping constant — weight = 1/(cte_coh + |res|/rms). "
                             "0.5 is recommended (same as Fortran). [default: 0.4]")
    parser.add_argument("--imref",       type=int,   default=1,    help="Ref image 1-based [1]")
    parser.add_argument("--dateslim",    default=None,  help="dmin,dmax e.g. 20141013,20220528")
    parser.add_argument("--nproc",       type=int,   default=4,    help="CPU cores [4]")
    parser.add_argument("--block_size",  type=int,   default=0,
                        help="Lines per parallel block. 0 = auto from RAM. [default: 0 = auto]")
    parser.add_argument("--ref_zone",    default=None,
                        help="Reference zone for APS std: l0,l1,c0,c1 (0-based). "
                             "Default: all pixels.")
    parser.add_argument("--rmspixel",    default=None,
                        help="Path to RMSpixel map (r4 or tiff). Pixels with "
                             "RMSpixel > rmsl are excluded from APS computation.")
    parser.add_argument("--rmsl",        type=float, default=10.0,
                        help="RMSpixel threshold for APS mask. [default: 10.0]")
    parser.add_argument("--plot",        default='no',  help="Show plots [no]")
    args = parser.parse_args()

    ts_dir = _resolve_ts_dir(args.track_dir)
    os.chdir(ts_dir)
    logger.info(f'Working directory: {ts_dir}')

    # ── auto-detect input files ──────────────────────────────────────────────
    cube_path = args.cube or _find(ts_dir, 'CNES_DTs_geo_*.tiff',
                                           'CNES_DTs_radar_*.tiff',
                                           'depl_cumule')
    list_path = args.list_images or _find(ts_dir, 'list_images.txt',
                                                   'images_retenues')
    rms_path = None if (args.rms and args.rms.lower() == 'none') \
               else (args.rms or _find(ts_dir, 'inrms.txt'))
    aps_path = None if (args.aps and args.aps.lower() == 'none') \
               else (args.aps or _find(ts_dir, 'inaps.txt'))

    for p, lbl in [(cube_path, 'cube'), (list_path, 'list_images')]:
        if not p or not os.path.exists(p):
            logger.critical(f'Required file not found: {lbl} ({p}). Exit.')
            sys.exit(1)

    logger.info(f'Cube        : {cube_path}')
    logger.info(f'List images : {list_path}')
    logger.info(f'RMS weights : {rms_path or "DISABLED (unit weights)"}')
    logger.info(f'APS weights : {aps_path or "DISABLED (unit weights)"}')

    # ── GeoTIFF projection (from cube if tiff, else from first tiff found) ──
    geotiff_ref = cube_path if cube_path.endswith(('.tif', '.tiff')) \
                  else _find(ts_dir, 'CNES_*.tiff')
    if not geotiff_ref or not os.path.exists(geotiff_ref):
        logger.critical('No GeoTIFF reference found for projection. Exit.')
        sys.exit(1)
    georef = gdal.Open(geotiff_ref)
    gt     = georef.GetGeoTransform()
    proj   = georef.GetProjection()
    driver = gdal.GetDriverByName('GTiff')
    logger.info(f'GeoTIFF projection from: {geotiff_ref}')

    # ── read list_images (cols 1=YYYYMMDD, 3=decimal date, 5=bperp) ─────────
    data   = np.loadtxt(list_path, comments='#', usecols=[1, 3, 5], unpack=True)
    idates = data[0].astype(int)
    dates  = data[1]

    # ── date limits ──────────────────────────────────────────────────────────
    if args.dateslim:
        dmin_s, dmax_s = args.dateslim.replace(',', ' ').split()
        datemin = int(date2dec(dmin_s)[0])
        datemax = int(date2dec(dmax_s)[0]) + 1
    else:
        datemin, datemax = int(np.min(dates)), int(np.max(dates)) + 1

    indexd = np.flatnonzero((dates >= datemin) & (dates <= datemax))
    idates, dates = idates[indexd], dates[indexd]
    N     = len(dates)
    imref = max(0, args.imref - 1)

    # ── load cube ────────────────────────────────────────────────────────────
    ds = gdal.Open(cube_path)
    if not ds:
        logger.critical(f'Cannot open cube: {cube_path}')
        sys.exit(1)
    ncol, nlines, N_cube = ds.RasterXSize, ds.RasterYSize, ds.RasterCount
    maps_temp = np.zeros((nlines, ncol, N_cube), dtype=np.float32)
    for bi in range(1, N_cube + 1):
        maps_temp[:, :, bi-1] = ds.GetRasterBand(bi).ReadAsArray()
    maps_temp[np.isin(maps_temp, [9990, 9999])] = np.nan
    logger.info(f'Cube shape: {maps_temp.shape}')

    # reference image subtraction
    cst = np.copy(maps_temp[:, :, imref])
    cst[np.isnan(cst)] = 0.
    for l in range(N_cube):
        maps_temp[:, :, l] -= cst
        if l != imref:
            maps_temp[:, :, l][maps_temp[:, :, l] == 0.] = np.nan

    # select date range
    maps = maps_temp[:, :, indexd].copy()
    del maps_temp
    new_lines, new_cols = maps.shape[0], maps.shape[1]
    std_maps = float(np.nanstd(maps))
    logger.info(f'Cube after date selection: {maps.shape}  std={std_maps:.4f}')

    # write lect_ts.in
    with open('lect_ts.in', 'w') as f:
        np.savetxt(f, (new_cols, new_lines, N), fmt='%6i', newline='\t')

    # ── reference zone for APS computation ───────────────────────────────────
    if args.ref_zone:
        l0, l1, c0, c1 = map(int, args.ref_zone.replace(',', ' ').split())
    else:
        l0, l1, c0, c1 = 0, new_lines, 0, new_cols
    logger.info(f'APS reference zone: lines {l0}:{l1}  cols {c0}:{c1}')

    # ── rmspixel mask (exclude bad pixels from APS computation) ──────────────
    # Auto-detect RMSpixel in TS/ or AUX/ if not explicitly provided
    rms_mask = np.ones((new_lines, new_cols), dtype=bool)  # True = use pixel
    rmspixel_path = args.rmspixel or _find(ts_dir, 'RMSpixel', 'RMSpixel.tif')
    if rmspixel_path and os.path.exists(rmspixel_path):
        try:
            from osgeo import gdal as _gdal
            ds = _gdal.Open(rmspixel_path)
            rms_px = ds.GetRasterBand(1).ReadAsArray().astype(np.float32)
            del ds
        except Exception:
            rms_px = np.fromfile(rmspixel_path, dtype=np.float32)
            rms_px = rms_px.reshape(new_lines, new_cols) if rms_px.size == new_lines * new_cols else None
        if rms_px is not None:
            rms_mask = rms_px <= args.rmsl
            logger.info(f'RMSpixel ({rmspixel_path}): {rms_mask.sum()} / {new_lines*new_cols} pixels used')
    else:
        logger.info('RMSpixel not found — using all pixels for APS computation')

    # ── build basis functions ─────────────────────────────────────────────────
    cos_times = (list(map(float, args.steps.replace(',', ' ').split()))
                 if args.steps else [])

    basis = [reference(name='reference', reduction='ref', date=datemin)]
    index = 1
    indexseas = None

    if args.linear == 'yes':
        basis.append(linear(name='linear', reduction='lin', date=datemin))
        index += 1

    if args.seasonal == 'yes':
        indexseas = index
        basis.append(cosvar(name='seasonal (cos)', reduction='cos', date=datemin))
        basis.append(sinvar(name='seasonal (sin)', reduction='sin', date=datemin))
        index += 2

    if args.semianual == 'yes':
        basis.append(cos2var(name='semi-annual (cos)', reduction='cosw2t', date=datemin))
        basis.append(sin2var(name='semi-annual (sin)', reduction='sinw2t', date=datemin))
        index += 2

    if args.bianual == 'yes':
        basis.append(cos5var(name='bi-annual (cos)', reduction='cos5wt', date=datemin))
        basis.append(sin5var(name='bi-annual (sin)', reduction='sin5wt', date=datemin))
        index += 2

    for i, t in enumerate(cos_times):
        basis.append(steps(name=f'step {i}', reduction=f'step{i}', date=t))
        index += 1

    Mbasis = len(basis)
    M      = Mbasis

    print(f'\nBasis functions ({Mbasis}):')
    for b in basis:
        b.info()
    print()

    for b in basis:
        b.m      = np.full((new_lines, new_cols), np.nan, dtype=np.float32)
        b.sigmam = np.full((new_lines, new_cols), np.nan, dtype=np.float32)

    # ── input weights ─────────────────────────────────────────────────────────
    def _load_weights(path):
        """
        Read a per-image weight file (1-column or last column of multi-column).
        Returns raw values without epsilon treatment.
        """
        if path and os.path.exists(path):
            raw = np.loadtxt(path, comments='#', dtype='f')
            w = raw[:, -1] if raw.ndim == 2 else raw.flatten()
            w = w[indexd] if len(w) > N else w[:N]
            return w.astype(np.float32)
        return np.ones(N, dtype=np.float32)

    # APS: σ_APS + ε  (paper: σ_APS(tk) + ε)
    _aps_raw = _load_weights(aps_path)
    in_aps  = _aps_raw + args.cte_coh       # initial APS + ε, kept fixed
    logger.info(f'APS: raw min={_aps_raw.min():.3f} max={_aps_raw.max():.3f} '
                f'→ in_aps (+ ε): min={in_aps.min():.3f}')

    # RMS: max(σm, ε)  (paper: max(σm(tk), ε))
    _rms_raw = _load_weights(rms_path)
    in_rms   = np.clip(_rms_raw, args.cte_coh, None)
    logger.info(f'RMS: raw min={_rms_raw.min():.3f} max={_rms_raw.max():.3f} '
                f'→ in_rms (max ε): min={in_rms.min():.3f}')

    # in_sigma = (σ_APS+ε) * max(σm,ε)  — pixel term (|r|+ε) handled in inner loop
    in_sigma = in_aps * in_rms  # initial uncertainty

    # ── save cube to memmap for parallel reads ────────────────────────────────
    mm = np.memmap('.tmp_meanmap_depl_cumule', dtype='float32', mode='w+',
                   shape=(new_lines, new_cols, N))
    mm[:] = maps[:]
    mm.flush()
    write_envi_hdr('.tmp_meanmap_depl_cumule', shape=(new_lines, new_cols, N))
    del mm, maps

    mm_models = np.memmap('.tmp_meanmap_disp_cumul_models', dtype='float32', mode='w+',
                           shape=(new_lines, new_cols, N))
    mm_models.flush()
    write_envi_hdr('.tmp_meanmap_disp_cumul_models', shape=(new_lines, new_cols, N))
    del mm_models

    # ── iteration loop ────────────────────────────────────────────────────────
    nproc = min(args.nproc, new_lines)
    if args.block_size > 0:
        block_size = min(args.block_size, new_lines)
        logger.info(f'Block size: {block_size} (manual)')
    else:
        # try:
        #     import psutil
        #     avail = psutil.virtual_memory().available
        # except ImportError:
        #     avail = 8 * 1024**3
        #     logger.warning('psutil not found, assuming 8 GB. pip install psutil')
        # bytes_per_line = new_cols * (N*4 + N*4 + M*4*3)
        # block_size = min(max(1, int(0.30 * avail / bytes_per_line)), new_lines)
        # logger.info(f'Block size: {block_size} (auto, RAM={avail/1024**3:.1f}GB)')
        block_size = 100
        logger.info(f'Block size: {block_size}')

    # avg(k): per-date mean on reference zone (Fortran: avg(:)=0)
    avg = np.zeros(N, dtype=np.float64)

    for ii in range(args.niter):
        print(f'{"─"*45}')
        print(f'  Iteration {ii+1}/{args.niter}')
        print(f'{"─"*45}')
        logger.info(f'Input uncertainties: {in_sigma}')

        with multiprocessing.Pool(
            processes=nproc,
            initializer=_init_workers,
            initargs=(N, M, dates, basis, Mbasis, args.cond, args.cte_coh, avg),
        ) as pool:
            for line in range(0, new_lines, block_size):
                end_line = min(line + block_size, new_lines)
                bsz = end_line - line
                logger.info(f'  line {line:4d}/{new_lines}  '
                             f'({time.time()-start_time:.1f}s)')

                cube_r = np.memmap('.tmp_meanmap_depl_cumule', dtype='float32',
                                   mode='r', shape=(new_lines, new_cols, N))
                block  = cube_r[line:end_line, :, :].copy()  # (bsz, ncol, N)
                del cube_r
                # subtract per-date reference mean (Fortran: tab_deplac=deplac-avg)
                block -= avg[np.newaxis, np.newaxis, :]

                ts_block = block.transpose(2, 0, 1).reshape(N, -1)  # (N, bsz*ncol)
                chunks   = np.array_split(ts_block, nproc, axis=1)
                results  = pool.map(_chunk_worker,
                                    [(c, in_sigma) for c in chunks])  # in_sigma = uncertainty

                m_all, sm_all, md_all = (
                    np.concatenate([r[i] for r in results], axis=0)
                    for i in range(3)
                )
                md_all = md_all.reshape(bsz, new_cols, N)

                # add avg back to models (Fortran: phapred fitted on deplac-avg)
                md_all += avg[np.newaxis, np.newaxis, :]
                mm_w = np.memmap('.tmp_meanmap_disp_cumul_models', dtype='float32',
                                 mode='r+', shape=(new_lines, new_cols, N))
                mm_w[line:end_line, :, :] = md_all
                mm_w.flush()
                del mm_w

                for idx in range(bsz * new_cols):
                    i = line + idx // new_cols
                    j = idx  % new_cols
                    for l in range(Mbasis):
                        basis[l].m[i, j]      = m_all[idx, l]
                        basis[l].sigmam[i, j]  = sm_all[idx, l]

                del m_all, sm_all, md_all, results, block
                gc.collect()

        # ── residuals → update in_sigma ───────────────────────────────────────
        cube_r  = np.memmap('.tmp_meanmap_depl_cumule',  dtype='float32', mode='r',
                            shape=(new_lines, new_cols, N))
        mod_r   = np.memmap('.tmp_meanmap_disp_cumul_models', dtype='float32', mode='r',
                            shape=(new_lines, new_cols, N))
        mod_c   = np.copy(mod_r)
        mod_c[np.abs(mod_c) > 9999] = 0.
        # Fortran: somme_rescoh(k) = std of residuals on reference zone
        #          = sqrt( mean(r²) - mean(r)² )  on masked stable pixels
        r = (np.nan_to_num(cube_r, nan=np.nan)
             - np.nan_to_num(mod_c, nan=np.nan))          # (new_lines, new_cols, N)
        del cube_r, mod_r, mod_c

        # apply ref_zone and rmspixel mask
        r_ref = r[l0:l1, c0:c1, :]                        # (zone_lines, zone_cols, N)
        mask3d = rms_mask[l0:l1, c0:c1, np.newaxis]       # broadcast over N
        r_ref  = np.where(mask3d, r_ref, np.nan)

        # Use median for robustness (Fortran uses mean on reference zone)
        median_r = np.nanmedian(r_ref, axis=(0, 1))         # (N,) robust reference
        mean_r2  = np.nanmean(r_ref**2, axis=(0, 1))        # (N,) for std
        mean_r   = np.nanmean(r_ref, axis=(0, 1))           # (N,) for std formula
        res = np.sqrt(np.clip(mean_r2 - mean_r**2, 0, None))  # std = sqrt(E[r²]-E[r]²)
        # Update avg with median (more robust than mean)
        avg += median_r
        del r, r_ref

        print('\n  Dates         APS_std    Median_ref')
        for l in range(N):
            print(f'  {idates[l]}    {res[l]:.4f}    {median_r[l]:.4f}')
        # save APS std
        np.savetxt(f'aps_{ii}.txt', res.T, fmt='%.6f')
        # save median on reference zone (col 0 = YYYYMMDD, col 1 = median, col 2 = APS_std)
        out = np.column_stack([idates.astype(np.float64), median_r, res])
        np.savetxt(f'ref_median_{ii}.txt', out,
                   fmt='%10.0f  %.6f  %.6f',
                   header='YYYYMMDD  median_ref(rad)  APS_std(rad)')
        #  W = 1/[(σ_APS+ε) * max(σm,ε) * (|r|+ε)]
        in_sigma = (res + args.cte_coh) * in_rms

    # ── save coefficient maps (GeoTIFF) ───────────────────────────────────────
    print('\nSaving GeoTIFF outputs …')
    for b in basis:
        save_tif(b.m,      f'{b.reduction}_coeff',    ts_dir, driver, gt, proj)
        save_tif(b.sigmam, f'{b.reduction}_sigcoeff', ts_dir, driver, gt, proj)

    # seasonal amplitude & phase
    if args.seasonal == 'yes' and indexseas is not None:
        cos_m = basis[indexseas].m
        sin_m = basis[indexseas+1].m
        cos_s = basis[indexseas].sigmam
        sin_s = basis[indexseas+1].sigmam
        amp    = np.sqrt(cos_m**2 + sin_m**2)
        phi    = np.arctan2(sin_m, cos_m)
        sigamp = np.sqrt(cos_s**2 + sin_s**2)
        sigphi = (cos_s*np.abs(sin_m) + sin_s*np.abs(cos_m)) / (cos_s**2 + sin_s**2 + 1e-12)
        for name, arr in [('ampwt_coeff',    amp),    ('ampwt_sigcoeff', sigamp),
                           ('phiwt_coeff',    phi),    ('phiwt_sigcoeff', sigphi)]:
            save_tif(arr, name, ts_dir, driver, gt, proj)

    # ── summary plot ─────────────────────────────────────────────────────────
    cmap = cm.jet; cmap.set_bad('white')
    fig, axes = plt.subplots(1, Mbasis, figsize=(4*Mbasis, 4))
    if Mbasis == 1:
        axes = [axes]
    for ax, b in zip(axes, basis):
        finite = b.m[np.isfinite(b.m)]
        vmax = np.percentile(np.abs(finite), 98) if len(finite) else 1.
        im = ax.imshow(b.m, cmap=cmap, vmax=vmax, vmin=-vmax, aspect='auto')
        ax.set_title(b.reduction, fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])
        fig.colorbar(im, ax=ax, shrink=0.3, orientation='vertical')
    plt.suptitle('Temporal decomposition — coefficients', fontsize=11)
    fig.tight_layout()
    eps_path = os.path.join(ts_dir, 'inversion.eps')
    fig.savefig(eps_path, format='EPS', dpi=150)
    logger.info(f'Saved: {eps_path}')
    if args.plot == 'yes':
        plt.show()
    plt.close('all')

    # # ── cleanup ───────────────────────────────────────────────────────────────
    # # ── write disp_cumul_flat = cube - avg ───────────────────────────────────
    # # Fortran equivalent: depl_cumule_ref(k) = deplac(k) - avg(k)
    # logger.info('Writing disp_cumul_flat …')
    # cube_f = np.memmap('.tmp_meanmap_depl_cumule',     dtype='float32', mode='r',
    #                    shape=(new_lines, new_cols, N))
    # flat_f = np.memmap('disp_cumul_flat', dtype='float32', mode='w+',
    #                    shape=(new_lines, new_cols, N))
    # cube_arr = np.array(cube_f)
    # flat_f[:] = np.where(np.isnan(cube_arr),
    #                      np.nan,
    #                      cube_arr - avg[np.newaxis, np.newaxis, :])
    # flat_f.flush()
    # write_envi_hdr('disp_cumul_flat', shape=(new_lines, new_cols, N))
    # del cube_f, flat_f, cube_arr

    # cleanup temporary working files (distinct names so we never touch
    # a real depl_cumule/disp_cumul_models the user may have as input)
    for tmp in ['.tmp_meanmap_depl_cumule', '.tmp_meanmap_depl_cumule.hdr',
                '.tmp_meanmap_disp_cumul_models', '.tmp_meanmap_disp_cumul_models.hdr']:
        if os.path.exists(tmp):
            os.remove(tmp)

    print(f'\nDone in {time.time()-start_time:.1f}s  —  outputs in {ts_dir}')


if __name__ == '__main__':
    print()
    print('# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #')
    print('#                                                             #')
    print('#      Temporal Inversion of FLATSIM InSAR time series       #')
    print('#      (simplified — no spatial iterations)                  #')
    print('#                                                             #')
    print('# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #')
    print()
    main()
