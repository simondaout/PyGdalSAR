#!/usr/bin/env python
# -*- coding: utf-8 -*-

def temporal_decomp(line):
  'Linear temporal parametric decomposition'

  # Passage global des variables et matrice
  global models, inaps
  global basis, kernels, Mbasis
  # shared index
  global M, new_cols, new_lines, N
  if ((line % 10) == 0):
    logger.info('Processing line: {} --- {} seconds ---'.format(line,time.time() - start_time))

  # Initialisationi of forward temporal displacements
  mdisp=np.ones((N*new_cols), dtype=np.float32)*float('NaN')

  # Inisilize model parameters to zero
  m = np.ones((M*new_cols), dtype=np.float32)*float('NaN')
  sigmam = np.ones((M*new_cols), dtype=np.float32)*float('NaN')

  # Initialisation of G
  G = np.zeros((N*new_cols,M*new_cols), dtype=np.float32)
  Gf = np.zeros((N*new_cols,M*new_cols), dtype=np.float32)
  d = np.zeros((N*new_cols), dtype=np.float32)
  sigmad = np.zeros((N*new_cols), dtype=np.float32)

  # build G
  for col in range(new_cols):
    col_beg = col*M; col_end = col*M + M
    line_beg = col*N; line_end = col*N + N

    disp = maps_flata[line,col,:]
    sigmad[line_beg:line_end] = inaps
    # do not take into account NaN data
    k = np.flatnonzero(~np.isnan(disp)) # invers of isnan
    kk = len(k); tabx = dates[k]; taby = disp[k]

    if kk > N/10:
      Gp = np.zeros((N,M)) # without data with NaN
      Gpf = np.zeros((N,M)) # full for forward model only
      # Build G family of function k1(t),k2(t),...,kn(t): #
      #                                                   #
      #           |k1(0) .. kM(0)|                        #
      # Gfamily = |k1(1) .. kM(1)|                        #
      #           |..    ..  ..  |                        #
      #           |k1(N) .. kM(N)|                        #
      #                                                   #
      for l in range((Mbasis)):
        Gp[:kk,l]=basis[l].g(tabx)
        Gpf[:,l]=basis[l].g(dates)
      for l in range((Mker)):
        Gp[:kk,Mbasis+l]=kernels[l].g(k)
        Gpf[:,Mbasis+l]=kernels[l].g(range(N))

      # mise à jour des matrices creuses G, Gf, d, et sigmad
      G[line_beg:line_end,col_beg:col_end] = Gp
      d[line_beg:line_end] = taby
      sigmad[line_beg:line_end] = inaps[k]
      Gf[line_beg:line_end,col_beg:col_end] = Gpf

  ##### Inversion #####
  m,sigmam = consInvert(G,d,sigmad,cond=arguments["--cond"],ineq=arguments["--ineq"],eguality=eguality)
  ####################
  # Forward model
  try:
    mdisp = np.dot(Gf,m)
  except:
    pass

  # filled shared matrix
  models[line,:,:] = mdisp.reshape(new_cols,N)
  for l in range((Mbasis)):
        basis[l].m[line,:] = m[l::M]
        basis[l].sigmam[line,:] = sigmam[l::M]
  for l in range((Mker)):
        kernels[l].m[line,:] = m[l+Mbasis::M]
        kernels[l].sigmam[line,:] = sigmam[l+Mbasis::M]

  del m, sigmam


