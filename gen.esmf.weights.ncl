load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

;; create ESMF weight files to regrid from grid to grid

;; Specify arguments to script via command-line:
;; ncl gen.esmf.weights.ncl sourcefile=\"$source\" targetfile=\"$target\" weightfile=\"$weight\"

opt = True
    
opt@InterpMethod   = "patch"
opt@SrcRegional    = True
opt@DstRegional    = True
opt@ForceOverwrite = True
opt@PrintTimings   = True
  
ESMF_regrid_gen_weights(sourcefile, targetfile, weightfile, opt)


;; Copyright 2016 Univ. Corp. for Atmos. Research
;; Author: Seth McGinnis, mcginnis@ucar.edu
