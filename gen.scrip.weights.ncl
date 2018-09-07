load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

;; create weight files to go from grid to grid

grids = (/"CRCM","ECP2","HRM3","MM5I","RCM3","WRFG","UDEL"/)
n = dimsizes(grids)

landmask = (/"grid","mask"/)

do L = 0,1
  land = landmask(L)
  do s=0,n-1
    source = grids(s)
    do t=0,n-1
      target = grids(t)
  
      if (source .ne. target) then

        print(land+" "+source+" "+target)
  
        sourcefile = "grids/scrip."+land+"."+source+".nc"
        targetfile = "grids/scrip."+land+"."+target+".nc"
        weightfile = "wgt."+source+"."+target+"."+land+".nc"
        
        opt = True
    
        opt@InterpMethod   = "patch"
        opt@SrcRegional    = True
        opt@DstRegional    = True
        opt@ForceOverwrite = True
        opt@PrintTimings   = True
  
        ESMF_regrid_gen_weights(sourcefile, targetfile, weightfile, opt)

      end if

    end do
  end do
end do    


;; Copyright 2013 Univ. Corp. for Atmos. Research
;; Author: Seth McGinnis, mcginnis@ucar.edu
