PRO SMOOTH3D, im, hd, $
              imout, hdout, $
              fbeam, psf_org=psf_org,$
              svel=svel, $
              mask=mask, scale=scale,$
              ifail=ifail,$
              keep0=keep0
;+
; NAME:
;   SMOOTH3D
;
; PURPOSE:
;   * convol FITS cube/image to the desired resolution
;   * smooth in velocity if requested
;
; INPUTS:
;   IM          data (either 2d or 3d)
;   HD          fits header
;   fbeam       desired beam(psf) size in arcsec (could be a scalar or 3-element vector)
;               e.g.  fbeam=5             5"X5" (0D)    2D Gaussian
;                     fbeam=[10,5,90]     10"x5"(+90d)  2D Gaussian
;   [psf_org]   original psf. if bmaj/bmin doesn't exist in the header, one could specify
;               a pesudo 2D gaussian psf (useful for 2d optical images)
;               if psf_org was not specified or the specified psf_org<0 then we will try
;               to get bmaj/bmin from the fits header.
;   [svel]      in km/s: for 3d case only, smooth spectra with a Gaussian function
;               if svel=0.0, no spectra smoothing is performed      
;   [mask]      *nonzero* pixels in mask array will be replaced by 0 before convolution.
;               note: missing (NaN) data in IM are treated in the same way.
;   [scale]     set the scaling factor. otherwise, the scaling factor is automatically
;               determined by pixel value units.
;   [ifail]     if smoothing failed, ifail!=0 
;   [keep0]     if set, any value=0 pixel in the orginal image will be ZERO in the
;               smoothed image. This may not be a good choice for smoothing a masked 
;               mom0 from spectral cube, because it may screen out some low-brightness 
;               emission which could be detected in a mom0 image from a smoothed cube.
;   
; OUTPUTS:
;   IMOUT       output data
;   HDOUT       output header
;   if the convolution fails, imout will the original data w/mask applied if applicable
;
; EXAMPLE:
;   1.  smooth3d,im,hd,imout,hdout,[5.,4.,30.],20.
;   2.  im=readfits('n0337hi.line.cm.fits',hd)
;       smooth3d,im,hd,imsmo,hdsmo,[20.,10.,40.],svel=30
;
; NOTES:
;   convolve.pro from astrolib works as a wrapper of convol() or an alternative
;   of convol_fft(); convol3d.pro (from F.Varosi) doesn't provide padding and may
;   require too much memory for large cubes, althought it could be faster on some
;   computer configurations.
;   
; HISTORY:
;
;   20110306  RX  introduced
;   20110328  RX  handle Nan & 3d cube
;   20110401  RX  merge <convol2beam.pro> into <smo3d.pro>
;   20130227  RX  change the name from smo3d.pro to smooth3d.pro
;                 add an external masking option
;                 use convol_fft rather than convolve.pro
;                 handle images with rotated pixels
;   20130410  RX  add option /err to calculate an error map after smoothing
;                 svel=0 will not trigger 3D smoothing
;   20130625  RX  improve the compatibility of the image in JY/PIXEL (e.g. deconvolution model) 
;   20150528  TW  streamlined for use in MOMMAPS package                 
;   20150529  RX  use convol() instead of convol_ff() for better performance in small-kernel cases.
;                 preformance improvement: x0.2-0.3(memory)/x0.5-1.5(speed)
;-

currentExcept = !Except
!Except = 0
void = Check_Math()

; DETERMINE PIXEL SIZE
; CALL GETROT FROM ASTROLIB TO GET IMAGE ROTATION
GETROT,hd,rotang,cdelt
psize=cdelt[0]
psize=abs(psize*60.*60.)                ; pixel size in arcsec
naxis1=abs(sxpar(hd,'naxis1'))
naxis2=abs(sxpar(hd,'naxis2'))

; DETERMINE KERNEL
fpsf=fbeam
if  n_elements(fpsf) eq 1 then fpsf=[fpsf,fpsf,0.0]
ipsf=[0.,0.,0.]
if  n_elements(psf_org) eq 0 then begin
    RADIOHEAD,hd,s=h
    ipsf[0]=h.bmaj  ; in arcsec
    ipsf[1]=h.bmin  ; in arcsec
    ipsf[2]=h.bpa   ; in degrees (astro convention)
endif
if  n_elements(psf_org) eq 1 then begin
    ipsf[0]=psf_org  
    ipsf[1]=psf_org 
    ipsf[2]=0.0
    if  psf_org[0] lt 0 then begin
        RADIOHEAD,hd,s=h
        ipsf[0]=h.bmaj  ; in arcsec
        ipsf[1]=h.bmin  ; in arcsec
        ipsf[2]=h.bpa   ; in degrees (astro convention)
    endif
endif
if  n_elements(psf_org) eq 3 then begin
    ipsf=psf_org
endif
GKERNEL,fpsf[0],fpsf[1],fpsf[2],ipsf[0],ipsf[1],ipsf[2],bmaj,bmin,bpa,ifail
im_bpa=bpa+rotang

; REPLACE MISSING DATA & MASKED PIXELS WITH ZERO 
impad=im
if  n_elements(mask) eq 0 then begin
    tagnan=where(finite(im,/NAN))
endif else begin
    tagnan=where(finite(im,/NAN) or mask ne 0.0)
endelse
if tagnan[0] ne -1 then impad[tagnan]=0.0
imout=impad
hdout=hd

if  ifail eq 0 then begin

    bmaj=bmaj/psize
    bmin=bmin/psize
    psfsize=ceil(bmaj)*6+1
    psfsize=min([psfsize,floor(naxis1/2)*2-1,floor(naxis2/2)*2-1])
  
    ; GENERATE SMOOTHING KERNEL
    psf=PSF_GAUSSIAN(npixel=[psfsize,psfsize],fwhm=[bmin,bmaj],/NORMALIZE,DOUBLE=0)
    psf=rot(psf,-im_bpa,1.0,(psfsize-1)/2,(psfsize-1)/2,/INTERP,missing=0.0,cubic=-0.5)
    psf=psf>0d
    psf=psf/total(psf)
    
    ; SPATIAL SMOOTHING
    ww=sqrt(n_elements(impad))
    kk=sqrt(n_elements(psf))/(!CPU.TPOOL_NTHREADS*0.8)
    ww=sqrt(8.*alog(ww)/alog(2.0))
    sz=size(impad)
    message,/info," data size: "+strjoin(strtrim(sz,2)," ")
    message,/info, "kernel  size  [pixel]: "+string(psfsize)
    if  kk gt ww then begin
        ; prefer convol_fft() when the kernel is large
        message,/info,' use convol_fft()'
        if  sz[0] eq 2 then begin
            imout=convol_fft(impad,psf)
        endif else begin
            for i=0,sz[3]-1 do begin
                imout[0,0,i]=convol_fft(impad[*,*,i],psf)
            endfor
        endelse
    endif else begin
        ; prefer convol() when the kernel is small (no padding overheads like fft)
        ; PSF is rotated by 180 degrees to produce a "true" convolution in convol()
        ; note: the options /edge*/normal/nan will slow convol() signifcantly.
        message,/info,' use convol()' 
        psf=rotate(psf,2)
        if  sz[0] eq 3 then psf=reform(psf,psfsize,psfsize,1)
        if  sz[0] eq 4 then psf=reform(psf,psfsize,psfsize,1,1)
        imout=convol(impad,psf,nan=0,normal=0,edge_wrap=0,TPOOL_NOTHREAD=1)
    endelse
    
    ; VELOCITY SMOOTHING
    if  n_elements(svel) eq 1 then begin
        if  svel ne 0.0 then begin
            csize=abs(sxpar(hd,'cdelt3')/1000.0)    ; channel size in km/s
            naxis3=abs(sxpar(hd,'naxis3'))
            fvel=svel/csize
            spesize=ceil(fvel)*6+1
            spesize=min([spesize,floor(naxis3/2)*2-1])
            message,/info, "channel width [km/s]   : "+string(csize)
            message,/info, "kernel  fwhm  [channel]: "+string(fvel)
            message,/info, "kernel  size  [channel]: "+string(spesize)
            lsf=psf_gaussian(npixel=spesize,fwhm=fvel,/NORMALIZE,DOUBLE=0,ndimen=1)
            for j=0,naxis2-1 do begin
                for i=0,naxis1-1 do begin
                    if total(imout[i,j,*],/nan) ne 0 then $  
                       imout[i,j,*]=convol(reform(imout[i,j,*]),lsf,edge_mirror=1)
                endfor
            endfor
        endif
    endif

    ; SCALE PIXEL VALUE IF UNITS IN JY/BEAM
    if  n_elements(scale) ne 1 then begin
        scale=1.0
        if  STRPOS(STRUPCASE(sxpar(hd, 'BUNIT')), 'JY/B') ne -1  then $
            scale=fpsf[0]*fpsf[1]/(ipsf[0]*ipsf[1])
        if  STRPOS(STRUPCASE(sxpar(hd, 'BUNIT')), 'JY/P') ne -1  then begin
            scale=abs((fpsf[0]*fpsf[1])/(psize^2.0)*2.*!dpi/(8.*alog(2.)))
            SXADDPAR, hdout, 'BUNIT','JY/BEAM'
        endif
    endif
    message,/info, "scale factor: "+strjoin(scale)
    imout=float(imout*scale)

    SXADDPAR, hdout, 'BMAJ', fpsf[0]/60./60.
    SXADDPAR, hdout, 'BMIN', fpsf[1]/60./60.
    SXADDPAR, hdout, 'BPA',  fpsf[2]

endif else begin
  
    message,/info,"**"
    message,/info,"the target beam is smaller than the original beam!"
    message,/info,"return original image (with optional externel masking) and header"
    message,/info,"**"
    scale=1.0

endelse

; RESTORE MISSING DATA & MASKED PIXELS
if  tagnan[0] ne -1 then imout[tagnan]=!VALUES.F_NAN
if  keyword_set(keep0) then imout[where(im eq 0.0,/null)]=0.0
SXADDPAR, hdout, 'DATAMAX', max(imout,/nan)
SXADDPAR, hdout, 'DATAMIN', min(imout,/nan)

floating_point_underflow=32
status = Check_Math()         ; Get status and reset accumulated math error register.
if  (status and not floating_point_underflow) ne 0 then $
    message, 'IDL Check_Math() error: ' + strtrim(status, 2)
!Except = currentExcept
  
END



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;   This is an older version of smooth3d.pro (to be removed)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


PRO SMOOTH3D_OLD, im, hd, $
    imout, hdout, $
    fbeam, psf_org=psf_org,$
    svel=svel, $
    padchan=padchan, $
    mask=mask, scale=scale,$
    ifail=ifail,$
    keep0=keep0


;+
; NAME:
;   SMOOTH3D
;
; PURPOSE:
;   * convol FITS cube/image to the desired resolution
;   * smooth in velocity if requested
;
; INPUTS:
;   IM          data (either 2d or 3d)
;   HD          fits header
;   fbeam       desired beam(psf) size in arcsec (could be a scalar or 3-element vector)
;               e.g.  fbeam=5             5"X5" (0D)    2D Gaussian
;                     fbeam=[10,5,90]     10"x5"(+90d)  2D Gaussian
;   [psf_org]   original psf. if bmaj/bmin doesn't exist in the header, one could specify
;               a pesudo 2D gaussian psf (useful for 2d optical images)
;               if psf_org was not specified or the specified psf_org<0 then we will try
;               to get bmaj/bmin from the fits header.
;   [svel]      in km/s: for 3d case only, smooth spectra with a Gaussian function
;               if svel=0.0, no spectra smoothing is performed
;   [padchan]   pad some empty channels at spectra edges for better convolution in 3d
;               default: 0
;               fft will generally assumes your signal is periodic beyond its edges,
;               which is not a bad choice (compared with zero padding) if one would
;               like to use edge channels to estimate the noise
;   [mask]      *nonzero* pixels in mask array will be replaced by 0 before convolution.
;               note: missing (NaN) data in IM are treated in the same way.
;   [scale]     set the scaling factor. otherwise, the scaling factor is automatically
;               determined by pixel value units.
;   [ifail]     if smoothing failed, ifail!=0
;   [keep0]     if set, any value=0 pixel in the orginal image will be ZERO in the
;               smoothed image. This may not be a good choice for smoothing a masked
;               mom0 from spectral cube, because it may screen out some low-brightness
;               emission which could be detected in a mom0 image from a smoothed cube.
;
; OUTPUTS:
;   IMOUT       output data
;   HDOUT       output header
;   if the convolution fails, imout will the original data w/mask applied if applicable
;
; EXAMPLE:
;   1.  smooth3d,im,hd,imout,hdout,[5.,4.,30.],20.
;   2.  co=readfits('/Users/Rui/Dropbox/test/n4254co.line.cm.fits',cohd)
;       smooth3d,co,cohd,cosmo,cohdsmo,[20.,10.,40.],svel=30
;
; HISTORY:
;
;   20110306  RX  original version
;   20110328  RX  handle Nan & 3d cube
;   20110401  RX  merge <convol2beam.pro> into <smo3d.pro>
;   20130227  RX  change the name from smo3d.pro to smooth3d.pro
;                 add an external masking option
;                 use convol_fft rather than convolve.pro
;                 handle images with rotated pixels
;   20130410  RX  add option /err to calculate an error map after smoothing
;                 svel=0 will not trigger 3D smoothing
;   20130625  RX  improve the compatibility of the image in JY/PIXEL (e.g. deconvolution model)
;   20150528  TW  Streamlined for use in MOMMAPS package
;
;-

currentExcept = !Except
!Except = 0
void = Check_Math()

; DETERMINE PIXEL SIZE
; CALL GETROT FROM ASTROLIB TO GET IMAGE ROTATION
GETROT,hd,rotang,cdelt
psize=cdelt[0]
psize=abs(psize*60.*60.)                ; pixel size in arcsec
naxis1=abs(sxpar(hd,'naxis1'))
naxis2=abs(sxpar(hd,'naxis2'))

; DETERMINE KERNEL
fpsf=fbeam
if  n_elements(fpsf) eq 1 then fpsf=[fpsf,fpsf,0.0]
ipsf=[0.,0.,0.]
if  n_elements(psf_org) eq 0 then begin
    RADIOHEAD, hd, s = h
    ipsf[0]=h.bmaj  ; in arcsec
    ipsf[1]=h.bmin  ; in arcsec
    ipsf[2]=h.bpa   ; in degrees (astro convention)
endif
if  n_elements(psf_org) eq 1 then begin
    ipsf[0]=psf_org
    ipsf[1]=psf_org
    ipsf[2]=0.0
endif
if  n_elements(psf_org) eq 3 then begin
    ipsf=psf_org
endif
GKERNEL,fpsf[0],fpsf[1],fpsf[2],ipsf[0],ipsf[1],ipsf[2],bmaj,bmin,bpa,ifail

im_bpa=bpa+rotang

; REPLACE MISSING DATA & MASKED PIXELS WITH ZERO
impad=im
if  n_elements(mask) eq 0 then begin
    immask=im & immask[*]=0.0
endif else begin
    immask=mask
endelse
tagnan=where(finite(im,/NAN) or immask ne 0.0)
if tagnan[0] ne -1 then impad[tagnan]=0.0
imout=impad
hdout=hd

if  ifail eq 0 then begin

    bmaj=bmaj/psize
    bmin=bmin/psize
    psfsize=ceil(bmaj)*6+1
    psfsize=min([psfsize,floor(naxis1/2)*2-1,floor(naxis2/2)*2-1])

    ; GENERATE SMOOTHING KERNEL
    psf=PSF_GAUSSIAN(npixel=[psfsize,psfsize],fwhm=[bmin,bmaj],/NORMALIZE,/DOUBLE)
    psf=rot(psf,-im_bpa,1.0,(psfsize-1)/2,(psfsize-1)/2,/INTERP,missing=0.0)
    psf=psf>0d
    psf=psf/total(psf)
    if n_elements(svel) eq 1 then begin
        if svel ne 0.0 then begin
            csize=abs(sxpar(hd,'cdelt3')/1000.0)    ; channel size in km/s
            naxis3=abs(sxpar(hd,'naxis3'))
            fvel=svel/csize
            spesize=ceil(fvel)*6+1
            spesize=min([spesize,floor(naxis3/2)*2-1])
            message,/info, "channel width [km/s]   : "+string(csize)
            message,/info, "kernel  fwhm  [channel]: "+string(fvel)
            psf=PSF_GAUSSIAN(npixel=[psfsize,psfsize,spesize],fwhm=[bmin,bmaj,fvel],/NORMALIZE,/DOUBLE,ndimen=3)
            for i=0,spesize-1 do begin
                psf[*,*,i]=rot(psf[*,*,i],-im_bpa,1.0,(psfsize-1)/2,(psfsize-1)/2,/INTERP,missing=0)
            endfor
            psf=psf>0.0
            psf=psf/total(psf)
        endif
    endif
    message,/info, "kernel size:"+strjoin(size(psf,/d),' ')

    ; SMOOTHING
    if (size(psf))[0] eq 3 then begin
        if  n_elements(padchan) eq 0 then padchan=0
        tmp=dblarr(naxis1,naxis2,naxis3+2*padchan)
        tmp[*,*,padchan:(naxis3+padchan-1)]=impad
        tmp=CONVOL3D(tmp,psf)
        imout=tmp[*,*,padchan:(naxis3+padchan-1)]
    endif
    if (size(psf))[0] eq 2 then begin
        if  (size(impad))[0] eq 3 then begin
            for i=0,(size(impad))[3]-1 do begin
                imout[*,*,i]=convol_fft(double(impad[*,*,i]), psf)
            endfor
        endif else begin
            imout=convol_fft(impad,psf)
        endelse
    endif

    ; SCALE PIXEL VALUE IF UNITS IN JY/BEAM
    if n_elements(scale) ne 1 then begin
        scale=1.0
        if  STRPOS(STRUPCASE(sxpar(hd, 'BUNIT')), 'JY/B') ne -1  then $
            scale=fpsf[0]*fpsf[1]/(ipsf[0]*ipsf[1])
        if  STRPOS(STRUPCASE(sxpar(hd, 'BUNIT')), 'JY/P') ne -1  then begin
            scale=abs((fpsf[0]*fpsf[1])/(psize^2.0)*2.*!dpi/(8.*alog(2.)))
            SXADDPAR, hdout, 'BUNIT','JY/BEAM'
        endif
    endif
    message,/info, "scale factor: "+strjoin(scale)
    imout=float(imout*scale)

    SXADDPAR, hdout, 'BMAJ', fpsf[0]/60./60.
    SXADDPAR, hdout, 'BMIN', fpsf[1]/60./60.
    SXADDPAR, hdout, 'BPA',  fpsf[2]

endif else begin

    message,/info,"**"
    message,/info,"the target beam is smaller than the original beam!"
    message,/info,"return original image (with optional externel masking) and header"
    message,/info,"**"
    scale=1.0

endelse

; RESTORE MISSING DATA & MASKED PIXELS
if tagnan[0] ne -1 then imout[tagnan]=!VALUES.F_NAN
if keyword_set(keep0) then begin
    imout[where(im eq 0.0,/null)]=0.0
endif
SXADDPAR, hdout, 'DATAMAX', max(imout,/nan)
SXADDPAR, hdout, 'DATAMIN', min(imout,/nan)

floating_point_underflow = 32
status = Check_Math()         ; Get status and reset accumulated math error register.
IF(status AND NOT floating_point_underflow) NE 0 THEN $
    Message, 'IDL Check_Math() error: ' + StrTrim(status, 2)

!Except = currentExcept

END
