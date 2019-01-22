FUNCTION PADDING,im,edge

; if edge>0 pad data with some zero-value pixels in each dimension 
; if edge<0 remove the padding pixels

if  edge gt 0 then begin
    struct=size(im, /STRUCTURE)
    imdims=size(im)
    tmp=MAKE_ARRAY(struct.dimensions[0:imdims[0]-1]+2*edge, TYPE=struct.TYPE)
    if imdims[0] eq 1 then tmp[edge]=im
    if imdims[0] eq 2 then tmp[edge,edge]=im
    if imdims[0] eq 3 then tmp[edge,edge,edge]=im
    if imdims[0] eq 4 then tmp[edge,edge,edge,edge]=im
endif
if  edge lt 0 then begin
    imdims=size(im)
    if imdims[0] eq 1 then tmp=im[-edge:(edge+imdims[1]-1)]
    if imdims[0] eq 2 then tmp=im[-edge:(edge+imdims[1]-1),-edge:(edge+imdims[2]-1)]
    if imdims[0] eq 3 then tmp=im[-edge:(edge+imdims[1]-1),-edge:(edge+imdims[2]-1),-edge:(edge+imdims[3]-1)]
    if imdims[0] eq 4 then tmp=im[-edge:(edge+imdims[1]-1),-edge:(edge+imdims[2]-1),-edge:(edge+imdims[3]-1),-edge:(edge+imdims[4]-1)]
endif
if  edge eq 0 then tmp=im

return,tmp
END


FUNCTION GENMASK, im, err=err, hd=hd,$
                  spar=spar,$
                  sig=sig,$
                  chmin=chmin,$
                  grow=grow,$
                  guard=guard,$
                  im_use=im_smo,$
                  err_use=sen_smo
;+
; NAME:
;   GENMASK
;
; PURPOSE:
;   generate a signal mask from a data cube using sigma clipping, with optional
;   dilated masking and pre-smoothing.
;
; INPUTS:
;   IM        data cube
;   ERR       error cube
;   [HD]      data header (required if spar[0]>0.0)
;   [SPAR]    2-element vector
;             We will degrade the cube angular resolution to <spar[0]> arcsec, and
;             smooth spectra with a Gaussian function of fwhm=<spar[1]> km/s for 
;             signal detection.
;   [SIG]     detection threshold, in units of cube/smoothed cube sigma
;   [CHMIN]   minimum number of velocity channels in initial mask
;   [GROW]    threshold for expanding core mask, in units of cube/smoothed cube sigma
;   [GUARD]   thickness of guard band in pixels
;
; OUTPUTS:
;   MASK      signal mask cube (1 as signal, 0 as noise, nan as missing data)
;   IM_USE    actual data cube used for SNR calculations
;   ERR_USE   actual error cube used for SNR calculations
;
; HISTORY:
;
;   20130624  RX  introduced (for replacing grmask.pro and gsmask.pro)
;   20150528  TW  removed eprop option.
;   20161112  TW  add CHMIN parameter
;   20170824  RX  add compatibility with 2D images
;   
;-

 
forward_function dilate_mask

if  not keyword_set(spar)   then spar=[0.0,0.0]
if  not keyword_set(sig)    then sig =0.0
if  not keyword_set(grow)   then grow=0.0
if  not keyword_set(chmin)  then chmin=1 
if  not keyword_set(guard)  then guard=0
im_smo=im
sen_smo=err

; SMOOTH DATA IF REQUESTED
if  spar[0] gt 0.0 then begin
    nchan=1
    if  size(im,/n_d) eq 3 then nchan=(size(im,/d))[2]
    if  n_elements(spar) eq 1 then spar=[spar,0.0]
    SMOOTH3D,im/err,hd,im_smo,hd_smo,[spar[0],spar[0],0.],svel=spar[1],ifail=ifail
    if  ifail ne 0 then begin
        errmsg = 'ERROR - the target beam is too small!'
        message,'ERROR - ' + errmsg,level=1
    endif
    sen_smo=ERR_CUBE(im_smo,planes=indgen(nchan))
    ; calculate smoothed rms after clipping based on MIRIAD sigest.for
    mask1 = abs(im_smo) lt 2.5*sqrt(!PI/2)*meanabsdev(im_smo,/nan)
    sen_smo=ERR_CUBE(im_smo,planes=indgen(nchan),mask=mask1)
endif

; GENERATE MASK BY SIGMA CLIPPING
if  sig gt 0.0 then begin
    mask = im*0.0
    mask[where(im_smo gt sen_smo*sig, ct, /null)]=1.0
;    print,'Before dilation:',size(mask),total(mask,/nan)
;    print,'Pixel count is',ct
    ; REQUIRE MINIMUM NUMBER OF CHANNELS
    if chmin gt 1 then begin
        mask = padding(mask,chmin-1)
        for j = 2, chmin do begin
            shmask = mask*0.0
            for k = 1, j-1 do begin
                shmask += shift(mask,0,0,k) + shift(mask,0,0,-k)
            endfor
            mask = mask * (shmask-j+2) gt 0
        endfor
    endif
    ; DILATE THE MASK IF REQUESTED
    if (total(mask,/nan) gt 0.0) and (grow gt 0.0) then begin
        ;print,'Excluding loners:',size(mask),total(mask,/nan)
        ; CONSTRAINT MASK
        growmask = im_smo gt grow*sen_smo
        ; REQUIRE MINIMUM NUMBER OF CHANNELS (CURRENTLY DISABLED)
        if chmin gt 1 then begin
            growmask = padding(growmask,chmin-1)
;            for j = 2, chmin do begin
;                shmask = growmask*0.0
;                for k = 1, j-1 do begin
;                    shmask += shift(growmask,0,0,k) + shift(growmask,0,0,-k)
;                endfor
;                growmask = growmask * (shmask-j+2) gt 0
;            endfor
        endif
        ;print,'Growth mask:',size(growmask),total(growmask,/nan)
        ; EXPAND MASK
        mask = float(dilate_mask(mask, constraint = growmask))
        ;print,'Expanded mask:',size(mask),total(mask,/nan)
    endif
    if chmin gt 1 then mask = padding(mask,-chmin+1)  
    ; ADD A GUARD BAND IF REQUESTED
    if total(guard) gt 0 then begin
        mask = MASKGUARD_3D(mask, guard = guard)
    endif
endif else begin
    mask = im*0.0 + 1.0
endelse

return,mask
  
END
