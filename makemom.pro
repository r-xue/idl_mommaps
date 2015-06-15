PRO MAKEMOM, filename, errfile=errfile, rmsest=rmsest, maskfile=maskfile, $
      vrange=vrange, xyrange=xyrange, dvmin=dvmin, baseroot=baseroot, $
      thresh=thresh, edge=edge, guard=guard, smopar=smopar, senmsk=senmsk, $
      replace0=replace0, mask0=mask0, kelvin=kelvin, dorms=dorms, pvmom0=pvmom0
      
;+
; NAME:
;   MAKEMOM
;
; PURPOSE:
;   produce moment maps with masking
;
; ARGUMENTS:
;   FILENAME  --  FITS data cube [no default].  3rd axis should be velocity (M/S).
;   ERRFILE   --  FITS error map or cube.  Note that any pixels blanked here will 
;                 result in the corresponding data pixels being dropped.
;                 default: unset => error is assumed constant with position
;   RMSEST    --  estimate of channel noise, in same units as cube.  If ERRFILE is
;                 given, this is taken to be the noise at the minimum (e.g. center).
;                 default: unset => use error map if given, or assume rms=1, or 
;                 get rms from data cube if /DORMS is requested.
;   MASKFILE  --  external masking cube (value=1 for valid data, 0 otherwise).
;                 Must match the coordinate grid of the input data cube.  This is 
;                 applied AFTER any other requested masking operations.
;                 default: unset => mask is generated by the program
;   VRANGE    --  choose the velocity range for moment calculations [km/s, km/s]
;                 default: unset => whole velocity range is used
;   XYRANGE   --  choose the xy-region for moment calculations [x0, x1, y0, y1]
;                 default: unset => whole image region is used
;   DVMIN     --  typical linewidth in km/s for calculating minimum value in emom0 map.
;                 default: unset => no floor applied in emom0 map 
;   BASEROOT  --  root name for output files
;                 default: based on FILENAME root
;   THRESH    --  starting intensity threshold for signal mask, in units of sigma
;                 default: unset => no masking applied
;   EDGE      --  final threshold for dilated mask, in units of sigma
;                 default: unset => mask is not dilated
;   GUARD     --  width of guard band in x, y, v dimensions [pix, pix, pix]
;                 The mask is extended by the specified no. of pixels in each direction
;                 if THRESH and/or EDGE are also specified.  If only one value is given
;                 it is used for v; if two values are given they are used for x, y.
;                 default: unset => mask is not extended in any direction
;   SMOPAR    --  for mask generation, first degrade the cube angular resolution to 
;                 <smopar[0]> arcsec and/or smooth spectra with a Gaussian function of 
;                 fwhm=<smopar[1]> km/s for signal detection.
;                 default: [0.,0.] => no pre-smoothing applied
;   SENMSK    --  regions where the noise level is higher than <senmsk>*min(errcube)
;                 (e.g. edges of field-of-view) will be masked out.
;                 default: unset => no additional masking applied
;   
; OPTION KEYWORDS:
;   REPLACE0  --  input cube values of 0.0 will be treated as missing data
;   MASK0     --  mask pixels where some channels are blanked (e.g. FOV varies w/channel)
;   KELVIN    --  force conversion from Jy/beam to Kelvin
;   DORMS     --  estimate the channel noise from the data.  This is used to scale
;                 ERRFILE if given.  Cannot be used if RMSEST is given.
;   PVMOM0    --  produce mom0 images by collapsing cube along x & y axes 
;                 (not constant RAs or Decs!)
;   
; OUTPUTS:
;   <BASEROOT>.*.fits
;
; HISTORY:
;
;   20150527  tw  initial version based on Rui Xue's idl_moments code
;   20150601  tw  implement guard parameter
;   20150610  tw  give estr format I0
;
;-

; DEFAULT SETTINGS AND FILE NAMES
if  n_elements(smopar) eq 0 then begin
    smopar = [0.,0.]
    smostr = ''
endif else if  n_elements(smopar) eq 1 then begin
    smopar = [smopar[0],0.]
    smostr = 'sm'+string(smopar[0],format='(I0)')+'_'
endif else begin
    smostr = 'sm'+string(smopar[0],format='(I0)')+'v'+$
        string(smopar[1],format='(I0)')+'_'
endelse
if  keyword_set(thresh) then begin
    tstr = 't'+string(thresh,format='(F3.1)')
endif else begin
    thresh = 0.
    tstr = ''
endelse
if  keyword_set(edge) then begin
    estr = 'e'+string(edge,format='(I0)')
endif else begin
    edge = 0.
    estr = ''
endelse
if  keyword_set(guard) then begin
    gstr = 'g'+string(guard[0],format='(I0)')
endif else begin
    guard = [0, 0, 0]
    gstr = ''
endelse
galname=file_basename(filename,'.fits')
if  keyword_set(baseroot) eq 0 then begin
    baseroot = galname + '.' + smostr + tstr + estr + gstr
    if strpos(baseroot,'.',/reverse_search) eq strlen(baseroot)-1 then $
      baseroot=strmid(baseroot,0,strlen(baseroot)-1)
endif
if  keyword_set(errfile) then errname=file_basename(errfile,'.fits')
if  keyword_set(maskfile) then mskname=file_basename(maskfile,'.fits')

; READ IN DATA, CONVERT UNITS
data = READFITS(filename, hd, /silent)
RADIOHEAD, hd, s = h
if  keyword_set(kelvin) and strpos(strupcase(sxpar(hd,'BUNIT')),'JY/B') ne -1 then begin
    data = data * h.jypb2k
    SXADDPAR, hd, 'BUNIT', 'K'
endif

; EXTRACT SUBREGION IF REQUESTED
if  n_elements(xyrange) eq 4 then begin
    HEXTRACT3D, data, hd, tmp, tmphd, xyrange
    SXADDPAR, tmphd, 'DATAMAX', max(tmp,/nan), before='HISTORY'
    SXADDPAR, tmphd, 'DATAMIN', min(tmp,/nan), before='HISTORY'
    ;WRITEFITS, galname+'.subreg.fits', float(tmp), tmphd
    data=tmp
    hd=tmphd
endif

; OUTPUT SOME INFO
print,replicate('-',35)
print,'spectral cube size: ',size(data,/d)
print,replicate('-',35)

; BLANKING
if  keyword_set(replace0) then begin
    data[where(data eq 0.0,/null)]=!values.f_nan
endif

; APPLY INPUT MASK AND VELOCITY WINDOW IF GIVEN
sz = size(data)
if  n_elements(maskfile) eq 0 then begin 
    exmask=data*0.0+1.0
endif else begin
    exmask=READFITS(maskfile, exmaskhd, /silent)
    if  n_elements(xyrange) eq 4 then begin
        HEXTRACT3D, exmask, exmaskhd, tmp, tmphd, xyrange
        SXADDPAR, tmphd, 'DATAMAX', max(tmp,/nan), before='HISTORY'
        SXADDPAR, tmphd, 'DATAMIN', min(tmp,/nan), before='HISTORY'
        ;WRITEFITS, mskname+'.subreg.fits', float(tmp), tmphd
        exmask=tmp
        exmaskhd=tmphd
    endif
endelse
if  n_elements(vrange) ne 0 then begin
    tag_outvrange = where(h.v lt vrange[0] or h.v gt vrange[1])
    if tag_outvrange[0] ne -1 then exmask[*,*,[tag_outvrange]]=0.0
endif else begin
    vrange = [min(h.v),max(h.v)]
    tag_outvrange = -1
endelse

; GENERATE ERROR CUBE AND 2D MAP IF NOT GIVEN
if  n_elements(errfile) eq 0 then begin 
    if  keyword_set(dorms) then begin
        ecube=ERR_CUBE(data)
    endif else begin
        if keyword_set(rmsest) then rms=rmsest else rms=1.
        ecube=data*0.0 + rms
    endelse
    emap = total(ecube, 3, /nan) / (total(ecube eq ecube, 3)>1)
    emap[where(emap eq 0.0,/null)]=!values.f_nan
; OR ELSE USE PROVIDED ERROR CUBE
endif else begin
    ecube = readfits(errfile, ehd, /silent)
    if  n_elements(xyrange) eq 4 then begin
        HEXTRACT3D,ecube,ehd,tmp,tmphd,xyrange
        SXADDPAR, tmphd, 'DATAMAX', max(tmp,/nan), before='HISTORY'
        SXADDPAR, tmphd, 'DATAMIN', min(tmp,/nan), before='HISTORY'
        ;WRITEFITS, errname+'.subreg.fits', float(tmp), tmphd
        ecube=tmp
        ehd=tmphd
    endif
    esz = size(ecube)
    if  esz[0] eq 2 then begin
        ecube0=ecube
        ecube=make_array(esz[1],esz[2],sz[3],/float,/nozero)
        for i=0,sz[3]-1 do ecube[0,0,i]=ecube0
        ;ecube[where(data ne data,/null)]=!values.f_nan
    endif
    ecube[where(ecube eq 0,/null)]=!values.f_nan
    if  keyword_set(dorms) then begin
        tmp=ERR_CUBE(data,pattern=ecube)
        ecube=tmp
    endif else begin
        if keyword_set(rmsest) then ecube=rmsest*ecube/min(ecube,/nan)
    endelse
    if  keyword_set(kelvin) and strpos(strupcase(sxpar(ehd,'BUNIT')),'JY/B') ne -1 then begin
        ecube=temporary(ecube) * h.jypb2k
        SXADDPAR, ehd, 'BUNIT', 'K'
    endif
    emap = total(ecube, 3, /nan) / (total(ecube eq ecube, 3)>1)
    emap[where(emap eq 0.0,/null)]=!values.f_nan
    data[where(ecube ne ecube,/null)]=!values.f_nan
endelse

; GENERATE MASKING CUBE (missing data locations are still kept in mask)
mask = GENMASK(data,err=ecube,hd=hd,spar=smopar,sig=thresh,grow=edge,guard=guard)
mask = mask*exmask

; MASK HIGH NOISE AT EDGES
if  n_elements(senmsk) ne 0 then begin
    mask1d=float((total(ecube, 3) le senmsk*min(total(ecube, 3),/nan)))
    mask1d[where(mask1d eq 0.0,/null)]=!values.f_nan
    for i = 0, sz[3]-1 do begin
        data[*,*,i]  = data[*,*,i]  + mask1d - mask1d
        ecube[*,*,i] = ecube[*,*,i] + mask1d - mask1d
        mask[*,*,i]  = mask[*,*,i]  + mask1d - mask1d
    endfor
endif

; CALCULATE MOMENTS
;   h.v: velocity in km/s
if  n_elements(dvmin) ne 0 then begin
    nchmin=ceil(dvmin/(abs(h.cdelt[2])/1.0e3))*1.0
endif else begin
    nchmin=0.
endelse
MASKMOMENT, data, mask, ecube, h.v, $
            mom0 = mom0, mom1 = mom1, mom2 = mom2, $
            emom0 = emom0, emom1 = emom1, emom2 = emom2, $
            peak = peak, snrpk=snrpk,$
            mask0=mask0, nchmin=nchmin
mhd = hd
;cmds = RECALL_COMMANDS()

histlabel = 'IDL_MOMMAPS: '
SXADDPAR, mhd, 'HISTORY', histlabel+systime()
SXADDPAR, mhd, 'HISTORY', histlabel+'filename='+filename
if n_elements(errfile) gt 0 then $
    SXADDPAR, mhd, 'HISTORY', histlabel+'errfile='+errfile
SXADDPAR, mhd, 'HISTORY', histlabel+'smopar=['+strcompress(smopar[0],/r)+','+$
    strcompress(smopar[1],/r)+']'
SXADDPAR, mhd, 'HISTORY', histlabel+'thresh='+strcompress(thresh,/r)
SXADDPAR, mhd, 'HISTORY', histlabel+'edge='+strcompress(edge,/r)
SXADDPAR, mhd, 'HISTORY', histlabel+'guard=['+strcompress(guard[0],/r)+$
    ','+strcompress(guard[1],/r)+','+strcompress(guard[2],/r)+']'
if n_elements(dvmin) gt 0 then $
    SXADDPAR, mhd, 'HISTORY', histlabel+'dvmin='+strcompress(dvmin,/r)

; OUTPUT MASK CUBE
SXADDPAR,mhd,'DATAMAX',2.0, before='HISTORY'
SXADDPAR,mhd,'DATAMIN',-1.0, before='HISTORY'
nan_tag=where(data ne data,nan_ct)
if  nan_ct ne 0 then mask[nan_tag]=!values.f_nan
if  thresh gt 0.0 then begin
    WRITEFITS,baseroot+'.mask.fits',float(mask),mhd
    SXADDPAR,mhd,'DATAMAX', max(mask*data,/nan), before='HISTORY'
    SXADDPAR,mhd,'DATAMIN', min(mask*data,/nan), before='HISTORY'
    WRITEFITS,baseroot+'.mskd.fits',mask*data,mhd
endif

; OUTPUT 2D ERROR MAP (SKIP IF ERROR MAP PROVIDED)
if  n_elements(errfile) eq 0 then begin 
    SXADDPAR, mhd, 'DATAMAX', max(emap,/nan), before='HISTORY'
    SXADDPAR, mhd, 'DATAMIN', min(emap,/nan), before='HISTORY'
    WRITEFITS, baseroot+'.rms.fits', float(emap), mhd
endif else if keyword_set(dorms) or keyword_set(rmsest) then begin
    SXADDPAR, mhd, 'DATAMAX', max(ecube,/nan), before='HISTORY'
    SXADDPAR, mhd, 'DATAMIN', min(ecube,/nan), before='HISTORY'
    WRITEFITS, baseroot+'.ecube.fits', float(ecube), mhd
endif

; PEAK INTENSITY
SXADDPAR, mhd, 'DATAMAX', max(peak,/nan), before='HISTORY'
SXADDPAR, mhd, 'DATAMIN', min(peak,/nan), before='HISTORY'
WRITEFITS, baseroot+'.peak.fits', float(peak), mhd
bunit = SXPAR(hd,'BUNIT')

; PEAK SNR
SXADDPAR, mhd, 'BUNIT', 'SNR', before='HISTORY'
SXADDPAR, mhd, 'DATAMAX', max(snrpk,/nan), before='HISTORY'
SXADDPAR, mhd, 'DATAMIN', min(snrpk,/nan), before='HISTORY'
WRITEFITS, baseroot+'.snrpk.fits', float(snrpk), mhd

; MOMENT 0 AND ERROR
SXADDPAR, mhd, 'BUNIT', strtrim(bunit,2)+'.KM/S', before='HISTORY'
SXDELPAR, mhd, 'CTYPE3'
SXDELPAR, mhd, 'CRVAL3'
SXDELPAR, mhd, 'CRPIX3'
SXDELPAR, mhd, 'CDELT3'
mom0gm = mom0 * abs(h.cdelt[2]) / 1.0e3
SXADDPAR, mhd, 'DATAMAX', max(mom0gm,/nan)
SXADDPAR, mhd, 'DATAMIN', 0, before='HISTORY'
WRITEFITS, baseroot+'.mom0.fits', float(mom0gm), mhd
emom0gm = emom0 * abs(h.cdelt[2]) / 1.0e3
SXADDPAR, mhd, 'DATAMAX', max(emom0gm,/nan), before='HISTORY'
SXADDPAR, mhd, 'DATAMIN', min(emom0gm,/nan), before='HISTORY'
WRITEFITS, baseroot+'.emom0.fits', float(emom0gm), mhd
print,'moment 0 flux: ',total(mom0,/nan),' ',strtrim(bunit,2)+'.KM/S.PIX'

; MOMENT 1 AND ERROR
SXADDPAR, mhd, 'BUNIT', 'KM/S', before='HISTORY'
SXADDPAR, mhd, 'DATAMAX', max(mom1,/nan), before='HISTORY'
SXADDPAR, mhd, 'DATAMIN', min(mom1,/nan), before='HISTORY'
WRITEFITS, baseroot+'.mom1.fits', float(mom1), mhd
SXADDPAR, mhd, 'DATAMAX', max(emom1,/nan), before='HISTORY'
SXADDPAR, mhd, 'DATAMIN', min(emom1,/nan), before='HISTORY'
WRITEFITS, baseroot+'.emom1.fits', float(emom1), mhd
PLTMOM, baseroot

; MOMENT 2 AND ERROR
SXADDPAR, mhd, 'DATAMAX', max(mom2,/nan), before='HISTORY'
SXADDPAR, mhd, 'DATAMIN', min(mom2,/nan), before='HISTORY'
WRITEFITS, baseroot+'.mom2.fits', float(mom2), mhd
SXADDPAR, mhd, 'DATAMAX', max(emom2,/nan), before='HISTORY'
SXADDPAR, mhd, 'DATAMIN', min(emom2,/nan), before='HISTORY'
WRITEFITS, baseroot+'.emom2.fits', float(emom2), mhd

; MOM0 in XV and YV (OPTIONAL)
if  keyword_set(pvmom0) then begin
    MASKMOMENT_PV,data,hd,mask,mom0xv,mom0vy,mom0xvhd=mom0xvhd,mom0vyhd=mom0vyhd,$
        vrange=vrange
    WRITEFITS, baseroot+'.mom0xv.fits',mom0xv,mom0xvhd
    WRITEFITS, baseroot+'.mom0vy.fits',mom0vy,mom0vyhd
    PLTMOM_PV, baseroot
endif

END
