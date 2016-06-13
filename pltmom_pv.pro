PRO PLTMOM_PV,prefix,label=label
;+
; NAME:
;   PLTMOM_PV
;
; PURPOSE:
;   summary plot for mom0 and mom0-pv maps
;
; INPUTS:
;   PREFIX    -- prefix for maps
;                e.g.'yourpath/n4254co.gm'
;
; OUTPUTS:
;   prefix.mom0pv.eps
;
; HISTORY:
;
;   20130214  RX  introduced
;   20150910  RX  merged from makemompv_pl.pro
;-

;   CHECK IMAGE FILES

im=prefix+'.mom0.fits'
eim=prefix+'.emom0.fits'
imsnr=prefix+'.snrpk.fits'
imxv=prefix+'.mom0xv.fits'
imvy=prefix+'.mom0vy.fits'

if  ~file_test(im) and $
    ~file_test(eim) and $
    ~file_test(imxv) and $
    ~file_test(imvy) then begin
    print,'no mom0/emom0/mom0xy/mom0vy image found'
    return
endif

im=READFITS(im,imhd,/SILENT)
eim=READFITS(eim,eimhd,/SILENT)
imxv=READFITS(imxv,imxvhd,/SILENT)
imvy=READFITS(imvy,imvyhd,/SILENT)

sz_xv=size(imxv,/d)
sz_vy=size(imvy,/d)
sz=size(im,/d)
scy=sz_xv[1]*1.0/sz[1]
scx=sz_vy[0]*1.0/sz[0]

nvels=(size(imvy))[1]
vels=sxpar(imvyhd,'CRVAL1')+(findgen(nvels)+1-sxpar(imvyhd,'CRPIX1'))*sxpar(imvyhd,'CDELT1')
vels=vels/1000.

szall=[sz[0]+sz_vy[0],sz[1]+sz_xv[1]]

posx=[0,0,sz[0],sz[0],szall[0],szall[0]]
posx[1:5]=posx[1:5]+szall[0]*0.15       ;   left-edge
posx[3:5]=posx[3:5]+max(szall)*0.01     ;   panel gap
posx[5]=posx[5]+szall[0]*0.04

posy=[0,0,sz[1],sz[1],szall[1],szall[1]]
posy[1:5]=posy[1:5]+szall[1]*0.10       ;   bottom-edge
posy[3:5]=posy[3:5]+max(szall)*0.01     ;   panel gap
posy[5]=posy[5]+szall[1]*0.04

posx=posx*1.0
posy=posy*1.0
posxn=posx/max(posx)
posyn=posy/max(posy)
epsxy=[posx[-1],posy[-1]]/max([posx[-1],posy[-1]])*11.0

set_plot, 'ps'
device, filename=prefix+'.mom0pv.eps', $
    bits_per_pixel=8,/encapsulated,$
    xsize=epsxy[0],ysize=epsxy[1],/inches,/col,xoffset=0,yoffset=0
;print,epsxy
!p.thick = 1.5
!x.thick = 1.5
!y.thick = 1.5
!z.thick = 1.5
!p.charsize=1.0
!p.charthick=1.5
xyouts,'!6'

; MOM-0 (XY) PLOT 
pos=[posxn[1],posyn[1],posxn[2],posyn[2]]
loadct,13,/silent
cgimage,im,pos=pos,stretch=1,/noe
imcontour,im,imhd,nlevels=10,$
    /noe,pos=pos,/nodata,color='red',AXISCOLOR='red',subtitle=' '
imcontour,im,imhd,nlevels=10,$
    /noe,pos=pos,/nodata,color='black',AXISCOLOR='black',subtitle=' ',ticklen=0
if  file_test(imsnr) then begin
    imsnr=readfits(imsnr)
    sz=size(imsnr,/d)
    bb=find_boundary(where(imsnr eq imsnr),xsize=sz[0],ysize=sz[1])
    oplot,bb[0,*],bb[1,*],color=cgcolor('yellow'),linestyle=2
endif

getrot,imhd,rotang,cdelt
imsz=size(im,/d)
rotang_s=rotang+45.0
ds=30./2.0/512*min(imsz)*sqrt(2)
xc=imsz[0]*0.92
yc=imsz[1]*0.08
one_arrow,xc+ds*sin(rotang_s/180*!dpi),yc-ds*cos(rotang_s/180*!dpi),+90+rotang,'N',color='yellow',/data,charsize=1.0
one_arrow,xc+ds*sin(rotang_s/180*!dpi),yc-ds*cos(rotang_s/180*!dpi),+180+rotang,'E',color='yellow',/data,charsize=1.0
    
RADIOHEAD,imhd,s=s
psize=abs(s.cdelt[0])*3600
tvellipse,s.bmaj/2.0/psize,s.bmin/2.0/psize,$
    sz[0]/10.0,sz[1]/10.0,$
    s.bpa-90.0+rotang,$
    /data,noclip=0,color=cgcolor('cyan'),/fill

if  keyword_set(label) then al_legend,label,/top,/right,textcolor='yellow',box=0

vinterval=ceil(abs(min(vels)-max(vels))/4.0)*1.0
    
subpos_xv=[posxn[1],posyn[3],posxn[2],posyn[4]]

loadct,13,/silent
cgimage,imxv,pos=subpos_xv,stretch=1,/noe
loadct,0,/silent

plot,[0,1],[min(vels),max(vels)],/nodata,/noe,pos=subpos_xv,$
    xstyle=5,ystyle=1,yTICKINTERVAL=vinterval,$
    ytitle='Velocity [km/s]',color=cgcolor('red')
plot,[0,1],[min(vels),max(vels)],/nodata,/noe,pos=subpos_xv,$
    xstyle=5,ystyle=1,yTICKINTERVAL=vinterval,$
    ytitle='Velocity [km/s]',ticklen=0.0,ytick_get=vticks   
imcontour,im,imhd,nlevels=10,$
    /noe,pos=subpos_xv,/nodata,color='red',AXISCOLOR='red',ystyle=5,xtickformat='(A1)',$
    subtitle=' ',xticklen=!p.ticklen/scy,xtitle=' ',ytitle=' '
    
subpos_yv=[posxn[3],posyn[1],posxn[4],posyn[2]]

loadct,13,/silent
cgimage,imvy,pos=subpos_yv,stretch=1,/noe;,minvalue=0.0
loadct,0,/silent
plot,[min(vels),max(vels)],[0,1],/nodata,/noe,pos=subpos_yv,$
    xstyle=1,ystyle=5,xticks=3,xtickinterval=vinterval,$
    xtitle='',color=cgcolor('red'),xtickformat='(A1)'
plot,[min(vels),max(vels)],[0,1],/nodata,/noe,pos=subpos_yv,$
    xstyle=1,ystyle=5,xticks=3,xtickinterval=vinterval,$
    xtitle='',ticklen=0.0,xtickformat='(A1)'
imcontour,im,imhd,nlevels=10,$
    /noe,pos=subpos_yv,/nodata,color='red',AXISCOLOR='red',xstyle=5,ytickformat='(A1)',$
    subtitle=' ',yticklen=!p.ticklen/scx,xtitle=' ',ytitle=' '

        
device, /close
set_plot,'X'

END
