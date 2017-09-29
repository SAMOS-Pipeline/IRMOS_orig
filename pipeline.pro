pro pipeline

while !d.window ne -1 do wdelete,!d.window
device,decomposed=0
loadct,5

;=====================================================================================
;          MMMM   MMMMMMMMMMD   DMMMMM    MMMMM8   ..MMMMMMMN.     MMMMMMMM .
;         MMMM   MMMMMMMMMMMM   MMMMMM.  MMMMMM.  MMMMMMMMMMMM   MMMM$~MMMMM.
;        MMMM$  ~MMMM.  OMMMM. MMM.MM~ .MMDMMM  .MMMM    MMMMM  MMMM   NMD?.
;        MMMM   MMMMOOMMMMMO  MMMMMMM.=MM MMM. .MMMM.    MMMM.. MMMMMMMM8.
;       MMMM   MMMMMMMMMM~.  .MMM MMM8MM.MMMM  MMMM     ~MMMM   .MMMMMMMMM
;      MMMM   MMMM8 MMMMM    MMM .MMMMM ,MMM  MMMM.    .MMMM .   .   MMMMM.
;     IMMMM. .MMMM  MMMMM   MMMN.MMMMM  MMM   MMMMM.  MMMMM. MMMM.  .MMMM
;    .MMMM   MMMM.  :MMMM. MMMM  MMMM. MMM$.  MMMMMMMMMMM    MMMMMMMMMM
;    MMMM.  MMMM..  .MMMM  MMM.  MMM  NMMM     .MMMMMM       .ZMMMMM+..
;=====================================================================================

;==============================
;the paragraph below containes the definination of the dataset we are using
;examples are contained in the files Run10/NightX_data.txt
;==============================

;------------------------------
;1° Orion field
;------------------------------
PREPATH='C:\IRMOS_pipeline\'
DATA_PATH='Run10\Night8\'
SOURCE_NAME='M42_Field18a_K3000_'
MODE='nodding' ;choose 'standard' or 'nodding'
FILTER='Kblue'
OBJTYPE='orion'
IMG_frames=[5]
SCI_A_frames=[6,8,10,12,14,16,18,20,22,24,26,28,30]   ;SCI_A and SCI_B refer to 2 nodding positions (shifting the stars left and right into the slit
SCI_B_frames=[7,9,11,13,15,17,19,21,23,25,27,29,31]
DK_frames=[32]  ;dark
DKFLT_frames=[33,34,35,36,37,38]  ;dark flat
FLT_frames=[39,40,41,42,43,44] ;flat
oldvectorshift=[-174,-200,-258,-65,0,-38,-110,-24,169]
rah=5 & ram=35 & ras=19.47   ;coordinates
deg=-5 & dem=27 & des=06.4



;-------------------------------------------------------------------------------------

mask=readfits('D:\science\IRMOS_DATA\Run10\newmaskcod.fits')

;-------------------------------------------------------------------------------------

flagnormy_sci=0      & stripey_sci=10       &  gapy_sci=251
flagnormy_dark=0     & stripey_dark=20      &  gapy_dark=30
flagnormy_darkflat=1 & stripey_darkflat=6  &  gapy_darkflat=1
flagnormy_flat=1     & stripey_flat=6      &  gapy_flat=1

flagnormx_sci=0      & stripex_sci=10       &  gapx_sci=1
flagnormx_dark=0     & stripex_dark=7      &  gapx_dark=1
flagnormx_darkflat=1 & stripex_darkflat=6 &  gapx_darkflat=1
flagnormx_flat=1     & stripex_flat=6    &  gapx_flat=1



shiftscreen=0
shiftnodding=15
flagoldvectorshift=0. ; & oldvectorshift=[269,288,-96,0,71,-52,1,-100,-117]

ohmethod='rousselot'
minOHflux=250.

;=====================================================================================
;=====================================================================================
;=====================================================================================
;-------------------------------------------------------------------------------------
;READING DATA
;-------------------------------------------------------------------------------------

PATH=PREPATH+DATA_PATH

if (mode eq 'standard') then begin
    dummy=size(SCI_frames)
    Nsource_frames=dummy[1]
    IRSPEC_load_random,PATH,SOURCE_NAME,SCI_frames,SOURCE_DCube   ;read sci frames
endif else if (mode eq 'nodding') then begin
    dummy=size(SCI_A_frames)
    Nsource_A_frames=dummy[1]
    dummy=size(SCI_B_frames)
    Nsource_B_frames=dummy[1]
    IRSPEC_load_random,PATH,SOURCE_NAME,SCI_A_frames,SOURCE_A_DCube   ;read si frames for nodding positions A and B
    IRSPEC_load_random,PATH,SOURCE_NAME,SCI_B_frames,SOURCE_B_DCube
endif
if (DK_frames eq ['none']) then begin
    DKSource_DCube=fltarr(1024,1024,1)
endif else begin
    IRSPEC_load_random,PATH,SOURCE_NAME,DK_frames,DKSource_DCube
endelse

IRSPEC_load_random,PATH,SOURCE_NAME,FLT_frames,Flat_DCube
IRSPEC_load_random,PATH,SOURCE_NAME,DKFLT_frames,DKFlat_DCube

if (objtype eq 'orion') then begin
    IRSPEC_load_random,PATH,SOURCE_NAME,IMG_frames,reference
    filenameout = PREPATH+DATA_PATH+'\red\'+Source_name+'reference.fits'
    reference2=transpose(reference)
    writefits,filenameout,reference2
    reference=0
    reference2=0
endif


;-------------------------------------------------------------------------------------
;NORMALIZING BIAS VARIATIONS (freaking crappy instrument had inconsistent bias level in the 4 quadrants of the detector)
;-------------------------------------------------------------------------------------

dummy=size(DK_frames) & numerodark=dummy[1]
dummy=size(DKFLT_frames) & numerodarkflat=dummy[1]
dummy=size(FLT_frames) & numeroflat=dummy[1]
levelDK=fltarr(numerodark)
levelDKFLT=fltarr(numerodarkflat)
levelFLT=fltarr(numeroflat)
if (mode eq 'standard') then begin
    dummy=size(SCI_frames) & numerosci=dummy[1]
    levelSCI=fltarr(numerosci)
endif else begin
    dummy=size(SCI_A_frames) & numerosciA=dummy[1]
    dummy=size(SCI_B_frames) & numerosciB=dummy[1]
    levelSCIA=fltarr(numerosciA)
    levelSCIB=fltarr(numerosciB)
endelse
for i=0,(numerodark-1) do begin
    image=DKSource_Dcube[*,*,i]
    if (flagnormx_dark eq 1) then  IRMOS_normchips_x,image,stripex_dark,gapx_dark
    if (flagnormy_dark eq 1) then  IRMOS_normchips_y,image,stripey_dark,gapy_dark
    DKSource_Dcube[*,*,i]=image
    levelDK[i]=median(image)
endfor
meanlevelDK=mean(levelDK)
for i=0,(numerodarkflat-1) do begin
    image=DKFlat_Dcube[*,*,i]
    if (flagnormx_darkflat eq 1) then  IRMOS_normchips_x,image,stripex_darkflat,gapx_darkflat
    if (flagnormy_darkflat eq 1) then  IRMOS_normchips_y,image,stripey_darkflat,gapy_darkflat
    DKFlat_Dcube[*,*,i]=image
    levelDKFLT[i]=median(image)
endfor
meanlevelDKFLT=mean(levelDKFLT)
for i=0,(numeroflat-1) do begin
    image=Flat_Dcube[*,*,i]
    if (flagnormx_flat eq 1) then IRMOS_normchips_x,image,stripex_flat,gapx_flat
    if (flagnormy_flat eq 1) then IRMOS_normchips_y,image,stripey_flat,gapy_flat
    Flat_Dcube[*,*,i]=image
    levelFLT[i]=median(image)
endfor
meanlevelFLT=mean(levelFLT)
if (mode eq 'standard') then begin
    for i=0,(numerosci-1) do begin
       image=Source_Dcube[*,*,i]
       if (flagnormx_sci eq 1) then IRMOS_normchips_x,image,stripex_sci,gapx_sci
       if (flagnormy_sci eq 1) then IRMOS_normchips_y,image,stripey_sci,gapy_sci
       Source_Dcube[*,*,i]=image
       levelSCI[i]=median(image)
    endfor
    meanlevelSCI=mean(levelSCI)
endif else begin
    for i=0,(numerosciA-1) do begin
       image=Source_A_Dcube[*,*,i]
       if (flagnormx_sci eq 1) then IRMOS_normchips_x,image,stripex_sci,gapx_sci
       if (flagnormy_sci eq 1) then IRMOS_normchips_y,image,stripey_sci,gapy_sci
       Source_A_Dcube[*,*,i]=image
       levelSCIA[i]=median(image)
    endfor
    for i=0,(numerosciB-1) do begin
       image=Source_B_Dcube[*,*,i]
       if (flagnormx_sci eq 1) then IRMOS_normchips_x,image,stripex_sci,gapx_sci
       if (flagnormy_sci eq 1) then IRMOS_normchips_y,image,stripey_sci,gapy_sci
       Source_B_Dcube[*,*,i]=image
       levelSCIB[i]=median(image)
    endfor
    meanlevelSCIA=mean(levelSCIA)
    meanlevelSCIB=mean(levelSCIB)
endelse
;applying relative bulk normalizaziont
for i=0,(numerodark-1) do begin
    DKSource_Dcube[*,*,i] = DKSource_Dcube[*,*,i] + meanlevelDK - levelDK[i]
endfor
for i=0,(numerodarkflat-1) do begin
    DKFlat_Dcube[*,*,i] = DKFlat_Dcube[*,*,i] + meanlevelDKFLT - levelDKFLT[i]
endfor
for i=0,(numeroflat-1) do begin
    Flat_Dcube[*,*,i] = Flat_Dcube[*,*,i] + meanlevelFLT - levelFLT[i]
endfor
if (mode eq 'standard') then begin
    for i=0,(numerosci-1) do begin
       Source_DCube[*,*,i] = Source_DCube[*,*,i] + meanlevelSCI - levelSCI[i]
    endfor
endif else begin
    for i=0,(numerosciA-1) do begin
       Source_A_DCube[*,*,i] = Source_A_DCube[*,*,i] + meanlevelSCIA - levelSCIA[i]
    endfor
    for i=0,(numerosciB-1) do begin
       Source_B_DCube[*,*,i] = Source_B_DCube[*,*,i] + meanlevelSCIB - levelSCIB[i]
    endfor
endelse


;-------------------------------------------------------------------------------------
;BLINKING raw images for visiual inspection first (optional)
; it prompts the option for the user
; moves to the next frame when the user clicks on the image
;-------------------------------------------------------------------------------------

read,wannablink,prompt='what do you want do to? 1=reduce spectra, 2=blink raw'

if wannablink eq 1 then begin
    goto, REDUCTION
endif
if wannablink eq 2 then begin
    goto, BLINKING
endif

BLINKING:
dummy=size(DK_frames) & numerodark=dummy[1]
dummy=size(DKFLT_frames) & numerodarkflat=dummy[1]
dummy=size(FLT_frames) & numeroflat=dummy[1]
if (mode eq 'standard') then begin
    dummy=size(SCI_frames) & numerosci=dummy[1]
endif
if (mode eq 'nodding') then begin
    dummy=size(SCI_A_frames) & numerosciA=dummy[1]
    dummy=size(SCI_B_frames) & numerosciB=dummy[1]
endif
for j = 0,numerodark-1 do begin
    window,21,xsize=512,ysize=512,title='DARK '+strtrim(string(j+1))+' of '+strtrim(string(numerodark))
    temp=rebin(DKSource_DCube[*,*,j],512,512)
    tvscl,plotsqrt(temp,99)
    cursor,xc,yc,/down,wait
endfor
for j = 0,numerodarkflat-1 do begin
    window,21,xsize=512,ysize=512,title='DARK of FLAT '+strtrim(string(j+1))+' of '+strtrim(string(numerodarkflat))
    temp=rebin(DKFlat_DCube[*,*,j],512,512)
    tvscl,plotsqrt(temp,99)
    cursor,xc,yc,/down,wait
endfor
for j = 0,numeroflat-1 do begin
    window,21,xsize=512,ysize=512,title='FLAT '+strtrim(string(j+1))+' of '+strtrim(string(numeroflat))
    temp=rebin(Flat_DCube[*,*,j],512,512)
    tvscl,plotsqrt(temp,99)
    cursor,xc,yc,/down,wait
endfor
if (mode eq 'standard') then begin
    for j=0,(numerosci-1) do begin
       window,21,xsize=512,ysize=512,title='SOURCE '+strtrim(string(j+1))+' of '+strtrim(string(numerosci))
        temp=rebin(Source_DCube[*,*,j],512,512)
        tvscl,plotsqrt(temp,99)
        cursor,xc,yc,/down,wait
    endfor
endif
if (mode eq 'nodding') then begin
    for j=0,(numerosciA-1) do begin
       window,21,xsize=512,ysize=512,title='SOURCE (A)'+strtrim(string(j+1))+' of '+strtrim(string(numerosciA))
        temp=rebin(Source_A_DCube[*,*,j],512,512)
        tvscl,plotsqrt(temp,99)
        cursor,xc,yc,/down,wait
    endfor
    for j=0,(numerosciB-1) do begin
       window,21,xsize=512,ysize=512,title='SOURCE (B)'+strtrim(string(j+1))+' of '+strtrim(string(numerosciB))
        temp=rebin(Source_B_DCube[*,*,j],512,512)
        tvscl,plotsqrt(temp,99)
        cursor,xc,yc,/down,wait
    endfor
endif

;;;;; blink images after subtracting dark
if (mode eq 'standard') then begin
    for j=0,(numerosci-1) do begin
       window,21,xsize=512,ysize=512,title='SOURCE dedarked'+strtrim(string(j+1))+' of '+strtrim(string(numerosci))
        temp=rebin(Source_DCube[*,*,j]-DKSource_DCube[*,*,0],512,512)
        tvscl,plotsqrt(temp,99)
        cursor,xc,yc,/down,wait
    endfor
endif
if (mode eq 'nodding') then begin
    for j=0,(numerosciA-1) do begin
       window,21,xsize=512,ysize=512,title='SOURCE dedarked (A)'+strtrim(string(j+1))+' of '+strtrim(string(numerosciA))
        temp=rebin(Source_A_DCube[*,*,j]-DKSource_DCube[*,*,0],512,512)
        tvscl,plotsqrt(temp,99)
        cursor,xc,yc,/down,wait
    endfor
    for j=0,(numerosciB-1) do begin
       window,21,xsize=512,ysize=512,title='SOURCE dedarked (B)'+strtrim(string(j+1))+' of '+strtrim(string(numerosciB))
        temp=rebin(Source_B_DCube[*,*,j]-DKSource_DCube[*,*,0],512,512)
        tvscl,plotsqrt(temp,99)
        cursor,xc,yc,/down,wait
    endfor
endif

;;;;;blink images after subtracting dark and correcting flat field
if (mode eq 'standard') then begin
    for j=0,(numerosci-1) do begin
       window,21,xsize=512,ysize=512,title='SOURCE dedarked and deflatted'+strtrim(string(j+1))+' of '+strtrim(string(numerosci))
        temp=rebin(Source_DCube[*,*,j]-DKSource_DCube[*,*,0]/(Flat_DCube[*,*,0]-DKFlat_DCube[*,*,0]),512,512)
        tvscl,plotsqrt(temp,99)
        cursor,xc,yc,/down,wait
    endfor
endif
if (mode eq 'nodding') then begin
    for j=0,(numerosciA-1) do begin
       window,21,xsize=512,ysize=512,title='SOURCE dedarked and deflatted (A)'+strtrim(string(j+1))+' of '+strtrim(string(numerosciA))
        temp=rebin(Source_A_DCube[*,*,j]-DKSource_DCube[*,*,0]/(Flat_DCube[*,*,0]-DKFlat_DCube[*,*,0]),512,512)
        tvscl,plotsqrt(temp,99)
        cursor,xc,yc,/down,wait
    endfor
    for j=0,(numerosciB-1) do begin
       window,21,xsize=512,ysize=512,title='SOURCE dedarked and deflatted (B)'+strtrim(string(j+1))+' of '+strtrim(string(numerosciB))
        temp=rebin(Source_B_DCube[*,*,j]-DKSource_DCube[*,*,0]/(Flat_DCube[*,*,0]-DKFlat_DCube[*,*,0]),512,512)
        tvscl,plotsqrt(temp,99)
        cursor,xc,yc,/down,wait
    endfor
endif

;-------------------------------------------------------------------------------------
;PREPARING FLATS and removing badpixels from flats
;-------------------------------------------------------------------------------------
REDUCTION:

medcomb,DKFlat_DCube,masterdarkflat  ;stacking together the cube of darks with a median filter
medcomb,Flat_DCube,masterflat
delvarx,DKFlat_DCube,Flat_DCube ;dealloc

;sigmarhp=8
r_in=2
r_out=5
masterflat_c=sigma_filter(masterflat,4)
masterdarkflat_c=sigma_filter(masterdarkflat,4)
delvarx,masterflat,masterdarkflat ;dealloc

flat0_c = masterflat_c-(masterdarkflat_c)
im=reverse(rotate(flat0_c,3),2)
flat0_cr=rot(im,1.21,/interp)   ;the spectra on the raw irmos data are rotated a bit from the xy axes of the detector

;-------------------------------------------------------------------------------------
;PREPARING DARKS (if we are in standard mode, nodding has no darks)
;-------------------------------------------------------------------------------------

    medcomb,DKSource_DCube,dark   ;stacking together darks
    delvarx,DKSource_DCube
    ;irmos_rhp,dark,dark0,sigmarhp,r_in,r_out
    dark0=sigma_filter(dark,4)
    delvarx,dark ;dealloc


;------------------------------------------------------------------------------------
;--------------------    FINDING SPECTRA  -------------------------------------------
; basically, rotated flats (where the spectra are aligned vertically) are used to
; select by eye on the image where each slit appears.
; It could be done automatically, but I didn't do it
;------------------------------------------------------------------------------------

modefinder = 2 ;forcing manual mode
if (modefinder eq 1) then begin
    threshold=20
    TARGET_FINDER_SIMPLE,flat0_cr,n_targets,threshold,x0_target,X1_target
endif
if (modefinder eq 2) then begin
    window,11,xsize=1024,ysize=200,ypos=200
    dummy=fltarr(1024)
    for i=0,1023,1 do begin
       dummy[i]=total(flat0_cr[i,*])
    endfor
    dummy2=fltarr(1024,200)
    for i=0,199,1 do begin
        dummy2[*,i]=dummy[*]
    endfor
    tvscl,dummy2
    if (OBJTYPE eq 'A0') then begin ;in this case we have only one slit
        howmanyslits = 1
        print,'point beginning and end of the slit'
    endif else begin
        read,howmanyslits,prompt='how many slits do you see?'
        print,'now for every slit point beginning and end'
    endelse
    x0_target=fltarr(howmanyslits)
    x1_target=fltarr(howmanyslits)
    for i=0,howmanyslits-1,1 do begin
       cursor,xc,yc,/device,/down
       x0_target[i]=xc
       print,i,' begin ',xc
       wait,0.2
       cursor,xc,yc,/device,/down
       x1_target[i]=xc
       print,i,' end ',xc
       wait,0.2
    endfor
    n_targets=howmanyslits
    wdelete,11
endif

print,x0_target
print,x1_target

;------------------------------------------------------------------------------------
; CREATING FLAT MASK (where in the image the slits are located
;------------------------------------------------------------------------------------

maskslittrak=fltarr(1024,1024)
for i=0,(n_targets-1) do begin
    maskslittrak[x0_target[i]:x1_target[i],*]=1
endfor
maskslittrak=rot(transpose(maskslittrak),1.21,/interp)


stop
;------------------------------------------------------------------------------------
; REMOVING OFFSET FROM FLATS
;------------------------------------------------------------------------------------


a=masterflat_c-(masterdarkflat_c)
sbatti_a_zero,a,flat0_c,maskslittrak
im=reverse(rotate(flat0_c,3),2)
flat0_cr=rot(im,1.21,/interp)

;------------------------------------------------------------------------------------
;COMPUTING FLAT FIELD NORMALIZATION
;------------------------------------------------------------------------------------

NormFF=fltarr(n_targets)
for ix = 0,N_targets-1 do begin
    NormFF[ix] = MEDIAN(flat0_cr[X0_Target[ix]:X1_Target[ix],*])
    NormFF[ix] = NormFF[ix]/max(NormFF[ix])
endfor

;-------------------------------------------------------------------------------------
;REMOVING BAXPIX from SCI image, FLATFIELD, DEDARK, and so on
;-------------------------------------------------------------------------------------

if (strpos(DATA_PATH,'Run5') gt -1) then begin
    shifts=[-3,-2,-1,0,1,2,3]
    if (mode eq 'standard') then begin
       for i=0,NSource_Frames-1 do begin
         sci_nohp=sigma_filter(Source_DCube[*,*,i],5)
         sci=(sci_nohp-dark0)/(flat0_c>0.0001) ;dedark and deflat cube
         maskcod_filter,sci,mask
         im = reverse(rotate(sci,3),2)
           Source_DCube[*,*,i] = rot(im,1.21,/interp) ;rotate whole cube
       endfor
       delvarx,im,sci,sci_nohp  ;dealloc
        v=median(Source_DCube,DIMENSION=2)
        shiftOKvec=fltarr(NSource_Frames)
        for ix=0,NSource_Frames-1 do begin
            result=c_correlate(v[300:700,0],v[300:700,ix],shifts)
            shiftOK=where(result EQ MAX(result))-3
            print,ix,result,shiftOK
            shiftOKvec[ix]=shiftOK
        endfor
        for ix=0,(Nsource_frames-1) do begin
            im=shift(Source_Dcube[*,*,ix],-shiftOKvec[ix],0)
            Source_Dcube[*,*,ix]=im
       endfor
        medcomb,Source_DCube,SPECTRA
        Source_DCube=0 ;dealloc
    endif else begin ;nodding case
        Source_A_DCube_Dedark=fltarr(1024,1024,NSource_A_Frames)
        for i=1,NSource_A_Frames-1 do begin ;removing hotpix
         sci_nohp=sigma_filter(Source_A_DCube[*,*,i],5)
         sci=sci_nohp/(flat0_c>0.0001)
         scidedark=(sci_nohp-dark0)/(flat0_c>0.0001)
         maskcod_filter,sci,mask
         maskcod_filter,scidedark,mask
         im = reverse(rotate(sci,3),2)
           Source_A_DCube[*,*,i] = rot(im,1.21,/interp) ;rotate whole cube
           im = reverse(rotate(scidedark,3),2)
           Source_A_DCube_Dedark[*,*,i] = rot(im,1.21,/interp)
       endfor
       for i=1,NSource_B_Frames-1 do begin ;removing hotpix
         sci_nohp=sigma_filter(Source_B_DCube[*,*,i],5)
         sci=sci_nohp/(flat0_c>0.0001)
         maskcod_filter,sci,mask
         im = reverse(rotate(sci,3),2)
           Source_B_DCube[*,*,i] = rot(im,1.21,/interp) ;rotate whole cube
       endfor
       delvarx,im,sci,sci_nohp ;dealloc
       vA=median(Source_A_DCube,DIMENSION=2)
       vB=median(Source_B_DCube,DIMENSION=2)
       shiftOKvecA=fltarr(NSource_A_Frames)
       shiftOKvecB=fltarr(NSource_B_Frames)
       for ix=0,NSource_A_Frames-1 do begin
            result=c_correlate(vA[300:700,0],vA[300:700,ix],shifts)
            shiftOK=where(result EQ MAX(result))-3
            print,ix,result,shiftOK
            shiftOKvecA[ix]=shiftOK
        endfor
        for ix=0,(Nsource_A_frames-1) do begin
            im=shift(Source_A_Dcube[*,*,ix],-shiftOKvecA[ix],0)
            Source_A_Dcube[*,*,ix]=im
       endfor
       ;
        for ix=0,NSource_B_Frames-1 do begin
            result=c_correlate(vB[300:700,0],vB[300:700,ix],shifts)
            shiftOK=where(result EQ MAX(result))-3
            print,ix,result,shiftOK
            shiftOKvecB[ix]=shiftOK
        endfor
        for ix=0,(Nsource_B_frames-1) do begin
            im=shift(Source_B_Dcube[*,*,ix],-shiftOKvecB[ix],0)
            Source_B_Dcube[*,*,ix]=im
       endfor
        medcomb,Source_A_DCube_Dedark,SPECTRA
        medcomb,Source_A_DCube,SPECTRAA
        medcomb,Source_B_DCube,SPECTRAB
        SPECTRA_difference_rot=SPECTRAA-SPECTRAB
       delvarx,Source_A_DCube_Dedark,Source_A_DCube,Source_B_DCube ;dealloc
    endelse

endif else begin ;i.e. if run is not 6

    if (mode eq 'standard') then begin
       Source_Dcube_c=fltarr(1024,1024,nsource_frames)
       for i=0,(nsource_frames-1) do begin
         Source_Dcube_c[*,*,i]=sigma_filter(Source_Dcube[*,*,i],4)
       endfor
       delvarx,Source_Dcube
       medcomb,Source_Dcube_c,sci_nohp
       delvarx,Source_Dcube_c
       a=(sci_nohp-dark0)
       maskcod_filter,a,mask
       sbatti_a_zero,a,sci0,maskslittrak
       sci=(sci0)/(flat0_c>1.)
       im = reverse(rotate(sci,3),2)
        SPECTRA = rot(im,1.21,/interp)
        delvarx,im,sci,a,sci0,sci_nohp
    endif else begin
       Source_A_Dcube_c=fltarr(1024,1024,nsource_A_frames)
       for i=0,(nsource_A_frames-1) do begin
         Source_A_Dcube_c[*,*,i]=sigma_filter(Source_A_Dcube[*,*,i],4)
       endfor
       delvarx,Source_A_Dcube
       Source_B_Dcube_c=fltarr(1024,1024,nsource_B_frames)
       for i=0,(nsource_B_frames-1) do begin
         Source_B_Dcube_c[*,*,i]=sigma_filter(Source_B_Dcube[*,*,i],4)
       endfor
       delvarx,Source_B_Dcube
       medcomb,Source_A_Dcube_c,sciA_nohp
       medcomb,Source_B_Dcube_c,sciB_nohp
       delvarx,Source_A_Dcube_c,Source_B_Dcube_c ;dealloc
         a=(sciA_nohp-dark0)
         maskcod_filter,a,mask
         sbatti_a_zero,a,out,maskslittrak
         sciA_dd=(out)/(flat0_c>1.)
         a=(sciB_nohp-dark0)
         maskcod_filter,a,mask
         sbatti_a_zero,a,out,maskslittrak
         sciB_dd=(out)/(flat0_c>1.)
         delvarx,a,out
       SPECTRA_difference=sciA_nohp-sciB_nohp
       maskcod_filter,SPECTRA_difference,mask
       sbatti_a_zero,SPECTRA_difference,im,maskslittrak
       SPECTRA_difference=im/(flat0_c>1)

       im = reverse(rotate(SPECTRA_difference,3),2)
       SPECTRA_difference_rot = rot(im,1.21,/interp)
       im = reverse(rotate(sciA_dd,3),2)
       SPECTRA=rot(im,1.21,/interp)
       delvarx,im
    endelse
endelse

;------------------------------------------------------------------------------------
; CREATING SMALLER IMAGE WITH JOINED SPECTRA and NORMALIZING
;------------------------------------------------------------------------------------

Xsize=X1_target-X0_Target+1
xsizetotal=total(xsize)
spectrajoined=fltarr(xsizetotal,1024)
if (mode eq 'nodding') then begin
    spectrajoined_diff=fltarr(xsizetotal,1024)
endif
newx0=fltarr(n_targets)
newx1=fltarr(n_targets)
bottom=0
beg=0
for is=0,n_targets-1 do begin
    spectrajoined[beg:beg+xsize[is]-1,*]=SPECTRA[x0_target[is]:x1_target[is],*];*NormFF[is]
    if (mode eq 'nodding') then begin
        spectrajoined_diff[beg:beg+xsize[is]-1,*]=SPECTRA_difference_rot[x0_target[is]:x1_target[is],*];*NormFF[is]
    endif
    beg=beg+xsize[is]
    newx0[is]=bottom
    newx1[is]=newx0[is]+xsize[is]-1
    bottom=newx1[is]+1
endfor
print,x0_target,x1_target,xsize

;filenameout = PREPATH+DATA_PATH+'\red\'+Source_name+'join.fits'
;writefits,filenameout,spectrajoined

;------------------------------------------------------------------------------------
; WARPING to straighten oblique sky lines
;------------------------------------------------------------------------------------
for is=0,n_targets-1 do begin
    prova=rebin(spectrajoined[newx0[is]:newx1[is],*],xsize[is]*4,1024*4)
    dummy=poly_2d(prova,[0,0,1,0],[0,1,0.0211,0])
    prova=rebin(dummy,xsize[is],1024)
    spectrajoined[newx0[is]:newx1[is],*]=prova[*,*]
    if (mode eq 'nodding') then begin
        prova=rebin(spectrajoined_diff[newx0[is]:newx1[is],*],xsize[is]*4,1024*4)
        dummy=poly_2d(prova,[0,0,1,0],[0,1,0.0211,0])
        prova=rebin(dummy,xsize[is],1024)
        spectrajoined_diff[newx0[is]:newx1[is],*]=prova[*,*]
    endif
endfor


;####################################################################################
; NOW BEGINNING WAVELENGHT CALIBRATION
;####################################################################################

;------------------------------------------------------------------------------------
; SELECTING SKY SPECTRUM
;------------------------------------------------------------------------------------

window,12,xpos=0,ypos=0,ysize=1024,xsize=xsizetotal,title='Select a good spectrum'
spectrajoined_show=sigrange(sqrt(spectrajoined>0.))
resistant_mean,spectrajoined_show,3,media
sigma=robust_sigma(spectrajoined_show)
minimo=(media-sigma)
massimo=(media+sigma)
minimoins=minimo
massimoins=massimo
for something=0,100 do begin
    TVSCL,spectrajoined_show>(minimo)<(massimo)
    print,'actual range: ',minimo,massimo
        read,minimoins,massimoins,prompt='insert range min,max for displaying the spectra on screen (0,0 if you are ok with the current display): '
    if ((minimoins eq 0.0) and (massimoins eq 0.0)) then begin
       break
    endif else begin
       minimo=minimoins
       massimo=massimoins
    endelse
endfor
if (n_targets eq 1) then begin
    selectedspectra=0 ;this means that I have only one spectrum, nothing to choose..
endif else begin
    print,'click to select reference spectrum. The others will be shifted to the same walenght scale.'
    cursor,xc,yc,/device
    for is=0,n_targets-1 do begin
       if (xc gt newx0[is] and xc lt newx1[is]) then selectedspectra=is
    endfor
endelse

;------------------------------------------------------------------------------------
; TRYING TO FIND AUTOMATIC SHIFT SOLUTION
;------------------------------------------------------------------------------------

fondo=fltarr(n_targets,1024)
vectorshift = fltarr(n_targets)
matrixout=fltarr(1024,3,n_targets)
shiftedspectra=fltarr(xsizetotal,2048)

for is=0,n_targets-1 do begin
    stripe=fltarr(xsize[is],1)
    for pix=0,1023 do begin
       stripe=spectrajoined[newx0[is]:newx1[is],pix]
       stripe2=stripe[sort(stripe)]
       quaranta=round(0.32*xsize[is])
       fondo[is,pix]=stripe2[quaranta]
    endfor
    d2spec,fondo[is,*],out,3,4
    out[0:100]=0
    out[923:1023]=0
    matrixout[*,1,is]=out
endfor

riferimento=matrixout[*,*,selectedspectra]

for i=0,n_targets-1 do begin
    imgtemp=matrixout[*,*,i]
    correl_optimize,riferimento,imgtemp,xa,ya
    vectorshift[i]=xa
    vectorshift[selectedspectra]=0
    if (flagoldvectorshift eq 1) then vectorshift[i] = oldvectorshift[i]
    shiftedspectra[newx0[i]:newx1[i],(512+vectorshift[i]):(1535+vectorshift[i])]=spectrajoined[newx0[i]:newx1[i],*]
endfor

wdelete,12
window,7,xpos=0,ypos=0,xsize=2048,ysize=xsizetotal,title='Shifted Spectra'
TVSCL,rotate(sigrange(sqrt(shiftedspectra>0.)),3)>(minimo)<(massimo)

;------------------------------------------------------------------------------------
; BY HAND CORRECTION
;------------------------------------------------------------------------------------

if (n_targets gt 1) then begin ;otherwise 1 single spectrum, and of course it IS ok..
read,isitok,prompt='is this ok, or do you want to correct the alignment by hand? 1=ok 2=no'
if (isitok eq 2) then begin
    read,howmany,prompt='how many spectra do you want to correct by hand?'
    for k=1,howmany do begin
        cursor,xc,yc,/device
        xc2=xsizetotal-yc
        for is=0,n_targets-1 do begin
          if (xc2 gt newx0[is] and xc2 lt newx1[is]) then guilty=is
        endfor
        print,'guilty spectra that we are correcting: ',guilty
        for ss=1,100 do begin
           print,'actual shift',vectorshift[guilty]
           read,newshift,prompt='insert newshift: (-1 will exit loop, keeping last value)'
           if (newshift eq -1) then begin
              print,'loop broken with newshift=',vectorshift[guilty],'. Next spectra'
              break
            endif
           vectorshift[guilty]=newshift
           shiftedspectra[newx0[guilty]:newx1[guilty],*]=0
           shiftedspectra[newx0[guilty]:newx1[guilty],(512+vectorshift[guilty]):(1535+vectorshift[guilty])]= $
                                          spectrajoined[newx0[guilty]:newx1[guilty],*]
          TVSCL,rotate(sigrange(sqrt(shiftedspectra>0.)),3)>(minimo)<(massimo)
        endfor
    endfor
endif
endif
filenameout = PREPATH+DATA_PATH+'\red\'+Source_name+'shifted.fits'
writefits,filenameout,shiftedspectra

print,'vectorshift= ',vectorshift

;------------------------------------------------------------------------------------
; WAVELENGHT CALIBRATION,
; this starts with a first order calibration done by eye (a catalog of sky lines is
; overplotted to the spectra, and the user fiddles with dlambda/dx and lambda0 until
; there is a decent overlap. Then the code refines the solution automatically
;------------------------------------------------------------------------------------

if (ohmethod eq 'oliva') then begin
    openr,/get_lun,u,(PREPATH+'OH_lines.dat')
    dummy = fltarr(500)
    k = 0
    while not eof(u) do begin
        readf,u,t
        dummy(k) = t
        k = k+1
    endwhile
    close,u
    OHlinesTOT=(dummy[where(dummy gt 0)])*1000
endif
if (ohmethod eq'rousselot') then begin
    readcol,PREPATH+'rousselot2000.dat',l,f
    l=l/10.
    good = where(f gt minOHflux)
    ohlinesTOT=l[good]
endif

lambdazero=2382 & disp=-0.285

CASE filter OF
   'Z':     BEGIN lambda_inf = 850 & lambda_sup = 1150
            END
   'J':     BEGIN lambda_inf = 1132 & lambda_sup = 1345 & lambdazero=1610 & disp=-0.3847
            END
   'Jblue': BEGIN lambda_inf = 1120 & lambda_sup = 1265
            END
   'Jred':  BEGIN lambda_inf = 1210 & lambda_sup = 1355
            END
   'H':     BEGIN lambda_inf = 1430 & lambda_sup = 1800 & lambdazero=2218 & disp=-0.612
            END
   'Hblue': BEGIN lambda_inf = 1440 & lambda_sup = 1625
            END
   'Hred':  BEGIN lambda_inf = 1585 & lambda_sup = 1835
            END
   'K':     BEGIN lambda_inf = 1900 & lambda_sup = 2460 & lambdazero=2382 & disp=-0.906
            END
   'Kblue': BEGIN lambda_inf = 1900 & lambda_sup = 2210 & lambdazero=2382 & disp=-0.285
            END
   'Kred':  BEGIN lambda_inf = 2100 & lambda_sup = 2460
            END
ENDCASE

bestlambdazero=lambdazero
bestdisp=disp

OHlines=OHlinesTOT[where((OHlinesTOT gt lambda_inf) and (OHlinesTOT lt lambda_sup))]


dummy=size(OHlines)
numberoflines=dummy[1]

print,'now you have to change dlambda/dx and lambda0 until a good match is found between the sky lines in the observed spectrum, and the catalog ones (vertical black lines)

for something=0,180 do begin
    xdilambda=(OHlines-bestlambdazero)/bestdisp
    TVSCL,rotate(sigrange(sqrt(shiftedspectra>0.)),3)>(minimo)<(massimo)
    for line=0,numberoflines-1 do begin
       arrow,xdilambda[line],0,xdilambda[line],1000,color=0
    endfor
       ;paschenbeta
       dummy=x_of_lambda(1290.2,bestdisp,bestlambdazero)
       arrow,dummy,0,dummy,1000,color=0,thick=3
    print,'current dlambda/dx and lambdazero:',disp,',',lambdazero
    bestlambdazero=lambdazero
    bestdisp=disp
    read,prompt='dlambda/dx,lambdazero (exit 0,0) ',disp,lambdazero
    if ((abs(disp) gt 0) and (abs(lambdazero) gt 0)) then begin
        bestdisp=disp
        bestlambdazero=lambdazero
    endif else begin
        break
    endelse
endfor

lambdazero=bestlambdazero
disp=bestdisp

fondo=fltarr(n_targets,2048)

bad_from_here=round(x_of_lambda(lambda_inf,bestdisp,bestlambdazero))<2047
bad_to_here=round(x_of_lambda(lambda_sup,bestdisp,bestlambdazero))>0

for is=0,n_targets-1 do begin & $
    shiftedspectra[newx0[is]:newx1[is],0:bad_to_here]=0
    shiftedspectra[newx0[is]:newx1[is],bad_from_here:2047]=0
    for pixel=0,2047 do begin & $
        stripe=shiftedspectra[newx0[is]:newx1[is],pixel] & $
        stripe2=stripe[sort(stripe)] & $
        quaranta=round(0.32*xsize[is]) & $
        fondo[is,pixel]=stripe[quaranta] & $
    endfor
endfor

TVSCL,rotate(sigrange(sqrt(shiftedspectra>0.)),3)>(minimo)<(massimo)

epsilon=7
read,epsilon,prompt='radius for automatic refinment of the solution? [pixel]: '
ohlinesvec=fltarr(n_targets,numberoflines)

for is=0,n_targets-1 do begin
    for riga=0,numberoflines-1 do begin
       in_inter=xdilambda[riga]-epsilon
       out_inter=xdilambda[riga]+epsilon
       print,'spettro, riga, inizio e fine intervallo', is, riga,in_inter,out_inter
       dummy=(where(fondo[is,in_inter:out_inter] EQ MAX(fondo[is,in_inter:out_inter]))+in_inter)
       print,size(dummy)
       dummy2=size(dummy)
       if (dummy2[1] gt 1) then begin
            dummy=0
       endif
       ohlinesvec[is,riga]=dummy
       plots,ohlinesvec[is,riga],xsizetotal-((newx0[is]+newx1[is])/2),psym=5,color=0,thick=3,/device
       plots,ohlinesvec[is,riga],xsizetotal-((newx0[is]+newx1[is])/2),psym=5,color=255,thick=2,/device
    endfor
    arrow,0,(xsizetotal-newx1[is]),4500,(xsizetotal-newx1[is]),color=255
endfor

dispvec=fltarr(n_targets)
lambdazerovec=fltarr(n_targets)

window,13,xsize=512,ysize=512
for is=0,n_targets-1 do begin
    pixels=ohlinesvec[is,*]
    lambdas=ohlines
         dummy2=where(pixels)
         pixels=pixels(dummy2)
         lambdas=lambdas(dummy2)
         plot,pixels,lambdas,psym=4
    dummy=robust_linefit(pixels,lambdas,lambdascorr)
    dispvec[is]=dummy[1]
    lambdazerovec[is]=dummy[0]
    oplot,pixels,lambdascorr,color=50
    wait,0.2
endfor
;wdelete,13
print,lambdazerovec,dispvec

;------------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;       PREVIOUS           ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;------------------------------------------------------------------------------------

;------------------------------------------------------------------------------------
; EXTRACTING SPECTRUM (STANDARD CASE)
;------------------------------------------------------------------------------------
if (mode eq 'standard') then begin
    bkg_spectrum=fltarr(2048,n_targets)
    bkg_spectrumA=fltarr(2048,n_targets)
    bkg_spectrumB=fltarr(2048,n_targets)
    target_spectrum=fltarr(2048,n_targets)
    lambdas_def=fltarr(2048,n_targets)
    iniziovec=fltarr(n_targets)
    finevec=fltarr(n_targets)
    pixel_inizio=fltarr(n_targets)
    pixel_fine=fltarr(n_targets)
    for is=0,n_targets-1 do begin
        proj=fltarr(xsize[is])
        dovevabene=total(shiftedspectra[newx0[is]:newx1[is],*])
        vabene=where(dovevabene gt 0)
        for pix=newx0[is],newx1[is] do begin
         indicepixel=pix-newx0[is]
         temp=(shiftedspectra[pix,*])
         temp_nezero=temp(where(temp ne 0))
         resistant_mean,temp_nezero,6,temp2
         proj[indicepixel]=temp2
       endfor
       maxSNR,proj,inizioskya,fineskya,inizio,fine,inizioskyb,fineskyb,1
       inizio=inizio+newx0[is]
       fine=fine+newx0[is]
       inizioskyA=inizioskyA+newx0[is]
       fineskyA=fineskyA+newx0[is]
       inizioskyB=inizioskyB+newx0[is]
       fineskyB=fineskyB+newx0[is]
       print,is,'these are the intervals',inizioskyA,fineskyA,inizio,fine,inizioskyB,fineskyB
       arrow,0,(xsizetotal-inizio),3000,(xsizetotal-inizio),color=0
       arrow,0,(xsizetotal-fine),3000,(xsizetotal-fine),color=0
       iniziovec[is]=inizio
       finevec[is]=fine
       ;
       ;BACKGROUND SUBTRACTION COMPUTATION
       ;
       width_target=fine-inizio+1
       width_bkgA=(fineskyA-inizioskyA+1)+(fineskyB-inizioskyB+1)
       for pixel=0,2047 do begin
         target_spectrum[pixel,is]=total(shiftedspectra[inizio:fine,pixel])/width_target
       endfor
       ;-------------------------------------
       if (inizioskyA lt fineskyA) then width_bkgA=(fineskyA-inizioskyA+1) else width_bkgA=0
       if (inizioskyB lt fineskyB) then width_bkgB=(fineskyB-inizioskyB+1) else width_bkgB=0
       width_bkg=(width_bkgA+width_bkgB)
       dummy=fltarr(width_bkg,2048)
       if (width_bkgA gt 0 and width_bkgB gt 0) then begin
              dummy[0:(width_bkgA-1),*]=shiftedspectra[inizioskyA:fineskyA,*]
              dummy[width_bkgA:(width_bkg-1),*]=shiftedspectra[inizioskyB:fineskyB,*]
       endif
       if (width_bkgA gt 0 and width_bkgB eq 0) then begin
              dummy=shiftedspectra[inizioskyA:fineskyA,*]
       endif
       if (width_bkgA eq 0 and width_bkgB gt 0) then begin
              dummy=shiftedspectra[inizioskyB:fineskyB,*]
       endif
       ;---------------------------------------------------------------
       for pixel=0,2047 do begin
            bkg_spectrum[pixel,is]=median(dummy[*,pixel])
       endfor
    endfor
    target_spectrum_2=(target_spectrum-bkg_spectrum)
     window,2,xsize=1100,ysize=800,title='spectrum'+string(is),xpos=(0+shiftscreen),ypos=0
    for is=0,n_targets-1 do begin
       for pixel=0,2047 do begin
         lambdas_def[pixel,is]=(dispvec[is]*float(pixel)+lambdazerovec[is])
       endfor
       pixel_inizio[is]=512+vectorshift[is]
       pixel_fine[is]=1536+vectorshift[is]
       plot,lambdas_def[pixel_inizio[is]:pixel_fine[is],is],target_spectrum[pixel_inizio[is]:pixel_fine[is],is],color=140
       oplot,lambdas_def[pixel_inizio[is]:pixel_fine[is],is],target_spectrum_2[pixel_inizio[is]:pixel_fine[is],is]
       bbb=lambdas_def[(1024+vectorshift[is]),is]
       arrow,bbb,0,bbb,10,/thick,/data,color=128
       cursor,pippo,pluto,/down,wait
    endfor
ENDIF

;------------------------------------------------------------------------------------
; EXTRACTING SPECTRUM (NODDING CASE)
;------------------------------------------------------------------------------------


if (mode eq 'nodding') then begin
    shiftedspectra_diff=fltarr(xsizetotal,2048)
    pixel_inizio=fltarr(n_targets)
    pixel_fine=fltarr(n_targets)
    for i=0,n_targets-1,1 do begin
        shiftedspectra_diff[newx0[i]:newx1[i],(512+vectorshift[i]):(1535+vectorshift[i])]=spectrajoined_diff[newx0[i]:newx1[i],*]
        shiftedspectra_diff[newx0[i]:newx1[i],0:bad_to_here]=0
        shiftedspectra_diff[newx0[i]:newx1[i],bad_from_here:2047]=0
        pixel_inizio[i]=512+vectorshift[i]
        pixel_fine[i]=1536+vectorshift[i]
    endfor
    resistant_mean,shiftedspectra_diff[*,700:1400],3,media2
    sigma2=robust_sigma(shiftedspectra_diff[*,700:1400])
        dummy=abs(rotate(shiftedspectra_diff,3))
        iniziovero=max([pixel_inizio[0],bad_to_here])
        finevero=min([pixel_fine[0],bad_from_here])
        spectrajoined_diff_trim=dummy[iniziovero:finevero,*]
        window,1,xsize=(finevero-iniziovero+1),ysize=xsizetotal,ypos=400,title='spectra'
        tvscl,spectrajoined_diff_trim<(median(spectrajoined_diff_trim)+stdev(spectrajoined_diff_trim))> $
                             (median(spectrajoined_diff_trim)-stdev(spectrajoined_diff_trim))
        target_spectrum_2=fltarr(2048,n_targets)
        target_spectrum=fltarr(2048,n_targets)
        stop
        lambdas_def=fltarr(2048,n_targets)
    for is=0,n_targets-1 do begin
       coaddedslit=fltarr(xsize[is],2048)
       coaddedslit[0:xsize[is]-1,*]=(shiftedspectra_diff[newx0[is]:newx1[is],*])
       proj2=fltarr(xsize[is])
       for riga=0,xsize[is]-1 do begin
            pippo=coaddedslit[riga,*]
            pippo2=pippo(where(pippo ne 0))
            ;resistant_mean,pippo2,3,mediapippo
         mediapippo=median(pippo2)
            proj2[riga]=mediapippo   ;this is a projection of the coadded nodded spectrum along the width of the slit
       endfor
       window,7,xsize=640,ysize=640,xpos=640
       plot,proj2
       arrow,0,0,40,0,/thick,/data,color=140
       ;cursor,xcc,ycc,/data,/down
       ;wait,0.2
       ;arrow,0,ycc,40,ycc,/thick,/data
       ycc=0.
       coaddedslit=coaddedslit-ycc
       for riga=0,xsize[is]-1 do begin
            pippo=(coaddedslit[riga,*])
            pippo2=pippo(where(pippo ne -ycc))
            ;resistant_mean,pippo2,3,mediapippo
         mediapippo=median(pippo2)
            proj2[riga]=abs(mediapippo)
       endfor
       print,'now click on the projected flux along the slit width 4 times:
       print,'to mark the start, end, start, end of the 2 broad peaks
       print,'the flux inside these ranges will be used to extract the spectrum'
       window,3,xsize=512,ysize=512,title='spectrum'+string(is)
       plot,proj2
        cursor,x1,dummy,/down,/DATA,wait & x1=max([round(x1),0])
        arrow,x1,0,x1,max(proj2),/data,color=50
        cursor,x2,dummy,/down,/DATA,wait & x2=round(x2)
        arrow,x2,0,x2,max(proj2),/data,color=150
        cursor,x3,dummy,/down,/DATA,wait & x3=round(x3)
        arrow,x3,0,x3,max(proj2),/data,color=50
        cursor,x4,dummy,/down,/DATA,wait & x4=min([round(x4),(xsize[is]-1)])
        arrow,x4,0,x4,max(proj2),/data,color=150
        wait,0.2
       for pix=0,2047 do begin
            lambdas_def[pix,is]=(dispvec[is]*float(pix)+lambdazerovec[is])
            target_spectrum_2[pix,is]=total(abs(coaddedslit[x1:x2,pix]))+ $
                                   total(abs(coaddedslit[x3:x4,pix]))
       endfor
       target_spectrum_2=target_spectrum_2/(x2+x4-x1-x3+2)
       target_spectrum=target_spectrum_2
       window,2,xsize=1100,ysize=800,title='spectrum'+string(is),xpos=(0+shiftscreen),ypos=0
       plot,lambdas_def[iniziovero:finevero,is],target_spectrum[iniziovero:finevero,is],color=255
       bbb=lambdas_def[(1024+vectorshift[is]),is]
       arrow,bbb,0,bbb,10,/thick,/data,color=128
       cursor,pippo,pluto,/down,wait
    endfor
    filenameout=PREPATH+DATA_PATH+'\red\'+Source_name+'shiftedspectradiff.fits'
    writefits,filenameout,shiftedspectra_diff
endif

;------------------------------------------------------------------------------------
; WRITING RESULTS
;------------------------------------------------------------------------------------

while !d.window ne -1 do wdelete,!d.window
pixels=indgen(2048)
read,doihavetowrite,prompt='do I have to write data on disk? (1=yes)'
if doihavetowrite eq 1 then begin
    for is=0,n_targets-1 do begin
       filenameout=PREPATH+DATA_PATH+'\red\'+Source_name+'spectra-source'+STRTRIM(STRING(is),2)+'.dat'
       lunghezza=pixel_fine[is]-pixel_inizio[is]+1
       openw,1,filenameout
       for pixel=0,2047 do begin
         printf,1,pixels[pixel],lambdas_def[(2047-pixel),is],' ',target_spectrum_2[(2047-pixel),is],' ',target_spectrum[(2047-pixel),is]
       endfor
       close,1
       filenameout=PREPATH+DATA_PATH+'\red\'+Source_name+'additional-data'+STRTRIM(STRING(is),2)+'.txt'
       openw,2,filenameout
         printf,2,n_targets
         printf,2,x0_target[is],' #x0_target'
         printf,2,x1_target[is],' #x1_target'
         printf,2,newx0[is],' #newx0'
         printf,2,newx1[is],' #newx1'
         ;printf,2,iniziovec[is],' #start point of target spectrum'
         ;printf,2,finevec[is],' #start point of target spectrum'
         printf,2,dispvec[is],' #dispersion dl/dx'
         printf,2,lambdazerovec[is],' #lambda of pixel 0'
         printf,2,pixel_inizio[is],' #first useful pixel'
         printf,2,pixel_fine[is],' #last useful pixel'
       close,2

    endfor
endif

stop


end







