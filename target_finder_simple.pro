PRO TARGET_FINDER_simple,flat,n_targets,threshold,x0_target,X1_target
Vtargets=fltarr(1024)
for i=0,1023 do begin
    Vtargets[i]=total(flat[i,*])/1024
    if (Vtargets[i] lt 0) then begin
        Vtargets[i] = 0
    endif
endfor
Vtargets=smooth(Vtargets,3)
;FINDER=(Vtargets-SHIFT(VTARGETS,1))
;FINDER2=smooth((FINDER-SHIFT(FINDER,-1)),3)

;n_targets=15
x0_target_ext=intarr(100)
x1_target_ext=intarr(100)

window,14,xsize=768,ysize=768
plot,Vtargets
read,threshold,prompt='INSERT THRESHOLD FOR TRACK DETECTION:'

arrow,0,threshold,4000,threshold,color=140,/data



s=0
flag=0



FOR I=10,1013 DO BEGIN
    print,i
    IF ((Vtargets[i]) GT threshold AND (Vtargets[i-1]) LT threshold) THEN BEGIN
        IF FLAG EQ 0 THEN BEGIN
          x0_target_EXT[s]=i
          flag=1
        ENDIF ELSE flag=0
    ENDIF
    IF ((Vtargets[i]) LT threshold AND (Vtargets[i-1]) GT threshold) THEN BEGIN
        IF Flag EQ 0 THEN BEGIN
          x1_target_EXT[s]=i
          s=s+1
        ENDIF
    ENDIF
;    print,i,s,flag,FINDER[i-1],FINDER[i],x0_target[(s-1)>0],x1_target[(s-1)>0]
ENDFOR
x1_target=fltarr(s)
x0_target=fltarr(s)
for i=0,s-1 do begin
    x0_target[i]=x0_target_EXT[i]
    x1_target[i]=x1_target_EXT[i]
endfor

n_targets=s

stop
end