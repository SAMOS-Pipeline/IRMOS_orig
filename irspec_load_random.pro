PRO IRSPEC_load_random,PATH,SOURCE_NAME,frames,name_DCube
dummy=size(frames)
n_frames=dummy[1]
name_DCube=fltarr(1024,1024,n_frames)
for i=0,n_frames-1 do begin
   fits_read,PATH+SOURCE_NAME+STRTRIM(STRING(frames[i]),2)+'.fit',ima
   name_DCube(*,*,i)=ima
   print,'reading file  '+PATH+SOURCE_NAME+STRTRIM(STRING(frames[i]),2)+'.fit'
endfor
end