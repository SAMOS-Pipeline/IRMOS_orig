pro medcomb,input,output

sizein = size(input)

if (sizein[0] eq 0 ) then begin
    print,'WARNING: input array not vadid'
endif else begin

xdim = sizein[1]
ydim = sizein[2]
zdim = sizein[3]

output=input[*,*,0]

for i=0,(xdim-1),1 do begin
    for j=0, (ydim-1),1 do begin
    vectortemp=input[i,j,*]
    ;output[i,j]=median(vectortemp,/EVEN)
    output[i,j]=mean(vectortemp)
    endfor
endfor

endelse


end