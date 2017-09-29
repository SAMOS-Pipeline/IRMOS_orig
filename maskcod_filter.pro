pro maskcod_filter,imagein,mask
	badpix=mask[*,100]
	dummy=size(badpix) & n=dummy[1]
	for i=0,(n-1) do begin
		good2=mask[i,0:99]
		good=good2[where(good2 gt 0)]
		imagein[badpix[i]]=median(imagein[good])
	endfor
end