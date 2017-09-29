PRO IRMOS_normchips_x,image,stripe,gap

;--------
; A | B |
;--------
; C | D |
;--------

;-----------------------------
;normalizing B on A and D on C
;-----------------------------
imageAB=image[*,512:1023]
imageCD=image[*,0:511]
projx,imageAB,projAB
projx,imageCD,projCD

valueA=median(projAB[(511-gap-stripe):(511-gap)],/EVEN)
valueB=median(projAB[(512+gap):(512+gap+stripe)],/EVEN)
valueC=median(projCD[(511-gap-stripe):(511-gap)],/EVEN)
valueD=median(projCD[(512+gap):(512+gap+stripe)],/EVEN)

proj_deriv_AB=projAB-shift(projAB,1)
proj_deriv_CD=projCD-shift(projCD,1)

derivA=median(proj_deriv_AB[(511-gap-stripe):(511-gap)],/EVEN)
derivB=median(proj_deriv_AB[(512+gap):(512+gap+stripe)],/EVEN)
derivC=median(proj_deriv_CD[(511-gap-stripe):(511-gap)],/EVEN)
derivD=median(proj_deriv_CD[(512+gap):(512+gap+stripe)],/EVEN)
derivAB=(derivA+derivB)/2
derivCD=(derivC+derivD)/2

jumpAB=valueB-valueA-(derivAB*(stripe+gap+gap))
jumpCD=valueD-valueC-(derivCD*(stripe+gap+gap))

image[512:1023,512:1023]=image[512:1023,512:1023]-jumpAB
image[512:1023,0:511]=image[512:1023,0:511]-jumpCD


END