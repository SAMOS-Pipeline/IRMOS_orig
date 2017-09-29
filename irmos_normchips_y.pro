PRO IRMOS_normchips_y,image,stripe,gap

;--------
; A | B |
;--------
; C | D |
;--------


;-----------------------------
;normalizing AB and CD
;-----------------------------

projy,image,projytot
valueDO=median(projytot[(511-gap-stripe):(511-gap)],/EVEN)
valueUP=median(projytot[(512+gap):(512+gap+stripe)],/EVEN)

proj_deriv=projytot-shift(projytot,1)

derivDO=median(proj_deriv[(511-gap-stripe):(511-gap)],/EVEN)
derivUP=median(proj_deriv[(512+gap):(512+gap+stripe)],/EVEN)

deriv=(derivUP+derivDO)/2

jump=valueUP-valueDO-(deriv*(stripe+gap+gap))
image[*,512:1023]=image[*,512:1023]-jump

END