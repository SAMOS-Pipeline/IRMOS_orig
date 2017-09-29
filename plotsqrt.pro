function plotsqrt,input,percent

dummy=sqrt(input>0)
dummy2=dummy[sort(dummy)]
_init=dummy2[((512.*512.)*((100.-percent)/200.))-1]
_end=dummy2[(512.*512.)-1-((512.*512.)*((100.-percent)/200.))]
dummy2=dummy>_init<_end

return,dummy2

end