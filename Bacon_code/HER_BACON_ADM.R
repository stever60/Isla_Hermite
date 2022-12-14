# Isla Hermite - age depth modelling 
setwd("/Users/Steve/Dropbox/BAS/Data/R/RBacon")
#check working directory
getwd()
library(rbacon)

#clear plot window
dev.off()

# Isla Hermite -----------------------------------------------------------------

# in order from north to south

#HER14-M5 USE = Hiatus of 8700 years centred on 75 cm depth (+/- 1 cm) - USE
Bacon("HER14L",depths.file=FALSE,thick=10,rotate.axes=TRUE,cc=3, postbomb=4,acc.mean=100,
      mem.mean=0.1,mem.strength=10,yr.max=20000,d.min=0,d.max=170,rounded=2,hiatus.depths=76,hiatus.max=8700,
      C14.border=rgb(0, 0, 1, 1),C14.col=rgb(0, 0, 1, 0.5),
      cal.col=rgb(0, 0.5, 0.5, 0.5), cal.border=rgb(0, 0.5, 0.5, 1),
      dark=0.7, height=2, title.location='topright',mgp=c(1.5, 0.7, 0),
      model.only=TRUE
      )

#HER14-M6 = Hiatus of 8700 years centred on 75 cm depth (+/- 1 cm) & DR of 700 applied to large Holocene age offsets
Bacon("HER14L",depths.file=TRUE,thick=10,rotate.axes=TRUE,cc=3, postbomb=4,acc.mean=100,rounded=2,
      mem.mean=0.1,mem.strength=20,yr.max=20000,d.min=0,d.max=160,hiatus.depths=c(75),hiatus.max=8700,
      C14.border=rgb(0, 0, 1, 1),C14.col=rgb(0, 0, 1, 0.5), 
      cal.col=rgb(0, 0.5, 0.5, 0.5), cal.border=rgb(0, 0.5, 0.5, 1),
      dark=0.7, height=20, title.location='topright',mgp=c(1.5, 0.7, 0),
      model.only=TRUE
      )

#HER24-M3 - revised strat depths 8/4/22
Bacon("HER24L",depths.file=FALSE,thick=5,rotate.axes=TRUE,cc=3, postbomb=4,acc.mean=100,rounded=2,
      mem.mean=0.4,mem.strength=20,yr.max=20000,d.min=0,d.max=100,
      C14.border=rgb(0, 0, 1, 1),C14.col=rgb(0, 0, 1, 0.5), 
      cal.col=rgb(0, 0.5, 0.5, 0.5), cal.border=rgb(0, 0.5, 0.5, 1),
      dark=0.7, height=0.5, title.location='topright',mgp=c(1.5, 0.7, 0),
      model.only=TRUE
      )

# HERPT1.1 - SS1 - HER24L
# HERPT1.1 peat record is above SS1 near HER24 and joins directly onto it at 105 cm depth - Unit4A is in both records
# Preliminary age-depth model is based on MS peak correlation of HER24 and HERPT1S2 peat record
# using 2 main HER24 MS peaks at c. 7 ka and 4 ka
Bacon("HERPT1.1_SS1_HER24",depths.file=TRUE,d.max=200,thick=10, cc=3, 
      postbomb=4,rotate.axes=TRUE,mem.mean=0.4,rounded=2,
      mem.strength=20,acc.mean=50, yr.max=20000,
      C14.border=rgb(0, 0, 1, 1), C14.col=rgb(0, 0, 1, 0.5), 
      cal.col=rgb(0, 0.5, 0.5, 0.5), cal.border=rgb(0, 0.5, 0.5, 1),
      dark=0.5, height=75, 
      title.location='topright', mgp=c(1.5, 0.7, 0),
      model.only=TRUE
)

#HER34-M1
Bacon("HER34L",depths.file=FALSE,thick=5,rotate.axes=TRUE,cc=3, postbomb=4,acc.mean=100,rounded=2,
      mem.mean=0.4,mem.strength=20,yr.max=5000,d.min=0,d.max=65,
      C14.border=rgb(0, 0, 1, 1),C14.col=rgb(0, 0, 1, 0.5), 
      cal.col=rgb(0, 0.5, 0.5, 0.5), cal.border=rgb(0, 0.5, 0.5, 1),
      dark=0.7, height=0.3, title.location='topright',mgp=c(1.5, 0.7, 0),
      model.only=TRUE
      )

#HER42L-M3
Bacon("HER42L",depths.file=TRUE,d.max=300,thick=5,cc=3,
      postbomb=4,rotate.axes=TRUE,rounded=2,
      mem.mean=0.5, mem.strength=20,
      boundary = 8, acc.mean=c(50,10),yr.max=5000,
      C14.border=rgb(0, 0, 1, 1), C14.col=rgb(0, 0, 1, 0.5), 
      cal.col=rgb(0, 0.5, 0.5, 0.5), cal.border=rgb(0, 0.5, 0.5, 1),
      dark=0.5, height=0.75, plot.pdf = TRUE,
      title.location='topright',
      model.only=TRUE
      )

#HER42PB-M2
Bacon("HER42PB",depths.file=TRUE,d.max=500,thick=10, cc=3, 
      postbomb=4,rotate.axes=TRUE,mem.mean=0.4,rounded=2,
      mem.strength=20,acc.mean=10, yr.max=15000,
      C14.border=rgb(0, 0, 1, 1), C14.col=rgb(0, 0, 1, 0.5), 
      cal.col=rgb(0, 0.5, 0.5, 0.5), cal.border=rgb(0, 0.5, 0.5, 1),
      dark=0.5, height=50, 
      title.location='topright', mgp=c(1.5, 0.7, 0),
      model.only=TRUE)

# HER42PB+SS2
# HER42PB peat record is above SS2
# SS2 is +815 cm below the base of HER42PB at 405 cm (590 cm) based on extrapolating C14 ages model
Bacon("HER42PB_PT2S1_SS2",depths.file=TRUE,d.max=1000,thick=10, cc=3, 
      postbomb=4,rotate.axes=TRUE,mem.mean=0.4,rounded=2,
      mem.strength=20,acc.mean=10, yr.max=20000,
      C14.border=rgb(0, 0, 1, 1), C14.col=rgb(0, 0, 1, 0.5), 
      cal.col=rgb(0, 0.5, 0.5, 0.5), cal.border=rgb(0, 0.5, 0.5, 1),
      dark=0.5, height=75, 
      title.location='topright', mgp=c(1.5, 0.7, 0),
      model.only=TRUE
      )

#HER42PB + SS1 (from above) + SS2 (from above) - FINAL FIGURE FOR PAPER
Bacon("HER42PB_SS2_SS1",depths.file=FALSE,d.max=1200,thick=10, cc=3, 
      postbomb=4,rotate.axes=TRUE,mem.mean=0.4,rounded=2,
      mem.strength=20,acc.mean=10, yr.max=20000,
      C14.border=rgb(0, 0, 1, 1), C14.col=rgb(0, 0, 1, 0.5), 
      cal.col=rgb(0, 0.5, 0.5, 0.5), cal.border=rgb(0, 0.5, 0.5, 1),
      dark=0.5, height=75, 
      title.location='topright', mgp=c(1.5, 0.7, 0),
      model.only=TRUE
      )

#HER44L-M2 - maximum age model
Bacon("HER44L",depths.file=TRUE,d.max=60,thick=4, cc=3,
      postbomb=4,rotate.axes=TRUE,rounded=2,
      mem.mean=0.2, mem.strength=20,acc.mean=200,yr.max=10000,
      C14.border=rgb(0, 0, 1, 1), C14.col=rgb(0, 0, 1, 0.5), 
      cal.col=rgb(0, 0.5, 0.5, 0.5), cal.border=rgb(0, 0.5, 0.5, 1),
      dark=0.5, height=0.2, 
      title.location='topright', mgp=c(1.5, 0.7, 0),
      model.only=TRUE
      )

#HER49L-M4 USE - dates 2 and 3  removed - large age reversal / reworking of older material 
#HER49L-M3  - all dates included - overlay plots
Bacon("HER49L",depths.file=FALSE,d.max=260,thick=5,cc=3,
      postbomb=4,rotate.axes=TRUE,rounded=2, 
      mem.mean=0.2, mem.strength=20,acc.mean=20,yr.max=10000,
      C14.border=rgb(0, 0, 1, 1), C14.col=rgb(0, 0, 1, 0.5), 
      cal.col=rgb(0, 0.5, 0.5, 0.5), cal.border=rgb(0, 0.5, 0.5, 1),
      dark=0.5,height=5, plot.pdf = TRUE,
      title.location='topright',
      model.only=TRUE)

