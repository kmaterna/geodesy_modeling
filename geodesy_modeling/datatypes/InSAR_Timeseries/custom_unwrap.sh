#!/bin/csh

# These functions are actually used in the ISCE/InSAR_Timeseries pipeline
# November 2020



set sharedir = `gmtsar_sharedir.csh`

gmt grdmath corr.grd 0.1 GE 0 NAN = mask2_patch.grd
gmt grdmath corr.grd 0. XOR 1. MIN  = corr_patch.grd
gmt grdmath mask2_patch.grd corr_patch.grd MUL = corr_tmp.grd 

gmt grd2xyz phase_interp.grd -ZTLf -do0 > phase.in
gmt grd2xyz corr_tmp.grd -ZTLf  -do0 > corr.in

snaphu phase.in `gmt grdinfo -C phase_interp.grd | cut -f 10` -f $sharedir/snaphu/config/snaphu.conf.brief -c corr.in -o unwrap.out -v -s

gmt xyz2grd unwrap.out -ZTLf `gmt grdinfo -I- phase_interp.grd` `gmt grdinfo -I phase_interp.grd` -Gtmp.grd 
mv tmp.grd unwrap.grd
gmt grdmath unwrap.grd phase_interp.grd OR = unwrap.grd
gmt grdgradient unwrap.grd -Nt.9 -A0. -Gunwrap_grad.grd
set tmp = `gmt grdinfo -C -L2 unwrap.grd`
set limitU = `echo $tmp | awk '{printf("%5.1f", $12+$13*2)}'`
set limitL = `echo $tmp | awk '{printf("%5.1f", $12-$13*2)}'`
set std = `echo $tmp | awk '{printf("%5.1f", $13)}'`
gmt makecpt -Cseis -I -Z -T"$limitL"/"$limitU"/1 -D > unwrap.cpt
set boundR = `gmt grdinfo unwrap.grd -C | awk '{print ($3-$2)/4}'`
set boundA = `gmt grdinfo unwrap.grd -C | awk '{print ($5-$4)/4}'`
gmt grdimage unwrap.grd -Iunwrap_grad.grd -Cunwrap.cpt -JX6.5i -B"$boundR":Range:/"$boundA":Azimuth:WSen -X1.3i -Y3i -P -K > unwrap.ps
gmt psscale -Dx3.3/-1.5+w5/0.2+h+e -Cunwrap.cpt -B"$std":"unwrapped phase, rad": -O >> unwrap.ps
