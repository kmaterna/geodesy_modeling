#!/bin/bash
# GMT Plot for T4/T5 intervals (fault slip)

slip_file1=$1
slip_file2=$2
slip_file3=$3
slip_file4=$4
slip_file5=$5
slip_file6=$6
slip_file7=$7
scale_low=$8
scale_high=$9
scale_int=${10}
title=${11}

range="-115.62/-115.46/32.97/33.08"
proj="M2.7i"
output=$title"_T4T5.ps"
field_file="../Mapping_Data/Fields_Boundaries.txt"
seismicity_T4A="../../Misc_Geophysics_Exps/QTM_exploring/Steps/T4A_depth_0_12_32.9_33.1_20120930_20130915/Brawley_QTM.txt"
seismicity_T4B="../../Misc_Geophysics_Exps/QTM_exploring/Steps/T4B_depth_0_12_32.9_33.1_20130915_20140610/Brawley_QTM.txt"
seismicity_T4C="../../Misc_Geophysics_Exps/QTM_exploring/Steps/T4C_depth_0_12_32.9_33.1_20140610_20141015/Brawley_QTM.txt"
seismicity_T4D="../../Misc_Geophysics_Exps/QTM_exploring/Steps/T4D_depth_0_12_32.9_33.1_20141015_20151130/Brawley_QTM.txt"
seismicity_T4E="../../Misc_Geophysics_Exps/QTM_exploring/Steps/T4E_depth_0_12_32.9_33.1_20151130_20161030/Brawley_QTM.txt"
seismicity_T4F="../../Misc_Geophysics_Exps/QTM_exploring/Steps/T4F_depth_0_12_32.9_33.1_20161030_20171101/Brawley_QTM.txt"
seismicity_T5="../../Misc_Geophysics_Exps/QTM_exploring/Steps/T5_depth_0_12_32.9_33.1_20171101_20181101/Brawley_QTM.txt"

gmt set MAP_FRAME_TYPE plain
gmt set FORMAT_GEO_MAP D

gmt makecpt -T$scale_low/$scale_high/$scale_int -Cpolar -D > mycpt.cpt
gmt makecpt -T0/12/0.1 -Cjet > depth.cpt

# Plotting the slip distribution
gmt pscoast -R$range -J$proj -Wthin,black -Di -N1 -N2 -B0.2WesN -Y5.5i -LjBR+c32+w5+o0.10i/0.25i -Gwhite -Slightgray -K > $output
gmt psxy $slip_file1 -R -J -Wthinnest,gray -Cmycpt.cpt -L -K -O >> $output
gmt psxy $field_file -R -J -Wthin,indianred -O -K >> $output
awk '{ print $2, $3, $4, $5*0.04}' $seismicity_T4A | gmt psxy -R -J -Sc -Cdepth.cpt -K -O >> $output
echo "-115.59 33.07 T4-A" | gmt pstext -F+f22p,Helvetica-bold -R -J -K -O >> $output

gmt pscoast -R$range -J$proj -Wthin,black -Di -N1 -N2 -B0.2wesn -X2.9i -LjBR+c32+w5+o0.10i/0.25i -Gwhite -Slightgray -K -O >> $output
gmt psxy $slip_file2 -R -J -Wthinnest,gray -Cmycpt.cpt -L -K -O >> $output
gmt psxy $field_file -R -J -Wthin,indianred -O -K >> $output
awk '{ print $2, $3, $4, $5*0.04}' $seismicity_T4B | gmt psxy -R -J -Sc -Cdepth.cpt -K -O >> $output
echo "-115.59 33.07 T4-B" | gmt pstext -F+f22p,Helvetica-bold -R -J -K -O >> $output

gmt pscoast -R$range -J$proj -Wthin,black -Di -N1 -N2 -B0.2wesN -X2.9i -LjBR+c32+w5+o0.10i/0.25i -Gwhite -Slightgray -K -O >> $output
gmt psxy $slip_file3 -R -J -Wthinnest,gray -Cmycpt.cpt -L -K -O >> $output
gmt psxy $field_file -R -J -Wthin,indianred -O -K >> $output
awk '{ print $2, $3, $4, $5*0.04}' $seismicity_T4C | gmt psxy -R -J -Sc -Cdepth.cpt -K -O >> $output
echo "-115.59 33.07 T4-C" | gmt pstext -F+f22p,Helvetica-bold -R -J -K -O >> $output

gmt pscoast -R$range -J$proj -Wthin,black -Di -N1 -N2 -B0.2Wesn -X-5.8i -Y-2.4i -LjBR+c32+w5+o0.10i/0.25i -Gwhite -Slightgray -K -O >> $output
gmt psxy $slip_file4 -R -J -Wthinnest,gray -Cmycpt.cpt -L -K -O >> $output
gmt psxy $field_file -R -J -Wthin,indianred -O -K >> $output
awk '{ print $2, $3, $4, $5*0.04}' $seismicity_T4D | gmt psxy -R -J -Sc -Cdepth.cpt -K -O >> $output
echo "-115.59 33.07 T4-D" | gmt pstext -F+f22p,Helvetica-bold -R -J -K -O >> $output

gmt pscoast -R$range -J$proj -Wthin,black -Di -N1 -N2 -B0.2wesn -X2.9i -LjBR+c32+w5+o0.10i/0.25i -Gwhite -Slightgray -K -O >> $output
gmt psxy $slip_file5 -R -J -Wthinnest,gray -Cmycpt.cpt -L -K -O >> $output
gmt psxy $field_file -R -J -Wthin,indianred -O -K >> $output
awk '{ print $2, $3, $4, $5*0.04}' $seismicity_T4E | gmt psxy -R -J -Sc -Cdepth.cpt -K -O >> $output
echo "-115.59 33.07 T4-E" | gmt pstext -F+f22p,Helvetica-bold -R -J -K -O >> $output

gmt pscoast -R$range -J$proj -Wthin,black -Di -N1 -N2 -B0.2wesn -X2.9i -LjBR+c32+w5+o0.10i/0.25i -Gwhite -Slightgray -K -O >> $output
gmt psxy $slip_file6 -R -J -Wthinnest,gray -Cmycpt.cpt -L -K -O >> $output
gmt psxy $field_file -R -J -Wthin,indianred -O -K >> $output
awk '{ print $2, $3, $4, $5*0.04}' $seismicity_T4F | gmt psxy -R -J -Sc -Cdepth.cpt -K -O >> $output
echo "-115.59 33.07 T4-F" | gmt pstext -F+f22p,Helvetica-bold -R -J -K -O >> $output

gmt pscoast -R$range -J$proj -Wthin,black -Di -N1 -N2 -B0.2weSn -X-5.8i -Y-2.4i -LjBR+c32+w5+o0.10i/0.25i -Gwhite -Slightgray -K -O >> $output
gmt psxy $slip_file7 -R -J -Wthinnest,gray -Cmycpt.cpt -L -K -O >> $output
gmt psxy $field_file -R -J -Wthin,indianred -O -K >> $output
awk '{ print $2, $3, $4, $5*0.04}' $seismicity_T5 | gmt psxy -R -J -Sc -Cdepth.cpt -K -O >> $output
echo "-115.59 33.07 T5" | gmt pstext -F+f22p,Helvetica-bold -R -J -K -O >> $output

gmt psscale -Cmycpt.cpt -Dx12.5c/4.2c+w8c/0.5c+jTC+h -Bxaf -By+l"Reverse slip(m)" -K -O >> $output
gmt psscale -Cdepth.cpt -Dx12.5c/2.1c+w8c/0.5c+jTC+h -Bxaf -By+l"Depth (km)" -O >> $output

gmt psconvert -Tg $output

# open $output