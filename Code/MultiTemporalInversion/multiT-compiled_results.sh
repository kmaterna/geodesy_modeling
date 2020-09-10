#!/bin/bash
# August 7, 2020
# An overview plot of seismicity, fault slip, and beachballs on T1/T2/T3

southrange="-115.68/-115.43/32.92/33.12"
proj="M2.4i"
faultmodelT1="../T1_MultiT/output/predicted_slip.txtthrust_gmt"
faultmodelT2="../T2_MultiT/output/predicted_slip.txtthrust_gmt"
faultmodelT3="../T3_MultiT/output/predicted_slip.txtthrust_FAULTORDERREV_gmt"
field_bounds="../../Misc_Geophysics_Exps/Injection_Data/Data/Fields_Boundaries.txt"
seismicity_T1="../../Misc_Geophysics_Exps/QTM_exploring/Steps/depth_0_12_32.9_33.1_20091015_20101215/Brawley_QTM.txt"
seismicity_T2="../../Misc_Geophysics_Exps/QTM_exploring/Steps/depth_0_12_32.9_33.1_20101215_20111120/Brawley_QTM.txt"
seismicity_T3="../../Misc_Geophysics_Exps/QTM_exploring/Steps/depth_0_12_32.9_33.1_20111120_20120930/Brawley_QTM.txt"
surface_rupture="../../Misc_Geophysics_Exps/Injection_Data/Data/M4p9_surface_rupture.txt"
FMfile="../../Misc_Geophysics_Exps/Injection_Data/Data/Wei_EPSL_2015_SupplementS1.txt"
output="T1T2T3.ps"


gmt set MAP_FRAME_TYPE plain
gmt set FORMAT_GEO_MAP D

gmt makecpt -T-.4/.4/0.001 -D -Cpolar > slip.cpt
gmt makecpt -T0/12/0.1 -Cjet > depth.cpt


# T1
gmt pscoast -J$proj -R$southrange -Gwhite -Slightgray -N1 -Wthin,black -B0.1WeSn -LjBR+c32+w5+o0.05i/0.25i -Y5.0 -Dh -K > $output
gmt psxy $faultmodelT3 -R$southrange -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
gmt psxy $faultmodelT1 -R$southrange -J$proj -W0.01p,gray -Cslip.cpt -L -K -O >> $output
awk '{ print $2, $3, $4, $5*0.04}' $seismicity_T1 | gmt psxy -R$southrange -J$proj -Sc -Cdepth.cpt -K -O >> $output
gmt psxy $field_bounds -J$proj -R$southrange -Wthick,indianred -K -O >> $output
echo "-115.59 33.105 A) T1: 2009-2010" | gmt pstext -R$southrange -J$proj -F+f15p,Helvetica-bold -K -O >> $output
echo "-115.61 32.935 M5.0 (aseismic)" | gmt pstext -R$southrange -J$proj -F+f12p,Helvetica-bold -K -O >> $output

# T2
gmt pscoast -J$proj -R$southrange -Gwhite -Slightgray -N1 -Wthin,black -B0.1weSn -LjBR+c32+w5+o0.05i/0.25i -X2.6i -Dh -K -O >> $output
gmt psxy $faultmodelT3 -R$southrange -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
gmt psxy $faultmodelT2 -R$southrange -J$proj -W0.01p,gray -Cslip.cpt -L -K -O >> $output
awk '{ print $2, $3, $4, $5*0.04}' $seismicity_T2 | gmt psxy -R$southrange -J$proj -Sc -Cdepth.cpt -K -O >> $output
gmt psxy $field_bounds -J$proj -R$southrange -Wthick,indianred -K -O >> $output
echo "-115.59 33.105 B) T2: 2010-2011" | gmt pstext -R$southrange -J$proj -F+f15p,Helvetica-bold -K -O >> $output
echo "-115.61 32.935 M5.1 (aseismic)" | gmt pstext -R$southrange -J$proj -F+f12p,Helvetica-bold -K -O >> $output

# T3
gmt pscoast -J$proj -R$southrange -Gwhite -Slightgray -N1 -Wthin,black -B0.1wESn -LjBR+c32+w5+o0.05i/0.25i -X2.6i -Dh -K -O >> $output
gmt psxy $faultmodelT3 -R$southrange -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
gmt psxy $faultmodelT3 -R$southrange -J$proj -W0.01p,gray -Cslip.cpt -L -K -O >> $output
awk '{ print $2, $3, $4, $5*0.015}' $seismicity_T3 | gmt psxy -R$southrange -J$proj -Sc -Cdepth.cpt -K -O >> $output
gmt psxy $field_bounds -J$proj -R$southrange -Wthick,indianred -K -O >> $output
gmt psxy $surface_rupture -J$proj -R$southrange -Wthick,black -K -O >> $output
# awk '{ print $3, $2, $4, $6, $7, $8, $5, 0, 0}' $FMfile | gmt psmeca -J$proj -R$southrange -Sa0.2 -Ggray40 -K -O >> $output
echo "-115.59 33.105 C) T3: 2011-2012" | gmt pstext -R$southrange -J$proj -F+f15p,Helvetica-bold -K -O >> $output
echo "-115.60 32.935 M5.8 (seis./aseis.)" | gmt pstext -R$southrange -J$proj -F+f12p,Helvetica-bold -K -O >> $output


gmt psscale -Cslip.cpt -Dx-0.25c/-1.3c+w7c/0.35c+jTC+h -Bxaf -By+l"Reverse Slip (m)" -O -K >> $output
gmt psscale -Cdepth.cpt -Dx-9.75c/-1.3c+w7c/0.35c+jTC+h -Bxaf -By+l"Depth (km)" -O >> $output

gmt psconvert $output -Tg

open $output

