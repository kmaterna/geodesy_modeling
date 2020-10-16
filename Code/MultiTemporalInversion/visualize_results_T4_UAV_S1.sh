#!/bin/bash

ifile=$1
proj="M1.2i"
horiz_scale=100
output="total_output_S1_UAV.ps"
files=`awk 'NR==21' "$1"`
ss_fault_file=`echo $files | awk '{print $1;}'`
field_file="../Mapping_Data/Fields_Boundaries.txt"

gmt set MAP_FRAME_TYPE plain
gmt set FORMAT_GEO_MAP D

gmt makecpt -T-10/10/0.5 -Cpolar -D > mycpt.cpt
gmt makecpt -T-30/30/0.2 -Cjet -D > los.cpt

# UAVSAR
uav_range="-115.60/-115.49/32.99/33.088"
files=`awk 'NR==13' "$1"`
obs_insar_file=`echo $files | awk '{print $1;}'`
model_insar_file=`echo $files | awk '{print $2;}'`
gmt pscoast -J$proj -R$uav_range -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -LjBL+c32+w3+o0.05i/0.19i -Y6i -X1.5i -Dh -K > $output
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output
awk '{print $1, $2, $3*1000}' $obs_insar_file | gmt psxy -Sc0.06 -Clos.cpt -R -J$proj -O -K >> $output
echo "-115.55 33.098 Data" | gmt pstext -R -J$projection -F+f18p,Helvetica-bold -N -K -O >> $output
echo "-115.545 33.075 AB / UAVSAR" | gmt pstext -R -J$projection -F+f12p,Helvetica-bold -K -O >> $output
gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output
gmt pscoast -J$proj -R -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X1.3i -Dh -K -O >> $output
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
awk '{print $1, $2, $3*1000}' $model_insar_file | gmt psxy -Sc0.06 -Clos.cpt -R -J$proj -O -K >> $output
gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output
echo "-115.55 33.098 Model" | gmt pstext -R -J$projection -F+f18p,Helvetica-bold -N -K -O >> $output


files=`awk 'NR==14' "$1"`
obs_insar_file=`echo $files | awk '{print $1;}'`
model_insar_file=`echo $files | awk '{print $2;}'`
gmt pscoast -J$proj -R$uav_range -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -LjBL+c32+w3+o0.05i/0.19i -X-1.3i -Y-1.5i -Dh -K -O >> $output
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output
awk '{print $1, $2, $3*1000}' $obs_insar_file | gmt psxy -Sc0.06 -Clos.cpt -R -J$proj -O -K >> $output
gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output
echo "-115.545 33.075 CDEF / UAVSAR" | gmt pstext -R -J$projection -F+f10p,Helvetica-bold -K -O >> $output
gmt pscoast -J$proj -R -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X1.3i -Dh -K -O >> $output
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
awk '{print $1, $2, $3*1000}' $model_insar_file | gmt psxy -Sc0.06 -Clos.cpt -R -J$proj -O -K >> $output
gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output

gmt psscale -Clos.cpt -Dx-0.5c/-1.0c+w5c/0.25c+jTC+h -Bxaf -By+l"mm" -K -O >> $output



# Sentinel-1
s1_range="-115.675/-115.435/32.91/33.125"
files=`awk 'NR==15' "$1"`
obs_insar_file=`echo $files | awk '{print $1;}'`
model_insar_file=`echo $files | awk '{print $2;}'`
gmt pscoast -J$proj -R$s1_range -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -LjBL+c32+w3+o0.15i/0.19i -X1.6i -Y1.5i -Dh -K -O >> $output
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output
awk '{print $1, $2, $3*1000}' $obs_insar_file | gmt psxy -Sc0.045 -Clos.cpt -R -J$proj -O -K >> $output
gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output
echo "-115.55 33.16 Data" | gmt pstext -R -J$projection -F+f18p,Helvetica-bold -N -K -O >> $output
echo "-115.55 33.105 D / S1_Asc" | gmt pstext -R -J$projection -F+f12p,Helvetica-bold -K -O >> $output
gmt pscoast -J$proj -R -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X1.3i -Dh -K -O >> $output
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
awk '{print $1, $2, $3*1000}' $model_insar_file | gmt psxy -Sc0.045 -Clos.cpt -R -J$proj -O -K >> $output
echo "-115.55 33.16 Model" | gmt pstext -R -J$projection -F+f18p,Helvetica-bold -N -K -O >> $output
gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output



files=`awk 'NR==16' "$1"`
obs_insar_file=`echo $files | awk '{print $1;}'`
model_insar_file=`echo $files | awk '{print $2;}'`
gmt pscoast -J$proj -R$s1_range -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -LjBL+c32+w3+o0.15i/0.19i -X-1.3i -Y-1.5i -Dh -K -O >> $output
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output
awk '{print $1, $2, $3*1000}' $obs_insar_file | gmt psxy -Sc0.045 -Clos.cpt -R -J$proj -O -K >> $output
gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output
echo "-115.55 33.105 E / S1_Asc" | gmt pstext -R -J$projection -F+f12p,Helvetica-bold -K -O >> $output
gmt pscoast -J$proj -R -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X1.3i -Dh -K -O >> $output
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
awk '{print $1, $2, $3*1000}' $model_insar_file | gmt psxy -Sc0.045 -Clos.cpt -R -J$proj -O -K >> $output
gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output


files=`awk 'NR==17' "$1"`
obs_insar_file=`echo $files | awk '{print $1;}'`
model_insar_file=`echo $files | awk '{print $2;}'`
gmt pscoast -J$proj -R$s1_range -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -LjBL+c32+w3+o0.15i/0.19i -X-1.3i -Y-1.5i -Dh -K -O >> $output
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output
awk '{print $1, $2, $3*1000}' $obs_insar_file | gmt psxy -Sc0.045 -Clos.cpt -R -J$proj -O -K >> $output
gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output
echo "-115.55 33.105 F / S1_Asc" | gmt pstext -R -J$projection -F+f12p,Helvetica-bold -K -O >> $output
gmt pscoast -J$proj -R -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X1.3i -Dh -K -O >> $output
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
awk '{print $1, $2, $3*1000}' $model_insar_file | gmt psxy -Sc0.045 -Clos.cpt -R -J$proj -O -K >> $output
gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output


files=`awk 'NR==18' "$1"`
obs_insar_file=`echo $files | awk '{print $1;}'`
model_insar_file=`echo $files | awk '{print $2;}'`
gmt pscoast -J$proj -R$s1_range -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -LjBL+c32+w3+o0.15i/0.19i -X1.6i -Y3.0i -Dh -K -O >> $output
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output
awk '{print $1, $2, $3*1000}' $obs_insar_file | gmt psxy -Sc0.045 -Clos.cpt -R -J$proj -O -K >> $output
gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output
echo "-115.55 33.16 Data" | gmt pstext -R -J$projection -F+f18p,Helvetica-bold -N -K -O >> $output
echo "-115.55 33.105 D / S1_Desc" | gmt pstext -R -J$projection -F+f12p,Helvetica-bold -K -O >> $output
gmt pscoast -J$proj -R -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X1.3i -Dh -K -O >> $output
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
awk '{print $1, $2, $3*1000}' $model_insar_file | gmt psxy -Sc0.045 -Clos.cpt -R -J$proj -O -K >> $output
echo "-115.55 33.16 Model" | gmt pstext -R -J$projection -F+f18p,Helvetica-bold -N -K -O >> $output
gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output




files=`awk 'NR==19' "$1"`
obs_insar_file=`echo $files | awk '{print $1;}'`
model_insar_file=`echo $files | awk '{print $2;}'`
gmt pscoast -J$proj -R$s1_range -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -LjBL+c32+w3+o0.15i/0.19i -X-1.3i -Y-1.5i -Dh -K -O >> $output
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output
awk '{print $1, $2, $3*1000}' $obs_insar_file | gmt psxy -Sc0.045 -Clos.cpt -R -J$proj -O -K >> $output
gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output
echo "-115.55 33.105 E / S1_Desc" | gmt pstext -R -J$projection -F+f12p,Helvetica-bold -K -O >> $output
gmt pscoast -J$proj -R -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X1.3i -Dh -K -O >> $output
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
awk '{print $1, $2, $3*1000}' $model_insar_file | gmt psxy -Sc0.045 -Clos.cpt -R -J$proj -O -K >> $output
gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output



files=`awk 'NR==20' "$1"`
obs_insar_file=`echo $files | awk '{print $1;}'`
model_insar_file=`echo $files | awk '{print $2;}'`
gmt pscoast -J$proj -R$s1_range -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -LjBL+c32+w3+o0.15i/0.19i -X-1.3i -Y-1.5i -Dh -K -O >> $output
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output
awk '{print $1, $2, $3*1000}' $obs_insar_file | gmt psxy -Sc0.045 -Clos.cpt -R -J$proj -O -K >> $output
gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output
echo "-115.55 33.105 F / S1_Desc" | gmt pstext -R -J$projection -F+f12p,Helvetica-bold -K -O >> $output
gmt pscoast -J$proj -R -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X1.3i -Dh -K -O >> $output
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
awk '{print $1, $2, $3*1000}' $model_insar_file | gmt psxy -Sc0.045 -Clos.cpt -R -J$proj -O -K >> $output
gmt psxy $field_file -R -J$proj -Wthin,indianred -O >> $output




rm gmt.history

gmt psconvert -Tg $output

open $output

