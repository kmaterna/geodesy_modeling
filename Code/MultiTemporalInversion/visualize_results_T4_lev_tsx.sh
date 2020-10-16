#!/bin/bash

ifile=$1
proj="M1.2i"
leveling_range="-115.60/-115.49/32.975/33.07"
gps_range="-115.7/-115.3867/32.8556/33.159"
horiz_scale=100
output="total_output_lev_tsx.ps"
files=`awk 'NR==21' "$1"`
ss_fault_file=`echo $files | awk '{print $1;}'`
field_file="../Mapping_Data/Fields_Boundaries.txt"

gmt set MAP_FRAME_TYPE plain
gmt set FORMAT_GEO_MAP D

gmt makecpt -T-10/10/0.5 -Cpolar -D > mycpt.cpt
gmt makecpt -T-30/30/1 -Cjet -D > los.cpt


## GPS PLOTS
#files=`awk 'NR==1' "$1"`
#observed_gps_file=`echo $files | awk '{print $1;}'`
#predicted_gps_file=`echo $files | awk '{print $2;}'`
#label=`echo $files | awk '{print $3, $4, $5;}'`
#gmt pscoast -J$proj -R$gps_range -Gwhite -Slightgray -N1 -Wthin,black -B0.2WesN -Y10i -X0.5i -Dh -K -P > $output
#awk '{print $1, $2, $5*1000}' $observed_gps_file | gmt psxy -Sh0.35 -Cmycpt.cpt -W0.1,black -R -J$proj -O -K >> $output
#awk '{print $1, $2, $5*1000}' $predicted_gps_file | gmt psxy -Sh0.2 -Cmycpt.cpt -W0.1,red -R -J$proj -O -K >> $output
#awk '{print $1, $2, $3, $4, $6, $7, 0}'  $observed_gps_file | gmt psvelo -A+e+gblack+pthickest -Wblack -Se$horiz_scale/0.68/10 -R -J$proj -O -K >> $output
#awk '{print $1, $2, $3, $4, $6, $7, 0}' $predicted_gps_file | gmt psvelo -A+e+gred+pthickest -Wred -Se$horiz_scale/0.68/10 -R -J$proj -O -K >> $output
#gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output
#echo $label | gmt pstext -R -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
#
#
#files=`awk 'NR==2' "$1"`
#observed_gps_file=`echo $files | awk '{print $1;}'`
#predicted_gps_file=`echo $files | awk '{print $2;}'`
#label=`echo $files | awk '{print $3, $4, $5;}'`
#gmt pscoast -J$proj -R$gps_range -Gwhite -Slightgray -N1 -Wthin,black -B0.2Wesn -Y-1.65i -Dh -K -O -P >> $output
#awk '{print $1, $2, $5*1000}' $observed_gps_file | gmt psxy -Sh0.35 -Cmycpt.cpt -W0.1,black -R -J$proj -O -K >> $output
#awk '{print $1, $2, $5*1000}' $predicted_gps_file | gmt psxy -Sh0.2 -Cmycpt.cpt -W0.1,red -R -J$proj -O -K >> $output
#awk '{print $1, $2, $3, $4, $6, $7, 0}'  $observed_gps_file | gmt psvelo -A+e+gblack+pthickest -Wblack -Se$horiz_scale/0.68/10 -R -J$proj -O -K >> $output
#awk '{print $1, $2, $3, $4, $6, $7, 0}' $predicted_gps_file | gmt psvelo -A+e+gred+pthickest -Wred -Se$horiz_scale/0.68/10 -R -J$proj -O -K >> $output
#gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output
#echo $label | gmt pstext -R -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
#
#
#files=`awk 'NR==3' "$1"`
#observed_gps_file=`echo $files | awk '{print $1;}'`
#predicted_gps_file=`echo $files | awk '{print $2;}'`
#label=`echo $files | awk '{print $3, $4, $5;}'`
#gmt pscoast -J$proj -R$gps_range -Gwhite -Slightgray -N1 -Wthin,black -B0.2Wesn -Y-1.65i -Dh -K -O -P >> $output
#awk '{print $1, $2, $5*1000}' $observed_gps_file | gmt psxy -Sh0.35 -Cmycpt.cpt -W0.1,black -R -J$proj -O -K >> $output
#awk '{print $1, $2, $5*1000}' $predicted_gps_file | gmt psxy -Sh0.2 -Cmycpt.cpt -W0.1,red -R -J$proj -O -K >> $output
#awk '{print $1, $2, $3, $4, $6, $7, 0}'  $observed_gps_file | gmt psvelo -A+e+gblack+pthickest -Wblack -Se$horiz_scale/0.68/10 -R -J$proj -O -K >> $output
#awk '{print $1, $2, $3, $4, $6, $7, 0}' $predicted_gps_file | gmt psvelo -A+e+gred+pthickest -Wred -Se$horiz_scale/0.68/10 -R -J$proj -O -K >> $output
#gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output
#echo $label | gmt pstext -R -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
#
#
#files=`awk 'NR==4' "$1"`
#observed_gps_file=`echo $files | awk '{print $1;}'`
#predicted_gps_file=`echo $files | awk '{print $2;}'`
#label=`echo $files | awk '{print $3, $4, $5;}'`
#gmt pscoast -J$proj -R$gps_range -Gwhite -Slightgray -N1 -Wthin,black -B0.2Wesn -Y-1.65i -Dh -K -O -P >> $output
#awk '{print $1, $2, $5*1000}' $observed_gps_file | gmt psxy -Sh0.35 -Cmycpt.cpt -W0.1,black -R -J$proj -O -K >> $output
#awk '{print $1, $2, $5*1000}' $predicted_gps_file | gmt psxy -Sh0.2 -Cmycpt.cpt -W0.1,red -R -J$proj -O -K >> $output
#awk '{print $1, $2, $3, $4, $6, $7, 0}'  $observed_gps_file | gmt psvelo -A+e+gblack+pthickest -Wblack -Se$horiz_scale/0.68/10 -R -J$proj -O -K >> $output
#awk '{print $1, $2, $3, $4, $6, $7, 0}' $predicted_gps_file | gmt psvelo -A+e+gred+pthickest -Wred -Se$horiz_scale/0.68/10 -R -J$proj -O -K >> $output
#gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output
#echo $label | gmt pstext -R -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
#
#
#files=`awk 'NR==5' "$1"`
#observed_gps_file=`echo $files | awk '{print $1;}'`
#predicted_gps_file=`echo $files | awk '{print $2;}'`
#label=`echo $files | awk '{print $3, $4, $5;}'`
#gmt pscoast -J$proj -R$gps_range -Gwhite -Slightgray -N1 -Wthin,black -B0.2Wesn -Y-1.65i -Dh -K -O -P >> $output
#awk '{print $1, $2, $5*1000}' $observed_gps_file | gmt psxy -Sh0.35 -Cmycpt.cpt -W0.1,black -R -J$proj -O -K >> $output
#awk '{print $1, $2, $5*1000}' $predicted_gps_file | gmt psxy -Sh0.2 -Cmycpt.cpt -W0.1,red -R -J$proj -O -K >> $output
#awk '{print $1, $2, $3, $4, $6, $7, 0}'  $observed_gps_file | gmt psvelo -A+e+gblack+pthickest -Wblack -Se$horiz_scale/0.68/10 -R -J$proj -O -K >> $output
#awk '{print $1, $2, $3, $4, $6, $7, 0}' $predicted_gps_file | gmt psvelo -A+e+gred+pthickest -Wred -Se$horiz_scale/0.68/10 -R -J$proj -O -K >> $output
#gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output
#echo $label | gmt pstext -R -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
#
#
#files=`awk 'NR==6' "$1"`
#observed_gps_file=`echo $files | awk '{print $1;}'`
#predicted_gps_file=`echo $files | awk '{print $2;}'`
#label=`echo $files | awk '{print $3, $4, $5;}'`
#gmt pscoast -J$proj -R$gps_range -Gwhite -Slightgray -N1 -Wthin,black -B0.2Wesn -Y-1.65i -Dh -K -O -P >> $output
#awk '{print $1, $2, $5*1000}' $observed_gps_file | gmt psxy -Sh0.35 -Cmycpt.cpt -W0.1,black -R -J$proj -O -K >> $output
#awk '{print $1, $2, $5*1000}' $predicted_gps_file | gmt psxy -Sh0.2 -Cmycpt.cpt -W0.1,red -R -J$proj -O -K >> $output
#awk '{print $1, $2, $3, $4, $6, $7, 0}'  $observed_gps_file | gmt psvelo -A+e+gblack+pthickest -Wblack -Se$horiz_scale/0.68/10 -R -J$proj -O -K >> $output
#awk '{print $1, $2, $3, $4, $6, $7, 0}' $predicted_gps_file | gmt psvelo -A+e+gred+pthickest -Wred -Se$horiz_scale/0.68/10 -R -J$proj -O -K >> $output
#gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output
#echo $label | gmt pstext -R -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
#


# TSX
#smallrange="-115.61/-115.47/32.96/33.08"  # wrong range
smallrange="-115.59/-115.48/32.975/33.07"
files=`awk 'NR==7' "$1"`
obs_insar_file=`echo $files | awk '{print $1;}'`
model_insar_file=`echo $files | awk '{print $2;}'`
gmt pscoast -J$proj -R$smallrange -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -LjBR+c32+w3+o0.05i/0.19i -Y8.25i -X1.5i -Dh -K -P > $output
#gmt pscoast -J$proj -R$smallrange -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X1.6i -Y0.25i -Dh -K -O -P >> $output
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output
awk '{print $1, $2, $3*1000}' $obs_insar_file | gmt psxy -Sc0.08 -Clos.cpt -R -J$proj -O -K >> $output
gmt psxy $field_file -R -J$proj -Wthinnest,indianred -O -K >> $output
echo "-115.54 33.08 Data" | gmt pstext -R -J$projection -F+f18p,Helvetica-bold -N -K -O >> $output
echo "-115.55 33.06 A / TSX" | gmt pstext -R -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
gmt pscoast -J$proj -R -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X1.3i -Dh -K -O -P >> $output
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
awk '{print $1, $2, $3*1000}' $model_insar_file | gmt psxy -Sc0.08 -Clos.cpt -R -J$proj -O -K >> $output
gmt psxy $field_file -R -J$proj -Wthinnest,indianred -O -K >> $output
echo "-115.54 33.08 Model" | gmt pstext -R -J$projection -F+f18p,Helvetica-bold -N -K -O >> $output


#LEVELING
files=`awk 'NR==8' "$1"`
obs_insar_file=`echo $files | awk '{print $1;}'`
model_insar_file=`echo $files | awk '{print $2;}'`
gmt pscoast -J$proj -R$leveling_range -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -LjBR+c32+w3+o0.05i/0.19i -X-1.3i -Y-1.65i -Dh -K -O -P >> $output
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output
awk '{print $1, $2, $3*1000}' $obs_insar_file | gmt psxy -Sc0.16 -Clos.cpt -R -J$proj -O -K >> $output
gmt psxy $field_file -R -J$proj -Wthinnest,indianred -O -K >> $output
echo '-115.57 33.06 A / Lev' | gmt pstext -R -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
gmt pscoast -J$proj -R -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X1.3i -Dh -K -O -P >> $output
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
awk '{print $1, $2, $3*1000}' $model_insar_file | gmt psxy -Sc0.16 -Clos.cpt -R -J$proj -O -K >> $output
gmt psxy $field_file -R -J$proj -Wthinnest,indianred -O -K >> $output



files=`awk 'NR==9' "$1"`
obs_insar_file=`echo $files | awk '{print $1;}'`
model_insar_file=`echo $files | awk '{print $2;}'`
gmt pscoast -J$proj -R$leveling_range -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -LjBR+c32+w3+o0.05i/0.19i -X-1.3i -Y-1.65i -Dh -K -O -P >> $output
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output
awk '{print $1, $2, $3*1000}' $obs_insar_file | gmt psxy -Sc0.16 -Clos.cpt -R -J$proj -O -K >> $output
gmt psxy $field_file -R -J$proj -Wthinnest,indianred -O -K >> $output
echo '-115.565 33.06 BC / Lev' | gmt pstext -R -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
gmt pscoast -J$proj -R -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X1.3i -Dh -K -O -P >> $output
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
awk '{print $1, $2, $3*1000}' $model_insar_file | gmt psxy -Sc0.16 -Clos.cpt -R -J$proj -O -K >> $output
gmt psxy $field_file -R -J$proj -Wthinnest,indianred -O -K >> $output

gmt psscale -Clos.cpt -Dx-0.5c/-0.5c+w5c/0.25c+jTC+h -Bxaf -By+l"mm" -K -O >> $output


files=`awk 'NR==10' "$1"`
obs_insar_file=`echo $files | awk '{print $1;}'`
model_insar_file=`echo $files | awk '{print $2;}'`
gmt pscoast -J$proj -R$leveling_range -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -LjBR+c32+w3+o0.05i/0.19i -X1.6i -Y3.3i -Dh -K -O -P >> $output
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output
awk '{print $1, $2, $3*1000}' $obs_insar_file | gmt psxy -Sc0.16 -Clos.cpt -R -J$proj -O -K >> $output
gmt psxy $field_file -R -J$proj -Wthinnest,indianred -O -K >> $output
echo '-115.57 33.06 D / Lev' | gmt pstext -R -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
echo "-115.55 33.08 Data" | gmt pstext -R -J$projection -F+f18p,Helvetica-bold -N -K -O >> $output
gmt pscoast -J$proj -R -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X1.3i -Dh -K -O -P >> $output
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
awk '{print $1, $2, $3*1000}' $model_insar_file | gmt psxy -Sc0.16 -Clos.cpt -R -J$proj -O -K >> $output
gmt psxy $field_file -R -J$proj -Wthinnest,indianred -O -K >> $output
echo "-115.55 33.08 Model" | gmt pstext -R -J$projection -F+f18p,Helvetica-bold -N -K -O >> $output


files=`awk 'NR==11' "$1"`
obs_insar_file=`echo $files | awk '{print $1;}'`
model_insar_file=`echo $files | awk '{print $2;}'`
gmt pscoast -J$proj -R$leveling_range -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -LjBR+c32+w3+o0.05i/0.19i -X-1.3i -Y-1.65i -Dh -K -O -P >> $output
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output
awk '{print $1, $2, $3*1000}' $obs_insar_file | gmt psxy -Sc0.16 -Clos.cpt -R -J$proj -O -K >> $output
gmt psxy $field_file -R -J$proj -Wthinnest,indianred -O -K >> $output
echo '-115.57 33.06 E / Lev' | gmt pstext -R -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
gmt pscoast -J$proj -R -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X1.3i -Dh -K -O -P >> $output
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
awk '{print $1, $2, $3*1000}' $model_insar_file | gmt psxy -Sc0.16 -Clos.cpt -R -J$proj -O -K >> $output
gmt psxy $field_file -R -J$proj -Wthinnest,indianred -O -K >> $output


files=`awk 'NR==12' "$1"`
obs_insar_file=`echo $files | awk '{print $1;}'`
model_insar_file=`echo $files | awk '{print $2;}'`
gmt pscoast -J$proj -R$leveling_range -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -LjBR+c32+w3+o0.05i/0.19i -X-1.3i -Y-1.65i -Dh -K -O -P >> $output
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output
awk '{print $1, $2, $3*1000}' $obs_insar_file | gmt psxy -Sc0.16 -Clos.cpt -R -J$proj -O -K >> $output
gmt psxy $field_file -R -J$proj -Wthinnest,indianred -O -K >> $output
echo '-115.57 33.06 F / Lev' | gmt pstext -R -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
gmt pscoast -J$proj -R -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X1.3i -Dh -K -O -P >> $output
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
awk '{print $1, $2, $3*1000}' $model_insar_file | gmt psxy -Sc0.16 -Clos.cpt -R -J$proj -O -K >> $output
gmt psxy $field_file -R -J$proj -Wthinnest,indianred -O -K >> $output

gmt psscale -Clos.cpt -Dx-0.5c/-0.5c+w5c/0.25c+jTC+h -Bxaf -By+l"mm" -O >> $output

gmt psconvert -Tg $output

rm gmt.history

open $output

