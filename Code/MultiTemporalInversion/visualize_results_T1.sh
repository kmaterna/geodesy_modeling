#!/bin/bash

ifile=$1
horiz_scale=$2
proj="M2.3i"
output="output/total_output.ps"

gmt set MAP_FRAME_TYPE plain
gmt set FORMAT_GEO_MAP D

gmt makecpt -T-10/10/0.5 -Cpolar -D > mycpt.cpt
gmt makecpt -T-40/40/1 -Cjet -D > los.cpt
gmt makecpt -T-.15/.15/0.001 -D -Cwysiwyg > slip.cpt

# Get files
files=`awk 'NR==5' "$1"`
ss_fault_file=`echo $files | awk '{print $1;}'`
thrust_file=`echo $files | awk '{print $2;}'`
field_file="../Mapping_Data/Fields_Boundaries.txt"

# FAULT PLOTS
range="-115.69/-115.4267/32.93/33.135"
gmt pscoast -J$proj -R$range -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X0.7i -Y9i -Dh -K -P > $output
gmt psxy $thrust_file -R -J$proj -W0.01p,gray -Cslip.cpt -L -K -O >> $output
echo "-115.55 33.12 Model Reverse Slip" | gmt pstext -R -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output
gmt pscoast -J$proj -R -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X2.35i -Y0i -Dh -K -O -P >> $output
gmt psxy $ss_fault_file -R -J$proj -W0.01p,gray -Cslip.cpt -L -K -O >> $output
echo "-115.55 33.12 Model Strike Slip (LL)" | gmt pstext -R -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output
gmt psscale -Cslip.cpt -Dx-0.5c/-0.5c+w5c/0.25c+jTC+h -Bxaf -By+l"m" -K -O >> $output



# GPS PLOTS
bigrange="-115.8/-115.35/32.84/33.185"
gmt pscoast -J$proj -R$bigrange -Gwhite -Slightgray -N1 -Wthin,black -B0.2wEsN -Y0i -X2.35i -Dh -K -O -P >> $output
files=`awk 'NR==1' "$1"`
observed_gps_file=`echo $files | awk '{print $1;}'`
predicted_gps_file=`echo $files | awk '{print $2;}'`
label=`echo $files | awk '{print $3, $4, $5;}'`
gmt psxy $thrust_file -R -J$proj -Wthinnest,gray -L -K -O >> $output
awk '{print $1, $2, $5*1000}' $observed_gps_file | gmt psxy -Sh0.35 -Cmycpt.cpt -W0.1,black -R -J$proj -O -K >> $output
awk '{print $1, $2, $5*1000}' $predicted_gps_file | gmt psxy -Sh0.2 -Cmycpt.cpt -W0.1,red -R -J$proj -O -K >> $output
awk '{print $1, $2, $3, $4, $6, $7, 0}'  $observed_gps_file | gmt psvelo -A+e+gblack+pthickest -Wblack -Se$horiz_scale/0.68/10 -R -J$proj -O -K >> $output
awk '{print $1, $2, $3, $4, $6, $7, 0}' $predicted_gps_file | gmt psvelo -A+e+gred+pthickest -Wred -Se$horiz_scale/0.68/10 -R -J$proj -O -K >> $output
echo "-115.47 33.16 GNSS" | gmt pstext -R -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
gmt psscale -Cmycpt.cpt -Dx3.0c/-0.5c+w5c/0.25c+jTC+h -Bxaf -By+l"mm" -K -O >> $output
gmt psvelo -A+e+gred+pthickest -Wred -Se$horiz_scale/0.68/10 -R -J$proj -O -K <<EOF>> $output
`echo -115.73 + 0.11 | bc` `echo 32.9256 + 0.00 | bc` .005 0 0 0 0 5mm model
EOF
gmt psvelo -A+e+gblack+pthickest -Wblack -Se$horiz_scale/0.68/10 -R -J$proj -O -K <<EOF>> $output
`echo -115.73 + 0.11 | bc` `echo 32.9256 + 0.02 | bc` .005 0 0.001 0.001 0 5mm obs
EOF
gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output



# LEVELING, ETC. 
smallrange="-115.63/-115.46/32.97/33.08"
files=`awk 'NR==2' "$1"`
obs_insar_file=`echo $files | awk '{print $1;}'`
model_insar_file=`echo $files | awk '{print $2;}'`
label=`echo $files | awk '{print $3, $4, $5;}'`
gmt pscoast -J$proj -R$smallrange -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X-4.7i -Y-2.9i -Dh -K -O -P >> $output
# echo "-115.51 33.07 Leveling Data" | gmt pstext -R -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
echo "-115.645 33.03 Leveling" | gmt pstext -R -J$projection -F+a90+f18p,Helvetica-bold -K -N -O -P >> $output
echo "-115.54 33.09 Observations" | gmt pstext -R -J$projection -F+f18p,Helvetica-bold -N -K -O -P >> $output
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
awk '{print $1, $2, $3*1000}' $obs_insar_file | gmt psxy -Sc0.16 -Clos.cpt -R -J$proj -O -K >> $output
gmt psxy $field_file -R -J$proj -Wthinnest,indianred -O -K >> $output
gmt pscoast -J$proj -R -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X2.35i -Dh -K -O -P >> $output
# echo "-115.51 33.07 Leveling Model" | gmt pstext -R -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
echo "-115.54 33.09 Models" | gmt pstext -R -J$projection -F+f18p,Helvetica-bold -N -K -O -P >> $output
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
awk '{print $1, $2, $3*1000}' $model_insar_file | gmt psxy -Sc0.16 -Clos.cpt -R -J$proj -O -K >> $output
paste $obs_insar_file $model_insar_file | awk '{print $1, $2, $3*1000, $10*1000}'  > temp_lev.txt
gmt psxy $field_file -R -J$proj -Wthinnest,indianred -O -K >> $output
gmt pscoast -J$proj -R -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X2.35i -Dh -K -O -P >> $output
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
awk '{print $1, $2, $3-$4}' temp_lev.txt | gmt psxy -Sc0.16 -Clos.cpt -R -J$proj -O -K >> $output
# echo "-115.51 33.07 Leveling Residual" | gmt pstext -R -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
gmt psxy $field_file -R -J$proj -Wthinnest,indianred -O -K >> $output
echo "-115.54 33.09 Residuals" | gmt pstext -R -J$projection -F+f18p,Helvetica-bold -N -K -O -P >> $output



# UAVSAR 26509
range="-115.73/-115.44/32.95/33.13"
files=`awk 'NR==3' "$1"`
obs_insar_file=`echo $files | awk '{print $1;}'`
model_insar_file=`echo $files | awk '{print $2;}'`
label=`echo $files | awk '{print $3, $4, $5;}'`
gmt pscoast -J$proj -R$range -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X-4.7i -Y-1.9i -Dh -K -O -P >> $output
echo "-115.78 33.04 UAVSAR" | gmt pstext -R -J$projection -F+a90+f17p,Helvetica-bold -K -N -O -P >> $output
echo "-115.75 33.04 T26509" | gmt pstext -R -J$projection -F+a90+f17p,Helvetica-bold -K -N -O -P >> $output
# echo "-115.51 33.14 UAVSAR Data" | gmt pstext -R -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
echo "-115.47 33.12 -0.99 -0.09 0 0 0" | gmt psvelo -R -J -Se1/0.95/0 -A0.3+a30+ea+gblack+p-2p,black -K -O >> $output  # Flight vector
echo "-115.47 33.12 0.05 -0.49 0 0 0" | gmt psvelo -R -J -Se1/0.95/0 -A0.5+a30+ea+gblack+p-3p,black -K -O >> $output  # LOS
echo "-115.47 33.092 LOS" | gmt pstext -R -J$projection -F+f11p,Helvetica-bold -K -O -P >> $output
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
awk '{print $1, $2, $3*1000}' $obs_insar_file | gmt psxy -Sc0.11 -Clos.cpt -R -J$proj -O -K >> $output
gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output
gmt pscoast -J$proj -R -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X2.35i -Dh -K -O -P >> $output
# echo "-115.52 33.14 UAVSAR Model" | gmt pstext -R -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
echo "-115.47 33.12 -0.99 -0.09 0 0 0" | gmt psvelo -R -J -Se1/0.95/0 -A0.3+a30+ea+gblack+p-2p,black -K -O >> $output  # Flight vector
echo "-115.47 33.12 0.05 -0.49 0 0 0" | gmt psvelo -R -J -Se1/0.95/0 -A0.5+a30+ea+gblack+p-3p,black -K -O >> $output  # LOS
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
awk '{print $1, $2, $3*1000}' $model_insar_file | gmt psxy -Sc0.11 -Clos.cpt -R -J$proj -O -K >> $output
gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output
paste $obs_insar_file $model_insar_file | awk '{print $1, $2, $3*1000, $10*1000}'  > temp_insar.txt
gmt pscoast -J$proj -R -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X2.35i -Dh -K -O -P >> $output
# echo "-115.53 33.14 UAVSAR Residual" | gmt pstext -R -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
echo "-115.47 33.12 -0.99 -0.09 0 0 0" | gmt psvelo -R -J -Se1/0.95/0 -A0.3+a30+ea+gblack+p-2p,black -K -O >> $output  # Flight vector
echo "-115.47 33.12 0.05 -0.49 0 0 0" | gmt psvelo -R -J -Se1/0.95/0 -A0.5+a30+ea+gblack+p-3p,black -K -O >> $output  # LOS
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
awk '{print $1, $2, $3-$4}' temp_insar.txt | gmt psxy -Sc0.11 -Clos.cpt -R -J$proj -O -K >> $output
gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output

# UAVSAR 08508
range="-115.71/-115.42/32.91/33.09"
files=`awk 'NR==4' "$1"`
obs_insar_file=`echo $files | awk '{print $1;}'`
model_insar_file=`echo $files | awk '{print $2;}'`
label=`echo $files | awk '{print $3, $4, $5;}'`
gmt pscoast -J$proj -R$range -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X-4.7i -Y-1.9i -Dh -K -O -P >> $output
# echo "-115.51 33.14 UAVSAR Data" | gmt pstext -R -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
echo "-115.76 33.01 UAVSAR" | gmt pstext -R -J$projection -F+a90+f17p,Helvetica-bold -K -N -O -P >> $output
echo "-115.73 33.01 T08508" | gmt pstext -R -J$projection -F+a90+f17p,Helvetica-bold -K -N -O -P >> $output
echo "-115.68 32.92 0.99 0.09 0 0 0" | gmt psvelo -R -J -Se1/0.95/0 -A0.3+a30+ea+gblack+p-2p,black -K -O >> $output  # Flight vector
echo "-115.68 32.92 -0.05 0.49 0 0 0" | gmt psvelo -R -J -Se1/0.95/0 -A0.5+a30+ea+gblack+p-3p,black -K -O >> $output  # LOS
echo "-115.68 32.95 LOS" | gmt pstext -R -J$projection -F+f11p,Helvetica-bold -K -O -P >> $output
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
awk '{print $1, $2, $3*1000}' $obs_insar_file | gmt psxy -Sc0.11 -Clos.cpt -R -J$proj -O -K >> $output
gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output
gmt pscoast -J$proj -R -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X2.35i -Dh -K -O -P >> $output
# echo "-115.52 33.14 UAVSAR Model" | gmt pstext -R -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
echo "-115.68 32.92 0.99 0.09 0 0 0" | gmt psvelo -R -J -Se1/0.95/0 -A0.3+a30+ea+gblack+p-2p,black -K -O >> $output  # Flight vector
echo "-115.68 32.92 -0.05 0.49 0 0 0" | gmt psvelo -R -J -Se1/0.95/0 -A0.5+a30+ea+gblack+p-3p,black -K -O >> $output  # LOS
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
awk '{print $1, $2, $3*1000}' $model_insar_file | gmt psxy -Sc0.11 -Clos.cpt -R -J$proj -O -K >> $output
gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output
paste $obs_insar_file $model_insar_file | awk '{print $1, $2, $3*1000, $10*1000}'  > temp_insar.txt
gmt pscoast -J$proj -R -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X2.35i -Dh -K -O -P >> $output
# echo "-115.53 33.14 UAVSAR Residual" | gmt pstext -R -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
echo "-115.68 32.92 0.99 0.09 0 0 0" | gmt psvelo -R -J -Se1/0.95/0 -A0.3+a30+ea+gblack+p-2p,black -K -O >> $output  # Flight vector
echo "-115.68 32.92 -0.05 0.49 0 0 0" | gmt psvelo -R -J -Se1/0.95/0 -A0.5+a30+ea+gblack+p-3p,black -K -O >> $output  # LOS
gmt psxy $ss_fault_file -R -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
awk '{print $1, $2, $3-$4}' temp_insar.txt | gmt psxy -Sc0.11 -Clos.cpt -R -J$proj -O -K >> $output
gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output

gmt psscale -Clos.cpt -Dx-3.2c/-0.5c+w5c/0.25c+jTC+h -Bxaf -By+l"mm" -O >> $output



rm gmt.history
rm temp_lev.txt
rm temp_insar.txt
rm *.cpt
rm gmt.conf

gmt psconvert -Tg $output

open $output

