#!/bin/bash

ifile=$1
ll_lon=$2
ll_lat=$3
ur_lon=$4
ur_lat=$5
proj="M2.3i"
horiz_scale=$6
output="output/total_output.ps"

gmt set MAP_FRAME_TYPE plain
gmt set FORMAT_GEO_MAP D

gmt makecpt -T-10/10/0.5 -Cpolar -D > mycpt.cpt
gmt makecpt -T-40/40/1 -Cjet -D > los.cpt
# gmt makecpt -T-24/24/1 -Cjet -D > los.cpt
gmt makecpt -T-.15/.15/0.001 -D -Cwysiwyg > slip.cpt

# Get files
files=`awk 'NR==5' "$1"`
ss_fault_file=`echo $files | awk '{print $1;}'`
thrust_file=`echo $files | awk '{print $2;}'`

# Set other ranges
range="$ll_lon/$ur_lon/$ll_lat/$ur_lat"
smallrange="-115.63/-115.45/32.96/33.08"
bigrange="-115.8/-115.35/32.84/33.185"


# GPS PLOTS
gmt pscoast -J$proj -R$bigrange -Gwhite -Slightgray -N1 -Wthin,black -B0.2WesN -Y9i -X0.7i -Dh -K -P > $output
files=`awk 'NR==1' "$1"`
observed_gps_file=`echo $files | awk '{print $1;}'`
predicted_gps_file=`echo $files | awk '{print $2;}'`
label=`echo $files | awk '{print $3, $4, $5;}'`
gmt psxy $thrust_file -R$bigrange -J$proj -Wthinnest,gray -L -K -O >> $output
awk '{print $1, $2, $5*1000}' $observed_gps_file | gmt psxy -Sh0.35 -Cmycpt.cpt -W0.1,black -R$bigrange -J$proj -O -K >> $output
awk '{print $1, $2, $5*1000}' $predicted_gps_file | gmt psxy -Sh0.2 -Cmycpt.cpt -W0.1,red -R$bigrange -J$proj -O -K >> $output
awk '{print $1, $2, $3, $4, $6, $7, 0}'  $observed_gps_file | gmt psvelo -A+e+gblack+pthickest -Wblack -Se$horiz_scale/0.68/10 -R$bigrange -J$proj -O -K >> $output
awk '{print $1, $2, $3, $4, $6, $7, 0}' $predicted_gps_file | gmt psvelo -A+e+gred+pthickest -Wred -Se$horiz_scale/0.68/10 -R$bigrange -J$proj -O -K >> $output
echo "-115.47 33.16 GPS" | gmt pstext -R$bigrange -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
gmt psscale -Cmycpt.cpt -Dx3.0c/-0.5c+w5c/0.25c+jTC+h -Bxaf -By+l"mm" -K -O >> $output
gmt psvelo -A+e+gred+pthickest -Wred -Se$horiz_scale/0.68/10 -R$bigrange -J$proj -O -K <<EOF>> $output
`echo $ll_lon + 0.11 | bc` `echo $ll_lat + 0.00 | bc` .005 0 0 0 0 5mm model
EOF
gmt psvelo -A+e+gblack+pthickest -Wblack -Se$horiz_scale/0.68/10 -R$bigrange -J$proj -O -K <<EOF>> $output
`echo $ll_lon + 0.11 | bc` `echo $ll_lat + 0.02 | bc` .005 0 0.001 0.001 0 5mm obs
EOF

# FAULT PLOTS
gmt pscoast -J$proj -R$range -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X2.35i -Y0i -Dh -K -O -P >> $output
gmt psxy $ss_fault_file -R$range -J$proj -W0.01p,gray -Cslip.cpt -L -K -O >> $output
echo "-115.52 33.14 Strike Slip (LL)" | gmt pstext -R$range -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
gmt pscoast -J$proj -R$range -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X2.35i -Y0i -Dh -K -O -P >> $output
gmt psxy $thrust_file -R$range -J$proj -W0.01p,gray -Cslip.cpt -L -K -O >> $output
echo "-115.52 33.14 Reverse Slip" | gmt pstext -R$range -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
gmt psscale -Cslip.cpt -Dx-0.5c/-0.5c+w5c/0.25c+jTC+h -Bxaf -By+l"m" -K -O >> $output




# LEVELING, ETC. 
files=`awk 'NR==2' "$1"`
obs_insar_file=`echo $files | awk '{print $1;}'`
model_insar_file=`echo $files | awk '{print $2;}'`
label=`echo $files | awk '{print $3, $4, $5;}'`
gmt pscoast -J$proj -R$smallrange -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X-4.7i -Y-2.6i -Dh -K -O -P >> $output
echo "-115.51 33.07 Leveling Data" | gmt pstext -R$smallrange -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
gmt psxy $ss_fault_file -R$smallrange -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
awk '{print $1, $2, $3*1000}' $obs_insar_file | gmt psxy -Sc0.12 -Clos.cpt -R$smallrange -J$proj -O -K >> $output
gmt pscoast -J$proj -R$smallrange -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X2.35i -Dh -K -O -P >> $output
echo "-115.51 33.07 Leveling Model" | gmt pstext -R$smallrange -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
gmt psxy $ss_fault_file -R$smallrange -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
awk '{print $1, $2, $3*1000}' $model_insar_file | gmt psxy -Sc0.12 -Clos.cpt -R$smallrange -J$proj -O -K >> $output
paste $obs_insar_file $model_insar_file | awk '{print $1, $2, $3*1000, $10*1000}'  > temp_lev.txt
gmt pscoast -J$proj -R$smallrange -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X2.35i -Dh -K -O -P >> $output
gmt psxy $ss_fault_file -R$smallrange -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
awk '{print $1, $2, $3-$4}' temp_lev.txt | gmt psxy -Sc0.12 -Clos.cpt -R$smallrange -J$proj -O -K >> $output
echo "-115.51 33.07 Leveling Residual" | gmt pstext -R$smallrange -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output



# UAVSAR
files=`awk 'NR==3' "$1"`
obs_insar_file=`echo $files | awk '{print $1;}'`
model_insar_file=`echo $files | awk '{print $2;}'`
label=`echo $files | awk '{print $3, $4, $5;}'`
gmt pscoast -J$proj -R$range -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X-4.7i -Y-2.3i -Dh -K -O -P >> $output
echo "-115.51 33.14 UAVSAR Data" | gmt pstext -R$range -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
gmt psxy $ss_fault_file -R$range -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
awk '{print $1, $2, $3*1000}' $obs_insar_file | gmt psxy -Sc0.07 -Clos.cpt -R$range -J$proj -O -K >> $output
gmt pscoast -J$proj -R$range -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X2.35i -Dh -K -O -P >> $output
echo "-115.52 33.14 UAVSAR Model" | gmt pstext -R$range -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
gmt psxy $ss_fault_file -R$range -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
awk '{print $1, $2, $3*1000}' $model_insar_file | gmt psxy -Sc0.07 -Clos.cpt -R$range -J$proj -O -K >> $output
paste $obs_insar_file $model_insar_file | awk '{print $1, $2, $3*1000, $10*1000}'  > temp_insar.txt
gmt pscoast -J$proj -R$range -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X2.35i -Dh -K -O -P >> $output
echo "-115.53 33.14 UAVSAR Residual" | gmt pstext -R$range -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
gmt psxy $ss_fault_file -R$range -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
awk '{print $1, $2, $3-$4}' temp_insar.txt | gmt psxy -Sc0.07 -Clos.cpt -R$range -J$proj -O -K >> $output


# UAVSAR
files=`awk 'NR==4' "$1"`
obs_insar_file=`echo $files | awk '{print $1;}'`
model_insar_file=`echo $files | awk '{print $2;}'`
label=`echo $files | awk '{print $3, $4, $5;}'`
gmt pscoast -J$proj -R$range -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X-4.7i -Y-2.3i -Dh -K -O -P >> $output
echo "-115.51 33.14 UAVSAR Data" | gmt pstext -R$range -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
gmt psxy $ss_fault_file -R$range -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
awk '{print $1, $2, $3*1000}' $obs_insar_file | gmt psxy -Sc0.07 -Clos.cpt -R$range -J$proj -O -K >> $output
gmt pscoast -J$proj -R$range -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X2.35i -Dh -K -O -P >> $output
echo "-115.52 33.14 UAVSAR Model" | gmt pstext -R$range -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
gmt psxy $ss_fault_file -R$range -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
awk '{print $1, $2, $3*1000}' $model_insar_file | gmt psxy -Sc0.07 -Clos.cpt -R$range -J$proj -O -K >> $output
paste $obs_insar_file $model_insar_file | awk '{print $1, $2, $3*1000, $10*1000}'  > temp_insar.txt
gmt pscoast -J$proj -R$range -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X2.35i -Dh -K -O -P >> $output
echo "-115.53 33.14 UAVSAR Residual" | gmt pstext -R$range -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output
gmt psxy $ss_fault_file -R$range -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks
awk '{print $1, $2, $3-$4}' temp_insar.txt | gmt psxy -Sc0.07 -Clos.cpt -R$range -J$proj -O -K >> $output

gmt psscale -Clos.cpt -Dx-3.2c/-0.5c+w5c/0.25c+jTC+h -Bxaf -By+l"mm" -K -O >> $output



rm gmt.history
rm temp_lev.txt
rm temp_insar.txt
rm *.cpt
rm gmt.conf

open $output

