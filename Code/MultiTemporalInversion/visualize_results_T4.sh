#!/bin/bash

ifile=$1
ll_lon=$2
ll_lat=$3
ur_lon=$4
ur_lat=$5
proj=$6
range="$ll_lon/$ur_lon/$ll_lat/$ur_lat"
horiz_scale=100
output="total_output.ps"

gmt set MAP_FRAME_TYPE plain
gmt set FORMAT_GEO_MAP D

gmt makecpt -T-10/10/0.5 -Cpolar -D > mycpt.cpt
gmt makecpt -T-30/30/1 -Cjet -D > los.cpt

# GPS PLOTS
gmt pscoast -J$proj -R$range -Gwhite -Slightgray -N1 -Wthin,black -B0.2WesN -Y9i -X0.5i -Dh -K -P > $output
files=`awk 'NR==1' "$1"`
observed_gps_file=`echo $files | awk '{print $1;}'`
predicted_gps_file=`echo $files | awk '{print $2;}'`
label=`echo $files | awk '{print $3, $4, $5;}'`
awk '{print $1, $2, $5*1000}' $observed_gps_file | gmt psxy -Sh0.35 -Cmycpt.cpt -W0.1,black -R$range -J$proj -O -K >> $output
awk '{print $1, $2, $5*1000}' $predicted_gps_file | gmt psxy -Sh0.2 -Cmycpt.cpt -W0.1,red -R$range -J$proj -O -K >> $output
awk '{print $1, $2, $3, $4, $6, $7, 0}'  $observed_gps_file | gmt psvelo -A+e+gblack+pthickest -Wblack -Se$horiz_scale/0.68/10 -R$range -J$proj -O -K >> $output
awk '{print $1, $2, $3, $4, $6, $7, 0}' $predicted_gps_file | gmt psvelo -A+e+gred+pthickest -Wred -Se$horiz_scale/0.68/10 -R$range -J$proj -O -K >> $output
echo $label | gmt pstext -R$range -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output

gmt pscoast -J$proj -R$range -Gwhite -Slightgray -N1 -Wthin,black -B0.2Wesn -Y-1.65i -Dh -K -O -P >> $output
files=`awk 'NR==2' "$1"`
observed_gps_file=`echo $files | awk '{print $1;}'`
predicted_gps_file=`echo $files | awk '{print $2;}'`
label=`echo $files | awk '{print $3, $4, $5;}'`
awk '{print $1, $2, $5*1000}' $observed_gps_file | gmt psxy -Sh0.35 -Cmycpt.cpt -W0.1,black -R$range -J$proj -O -K >> $output
awk '{print $1, $2, $5*1000}' $predicted_gps_file | gmt psxy -Sh0.2 -Cmycpt.cpt -W0.1,red -R$range -J$proj -O -K >> $output
awk '{print $1, $2, $3, $4, $6, $7, 0}'  $observed_gps_file | gmt psvelo -A+e+gblack+pthickest -Wblack -Se$horiz_scale/0.68/10 -R$range -J$proj -O -K >> $output
awk '{print $1, $2, $3, $4, $6, $7, 0}' $predicted_gps_file | gmt psvelo -A+e+gred+pthickest -Wred -Se$horiz_scale/0.68/10 -R$range -J$proj -O -K >> $output
echo $label | gmt pstext -R$range -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output

gmt pscoast -J$proj -R$range -Gwhite -Slightgray -N1 -Wthin,black -B0.2Wesn -Y-1.65i -Dh -K -O -P >> $output
files=`awk 'NR==3' "$1"`
observed_gps_file=`echo $files | awk '{print $1;}'`
predicted_gps_file=`echo $files | awk '{print $2;}'`
label=`echo $files | awk '{print $3, $4, $5;}'`
awk '{print $1, $2, $5*1000}' $observed_gps_file | gmt psxy -Sh0.35 -Cmycpt.cpt -W0.1,black -R$range -J$proj -O -K >> $output
awk '{print $1, $2, $5*1000}' $predicted_gps_file | gmt psxy -Sh0.2 -Cmycpt.cpt -W0.1,red -R$range -J$proj -O -K >> $output
awk '{print $1, $2, $3, $4, $6, $7, 0}'  $observed_gps_file | gmt psvelo -A+e+gblack+pthickest -Wblack -Se$horiz_scale/0.68/10 -R$range -J$proj -O -K >> $output
awk '{print $1, $2, $3, $4, $6, $7, 0}' $predicted_gps_file | gmt psvelo -A+e+gred+pthickest -Wred -Se$horiz_scale/0.68/10 -R$range -J$proj -O -K >> $output
echo $label | gmt pstext -R$range -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output

gmt pscoast -J$proj -R$range -Gwhite -Slightgray -N1 -Wthin,black -B0.2Wesn -Y-1.65i -Dh -K -O -P >> $output
files=`awk 'NR==4' "$1"`
observed_gps_file=`echo $files | awk '{print $1;}'`
predicted_gps_file=`echo $files | awk '{print $2;}'`
label=`echo $files | awk '{print $3, $4, $5;}'`
awk '{print $1, $2, $5*1000}' $observed_gps_file | gmt psxy -Sh0.35 -Cmycpt.cpt -W0.1,black -R$range -J$proj -O -K >> $output
awk '{print $1, $2, $5*1000}' $predicted_gps_file | gmt psxy -Sh0.2 -Cmycpt.cpt -W0.1,red -R$range -J$proj -O -K >> $output
awk '{print $1, $2, $3, $4, $6, $7, 0}'  $observed_gps_file | gmt psvelo -A+e+gblack+pthickest -Wblack -Se$horiz_scale/0.68/10 -R$range -J$proj -O -K >> $output
awk '{print $1, $2, $3, $4, $6, $7, 0}' $predicted_gps_file | gmt psvelo -A+e+gred+pthickest -Wred -Se$horiz_scale/0.68/10 -R$range -J$proj -O -K >> $output
echo $label | gmt pstext -R$range -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output

gmt pscoast -J$proj -R$range -Gwhite -Slightgray -N1 -Wthin,black -B0.2Wesn -Y-1.65i -Dh -K -O -P >> $output
files=`awk 'NR==5' "$1"`
observed_gps_file=`echo $files | awk '{print $1;}'`
predicted_gps_file=`echo $files | awk '{print $2;}'`
label=`echo $files | awk '{print $3, $4, $5;}'`
awk '{print $1, $2, $5*1000}' $observed_gps_file | gmt psxy -Sh0.35 -Cmycpt.cpt -W0.1,black -R$range -J$proj -O -K >> $output
awk '{print $1, $2, $5*1000}' $predicted_gps_file | gmt psxy -Sh0.2 -Cmycpt.cpt -W0.1,red -R$range -J$proj -O -K >> $output
awk '{print $1, $2, $3, $4, $6, $7, 0}'  $observed_gps_file | gmt psvelo -A+e+gblack+pthickest -Wblack -Se$horiz_scale/0.68/10 -R$range -J$proj -O -K >> $output
awk '{print $1, $2, $3, $4, $6, $7, 0}' $predicted_gps_file | gmt psvelo -A+e+gred+pthickest -Wred -Se$horiz_scale/0.68/10 -R$range -J$proj -O -K >> $output
echo $label | gmt pstext -R$range -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output

gmt pscoast -J$proj -R$range -Gwhite -Slightgray -N1 -Wthin,black -B0.2Wesn -Y-1.65i -Dh -K -O -P >> $output
files=`awk 'NR==6' "$1"`
observed_gps_file=`echo $files | awk '{print $1;}'`
predicted_gps_file=`echo $files | awk '{print $2;}'`
label=`echo $files | awk '{print $3, $4, $5;}'`
awk '{print $1, $2, $5*1000}' $observed_gps_file | gmt psxy -Sh0.35 -Cmycpt.cpt -W0.1,black -R$range -J$proj -O -K >> $output
awk '{print $1, $2, $5*1000}' $predicted_gps_file | gmt psxy -Sh0.2 -Cmycpt.cpt -W0.1,red -R$range -J$proj -O -K >> $output
awk '{print $1, $2, $3, $4, $6, $7, 0}'  $observed_gps_file | gmt psvelo -A+e+gblack+pthickest -Wblack -Se$horiz_scale/0.68/10 -R$range -J$proj -O -K >> $output
awk '{print $1, $2, $3, $4, $6, $7, 0}' $predicted_gps_file | gmt psvelo -A+e+gred+pthickest -Wred -Se$horiz_scale/0.68/10 -R$range -J$proj -O -K >> $output
echo $label | gmt pstext -R$range -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output



# TSX, LEVELING, ETC. 
smallrange="-115.63/-115.45/32.94/33.11"
files=`awk 'NR==7' "$1"`
obs_insar_file=`echo $files | awk '{print $1;}'`
model_insar_file=`echo $files | awk '{print $2;}'`
label=`echo $files | awk '{print $3, $4, $5;}'`
gmt pscoast -J$proj -R$smallrange -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X1.8i -Y8.25i -Dh -K -O -P >> $output
awk '{print $1, $2, $3*1000}' $obs_insar_file | gmt psxy -Sc0.04 -Clos.cpt -R$smallrange -J$proj -O -K >> $output
gmt pscoast -J$proj -R$smallrange -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X1.3i -Dh -K -O -P >> $output
awk '{print $1, $2, $3*1000}' $model_insar_file | gmt psxy -Sc0.04 -Clos.cpt -R$smallrange -J$proj -O -K >> $output
echo $label | gmt pstext -R$smallrange -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output


files=`awk 'NR==8' "$1"`
obs_insar_file=`echo $files | awk '{print $1;}'`
model_insar_file=`echo $files | awk '{print $2;}'`
label=`echo $files | awk '{print $3, $4, $5;}'`
gmt pscoast -J$proj -R$smallrange -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X-1.3i -Y-1.65i -Dh -K -O -P >> $output
awk '{print $1, $2, $3*1000}' $obs_insar_file | gmt psxy -Sc0.09 -Clos.cpt -R$smallrange -J$proj -O -K >> $output
gmt pscoast -J$proj -R$smallrange -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X1.3i -Dh -K -O -P >> $output
awk '{print $1, $2, $3*1000}' $model_insar_file | gmt psxy -Sc0.09 -Clos.cpt -R$smallrange -J$proj -O -K >> $output
echo $label | gmt pstext -R$smallrange -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output


files=`awk 'NR==9' "$1"`
obs_insar_file=`echo $files | awk '{print $1;}'`
model_insar_file=`echo $files | awk '{print $2;}'`
label=`echo $files | awk '{print $3, $4, $5;}'`
gmt pscoast -J$proj -R$smallrange -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X-1.3i -Y-1.65i -Dh -K -O -P >> $output
awk '{print $1, $2, $3*1000}' $obs_insar_file | gmt psxy -Sc0.09 -Clos.cpt -R$smallrange -J$proj -O -K >> $output
gmt pscoast -J$proj -R$smallrange -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X1.3i -Dh -K -O -P >> $output
awk '{print $1, $2, $3*1000}' $model_insar_file | gmt psxy -Sc0.09 -Clos.cpt -R$smallrange -J$proj -O -K >> $output
echo $label | gmt pstext -R$smallrange -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output


files=`awk 'NR==10' "$1"`
obs_insar_file=`echo $files | awk '{print $1;}'`
model_insar_file=`echo $files | awk '{print $2;}'`
label=`echo $files | awk '{print $3, $4, $5;}'`
gmt pscoast -J$proj -R$smallrange -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X-1.3i -Y-1.65i -Dh -K -O -P >> $output
awk '{print $1, $2, $3*1000}' $obs_insar_file | gmt psxy -Sc0.09 -Clos.cpt -R$smallrange -J$proj -O -K >> $output
gmt pscoast -J$proj -R$smallrange -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X1.3i -Dh -K -O -P >> $output
awk '{print $1, $2, $3*1000}' $model_insar_file | gmt psxy -Sc0.09 -Clos.cpt -R$smallrange -J$proj -O -K >> $output
echo $label | gmt pstext -R$smallrange -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output


files=`awk 'NR==11' "$1"`
obs_insar_file=`echo $files | awk '{print $1;}'`
model_insar_file=`echo $files | awk '{print $2;}'`
label=`echo $files | awk '{print $3, $4, $5;}'`
gmt pscoast -J$proj -R$smallrange -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X-1.3i -Y-1.65i -Dh -K -O -P >> $output
awk '{print $1, $2, $3*1000}' $obs_insar_file | gmt psxy -Sc0.09 -Clos.cpt -R$smallrange -J$proj -O -K >> $output
gmt pscoast -J$proj -R$smallrange -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X1.3i -Dh -K -O -P >> $output
awk '{print $1, $2, $3*1000}' $model_insar_file | gmt psxy -Sc0.09 -Clos.cpt -R$smallrange -J$proj -O -K >> $output
echo $label | gmt pstext -R$smallrange -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output


files=`awk 'NR==12' "$1"`
obs_insar_file=`echo $files | awk '{print $1;}'`
model_insar_file=`echo $files | awk '{print $2;}'`
label=`echo $files | awk '{print $3, $4, $5;}'`
gmt pscoast -J$proj -R$smallrange -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X-1.3i -Y-1.65i -Dh -K -O -P >> $output
awk '{print $1, $2, $3*1000}' $obs_insar_file | gmt psxy -Sc0.09 -Clos.cpt -R$smallrange -J$proj -O -K >> $output
gmt pscoast -J$proj -R$smallrange -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X1.3i -Dh -K -O -P >> $output
awk '{print $1, $2, $3*1000}' $model_insar_file | gmt psxy -Sc0.09 -Clos.cpt -R$smallrange -J$proj -O -K >> $output
echo $label | gmt pstext -R$smallrange -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output

gmt psscale -Clos.cpt -Dx-0.5c/-0.5c+w5c/0.25c+jTC+h -Bxaf -By+l"mm" -K -O >> $output



# UAVSAR
files=`awk 'NR==13' "$1"`
obs_insar_file=`echo $files | awk '{print $1;}'`
model_insar_file=`echo $files | awk '{print $2;}'`
label=`echo $files | awk '{print $3, $4, $5;}'`
gmt pscoast -J$proj -R$range -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X1.8i -Y8.25i -Dh -K -O -P >> $output
awk '{print $1, $2, $3*1000}' $obs_insar_file | gmt psxy -Sc0.03 -Clos.cpt -R$range -J$proj -O -K >> $output
gmt pscoast -J$proj -R$range -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X1.3i -Dh -K -O -P >> $output
awk '{print $1, $2, $3*1000}' $model_insar_file | gmt psxy -Sc0.03 -Clos.cpt -R$range -J$proj -O -K >> $output
echo $label | gmt pstext -R$range -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output


files=`awk 'NR==14' "$1"`
obs_insar_file=`echo $files | awk '{print $1;}'`
model_insar_file=`echo $files | awk '{print $2;}'`
label=`echo $files | awk '{print $3, $4, $5;}'`
gmt pscoast -J$proj -R$range -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X-1.3i -Y-1.65i -Dh -K -O -P >> $output
awk '{print $1, $2, $3*1000}' $obs_insar_file | gmt psxy -Sc0.06 -Clos.cpt -R$range -J$proj -O -K >> $output
gmt pscoast -J$proj -R$range -Gwhite -Slightgray -N1 -Wthin,black -B0.2wesn -X1.3i -Dh -K -O -P >> $output
awk '{print $1, $2, $3*1000}' $model_insar_file | gmt psxy -Sc0.06 -Clos.cpt -R$range -J$proj -O -K >> $output
echo $label | gmt pstext -R$range -J$projection -F+f12p,Helvetica-bold -K -O -P >> $output





rm gmt.history

open $output

