#!/bin/bash
# GMT Plot for A/B/C time intervals (fault slip)

slip_file1=$1
slip_file2=$2
slip_file3=$3
slip_file4=$4
slip_file5=$5
slip_file6=$6
scale_low=$7
scale_high=$8
scale_int=$9
outdir=${10}
title=${11}

range="-115.62/-115.46/32.97/33.08"
proj="M2.7i"
output=$outdir/$title"_ABCDEF.ps"
field_file="../Mapping_Data/Fields_Boundaries.txt"

gmt makecpt -T$scale_low/$scale_high/$scale_int -Cpolar -D > mycpt.cpt

# Plotting the slip distribution
gmt pscoast -R$range -J$proj -Wthin,black -Di -N1 -N2 -B0.2WesN -Y4.5i -Gwhite -Slightgray -K > $output
echo "-115.59 33.07 T4-A" | gmt pstext -F+f18p,Helvetica-bold -R -J -K -O >> $output
gmt psxy $slip_file1 -R -J -Wthinnest,gray -Cmycpt.cpt -L -K -O >> $output
gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output

gmt pscoast -R$range -J$proj -Wthin,black -Di -N1 -N2 -B0.2wesn -X3.0i -Gwhite -Slightgray -K -O >> $output
gmt psxy $slip_file2 -R -J -Wthinnest,gray -Cmycpt.cpt -L -K -O >> $output
echo "-115.59 33.07 T4-B" | gmt pstext -F+f18p,Helvetica-bold -R -J -K -O >> $output
gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output

gmt pscoast -R$range -J$proj -Wthin,black -Di -N1 -N2 -B0.2wesN -X3.0i -Gwhite -Slightgray -K -O >> $output
echo "-115.59 33.07 T4-C" | gmt pstext -F+f18p,Helvetica-bold -R -J -K -O >> $output
gmt psxy $slip_file3 -R -J -Wthinnest,gray -Cmycpt.cpt -L -K -O >> $output
gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output

gmt pscoast -R$range -J$proj -Wthin,black -Di -N1 -N2 -B0.2WeSn -X-6.0i -Y-2.5i -Gwhite -Slightgray -K -O >> $output
echo "-115.59 33.07 T4-D" | gmt pstext -F+f18p,Helvetica-bold -R -J -K -O >> $output
gmt psxy $slip_file4 -R -J -Wthinnest,gray -Cmycpt.cpt -L -K -O >> $output
gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output

gmt pscoast -R$range -J$proj -Wthin,black -Di -N1 -N2 -B0.2weSn -X3.0i -Gwhite -Slightgray -K -O >> $output
echo "-115.59 33.07 T4-E" | gmt pstext -F+f18p,Helvetica-bold -R -J -K -O >> $output
gmt psxy $slip_file5 -R -J -Wthinnest,gray -Cmycpt.cpt -L -K -O >> $output
gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output

gmt pscoast -R$range -J$proj -Wthin,black -Di -N1 -N2 -B0.2weSn -X3.0i -Gwhite -Slightgray -K -O >> $output
echo "-115.59 33.07 T4-F" | gmt pstext -F+f18p,Helvetica-bold -R -J -K -O >> $output
gmt psxy $slip_file6 -R -J -Wthinnest,gray -Cmycpt.cpt -L -K -O >> $output
gmt psxy $field_file -R -J$proj -Wthin,indianred -O -K >> $output

gmt psscale -Cmycpt.cpt -Dx-4c/-0.9c+w12c/0.5c+jTC+h -Bxaf -By+l"slip(m)" -K -O >> $output

# open $output