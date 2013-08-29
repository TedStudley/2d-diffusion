./../build/2d-diffusion

echo 'set terminal png' >> my_plot.gps

for i in *.dat
do
    filename=`echo $i | sed 's/dat/png/'`

    echo "set output '$filename'" >> my_plot.gps
    echo "set xrange [0:32]" >> my_plot.gps
    echo "set yrange [0:32]" >> my_plot.gps
    echo "set zrange [-0.1:1.1]" >> my_plot.gps
    echo "set hidden3d front" >> my_plot.gps
    echo "splot '$i' matrix with lines" >> my_plot.gps
    echo "reset" >> my_plot.gps
done

gnuplot my_plot.gps

rm my_plot.gps

convert -delay 50 -loop 0 *.png ani.gif
