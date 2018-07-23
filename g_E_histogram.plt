# file: ising_model_WL_histogram.plt
#
# gnuplot plotfile for Wang-Landau Sampling of 2D Ising model
#  
#  Programmer:  Tony Menzo
# 
#

# record the time and date the graph was generated
set timestamp

# titles and labels
set title "Density of States 2D Ising Model:Histogram"
set xlabel "# Sampled"
set ylabel "H[E]"

set yrange [0:170000]

# set the terminal type to be the screen (which is x11 here)
set term x11 

# plot


plot \
     "g_E.dat" using 1:3 title 'H(E)'


# output the plot to the file g_E_histogram_plt.ps   
set out "g_E_histogram_plt.ps"
set term postscript color 
replot

reset
set term x11
