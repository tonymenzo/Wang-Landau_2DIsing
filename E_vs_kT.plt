# file: E_vs_kT.plt
#
# gnuplot plotfile for Wang-Landau Sampling of 2D Ising model
#  
#  Programmer:  Tony Menzo
# 
#

# record the time and date the graph was generated
set timestamp

# titles and labels
set title "Density of States 2D Ising Model: Average Energy vs kT"
set xlabel "kT"
set ylabel "Energy"



# set the terminal type to be the screen (which is x11 here)
set term x11 

# plot


plot \
     "g_E_fit_Z.dat" using 1:2 with linespoints title 'Energy'
     


# output the plot to the file E_kT_plt.ps   
set out "E_vs_kT_plt.ps"
set term postscript color 
replot

reset
set term x11
