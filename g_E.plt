# file: ising_model_WL.plt
#
# gnuplot plotfile for Wang-Landau Sampling of 2D Ising model
#  
#  Programmer:  Tony Menzo
# 
#

# record the time and date the graph was generated
set timestamp

# titles and labels
set title "Density of States 2D Ising Model"
set xlabel "energy E"
set ylabel "g(E)"



f(x) = a*x**2 + b*x + c
fit f(x)  "g_E.dat" using ($1):($2) via a,b,c
fit_title1 = sprintf("%-+4.1f*x %-+4.1f",a,b,c)


# set the terminal type to be the screen (which is x11 here)
set term x11 

# plot


plot \
     "g_E.dat" using 1:2 title 'g(E)', f(x) title 'gnuplot fit'
     #"g_E_fit.dat" using 1:2 with linespoints title 'cubic spline fit', \


# output the plot to the file g_E__1plt.ps   
set out "g_E_plt.ps"
set term postscript color 
replot

reset
set term x11
