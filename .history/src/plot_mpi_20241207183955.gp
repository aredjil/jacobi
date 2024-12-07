
# set terminal pngcairo size 800,600 enhanced font 'Verdana,10'
# set output 'rowstacked_histogram.png'

# # Set titles and labels
# set title "Communication vs Computation Times, for a non-blocking mpi jacobi solver of 12K mesh size for 10 iterations"
# set xlabel "Number of Processes"
# set ylabel "Time (ms)"

# # Data is in 'data_mpi.csv'
# set datafile separator ','

# # Configure histogram and bar styles
# set style data histogram
# set style histogram rowstacked
# set style fill solid border -1
# set boxwidth 0.8

# # Configure key (legend) in the top right corner
# set key top right

# # Plot the data (reverse the order so communication is on top)
# plot 'data_mpi.csv' using 3:xtic(1) title "Computation Time" lc rgb "blue", \
#      '' using 2 title "Communication Time" lc rgb "red"

