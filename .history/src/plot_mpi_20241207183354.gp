set terminal pngcairo size 800,600 enhanced font 'Verdana,10'
set output 'comparison_histogram.png'

# Set titles and labels
set title "Communication vs Computation Times: MPI vs OpenACC for a 12K Mesh Size (10 Iterations)"
set xlabel "Number of Processes"
set ylabel "Time (ms)"

# Data is in 'data_mpi.csv' and 'data_openacc.csv'
set datafile separator ','

# Configure histogram and bar styles
set style data histogram
set style histogram rowstacked
set style fill solid border -1
set boxwidth 0.4

# Configure key (legend) in the top right corner
set key top right

# Define the two datasets as blocks for plotting side by side
plot 'data_mpi.csv' using (column(0)*2):3 title "MPI Computation Time" lc rgb "blue", \
     '' using (column(0)*2):2 title "MPI Communication Time" lc rgb "red", \
     'data_openacc.csv' using (column(0)*2+1):3 title "OpenACC Computation Time" lc rgb "green", \
     '' using (column(0)*2+1):2 title "OpenACC Communication Time" lc rgb "orange"

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

