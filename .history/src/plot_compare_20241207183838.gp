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

# Define manual x-offsets for side-by-side plotting
plot 'data_mpi.csv' using (2*${0}-0.4):3 title "MPI Computation Time" lc rgb "blue", \
     '' using (2*${0}-0.4):2 title "MPI Communication Time" lc rgb "red", \
     'data_openacc.csv' using (2*$0+0.4):3 title "OpenACC Computation Time" lc rgb "green", \
     '' using (2*${0}+0.4):2 title "OpenACC Communication Time" lc rgb "orange"
