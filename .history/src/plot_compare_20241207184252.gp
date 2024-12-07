set terminal pngcairo size 800,600 enhanced font 'Verdana,10'
set output 'side_by_side_histogram.png'

# Titles and labels
set title "MPI vs OpenACC Communication and Computation Times"
set xlabel "Number of Processes"
set ylabel "Time (ms)"
set datafile separator ','

# Configure histogram styles
set style data histogram
set style histogram cluster gap 1
set style fill solid border -1
set boxwidth 0.6
set key top right

# Plot data side by side
plot 'data_mpi.csv' using 2:xtic(1) title "MPI Communication Time" lc rgb "red", \
     '' using 3 title "MPI Computation Time" lc rgb "blue", \
     'data_openacc.csv' using 2 title "OpenACC Communication Time" lc rgb "orange", \
     '' using 3 title "OpenACC Computation Time" lc rgb "green"
