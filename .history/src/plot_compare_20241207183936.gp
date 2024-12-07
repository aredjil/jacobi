set terminal pngcairo size 800,600 enhanced font 'Verdana,10'
set output 'simple_comparison.png'

# Titles and labels
set title "MPI vs OpenACC Performance Comparison"
set xlabel "Number of Processes"
set ylabel "Time (ms)"
set datafile separator ','

# Histogram configuration
set style data histogram
set style histogram rowstacked
set style fill solid border -1
set boxwidth 0.4
set key top right

# Plot data with simple x-offsets
plot 'data_mpi.csv' using ($0-0.2):3 title "MPI Computation" lc rgb "blue", \
     '' using ($0-0.2):2 title "MPI Communication" lc rgb "red", \
     'data_openacc.csv' using ($0+0.2):3 title "OpenACC Computation" lc rgb "green", \
     '' using ($0+0.2):2 title "OpenACC Communication" lc rgb "orange"
