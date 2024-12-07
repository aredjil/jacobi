# plot_timing_histogram.gp
# Gnuplot script to create a histogram of Communication vs Computation Time

# Set the input file directly in the plot command
set datafile separator ","

# Set the output to an image file
set terminal pngcairo size 800,600 enhanced font "Arial,12"
set output "communication_vs_computation_histogram.png"

# Set titles and labels
set title "Communication vs Computation Time for Different Process Counts"
set xlabel "Process Count"
set ylabel "Time (ms)"

# Set the style for the histogram
set style data histograms
set style histogram cluster gap 1
set style fill solid border -1

# Define colors for communication and computation bars
set boxwidth 0.8 relative
set key outside top center

# Plot the data from timing_results.csv
plot "timing_results.csv" using 2:xtic(1) title "Computation Time (ms)", \
     "timing_results.csv" using 3 title "Communication Time (ms)"

