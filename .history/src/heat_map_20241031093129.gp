set terminal gif animate delay 5  loop 0 # Set the terminal to animated GIF with 5 ms delay between frames
set output '../gif/jacobi.gif'           # Set the output file for the GIF
set title "Laplace Equation solution using Jacobi Iteration Technique"
set size square                   # Ensure square aspect ratio for the plots
unset colorbox                    # Hide the color bar
set palette rgb 33,13,10          # Use a custom color palette
set yrange [] reverse             # Reverse the y-axis so that the image is not flipped

# Loop over the files and plot each one
do for [i=1:1000] {                # Assuming 100 files named like file_001.txt, file_002.txt, etc.
    file = sprintf('../data/file_%03d.dat', i)  # Format the file name with leading zeros
    plot file matrix with image   # Plot the file as a matrix with an image
}

reset                             # Reset Gnuplot settings to defaults
set output                        # Close the output GIF
set terminal qt                   # Set terminal back to normal (e.g., qt or appropriate for your system)
