# plot_script.gp

set datafile separator ','  # Specify that the data file uses commas as separators
set grid  # Enable grid for better readability
set term pngcairo size 1600, 900   # Set output terminal to PNG
set output 'plot.png'  # Output file

# Constant for conversion from seconds to years
seconds_in_year = 365 * 24 * 60 * 60


set multiplot layout 2,2  
set logscale y  # Set both axes to logarithmic scale
set format y "10^{%L}"  # Format y-axis labels for better readability in log scale
set xlabel 't (years)'  # Label for x-axis
set ylabel 'V (m/s)'  # Label for y-axis
plot 'out.txt' using ($1/seconds_in_year):4 with line lw 2 notitle
unset logscale
unset format y

# Second subplot (psi vs t)
set xlabel 't (years)'
set ylabel 'psi'
plot 'out.txt' using ($1/seconds_in_year):3 with lines lw 2 notitle

# Third subplot (D vs t)
set xlabel 't (years)'
set ylabel 'D (m)'
plot 'out.txt' using ($1/seconds_in_year):2 with lines lw 2 notitle

# Fourth subplot (tau vs t)
set xlabel 't (years)'
set ylabel 'tau (MPa)'
plot 'out.txt' using ($1/seconds_in_year):($5/1000000) with lines lw 2 notitle

unset multiplot
