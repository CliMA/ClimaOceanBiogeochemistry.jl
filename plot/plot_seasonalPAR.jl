
function seasonal_PAR(y, z, t)
    # Constants
    day_in_year = 365.25
    axial_tilt = 23.5  # Earth's axial tilt in degrees
    latitude = (y/Ly - 0.5) * 180.0  # Convert y to latitude in degrees
    
    # Solar declination (varies over the year due to Earth's tilt)
    solar_declination = axial_tilt * sinpi(2 * t / day_in_year)
    
    # Solar angle of incidence (latitude and solar declination combined)
    solar_angle = latitude - solar_declination
    
    # PAR depends on the cosine of the solar angle 
    angle_effect = cosd(solar_angle)
    angle_effect = max(0.0, angle_effect)  # Ensure no negative values
    
    # Scale PAR by angle effect
    PAR = 700 * angle_effect
    return PAR
end

time_duration = 0:365
y_positions = range(0, Ly, length=20)

# Create data arrays for different plots
time_series_data = [seasonal_PAR(y, 0, t) for y in y_positions, t in time_duration]

# Create figure with multiple subplots
fig = Figure(size=(1000, 800))

# Plot 1: Heatmap
ax1 = Axis(fig[2, 1], 
    title="PAR Distribution Over Space and Time",
    xlabel="Day of Year",
    ylabel="y position (0=South, 1=North)")
    
hm = heatmap!(ax1, time_duration, y_positions, time_series_data',
    colorrange = (0,750), colormap=:thermal, 
    interpolate=true)
Colorbar(fig[2, 2], hm, label="PAR (W/m²)")

# Plot 2: Line plot for different positions
ax2 = Axis(fig[1, 3],
    title="PAR Over Time at Different Positions",
    xlabel="Day of Year",
    ylabel="PAR (W/m²)")
selected_y = [0.0, 0.25, 0.5, 0.75, 1.0]*Ly
for y in selected_y
    lines!(ax2, time_duration, [seasonal_PAR(y, 0, t) for t in time_duration],
           label="y = $y")
end
ylims!(ax2,0,750)
axislegend(ax2, position=:rb)

# Plot 3: Line plot for different days
ax3 = Axis(fig[2, 3],
    title="PAR Distribution Along y at Different Days",
    xlabel="y position",
    ylabel="PAR (W/m²)")

selected_days = [10, 91, 182, 273, 355]
for day in selected_days
    lines!(ax3, y_positions, [seasonal_PAR(y, 0, day) for y in y_positions],
           label="Day $day")
end
ylims!(ax3,0,750)
axislegend(ax3, position=:rt)

display(fig)