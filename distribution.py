import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Data
radius = np.array([0.17, 0.232, 0.294, 0.355, 0.417, 0.479, 0.541, 0.603, 0.665, 0.726, 0.788, 0.85])
chord = np.array([0.087, 0.079, 0.07, 0.062, 0.054, 0.048, 0.043, 0.037, 0.032, 0.026, 0.019, 0.017])
twist = np.array([31.294, 23.542, 17.914, 14.741, 13.178, 11.567, 9.516, 7.876, 4.101, 3.27, 2.468, 1.671])
airfoils = ['NACA 4412', 'NACA 4412', 'NACA 4412', 'NACA 4412', 'NACA 2412', 'NACA 2412', 
            'NACA 2412', 'NACA 2413', 'NACA 63-415', 'NACA 63-415', 'NACA 63-415', 'NACA 63-416']

# Create figure and axis objects with a single subplot
fig, ax1 = plt.subplots(figsize=(12, 8))

# Style
ax1.set_facecolor('#FBFBF9')

# Plot chord length on primary y-axis
color1 = '#1f77b4'  # Blue
ax1.set_xlabel('Radius (m)')
ax1.set_ylabel('Chord Length (m)', color=color1)
line1 = ax1.plot(radius, chord, color=color1, marker='o', label='Chord Length')
ax1.tick_params(axis='y', labelcolor=color1)

# Create secondary y-axis and plot twist angle
ax2 = ax1.twinx()
color2 = '#000000'  # Black
ax2.set_ylabel('Twist Angle (degrees)', color=color2)
line2 = ax2.plot(radius, twist, color=color2, marker='s', label='Twist Angle')
ax2.tick_params(axis='y', labelcolor=color2)

# Add title
plt.title('Wind Turbine Blade Characteristics', pad=20)

# Combine lines for legend
lines = line1 + line2
labels = [l.get_label() for l in lines]
ax1.legend(lines, labels, loc='lower left')

# Add airfoil transition markers
prev_airfoil = airfoils[0]
transition_points = []
for i, airfoil in enumerate(airfoils[1:], 1):
    if airfoil != prev_airfoil:
        transition_points.append(radius[i])
        # Add vertical line at transition point
        plt.axvline(x=radius[i], color='gray', linestyle='--', alpha=0.5)
        # Add annotation
        plt.text(radius[i], plt.ylim()[1], f'\n{airfoil}', 
                rotation=90, verticalalignment='top', horizontalalignment='right')
    prev_airfoil = airfoil

# Add initial airfoil annotation
plt.text(radius[0], plt.ylim()[1], f'\n{airfoils[0]}', 
         rotation=90, verticalalignment='top', horizontalalignment='right')
         
# Show plot
plt.tight_layout()

ax1.grid(True, alpha=0.3)

plt.show()
