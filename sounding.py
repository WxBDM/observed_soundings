# -*- coding: utf-8 -*-
from datetime import datetime
import matplotlib.gridspec as gridspec

from time import gmtime, strftime
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import metpy.calc as mpcalc

from metpy.calc import resample_nn_1d
from metpy.io import get_upper_air_data
from metpy.plots import Hodograph, SkewT
from metpy.units import units


def write_text(skew, x, y, text, halign, valign, fontsize = 'x-small', color = 'k', bbox = False):
    if bbox is True:
        bbox_props = dict(boxstyle="square,pad=0.3", fc="white", ec="black", lw=1)
        skew.ax.text(x, y, text, horizontalalignment = halign, verticalalignment = valign,
            transform=skew.ax.transAxes, fontsize = fontsize, color = color, bbox = bbox_props)
    else:
        skew.ax.text(x, y, text, horizontalalignment = halign, verticalalignment = valign,
            transform=skew.ax.transAxes, fontsize = fontsize, color = color)

wfos = {
 'ALB': ['Albany', 'NY'],
 'BUF': ['Buffalo', 'NY'],
 'CAR': ['Caribou', 'ME'],
 'GYX': ['Portland', 'ME'],
 'IAD': ['Sterling', 'VA'],
 'ILN': ['Cincinnati', 'OH'],
 'OKX': ['NYC', 'NY'],
 'PIT': ['Pittsburgh', 'PA'],
 'LBF' : ['North Platte', 'NE'],
 'DDC' : ['Dodge City', 'KS']
}

year = int(input("What year? "))
month = int(input("What month? "))
day = int(input("What day? "))
zulu = int(input("Time (zulu)? "))

station = input('Input the station (ex. DDC, GYX, IAD): ').upper()

plt.close('all')

try:
    dataset = get_upper_air_data(datetime(year, month, day, zulu), station)
except ValueError:
    print('Trying iowa state data.')
    try:
        dataset = get_upper_air_data(datetime(year, month, day, zulu), station, source = 'iastate')
    except ValueError:
        print('No data available for {}z'.format(zulu))
        print('Moving onto next station.')
        
p = dataset.variables['pressure'][:]
T = dataset.variables['temperature'][:]
Td = dataset.variables['dewpoint'][:]
u = dataset.variables['u_wind'][:]
v = dataset.variables['v_wind'][:]
hgt = dataset.variables['height'][:]

# Create a new figure.
fig = plt.figure(figsize=(15, 9))

#create a new axis to put potential temperature
gs = gridspec.GridSpec(20,10)
pot_temp = fig.add_subplot(gs[-7:-1, -2:]) #for equiv. pot. temp

# Grid for plots
skew = SkewT(fig, subplot = gs[:, :-3])

#grid for required required levels
req_levels = fig.add_subplot(gs[:2, -2:])
req_levels.axis('off')

#Line at constant T
skew.ax.axvline(0, color='c', linestyle='solid', linewidth=1, alpha = 0.8)
skew.ax.axvline(-10, color='c', linestyle='solid', linewidth=1, alpha = 0.8)

# Plot the data using normal plotting functions, in this case using
# log scaling in Y, as dictated by the typical meteorological plot
skew.plot(p, T, 'r')
skew.ax.set_xlabel('Temperature (Celsius)')
skew.plot(p, Td, 'g')

#plots the barbs
my_interval = np.arange(100, 1000, 20) * units('mbar')
ix = resample_nn_1d(p, my_interval)

skew.plot_barbs(p[ix], u[ix], v[ix])
skew.ax.set_ylim(1075, 100)
skew.ax.set_ylabel('Pressure (hPa)')

lcl_pressure, lcl_temperature = mpcalc.lcl(p[0], T[0], Td[0]) #LCL
pwat = mpcalc.precipitable_water(Td, p, 500 * units.hectopascal).to('in') #PWAT
cape, cin = mpcalc.most_unstable_cape_cin(p[:], T[:], Td[:]) #MUCAPE
cape_sfc, cin_sfc = mpcalc.surface_based_cape_cin(p, T, Td) #SBCAPE
prof = mpcalc.parcel_profile(p, T[0], Td[0]).to('degC') #parcel profile
equiv_pot_temp = mpcalc.equivalent_potential_temperature(p, T, Td) #equivalent potential temperature
el_pressure, el_temperature = mpcalc.el(p, T, Td) #elevated level
lfc_pressure, lfc_temperature = mpcalc.lfc(p, T, Td) #LFC

#calculates shear
u_threekm_bulk_shear, v_threekm_bulk_shear = mpcalc.bulk_shear(p, u, v, hgt, bottom = min(hgt), depth = 3000 * units.meter)
threekm_bulk_shear = mpcalc.get_wind_speed(u_threekm_bulk_shear, v_threekm_bulk_shear)
u_onekm_bulk_shear, v_onekm_bulk_shear = mpcalc.bulk_shear(p, u, v, hgt, bottom = min(hgt), depth = 1000 * units.meter)
onekm_bulk_shear = mpcalc.get_wind_speed(u_onekm_bulk_shear, v_onekm_bulk_shear)

#shows the level of the LCL, LFC, and EL.
skew.ax.text(T[0].magnitude, p[0].magnitude + 5, str(int(np.round(T[0].to('degF').magnitude))), fontsize = 'medium', horizontalalignment = 'left', verticalalignment = 'top', color = 'red')
skew.ax.text(Td[0].magnitude, p[0].magnitude + 5, str(int(np.round(Td[0].to('degF').magnitude))), fontsize = 'medium', horizontalalignment = 'right', verticalalignment = 'top', color = 'green')
skew.ax.text(lcl_temperature.magnitude + 5, lcl_pressure.magnitude, "---- LCL", fontsize = 'medium', verticalalignment = 'center')
skew.ax.text(Td[0].magnitude - 10, p[0].magnitude, 'SFC: {}hPa ----'.format(p[0].magnitude), fontsize = 'medium', horizontalalignment = 'right', verticalalignment = 'center', color = 'black')

if str(lfc_temperature.magnitude) != 'nan': #checks to see if LFC/EL exists. If not, skip.
    skew.ax.text(lfc_temperature.magnitude + 5, lfc_pressure.magnitude, "---- LFC", fontsize = 'medium', verticalalignment = 'center')
    skew.ax.text(el_temperature.magnitude + 5, el_pressure.magnitude, "---- EL", fontsize = 'medium', verticalalignment = 'center')

skew.plot(p, prof, 'k-', linewidth=1) #plots parcel profile
skew.shade_cape(p, T, prof) #shades cape

#plots potential temperature
pot_temp.plot(equiv_pot_temp, p)
pot_temp.invert_yaxis()
pot_temp.set_title('Equivalent Potential Temperature', fontsize = 'small')
pot_temp.set_ylabel('')
pot_temp.set_xlabel('')
pot_temp.set_ylim(1075, 525, 200)
pot_temp.set_xlim(min(equiv_pot_temp.magnitude) - 20, 390, 20)

#shows variables
fig.text(0.70, 0.42, "MUCAPE: %0.2f J/kg" % cape.magnitude, fontsize = 'large', color = '#ff9933')
fig.text(0.70, 0.45, "SBCAPE: %0.2f J/kg" % cape_sfc.magnitude, fontsize = 'large', color = '#ff9933')
fig.text(0.70, 0.48, "PWAT: %0.2f .in" % pwat.magnitude, fontsize = 'large', color = '#3399ff')
fig.text(0.70, 0.51, "0-3km Shear: %0.0d knots" % np.round(threekm_bulk_shear.magnitude), fontsize = 'large', color = 'k')
fig.text(0.70, 0.54, "0-1km Shear: %0.0d knots" % np.round(onekm_bulk_shear.magnitude), fontsize = 'large', color = 'k')

# Add the relevant special lines
skew.plot_dry_adiabats(alpha = 0.1, linestyle = '-')
skew.plot_moist_adiabats(alpha = 0.1, linestyle = '-')
skew.plot_mixing_lines(alpha = 0.1, linestyle = '-')

skew.ax.set_xlim(-50, 60) #sets x axis

# Create a hodograph
ax_hod = inset_axes(skew.ax, '40%', '40%', loc=1, borderpad = 1.0)
h = Hodograph(ax_hod, component_range=85.)
h.add_grid(increment=20)
h.plot_colormapped(u, v, np.hypot(u, v))

pot_temp.text(0.5, -0.28, '@WxBDM', horizontalalignment = 'center', verticalalignment = 'bottom',
        transform=pot_temp.transAxes, fontsize = 'large', color = '#6b6964', zorder = 10)
write_text(skew, 0.002, 1, '{0} ({4}, {5}) Observed Sounding for {1}/{2}/{3}: %0.2dz'.format(station, str(month), str(day), str(year), wfos[station][0], wfos[station][1]) % zulu,
    'left', 'bottom', fontsize = 'x-large')
