import matplotlib.pyplot as plt
import numpy as np
import esutil as es
from astropy.table import Table
from mpl_toolkits.basemap import Basemap
import astropy.coordinates as coo

htm_grid = es.htm.HTM(4)

centres_tab = Table.read("level4_grid.txt",format="ascii")
theta = centres_tab["theta"]
phi = centres_tab["phi"]

plot_theta = 90 - theta

new_map = Basemap(projection="moll",resolution="l", lat_0=45, lon_0=0)

points_done = []

for i in range(2048,2048*2):
    v1, v2, v3 = htm_grid.get_vertices(i)

    theta_1 = np.rad2deg(np.arccos(v1[2]))
    phi_1 = np.rad2deg(np.arctan2(v1[0],v1[1]))
    theta_2 = np.rad2deg(np.arccos(v2[2]))
    phi_2 = np.rad2deg(np.arctan2(v2[0],v2[1]))
    theta_3 = np.rad2deg(np.arccos(v3[2]))
    phi_3 = np.rad2deg(np.arctan2(v3[0],v3[1]))

    if (phi_1 < 0):
        phi_1 = 360.0 + phi_1
    if (phi_2 < 0):
        phi_2 = 360.0 + phi_2
    if (phi_3 < 0):
        phi_3 = 360.0 + phi_3


    xpt1, ypt1 = new_map(phi_1, 90 - theta_1)
    xpt2, ypt2 = new_map(phi_2, 90 - theta_2)
    xpt3, ypt3 = new_map(phi_3, 90 - theta_3)

    new_map.plot(xpt1, ypt1, "bo", ms=1)
    new_map.plot(xpt2, ypt2, "bo", ms=1)
    new_map.plot(xpt3, ypt3, "bo", ms=1)

    

    new_map.drawgreatcircle(phi_1,90-theta_1,phi_2,90-theta_2,linewidth=0.7,color="gray")
    new_map.drawgreatcircle(phi_2,90-theta_2,phi_3,90-theta_3,linewidth=0.7,color="gray")
    new_map.drawgreatcircle(phi_1,90-theta_1,phi_3,90-theta_3,linewidth=0.7,color="gray")
    
    

plt.show()
