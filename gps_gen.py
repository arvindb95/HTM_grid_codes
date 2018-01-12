"""
    gps_gen.py

    Aim : Generate gps.mac and other accessory files to be used for executing by GEANT4 on a cluster


"""
import numpy as np
import esutil as es
from astropy.table import Table
import os
import astropy.coordinates as coo
import matplotlib.pyplot as plt

def gen_gps(theta,phi,energy,nevt=6.97e6,satIllum=1,R=500.0,path=''):
    """
    Returns a string which should go into the gps file

    Required inputs:
    theta = theta of the point (deg)
    phi = phi of the point (deg)
    energy = energy to be run (keV) 
    Optional inputs:
    nevt = No. of photons to be shined 
    satIllum = Boolean to decide wheteher we need to illuminate only czti(0) or whole satellite(1), default=1
    R = radius of the celestial sphere (a large number), default = 500.0
    path = path where the gps file needs to be saved

    Output:
    gps output file Tddd.dd_Pddd.dd_Edddd.d_gps.mac
    """
    
    # Store energy in E, theta and phi
    E = energy
    th = np.deg2rad(theta)
    ph = np.deg2rad(phi)

    # Calculating the vector from source plane to czti (minus sign because of opposite direction)
    vx = -np.sin(th)*np.sin(ph)
    vy = -np.cos(th)
    vz = -np.sin(th)*np.cos(ph)
    
    # Calculating the position of the source plane
    if (satIllum == 0):
        posx = -R*vx
        posy = -R*vy
        posz = -R*vz
        radius = 50.0
    else :
        posx = -R*vx + 50.0
        posy = -R*vy
        posz = -R*vz + 60.0
        radius = 200.0
        

    # Calculating the rotation of the source plane
    rot1 = np.zeros(3)
    rot2 = np.zeros(3)

    rot1[0] = -vz/np.sqrt(1.0 - vy**2.0)
    rot1[1] = 0
    rot1[2] = vx/np.sqrt(1.0 - vy**2.0)

    rot2[0] = vy*vx/np.sqrt(1.0 - vy**2.0)
    rot2[1] = -np.sqrt(vx**2.0 + vz**2.0)
    rot2[2] = vy*vz/(1.0 - vy**2)

    # Writing to a file
    
    gps_file = open(path+"/T{th:06.2f}_P{ph:06.2f}_E{e:07.2f}_gps.mac".format(th=theta,ph=phi,e=energy),"w")
    
    gps_file.write("""/gps/particle gamma
/gps/ene/type Mono
/gps/ene/mono {} keV
/gps/pos/type Plane
/gps/pos/shape Circle
/gps/pos/centre {} {} {} cm
/gps/direction {} {} {}
/gps/pos/rot1 {} {} {}
/gps/pos/rot2 {} {} {}
/gps/pos/radius {} cm
/tracking/verbose {}
/gps/verbose {}
/run/beamOn {}
            """.format(energy,posx,posy,posz,vx,vy,vz,rot1[0],rot1[1],rot1[2],rot2[0],rot2[1],rot2[2],radius,0,0,int(nevt)))
    gps_file.close()

    gps_filename = "T{th:06.2f}_P{ph:06.2f}_E{e:07.2f}_gps.mac".format(th=theta,ph=phi,e=energy)
    return gps_filename


def get_centers(depth,theta_filename="theta.txt",phi_filename="phi.txt",R=500.0):
    """
    Returns the coordinates of the centers of the trixel of given depth

    Required inputs:
    depth = level of the HTM grid (int from 0 to 12)
    Optional inputs:
    theta_file = Name of file in which calculated theta values are to be put
    phi_file = Name of file in which calculated values of phi are to be put
    R = Radius of the celestial sphere (a large value), default = 500.0

    Returns:
    Arrays of theta and phi values of the centers of the trixels
    """
    htm_map = es.htm.HTM(depth)
    no_of_trixels = 8*4**(depth)
    
    # To get a list of all possible number of trixels for all depths from 0 to 13
    n_of_t = np.zeros(10)
    for i in range(0,len(n_of_t)):
        n_of_t[i] = 4**(i)*8

    n_id_arr = np.zeros(len(n_of_t)+1)
    n_id_arr[0] = 1
    n_id_arr[1:] = n_of_t
    n_id_range = np.cumsum(n_id_arr)
    
    n_id_low =  no_of_trixels#int(n_id_range[depth]) # The lower limit of the number id for the given depth
    n_id_high = no_of_trixels*2#int(n_id_range[depth+1]) # The higher limit of the number id for the given depth
   

    print "N_id low :",n_id_low 
    print "N_id high : ",n_id_high

    theta_arr = []
    phi_arr = []

    for nid in range(n_id_low,n_id_high):
        v1,v2,v3 = htm_map.get_vertices(nid)
        center = (v1+v2+v3)/3.0
        center = center/np.sqrt(center[0]**2.0 + + center[1]**2.0 + center[2]**2.0)
        theta =  np.rad2deg(np.arccos(center[2]))
        phi = np.rad2deg(np.arctan2(center[0],center[1]))
        if (phi < 0):
            phi = 360.0 + phi

        theta_arr.append(theta)
        phi_arr.append(phi)
    
    theta_file = open(theta_filename,"w")
    phi_file = open(phi_filename,"w")
    theta_table = Table([np.array(theta_arr)],names=['theta'])
    theta_table.write(theta_filename,format='ascii',overwrite=True)
    
    phi_table = Table([np.array(phi_arr)],names=['phi'])
    phi_table.write(phi_filename,format='ascii',overwrite=True)

    theta_file.close()
    phi_file.close()
    return theta_arr,phi_arr

# Main function begins

if __name__=='__main__':
    
    # Calling get_centers to write the theta and phi values for the centers of the trixels
    depth = 5

    #get_centers(depth)
    

    lvl5_minus_lvl4_file = Table.read("lvl5_minus_lvl4_grid.txt",format='ascii')
    full_theta = lvl5_minus_lvl4_file['theta'].data
    full_phi = lvl5_minus_lvl4_file['phi'].data

    grb_theta = np.array([60.83,106.11,116.86,0.66,105.66,138.85,140.52,10.15,52.95,156.18,65.32,55.30,65.53])
    grb_phi = np.array([67.57,255.69,184.86,159.48,85.52,315.78,118.09,95.08,273.18,59.18,332.32,35.19,16.19])
    grb_names = np.array(['GRB151006A','GRB160106A','GRB160131A','GRB160325A','GRB160509A','GRB160607A','GRB160623A','GRB160703','GRB160802A','GRB160821A','GRB160910A','GRB171010A','GRB171103A'])

#    for i in range(len(grb_theta)):
#        print "Start %s"%i
#        trans_theta = grb_theta[i]
#        trans_phi = grb_phi[i]
#        trans_pos = coo.SkyCoord(trans_phi,trans_theta-90,unit="deg")
#        neigh_theta_lst = []
#        neigh_phi_lst = []
#        for j in range(len(full_theta)):
#            point_pos = coo.SkyCoord(full_phi[j],full_theta[j]-90,unit="deg")
#            sep = abs(trans_pos.separation(point_pos).value)
#            if (sep <= 10):
#                neigh_theta_lst.append(full_theta[j])
#                neigh_phi_lst.append(full_phi[j])
#        sel_theta_arr = np.array(neigh_theta_lst)
#        sel_phi_arr = np.array(neigh_phi_lst)
#        print "Neighbour theta : ",sel_theta_arr
#        print "Neighbour phi : ",sel_phi_arr
#
#        neigh_file = open(grb_names[i]+"_neigh_grid.txt","w")
#        neigh_table = Table([sel_theta_arr,sel_phi_arr],names=["theta","phi"])
#        neigh_table.write(neigh_file,format='ascii')
#        neigh_file.close()
#        print "============================="+grb_names[i]+"_neigh_grid.txt created ================================"

    final_tab = Table.read("final_finer_grid.txt",format="ascii")
    theta_arr = final_tab["theta"].data
    phi_arr = final_tab["phi"].data

    # Open a file to note down the location of the gps files
    gps_file_list = open("gps_file_list.txt","w")
    
    path_file = open("path.txt","r")
    master_path = path_file.read()[:-1]

    no_of_files = 0
    
    for i in range(len(theta_arr)):
        path = master_path+"/T{th:06.2f}_P{ph:06.2f}".format(th=theta_arr[i],ph=phi_arr[i])
        if not os.path.exists(path):
            os.makedirs(path)
            energy_table = Table.read("grid_energy.txt",format='ascii')
            energy = energy_table['energy'].data
            for ene in energy:
                gps_filename = gen_gps(theta_arr[i],phi_arr[i],ene,path=path)
                gps_file_list.write(path+"/"+gps_filename+"\n")
                no_of_files += 1
    gps_file_list.close()
    print no_of_files," gps.mac files created !!!"


