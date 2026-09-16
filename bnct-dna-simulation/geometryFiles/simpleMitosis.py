"""

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from fractaldna.structure_models import voxelisation as v

def get_cube():   
    phi = np.arange(1,10,2)*np.pi/4
    Phi, Theta = np.meshgrid(phi, phi)

    x = np.cos(Phi)*np.sin(Theta)
    y = np.sin(Phi)*np.sin(Theta)
    z = np.cos(Theta)/np.sqrt(2)
    return x,y,z

def path_M(densityAim, volume, densityMAim,voxelheightStraight, voxelheightTurn):

    volumeChromosome = (densityAim*volume) /densityMAim

   
    chromosomeZ = voxelheightStraight*6 + voxelheightTurn*2
    chromosomeY = voxelheightTurn*6
    # chromosomeX = voxelheightTurn*14
    chromosomeX = volumeChromosome/chromosomeZ/chromosomeY + voxelheightTurn*3

    minLim = np.array([-chromosomeX/2+45/2,-chromosomeY/2+45/2,-chromosomeZ/2+45/2])
    maxLim = np.array([chromosomeX/2-45/2,chromosomeY/2-45/2,chromosomeZ/2-45/2])


    nPoints =10000 
    points = np.zeros((nPoints,3))
    turns = np.zeros(nPoints)


    i = 0
    direction = [0,0,1]
    sign = 1
    turn = 1

    while max(points[:,0])<chromosomeX/2 - voxelheightStraight/2:

        direction = direction/(np.linalg.norm(direction))

        if i==0:
            nextPoint = minLim 
            direction = [0,0,1]
            points[i,:] = nextPoint
            turns[i] = 0
            i+=1

        elif turns[i-1]==1:
            # add straight segment after turn
            nextPoint = direction*(voxelheightStraight/2+ voxelheightTurn/2)
            points[i,:] = points[i-1,:] + nextPoint
            turns[i] = 0

            i+=1
        # will a straight and turn section fit in the direction of travel? if so add if not turn
        elif (points[i-1,:] + direction*(1.5*voxelheightStraight+voxelheightTurn/2)) [np.where(abs(direction)>0)[0][0]] < maxLim[np.where(abs(direction)>0)[0][0]] and (points[i-1,:] + direction*(1.5*voxelheightStraight+voxelheightTurn/2)) [np.where(abs(direction)>0)[0][0]] > minLim[np.where(abs(direction)>0)[0][0]]:
            # add straight segment
            nextPoint = direction*voxelheightStraight
            points[i,:] = points[i-1,:] + nextPoint
            turns[i] = 0

            i+=1

        else:
            # add double turn segment 
            nextPoint = direction*(voxelheightTurn/2+voxelheightStraight/2)
            points[i,:] = points[i-1,:] + nextPoint
            turns[i] = 1

            if ((points[i-1,:] + direction*voxelheightTurn)[1] < maxLim[1]) and ((points[i-1,:] + direction*voxelheightTurn)[1] >= minLim[1]) and turn>0:
                direction = np.array([0,1,0])
            elif ((points[i-1,:] + direction*voxelheightTurn)[1] > minLim[1]) and turn<0:
                direction = np.array([0,-1,0])
            else:
                direction = np.array([1,0,0])
                sign *=-1
                turn *= -1

            i+=1

            nextPoint = direction*(voxelheightTurn)
            points[i,:] = points[i-1,:] + nextPoint
            turns[i] = 1

            if points[i,2] >0 :
                direction = [0,0,-1]
            else:
                direction = [0,0,1]

            i+=1

    points = points[:i-1]
    turns = turns[:len(points)]

    # centre points
    minvals = np.array([np.inf, np.inf, np.inf])
    maxvals = np.array([-np.inf, -np.inf, -np.inf])

    # identify max/min values
    for point in points:
        for (ii, (minv, v, maxv)) in enumerate(zip(minvals, point, maxvals)):
            if v < minv:
                minvals[ii] = v
            elif v > maxv:
                maxvals[ii] = v

    # transform
    transform = -minvals - (maxvals - minvals) / 2.0
    for ind, point in enumerate(points):
        oldpos = point
        points[ind] = oldpos + transform

    return points, turns

if __name__ == "__main__":
    densityAim = 0.0124
    volume = 750**3
    densityMAim = 0.07
    voxelheightStraight = 75
    voxelheightTurn = 47.5

    positions, turns = path_M(densityAim, volume, densityMAim, voxelheightStraight, voxelheightTurn)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    x,y,z = get_cube()

    # for plot and estimated values
    volumeChromosome = (densityAim*volume) /densityMAim

    chromosomeZ = voxelheightStraight*7
    chromosomeY = voxelheightStraight*7
    chromosomeX = volumeChromosome/chromosomeZ/chromosomeY

    bpStraight = 27770/2 #70 histones
    bpTurn = 13912/2 #50 histones

    ax.plot_surface(x*chromosomeX, y*chromosomeY, z*chromosomeZ, color='b', alpha=0.1)

    ax.plot(positions[:,0], positions[:,1], positions[:,2], ".-k")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    bp = (np.sum([1 for a in turns if a == 0])*bpStraight + np.sum([1 for a in turns if a == 1])*bpTurn)

    estDensity = bp/(chromosomeZ*chromosomeX*chromosomeY)
    estNucleusDensity =bp/volume

    print(f"estimated chromosome density = {estDensity} bp/nm3")
    print(f"estimated nucleus density = {estNucleusDensity} bp/nm3")
    print(f"estimated total bp  = {bp}")
    print(f"estimated linear packing density  = {bp/1e6/(chromosomeX/1e3)} Mbp/um") # to achieve more realistic need y/z to be double (radius ~520 nm) but x would be only one layer with so little DNA. 

    plt.show()