'''
Functions for creating complex DNA structure using fractalDNA, convert to input format for RBE simulation.
'''

from fractaldna.dna_models import dnachain as dna
from mayavi import mlab
import numpy as np
from fractaldna.structure_models import voxelisation as v
from fractaldna.structure_models import hilbert as h
import pandas as pd
import matplotlib.pyplot as plt

def chromatinFibrePosition(solenoidType, voxelheight = 75, radius = 11, nhistones = 33, histone_angle = 50):
    voxelheight = voxelheight *10 # convert to A
    radius = radius *10 # convert to A
    Angstrum2nm = 1e-1 #A positions convert to nm

    if solenoidType == "Straight":
        chain = dna.Solenoid(
            voxelheight = voxelheight, #A
            radius = radius, #A
            nhistones = nhistones,
            histone_angle = histone_angle,
            twist = False,
            chain = 0, # not used - provides chain number for molecularDNA
        )
    elif solenoidType == "Turntwist":
        chain = dna.TurnedSolenoid(
            voxelheight = voxelheight, #A
            radius = radius, #A
            nhistones = nhistones,
            histone_angle = histone_angle,
            twist = True,
            chain = 0,
        )
    if solenoidType == "Turn":
        chain = dna.TurnedSolenoid(
            voxelheight = voxelheight, #A
            radius = radius, #A
            nhistones = nhistones,
            histone_angle = histone_angle,
            twist = False,
            chain = 0,
        )

    chain.translate([0, 0, -voxelheight / 2.0])
    # MayaVI plots are best for visualisation here
    # plot = chain.to_strand_plot()

    # mlab.show()
    # plot.scene.save_jpg("single_solenoid_strand_plot.jpg")

    # plot = chain.to_line_plot()
    # mlab.show()

    # plot.scene.save_jpg("single_solenoid_line_plot.jpg")
    sugar = chain.to_frame()
    hist = chain.histones_to_frame()

    positionSugars = np.zeros((max(sugar.bp_idx)+1, 12))
    positionHistones = np.zeros((max(hist.histone_idx)+1, 6))

    for i in range(max(sugar.bp_idx)+1):
        atThatIdx = sugar[sugar.bp_idx==i]

        base0 = atThatIdx[atThatIdx.strand_idx==0].set_index("#name").drop(index = "Sugar").drop(index = "Phosphate")
        base1 = atThatIdx[atThatIdx.strand_idx==1].set_index("#name").drop(index = "Sugar").drop(index = "Phosphate")
        sugar0 = atThatIdx.loc[(atThatIdx["#name"]=="Sugar") & (atThatIdx["strand_idx"] ==0)]
        sugar1 = atThatIdx.loc[(atThatIdx["#name"]=="Sugar") & (atThatIdx["strand_idx"] ==1)]

        positionSugars[i][0:3] = [sugar0["pos_x"].values[0]*Angstrum2nm,sugar0["pos_y"].values[0]*Angstrum2nm,sugar0["pos_z"].values[0]*Angstrum2nm]
        positionSugars[i][3:6] = [sugar1["pos_x"].values[0]*Angstrum2nm,sugar1["pos_y"].values[0]*Angstrum2nm,sugar1["pos_z"].values[0]*Angstrum2nm]
        positionSugars[i][6:9] = [base0["pos_x"].values[0]*Angstrum2nm,base0["pos_y"].values[0]*Angstrum2nm,base0["pos_z"].values[0]*Angstrum2nm]
        positionSugars[i][9:12] = [base1["pos_x"].values[0]*Angstrum2nm,base1["pos_y"].values[0]*Angstrum2nm,base1["pos_z"].values[0]*Angstrum2nm]

    for i in range(max(hist.histone_idx)+1):
        histone = hist.loc[(hist["histone_idx"]==i)]

        positionHistones[i][:] = [histone.pos_x.values[0]*Angstrum2nm, histone.pos_y.values[0]*Angstrum2nm, histone.pos_z.values[0]*Angstrum2nm, 
            histone.rot_x.values[0], histone.rot_y.values[0], histone.rot_z.values[0]]

    return positionSugars, positionHistones

def complexGeometryFractal(voxelheight, nn=2):
    nanometer = 1 # save in nm

    initial_string = "X"
    # Iterate it as required (here, nn=2)
    # for nn = 8, this takes about two hours on my 16GB RAM Mac
    # nn = 2
    iterated_lstring = h.iterate_lstring(initial_string)
    for _ in range(nn - 1):
        iterated_lstring = h.iterate_lstring(iterated_lstring)

    vf = v.VoxelisedFractal.fromLString(iterated_lstring, pbar=True)
    vf.center_fractal()
    # fig = vf.to_plot()
    # fig.savefig('plot.png')

    # fig = vf.to_pretty_plot()
    # fig.savefig('fractalprettyplot.png')

    simplePos = vf.to_frame()

    # scale to nm, based on one step being one chromatin fibre length which is 75 nm
    positions = simplePos.copy()
    positions.POS_X = simplePos.POS_X*voxelheight*nanometer
    positions.POS_Y = simplePos.POS_Y*voxelheight*nanometer
    positions.POS_Z = simplePos.POS_Z*voxelheight*nanometer

    boxSize = max(positions["POS_X"])*2/nanometer + voxelheight #size of box return in nm (always cube)

    return positions,  boxSize

def complexGeometryFromFile(voxelheight, filename):
    nanometer = 1 # save nm

    df = pd.read_csv(filename, sep="\t")

    vf = v.VoxelisedFractal.from_path(df.values, pbar=False)

    vf.center_fractal()
    # fig = vf.to_plot()
    # fig.savefig('plot.png')

    # fig = vf.to_pretty_plot()
    # fig.savefig('fractalprettyplot.png')

    simplePos = vf.to_frame()

    # scale to nm, based on one step being one chromatin fibre length which is 75 nm
    positions = simplePos.copy()
    positions.POS_X = simplePos.POS_X*voxelheight*nanometer
    positions.POS_Y = simplePos.POS_Y*voxelheight*nanometer
    positions.POS_Z = simplePos.POS_Z*voxelheight*nanometer

    boxSize = max(positions["POS_X"])*2/nanometer + voxelheight #size of box return in nm (always cube)

    return positions,  boxSize

def complexGeometryFromArray(positions, turn):
    """ Create a dataframe matching the fractaldna structure with positions and angles for each solenoid voxel. First segment must be straight

    Args:
        positions (array): x,y,z positions
        turn (array): whether point is turn or not
    """
    defaultPrincipal = np.array([1, 0, 0])
    defaultHeading = np.array([0, 0, 1])
    defaultOrtho = np.cross(defaultPrincipal, defaultHeading)
    defaultAxis = np.array([defaultPrincipal, -defaultOrtho, defaultHeading])
    defaultAxis = np.transpose(defaultAxis)
    defaultAxisInv = np.linalg.inv(defaultAxis)
    inPrincipal = np.zeros((len(positions),3))
    inHeading =  np.zeros((len(positions),3))
    outPrincipal =  np.zeros((len(positions),3))
    outHeading =  np.zeros((len(positions),3))

    d = [] #for final data frame
    count = 0

    for ind in range(len(positions)-1):
        if ind == 0: # first line is straight so calculate from points 0 and 1
            currpos = positions[ind,:]
            nextpos = positions[ind+1,:]
            direction = nextpos-currpos
            direction  = direction/np.linalg.norm(direction)
            inHeading[ind,:] = direction
            inPrincipal[ind,:] = np.array([1, 0, 0])
            if (direction == inPrincipal[ind,:]).all():
                inPrincipal[ind,:] = np.array([0, 1, 0])
            outHeading[ind,:] = direction
            outPrincipal[ind,:] = inPrincipal[ind,:]
            continue
        else:

            currpos = positions[ind,:]
            nextpos = positions[ind+1,:]
            prevpos = positions[ind-1,:]

            firstChange = currpos - prevpos
            secondChange = nextpos - currpos
            firstChange = firstChange /np.linalg.norm(firstChange)
            secondChange = secondChange/ np.linalg.norm(secondChange)
            perp = np.cross(firstChange, secondChange)

            if turn[ind]==1:
                inHeading[ind,:] = outHeading[ind-1,:]
                inPrincipal[ind,:] = outPrincipal[ind-1,:]
                outHeading[ind,:] = secondChange[:]
                outPrincipal[ind,:] = perp
                paxis = outHeading[ind,:]
            else:
                inHeading[ind,:] = outHeading[ind-1,:]
                inPrincipal[ind,:] = outPrincipal[ind-1,:]
                outHeading[ind,:] = outHeading[ind-1,:]
                outPrincipal[ind,:] = outPrincipal[ind-1,:]
                paxis = inPrincipal[ind,:]
                

            orth = np.cross(paxis,inHeading[ind,:])

            axis = np.array([paxis, -orth, inHeading[ind,:]])
            axis = np.transpose(axis)

            direction = np.dot(axis, defaultAxisInv)

            psi, theta, phi = v.getEulerAngles(direction)

            if turn[ind]==0:
                turnType ="straight"
            elif turn[ind]==1:
                turnType ="turn"   
            else:
                raise ValueError
            
            d.append([count, turnType, positions[ind,0], positions[ind,1], positions[ind,2], psi, theta, phi])

            count += 1

    chromatinFibrePosition = pd.DataFrame(d,columns=["#IDX","KIND", "POS_X", "POS_Y", "POS_Z", "EUL_PSI", "EUL_THETA", "EUL_PHI"])

    return chromatinFibrePosition

def transformSugarPositions(positionSugars, row, offset=[0,0,0]):
    #  add rotation
    new = positionSugars.copy()

    alpha = row.EUL_PSI
    beta = row.EUL_THETA
    gamma = row.EUL_PHI 

    rotMatrix = np.array([
    [np.cos(beta)*np.cos(gamma),np.cos(gamma)*np.sin(beta)*np.sin(alpha)-np.sin(gamma)*np.cos(alpha),np.cos(alpha)*np.sin(beta)*np.cos(gamma)+np.sin(alpha)*np.sin(gamma)],
    [np.sin(gamma)*np.cos(beta),np.sin(alpha)*np.sin(beta)*np.sin(gamma)+np.cos(alpha)*np.cos(gamma),np.sin(gamma)*np.sin(beta)*np.cos(alpha)-np.cos(gamma)*np.sin(alpha)],
    [-1*np.sin(beta), np.cos(beta)*np.sin(alpha), np.cos(beta)*np.cos(alpha)]]).T

    # angle transform 
    for i in range(new.shape[0]):
        new[i,0:3] = np.matmul([positionSugars[i,0],positionSugars[i,1],positionSugars[i,2]],rotMatrix)
        new[i,3:6] = np.matmul([positionSugars[i,3],positionSugars[i,4],positionSugars[i,5]],rotMatrix)
        new[i,6:9] = np.matmul([positionSugars[i,6],positionSugars[i,7],positionSugars[i,8]],rotMatrix)
        new[i,9:12] = np.matmul([positionSugars[i,9],positionSugars[i,10],positionSugars[i,11]],rotMatrix)

    # translation
    new[:,0] += row.POS_X + offset[0] # sugar x
    new[:,1] += row.POS_Y + offset[1]  # sugar y
    new[:,2] += row.POS_Z + offset[2]  # sugar z
    new[:,3] += row.POS_X + offset[0]  # sugar x
    new[:,4] += row.POS_Y + offset[1]  # sugar y
    new[:,5] += row.POS_Z + offset[2]  # sugar z
    new[:,6] += row.POS_X + offset[0]  # base x
    new[:,7] += row.POS_Y + offset[1]  # base y
    new[:,8] += row.POS_Z + offset[2]  # base z
    new[:,9] += row.POS_X + offset[0]  # base x
    new[:,10] += row.POS_Y + offset[1]  # base y
    new[:,11] += row.POS_Z + offset[2]  # base z

    return new

def transformHistonePositions(positionHistone, row, offset=[0,0,0]):
    #  add rotation for more complex geometry
    new = positionHistone.copy()

    alpha = row.EUL_PSI
    beta = row.EUL_THETA
    gamma = row.EUL_PHI 

    rotMatrix = np.array([
    [np.cos(beta)*np.cos(gamma),np.cos(gamma)*np.sin(beta)*np.sin(alpha)-np.sin(gamma)*np.cos(alpha),np.cos(alpha)*np.sin(beta)*np.cos(gamma)+np.sin(alpha)*np.sin(gamma)],
    [np.sin(gamma)*np.cos(beta),np.sin(alpha)*np.sin(beta)*np.sin(gamma)+np.cos(alpha)*np.cos(gamma),np.sin(gamma)*np.sin(beta)*np.cos(alpha)-np.cos(gamma)*np.sin(alpha)],
    [-1*np.sin(beta), np.cos(beta)*np.sin(alpha), np.cos(beta)*np.cos(alpha)]]).T

    # angle transform 
    for i in range(new.shape[0]):
        new[i,0:3] = np.matmul([new[i,0],new[i,1],new[i,2]],rotMatrix)

    new[:,0] += row.POS_X + offset[0] # histone x
    new[:,1] += row.POS_Y + offset[1]  # histone y
    new[:,2] += row.POS_Z + offset[2]  # histone z

    # y and z reversed in RBE - rotation not relevant now sphere
    new[:,3] = positionHistone[:,3]
    new[:,5] = positionHistone[:,4] 
    new[:,4] = positionHistone[:,5]

    return new

def createFiles(chromatinFibrePositions, positionStraight, positionTurn, positionTurnTwist, width, voxelheightStraight, voxelheightTurn, f, offset = [0,0,0], plot = False, write = True, mitosis = False, densityAim = 0.012, densityMAim = 0.07):
    # chromatinFibrePositions, width = complexGeometry(voxelheight, n)

    numSugars = 0
    numHistones = 0

    straight = 0
    turn = 0
    turntwist = 0

    s = []
    h = []
    for ind, row in enumerate(chromatinFibrePositions.itertuples(index=True, name='Pandas')):
        if mitosis:
            density =  numSugars/((densityAim*(width)**3) /densityMAim)
            if density > densityMAim:
                continue # place segments until the required chromosome density is achieved 
        if row.KIND == "straight":
            newSugarPositions = transformSugarPositions(positionStraight[0], row, offset)
            newHistonePositions = transformHistonePositions(positionStraight[1], row, offset)
            straight+=1
        elif row.KIND == "turn":
            newSugarPositions = transformSugarPositions(positionTurn[0], row, offset)
            newHistonePositions = transformHistonePositions(positionTurn[1], row, offset)
            turn +=1
        elif row.KIND == "turntwist":
            newSugarPositions = transformSugarPositions(positionTurnTwist[0], row, offset)
            newHistonePositions = transformHistonePositions(positionTurnTwist[1], row, offset)
            turntwist+=1
        else:
            raise ValueError()

        for i in range(len(newSugarPositions)):
            if write:
                s.extend(newSugarPositions[i])
            numSugars+=1

        for i in range(len(newHistonePositions)):
            if write:
                h.extend(newHistonePositions[i])
            numHistones+=1
            if plot:
                plot.plot(newHistonePositions[i][0], newHistonePositions[i][1], newHistonePositions[i][2], ".-")

    print("density = ", numSugars/(width)**3, file=f)
    if mitosis:

        # calculate chromosomeVolume based on sugar positions. Add 5 nm margin to each point 
        chromosomeVolume = (-1*min(s[0::12])+max(s[0::12])+10)*(-1*min(s[1::12])+max(s[1::12])+10)*(-1*min(s[2::12])+max(s[2::12])+10)
        print("chromosome volume = ", chromosomeVolume , file=f)
        print("chromosome density = ", numSugars/chromosomeVolume, file=f)

    print("box size = ", (width), "nm", file=f)
    print("straight segments = ", straight, "turn segments = ", turn, "turn twist segments = ", turntwist, file=f)

    print("number histones = ", numHistones, file=f)
    effectiveLength = straight*voxelheightStraight + turn*voxelheightTurn*np.pi/4 + turntwist*voxelheightTurn*np.pi/4 #estimate length of chromatin, turn have a radius of voxel size/2 and are 1/4 of a circle
    print("packing fraction = ", numHistones*11/effectiveLength, " nucleosomes per 11 nm", file=f)
    print("packing fraction (straight segments)= ", positionStraight[1].shape[0]*11/voxelheightStraight, " nucleosomes per 11 nm", file=f)
    print("packing fraction (turn segments)= ", positionTurn[1].shape[0]*11/(voxelheightTurn*np.pi/4), " nucleosomes per 11 nm", file=f)


    print("Number of base pairs", numSugars, file = f)

    return s,h
# DNA structure parameters

def loadOrCreateSolenoids(voxelheight, nhistones, radius, histone_angle, voxelType):
    """ Load previously saved solenoids or create if not found. Solenoids saved to save time.

    Args:
        voxelheight (float): voxel height (cube)
        nhistones (int): number of histones in solenoid
        radius (float): radius of solenoid
        histone_angle (float): angle of histones
        voxelType (str): "Straight", "Turn", "Turntwist"
    """

    try:
        positionSugars = np.load(f"posSugars{voxelType}_H{voxelheight}nm_R{radius}nm_{nhistones}hist_{histone_angle}deg.npy")
        positionHistones = np.load(f"posHistones{voxelType}_H{voxelheight}nm_R{radius}nm_{nhistones}hist_{histone_angle}deg.npy")
    except:
        print(f"creating {voxelType} solenoid")
        positionSugars, positionHistones = chromatinFibrePosition(voxelType, voxelheight, radius, nhistones, histone_angle)

        np.save(f"posSugars{voxelType}_H{voxelheight}nm_R{radius}nm_{nhistones}hist_{histone_angle}deg", positionSugars)
        np.save(f"posHistones{voxelType}_H{voxelheight}nm_R{radius}nm_{nhistones}hist_{histone_angle}deg", positionHistones)

    return positionSugars, positionHistones