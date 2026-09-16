'''
Create more complex DNA structure using fractalDNA, convert to input format for RBE simulation.
'''
from continuousGeometry import *
import struct 
from simpleMitosis import path_M
# Model type: fractal, file, cellCycle
model = "fractal"

# If cell cycle type mitosis or interphase
cellCycleType = "mitosis"

# DNA structure parameters
n = 2
voxelheightStraight = 75 #nm
voxelheightTurn = 75 #nm
radius = 11 #nm  chromatin fibre radius
nucleusSize = 300 #nm - only for cellCycle
nhistonesStraight = 33
nhistonesTurn = 33  #actual number in segment will be less because of turn
histone_angle = 50
densityAim = 0.013 #bp/nm3
densityMAim = 0.07 #bp/nm3

inputFilename = "exampleInput.csv"
outputFilename = "output"

f = open(f"{outputFilename}.txt", "w") #write important variables to file

print("Parameters:", file=f)
print("", file=f)

if model == "fractal":

    outputFilename = f"n{n}_{nucleusSize}nm_H{voxelheightStraight}nm_R{radius}nm_{nhistonesStraight}hist_{histone_angle}deg"

    voxelheightTurn = voxelheightStraight
    nhistonesTurn = nhistonesStraight
    print("straight and turn voxels must be the same size, using straight voxel parameters")

    print("model type = ",model, file=f)
    print("n = ",n, file=f)
    print("voxelheightStraight = ",voxelheightStraight, file=f)
    print("nhistonesStraight = ",nhistonesStraight, file=f)
    print("voxelheightTurn = ",voxelheightTurn, file=f)
    print("nhistonesTurn = ",nhistonesTurn, file=f)
    print("radius = ",radius, file=f)
    print("histone_angle = ",histone_angle, file=f)
    print("", file=f)

    positionStraight = loadOrCreateSolenoids(voxelheightStraight, nhistonesStraight, radius, histone_angle, "Straight")
    positionTurn = loadOrCreateSolenoids(voxelheightTurn, nhistonesTurn, radius, histone_angle, "Turn")
    positionTurnTwist = loadOrCreateSolenoids(voxelheightTurn, nhistonesTurn, radius, histone_angle, "Turntwist")

    chromatinFibrePositions, width = complexGeometryFractal(voxelheightStraight, n)
    output_file_sugar = open(f"sugarPos_{outputFilename}.bin", "wb")
    output_file_hist = open(f"histonePos_{outputFilename}.bin", "wb")

    sugar,histone = createFiles(chromatinFibrePositions, positionStraight, positionTurn, positionTurnTwist, width, voxelheightStraight, voxelheightTurn, f)

    s = struct.pack('f'*len(sugar), *sugar)
    output_file_sugar.write(s)
    output_file_sugar.close()

    s = struct.pack('f'*len(histone), *histone)

    output_file_hist.write(s)
    output_file_hist.close()

elif model == "file":
    if voxelheightTurn != voxelheightStraight:
        raise ValueError("straight and turn voxels must be the same size, using straight voxel parameters")

    print("model type = ",model, file=f)
    print("inputFilename = ",inputFilename, file=f)
    print("voxelheightStraight = ",voxelheightStraight, file=f)
    print("nhistonesStraight = ",nhistonesStraight, file=f)
    print("voxelheightTurn = ",voxelheightTurn, file=f)
    print("nhistonesTurn = ",nhistonesTurn, file=f)
    print("radius = ",radius, file=f)
    print("histone_angle = ",histone_angle, file=f)

    print("", file=f)

    positionStraight = loadOrCreateSolenoids(voxelheightStraight, nhistonesStraight, radius, histone_angle, "Straight")
    positionTurn = loadOrCreateSolenoids(voxelheightTurn, nhistonesTurn, radius, histone_angle, "Turn")
    positionTurnTwist = loadOrCreateSolenoids(voxelheightTurn, nhistonesTurn, radius, histone_angle, "Turntwist")

    chromatinFibrePositions, width = complexGeometryFromFile(voxelheightStraight, inputFilename)
    output_file_sugar = open(f"sugarPos_{outputFilename}.bin", "wb")
    output_file_hist = open(f"histonePos_{outputFilename}.bin", "wb")

    sugar,histone = createFiles(chromatinFibrePositions, positionStraight, positionTurn, positionTurnTwist, width, voxelheightStraight, voxelheightTurn, f)

    s = struct.pack('f'*len(sugar), *sugar)
    output_file_sugar.write(s)
    output_file_sugar.close()

    s = struct.pack('f'*len(histone), *histone)

    output_file_hist.write(s)
    output_file_hist.close()

elif model == "cellCycle":
    if cellCycleType == "interphase":
        # create 2 geometry files. One with single set of interphase DNA, one with duplicate interphase DNA

        # 1 set
        voxelheightTurn = voxelheightStraight
        nhistonesTurn = nhistonesStraight
        print("straight and turn voxels must be the same size, using straight voxel parameters")

        print("model type = ",model, file=f)
        print("cell cycle type = ", cellCycleType, file=f)
        print("inputFilename = ",inputFilename, file=f)
        print("voxelheightStraight = ",voxelheightStraight, file=f)
        print("nhistonesStraight = ",nhistonesStraight, file=f)
        print("voxelheightTurn = ",voxelheightTurn, file=f)
        print("nhistonesTurn = ",nhistonesTurn, file=f)
        print("histone_angle = ",histone_angle, file=f)
        print("radius = ",radius, file=f)

        print("", file=f)

        positionStraight = loadOrCreateSolenoids(voxelheightStraight, nhistonesStraight, radius, histone_angle, "Straight")
        positionTurn = loadOrCreateSolenoids(voxelheightTurn, nhistonesTurn, radius, histone_angle, "Turn")
        positionTurnTwist = loadOrCreateSolenoids(voxelheightTurn, nhistonesTurn, radius, histone_angle, "Turntwist")

        chromatinFibrePositions, width = complexGeometryFromFile(voxelheightStraight, inputFilename)
        output_file_sugar = open(f"sugarPos_{outputFilename}.bin", "wb")
        output_file_hist = open(f"histonePos_{outputFilename}.bin", "wb")

        print("", file=f)
        print("Single DNA file", file=f)

        sugar,histone = createFiles(chromatinFibrePositions, positionStraight, positionTurn, positionTurnTwist, nucleusSize, voxelheightStraight, voxelheightTurn, f)

        s = struct.pack('f'*len(sugar), *sugar)
        output_file_sugar.write(s)
        output_file_sugar.close()

        s = struct.pack('f'*len(histone), *histone)

        output_file_hist.write(s)
        output_file_hist.close()

        # duplicate set

        output_file_sugar = open(f"sugarPos_{outputFilename}_Duplicate.bin", "wb")
        output_file_hist = open(f"histonePos_{outputFilename}_Duplicate.bin", "wb")


        # create 2 copies of the DNA offset so they do not overlap - offset must be more than 60 nm so that fibres do not overlap

        print("", file=f)
        print("Duplicate DNA file, first set of DNA", file=f)

        sugar1,histone1 = createFiles(chromatinFibrePositions, positionStraight, positionTurn, positionTurnTwist, nucleusSize, voxelheightStraight, voxelheightTurn, f, offset = [-0.25*voxelheightStraight]*3)

        print("", file=f)
        print("Duplicate DNA file, second set of DNA", file=f)

        sugar2,histone2 = createFiles(chromatinFibrePositions, positionStraight, positionTurn, positionTurnTwist, nucleusSize, voxelheightStraight, voxelheightTurn, f, offset = [0.25*voxelheightStraight]*3) 

        sugar1.extend(sugar2)
        histone1.extend(histone2)

        s = struct.pack('f'*len(sugar1), *sugar1)
        output_file_sugar.write(s)
        output_file_sugar.close()

        s = struct.pack('f'*len(histone1), *histone1)

        output_file_hist.write(s)
        output_file_hist.close()

    elif cellCycleType == "mitosis":
        # Very compact chromatin density structure with duplicate sets of DNA

        print("model type = ",model, file=f)
        print("cell cycle type = ", cellCycleType, file=f)
        print("simpleMitosis.py", file=f)
        print("voxelheightStraight = ",voxelheightStraight, file=f)
        print("nhistonesStraight = ",nhistonesStraight, file=f)
        print("voxelheightTurn = ",voxelheightTurn, file=f)
        print("nhistonesTurn = ",nhistonesTurn, file=f)
        print("radius = ",radius, file=f)
        print("histone_angle = ",histone_angle, file=f)
        print("densityAim = ",densityAim, file=f)
        print("densityMAim = ",densityMAim, file=f)

        print("", file=f)

        print("Results:\n", file=f)
    
        positions, turns = path_M(densityAim, nucleusSize**3, densityMAim, voxelheightStraight, voxelheightTurn)
        chromatinFibrePositions = complexGeometryFromArray(positions, turns)

        positionStraight = loadOrCreateSolenoids(voxelheightStraight, nhistonesStraight, radius, histone_angle, "Straight")
        positionTurn = loadOrCreateSolenoids(voxelheightTurn, nhistonesTurn, radius, histone_angle, "Turn")
        positionTurnTwist = loadOrCreateSolenoids(voxelheightTurn, nhistonesTurn, radius, histone_angle, "Turntwist")

        print("Mitosis first set of DNA", file=f)
        print("", file=f)

        # fig = plt.figure(figsize=(5, 5))
        # ax = fig.add_subplot(111, projection='3d')
        sugar1,histone1 = createFiles(chromatinFibrePositions, positionStraight, positionTurn, positionTurnTwist, nucleusSize, voxelheightStraight, voxelheightTurn, f, offset = [0,-nucleusSize/4,-nucleusSize/8], mitosis=True, densityAim = densityAim, densityMAim = densityMAim, plot=False) # condensed with 2x DNA

        print("", file=f)
        print("Mitosis second set of DNA", file=f)
        print("", file=f)

        sugar2,histone2 = createFiles(chromatinFibrePositions, positionStraight, positionTurn, positionTurnTwist, nucleusSize, voxelheightStraight, voxelheightTurn, f, offset = [0,+nucleusSize/4,+nucleusSize/8], mitosis=True, densityAim = densityAim, densityMAim = densityMAim) # 2nd set of DNA

        output_file_sugar = open(f"sugarPos_{outputFilename}.bin", "wb")
        output_file_hist = open(f"histonePos_{outputFilename}.bin", "wb")


        sugar1.extend(sugar2)
        histone1.extend(histone2)

        s = struct.pack('f'*len(sugar1), *sugar1)
        output_file_sugar.write(s)
        output_file_sugar.close()

        s = struct.pack('f'*len(histone1), *histone1)
        output_file_hist.write(s)
        output_file_hist.close()
    else:
        raise ValueError("Incorrect cell cycle type")

else:
    raise ValueError("Incorrect model type")

f.close()