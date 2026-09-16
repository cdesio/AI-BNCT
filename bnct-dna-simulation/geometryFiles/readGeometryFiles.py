import numpy as np

def readHistones(path):
    try:
        with open(path, "rb") as f:
            numpy_data = np.fromfile(f,np.float32)
    except IOError:
        print('Error While Opening the file!')  

    data = numpy_data.reshape(int(len(numpy_data)/6),6)

    return data

def readSugar(path):
    try:
        with open(path, "rb") as f:
            numpy_data = np.fromfile(f,np.float32)
    except IOError:
        print('Error While Opening the file!')  

    data = numpy_data.reshape(int(len(numpy_data)/12),12)

    return data
