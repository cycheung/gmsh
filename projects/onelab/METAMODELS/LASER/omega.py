import sys, string, numpy

##import functions as gf

def DeleteCommentsInFile(listOfLines):
    """Search for lines, starting with a comment and delete this lines
    Needs the list of lines"""
    cleanedLines = []
    for line in listOfLines:
        if '#' not in line.strip().split()[0]:
            cleanedLines.append(line)
    return cleanedLines

def ReadFile(fileName):
    """Returns the content of the file filename line wise"""
    file=open(fileName)
    lines=file.readlines()
    file.close()
    return lines

def ReadArrayFromFile(fileName):
    """Read an numpy array from a file (floats), seperated by white spaces"""
    linesInFile = ReadFile(fileName)
    # Clean comments
    linesInFile = DeleteCommentsInFile(linesInFile)
    numOfColumns = linesInFile[-1].split()
    numpyArray = numpy.zeros((len(linesInFile),len(numOfColumns)),dtype=float)
    for i in range(len(linesInFile)):
        for k in range(len(numOfColumns)):
            numpyArray[i,k] = float(linesInFile[i].split()[k])
    return numpyArray


if len(sys.argv) < 3:
    print "\nUsage: python omega.py filename threshold\n"
    raise RuntimeError, 'Args missing'
filename = sys.argv[1]
threshold = float(sys.argv[2])

R = 8.3143
Ea = 7.82e5 #6.28e5
H = 2.185e124#3.1e98
RCEM = 0.5

Array = ReadArrayFromFile(filename)
N = len(Array)

filename=open('omega.txt','w')
filename2=open('cem.txt','w')

Columns = [2, 11]
Omega = [0, 0]
Cem = [0, 0]
for j in range(1,N):
    filename.write(str(Array[j,0]) + ' ')
    filename2.write(str(Array[j,0]) + ' ')
    dt = Array[j,0] - Array[j-1,0]
    for k in range(0,len(Columns)):
        T = Array[j,Columns[k]-1]
        if(T >= threshold):
            myexp = numpy.exp(-Ea/(R*(T+273.15)))
            omega_dt = H * myexp * dt
            Omega[k] += H * myexp * dt
            Cem[k]   += RCEM**(threshold-T)*dt  
        filename.write(str(T) + ' ' + str(Omega[k]) + ' ')
        filename2.write(str(T) + ' ' + str(Cem[k]) + ' ')
    filename.write('\n')
    filename2.write('\n')

print 'Omega =', Omega[0], Omega[1]
print 'Cem =', Cem[0], Cem[1]




