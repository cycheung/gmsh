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

R = 8.314
E = 1000

Array = ReadArrayFromFile(filename)
N = len(Array)


filename=open('omega.txt','w')

Columns = [2, 11]
Omega = [0, 0]
for j in range(1,N):
    filename.write(str(Array[j,0]) + ' ')
    dt = Array[j,0] - Array[j-1,0]
    for k in range(0,len(Columns)):
        T = Array[j,Columns[k]-1]
        if(T >= threshold):
            Omega[k] += numpy.exp(-E/R/(T+273.15) * dt)
        filename.write(str(T) + ' ' + str(Omega[k]) + ' ')
    filename.write('\n')




