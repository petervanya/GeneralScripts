import subprocess
import os
import numpy
import math
import matplotlib.pyplot as plt

def plotPDB(_pdbfile):
    def PATHlist():
    # return a list of directories in $PATH
        # fetch current $PATH
        process = subprocess.Popen(['echo $PATH'],stdout=subprocess.PIPE,shell=True)
        (strtmp,error) = process.communicate()

        # split directories into elements of a list
        pathlist = strtmp.split(":")

        # remove '\n' character from last list element
        strtmp = pathlist[len(pathlist)-1].split("\n")
        pathlist[len(pathlist)-1] = strtmp[0]

        return pathlist
    ###################

    def searchpathlist(pathlist,_file):
        for _dir in pathlist:
            if os.path.isdir(_dir):
                files = os.listdir(_dir)
                if _file in files:
                    return True    
        return False
    ################

    def searchplottingtools(pathlist):
        ##### search for vesta #####
        ############################
        
        foundvesta = False
        binaries = ['vesta','Vesta','VESTA']
        binvesta = 'null'

        for _bin in binaries:
            if searchpathlist(pathlist,_bin):
                foundvesta = True
                binvesta = _bin
        
        ###### search for vmd #####
        ############################

        foundvmd = False
        binaries = ['vmd','Vmd','VMD']
        binvmd = 'null'

        for _bin in binaries:
            if searchpathlist(pathlist,_bin):
                foundvmd = True
                binvmd = _bin

        ###### search for pymol #####
        #############################

        foundpymol = False
        binaries = ['pymol','Pymol','PYMOL']
        binpymol = 'null'

        for _bin in binaries:
            if searchpathlist(pathlist,_bin):
                foundpymol = True
                binpymol = _bin
        return [foundvesta,foundvmd,foundpymol] , [binvesta,binvmd,binpymol]
    ########################################################################

    def twodmatplotlib(_pdbfile):
        with open(_pdbfile,'r') as f:
            flines = f.readlines()
        x = []
        y = []
        for line in flines:
            if 'ATOM' in line:
              strtmp = line.split()
              x.append(float(strtmp[3]))
              y.append(float(strtmp[4]))
        plt.plot(x,y,marker='o',linestyle='none',color='black')
        plt.xlabel('x /(\AA)')
        plt.ylabel('y /(\AA)')
        plt.show()
        plt.close()
        plt.cla()
        plt.clf()
        return
    ##########

    
    # get path list from user's current environment
    pathlist = PATHlist()

    # see which visualisation tools are in user's $PATH
    existingprogs,progbinaries = searchplottingtools(pathlist)

    usepython = True    
    # run left-most existing binary first #
    for ia in xrange(len(existingprogs)):
        if existingprogs[ia]:
            os.system(progbinaries[ia]+' '+_pdbfile)
            usepython = False
            break

    # if nothing fancy can be found, use matplotlib
    if usepython:
        twodmatplotlib(_pdbfile)

    return
##########

def createPDB(x,y,z,cell,atomtype,_pdbfile):
    def abcalphabetagamma(cell):
        a = math.sqrt(cell[0]**2+cell[1]**2+cell[2]**2)
        b = math.sqrt(cell[3]**2+cell[4]**2+cell[5]**2)
        c = math.sqrt(cell[6]**2+cell[7]**2+cell[8]**2)

        def angle(vect1,vect2):
            v1 = 0.0
            v2 = 0.0
            for ia in xrange(3):
                v1 = v1 + vect1[ia]**2
                v2 = v2 + vect2[ia]**2
            v1 = math.sqrt(v1)
            v2 = math.sqrt(v2)

            return 180.0/math.pi*math.acos(numpy.dot(vect1,vect2)/(v1*v2))
        ####################################################

        alpha   = angle(cell[3:6],cell[6:9])
        beta    = angle(cell[0:3],cell[6:9])
        gamma   = angle(cell[0:3],cell[3:6])

        return a,b,c,alpha,beta,gamma
    #################################

    # calculate unit cell vector magnitudes and angular displacements
    a,b,c,alpha,beta,gamma = abcalphabetagamma(cell)

    ##### header #####
    ##################
   
    header = "CRYST1"
    header = header+'{:9.3f}'.format(a)+'{:9.3f}'.format(b)+'{:9.3f}'.format(c)
    header = header+'{:7.2f}'.format(alpha)+'{:7.2f}'.format(beta)+'{:7.2f}'.format(gamma)+'\n'

    flines = [header]

    ##### atoms #####
    #################
    
    for ia in xrange(len(x)):
        strtmp=''
        for ib in xrange(4-len(atomtype[ia])):
            strtmp=strtmp+' '
        strtmp2=''
        for ib in xrange(15):
            strtmp2=strtmp2+' '
        flines.append('ATOM  '+'{:5d}'.format(ia+1)+strtmp+atomtype[ia]+strtmp2+\
        '{:8.3f}'.format(x[ia])+'{:8.3f}'.format(y[ia])+'{:8.3f}'.format(z[ia])+'\n')
        
    f = open(_pdbfile,'w')
    f.writelines(flines)
    f.close()
    del flines
    return
##########


#pathlist = PATHlist()
#searchplottingtools(pathlist)
#cell = [1.0 , 0.0 , 0.0 , 1.0 , 1.0 , 0.0 , 0.0 ,0.0 ,5.0]
#createPDB([0.0,0.0],[1.0,1.0],[2.0,2.0],cell,['C','C1'],'scratch.pdb')
#os.system('vesta scratch.pdb')
plotPDB('scratch.pdb')
