"""Convert a madx beamline to G4Beamline format

This scirpt converts MADX outputs, to a .g4bl file.
Currently only supports markers, drifts, sector dipoles and quads,
not for higher order magnets. They will be added later.

The MADsequenceFile is the sequence file saved by the "save, sequence"
command. It requires that only one sequence is saved in the save command.

The TWISSfile is the TWISS file saved by the "twiss, save, file=TWISSfile"
command. It contains the initial parameters of the beam.

The output file will be saved in the current working directory.
The file name will use the same name of the sequence and Twiss files if
their names (ignoring the last extension) are the same. If not, the name (
still ignoring the last ext.) of the sequence file will be used.

Then follow the prompt instructions.

Copyright: Ao Liu (frankliuao), 2015
"""

from __future__ import division
import math
import os
import sys
import time
import string
import re
import copy
import urllib2


__version__ = "1.01"
__author__ = "frankliuao"
__all__ = ['mad2g4bl']
__history__ = {"07/08/15": "Modified the script to make it work as a module.",
               "09/10/15": "Changed the conversion of a dipole, "
                           "added cornerarc; Followed PEP8 more closely"}
__help__ = """

The MADsequenceFile is the sequence file saved by the "save, sequence"
command. It requires that only one sequence is saved in the save command.

The TWISSfile is the TWISS file saved by the "twiss, save, file=TWISSfile"
command. It contains the initial parameters of the beam.

The output file will be saved in the current working directory.
The file name will be the same of the sequence file with the last extension
name replaced by ".g4bl". e.g. FOO.sequence --> FOO.g4bl

Then follow the prompt instructions.

Copyright: Ao Liu (frankliuao), 2015
"""


######
# Define function g4blwrite:
def g4blWrite(Cont, kind, currPos, tempPos):
    global elementIndex

    # Evaluate all the variables:
    for mm in range(0,len(SeqCont),1):
        SeqContCp = re.sub('[; ]','',SeqCont[mm])
        SeqContCp = re.split(':=|=',SeqContCp)
        if len(SeqContCp)==2 and len(re.split('\,',re.sub('[; ]','',SeqCont[mm])))==1:
            exec string.join(SeqContCp,'=')
    if kind == 'b':
        ContPrime = Cont
        Cont = re.sub('[ {};]','',Cont)
        ContSplit = string.split(Cont,',')
        name = Cont.split(':')[0]
        length, field, quadfield, tilt = (0,)*4
        e1flag, e2flag, apertFlag = (0,)*3
        Ax, Ay = (0,)*2
        for ll in range(1,len(ContSplit),1):
            # param is the parameter name of this magnet:
            param = re.split(':=|=',ContSplit[ll])[0].replace(' ','')
            # param2 is the value for param:
            param2 = re.split(':=|=',ContSplit[ll])[1].replace(' ','')
            # If param2 is not a number, it's another independent variable, need to find it;
            try:
                param2 = float(param2)
            except:
                # If param2 cannot be converted to a float, either it is the apertype or aperture,
                # Or it is another variable name;
                if param == 'apertype':
                    if param2 == 'circle':
                        apertFlag = 1
                    elif param2 == 'ellipse':
                        apertFlag = 2
                    elif param2 == 'rectangle':
                        apertFlag = 3
                    else:
                        print 'Error:Non-defined apertype %s'%param2
                        sys.exit()
                        
                elif param == 'aperture':
                    param2 = param2.split(',')
                    apertSize = param2
                else:
                    try:
                        param2 = eval(param2)
                    except:
                        print 'Error: the variable %s'%param2+' is not defined!'
                        sys.exit()
                    
            # Set the direct parameters:
            if param == 'l':
                length = param2*1000
            elif param == 'angle':
                angle = param2
            elif param == 'k1':
                quadfield = param2*rigidity
            elif param == 'tilt':
                tilt = param2
            elif param == 'e1':
                if param2 != 0:
                    e1Flag = 1
                    if angle >= 0:
                        e1 = param2*180/math.pi
                    else:
                        e1 = -1*param2*180/math.pi
            elif param == 'e2':
                if param2 != 0:
                    e2Flag = 1
                    if angle >= 0:
                        e2 = param2*180/math.pi
                    else:
                        e2 = -1*param2*180/math.pi
            elif param!='apertype' and param!='aperture':
                print 'Error: attribute %s of %s'%(param,name)+' is not supported!'
                print 'Contact the author (www.frankliuao.com) to add.'
                sys.exit()
            
        # Set the other parameters:
        field = angle*rigidity/(length/1000)
        centerR = math.fabs(length/angle)
        if apertFlag == 1:
            # 1 corresponds to a circle, only one parameter is needed;
            Ax, Ay = (float(apertSize[0])*1000,)*2
        elif apertFlag == 2:
            # 2 corresponds to an ellipse:
            Ax = float(apertSize[0])*1000
            Ay = float(apertSize[1])*1000
        elif apertFlag == 3:
            # 3 corresponds to a rectangle:
            Ax = float(apertSize[0])*1000
            Ay = float(apertSize[1])*1000
        
        # Write the informaiton in G4bl format:
        # First the drift before this dipole:
        if math.fabs(currPos+length/2-tempPos)>1:
            fileG4bl.write('# The drift before %s'%ContPrime+'\n')
            fileG4bl.write('param LO_to_%s'%name+'_%s'%(str(elementIndex))+'=%s'%(str(math.fabs(currPos+length/2-tempPos)))+'\n')
            fileG4bl.write('tubs O_to_%s'%name+'_%s'%(str(elementIndex))+' innerRadius=$IR outerRadius=$IR*1.2 length=$'+'LO_to_%s'%name+'_%s'%(str(elementIndex)))
            fileG4bl.write(' material=nusaver\\'+'\n'+'color=1,1,1,0.5\n')
            fileG4bl.write('place '+'O_to_%s'%name+'_%s'%(str(elementIndex))+' z=$Z front=1\n')
            fileG4bl.write('param Z=$Z+$LO_to_%s'%name+'_%s'%(str(elementIndex))+'\n\n')
            elementIndex += 1
        
        # If aperture set, limit the aperture of this dipole;
        if Ax!=0 and Ay!=0:
            Amin = min(Ax,Ay)
            fileG4bl.write('param IR=%s'%(str(Amin))+'\n')
        
        # Now the dipole itself:
        fileG4bl.write('# %s'%ContPrime+'\n' +\
            'param L%s'%name+'_%s'%(str(elementIndex)) + '=%s'%str(length) + '\n' + \
            'param B%s'%name+'_%s'%(str(elementIndex)) + '=%s'%str(field) + '*%s'%str(particleSign*-1) + '\n' + \
            'param Bend%s'%name + '_%s'%(str(elementIndex)) + '=%s'%str(angle*180/math.pi) + '\n')
        fileG4bl.write('param centerR%s'%name + '_%s'%(str(elementIndex)) + '=%s'%str(centerR) + '\n' + \
            'param IR%s'%name + '_%s'%(str(elementIndex)) + '=$centerR%s'%name + '_%s'%(str(elementIndex)) + '-$IR ' + \
            'OR%s'%name + '_%s'%(str(elementIndex)) + '=$centerR%s'%name + '_%s'%(str(elementIndex)) + '+$IR \n')
        fileG4bl.write('idealsectorbend %s'%name + '_%s'%(str(elementIndex)) + ' fieldCenterRadius=$centerR%s'%name + '_%s'%(str(elementIndex)) + ' ' + \
            'fieldInnerRadius=$IR%s'%name + '_%s'%(str(elementIndex)) + ' ' + \
            'fieldOuterRadius=$OR%s'%name + '_%s'%(str(elementIndex)) + ' angle=$Bend%s'%name + '_%s'%(str(elementIndex)) + '\\' + '\n' + \
            'fieldHeight=$IR*2' + ' ironInnerRadius=$IR%s'%name + '_%s'%(str(elementIndex)) + \
            '-0.2*$IR ' + 'ironOuterRadius=$OR%s'%name+'_%s'%(str(elementIndex)) + '+0.2*$OR' + ' ironHeight=$IR*2.4' + '\\' + '\n' + \
            'By=$B%s'%name+'_%s'%(str(elementIndex)) + ' ironColor=0,0,1,0.5 ironMaterial=nusaver\n')
        fileG4bl.write('place %s'%name + '_%s'%(str(elementIndex)) + ' z=$Z\n' + \
            'cornerarc z=$Z angle=$Bend%s'%name + '_%s'%(str(elementIndex)) + ' centerRadius=$centerR%s'%name + '_%s\n'%(str(elementIndex)) +\
            'param Z=$Z+$L%s'%name + '_%s'%(str(elementIndex)) + '\n\n')
        
        
        # Return the half length of this dipole so that the current position can be set:
        return length/2
        
    if kind == 'q':
        ContPrime = Cont
        Cont = re.sub('[ {};]','',Cont)
        ContSplit = string.split(Cont,',')
        name = Cont.split(':')[0]
        length, gradient, apertFlag, tilt = (0,)*4
        Ax, Ay = (0,)*2
        for ll in range(1,len(ContSplit),1):
            # param is the parameter name of this magnet:
            param = re.split(':=|=',ContSplit[ll])[0].replace(' ','')
            # param2 is the value for param:
            param2 = re.split(':=|=',ContSplit[ll])[1].replace(' ','')
            # If param2 is not a number, it's another independent variable, need to find it;
            try:
                param2 = float(param2)
            except:
                # If param2 cannot be converted to a float, either it is the apertype or aperture,
                # Or it is another variable name;
                if param == 'apertype':
                    if param2 == 'circle':
                        apertFlag = 1
                    elif param2 == 'ellipse':
                        apertFlag = 2
                    elif param2 == 'rectangle':
                        apertFlag = 3
                    else:
                        print 'Error:Non-defined apertype %s'%param2
                        sys.exit()
                elif param == 'aperture':
                    param2 = param2.split(',')
                    apertSize = param2
                else:
                    try:
                        param2 = eval(param2)
                    except:
                        print 'Error: the variable %s'%param2+' is not defined!'
                        sys.exit()
                        
            # Set the direct parameters:
            if param == 'l':
                length = param2*1000
            elif param == 'k1':
                gradient = param2*rigidity
            elif param == 'tilt':
                tilt = param2
            elif param!='apertype' and param!='aperture':
                print 'Error: attribute %s of %s'%(param,name)+' is not supported!'
                print 'Contact the author (www.frankliuao.com) to add.'
                sys.exit()
            
            # Set the other parameters:
            if apertFlag == 1:
                # 1 corresponds to a circle, only one parameter is needed;
                Ax, Ay = (float(apertSize[0])*1000,)*2
            elif apertFlag == 2:
                # 2 corresponds to an ellipse:
                Ax = float(apertSize[0])*1000
                Ay = float(apertSize[1])*1000
            elif apertFlag == 3:
                # 3 corresponds to a rectangle:
                Ax = float(apertSize[0])*1000
                Ay = float(apertSize[1])*1000
            
        
        # Write the informaiton in G4bl format:
        # First the drift:
        if math.fabs(currPos+length/2-tempPos)>1:
            fileG4bl.write('# The drift before %s'%name+'\n')
            fileG4bl.write('param LO_to_%s'%name+'_%s'%(str(elementIndex))+'=%s'%(str(math.fabs(currPos+length/2-tempPos)))+'\n')
            fileG4bl.write('tubs O_to_%s'%name+'_%s'%(str(elementIndex))+' innerRadius=$IR outerRadius=$IR*1.2 length=$'+'LO_to_%s'%name+'_%s'%(str(elementIndex)))
            fileG4bl.write(' material=nusaver color=1,1,1,0.5\n')
            fileG4bl.write('place '+'O_to_%s'%name+'_%s'%(str(elementIndex))+' z=$Z front=1\n')
            fileG4bl.write('param Z=$Z+$LO_to_%s'%name+'_%s'%(str(elementIndex))+'\n\n')
            elementIndex += 1
            
        # If aperture set, limit the aperture of this quadrupole;
        if Ax!=0 and Ay!=0:
            Amin = min(Ax,Ay)
            fileG4bl.write('param IR=%s'%(str(Amin))+'\n')
            
        # Then the quadrupole:
        fileG4bl.write('# %s\n'%ContPrime)
        fileG4bl.write('param L%s'%name+'_%s'%(str(elementIndex)) + '=%s'%(str(length)) + ' G%s'%name + '_%s'%(str(elementIndex)) + '=%s*%s \n'%(str(gradient),str(particleSign)))
        fileG4bl.write('genericquad %s'%name + '_%s'%str(elementIndex) + ' apertureRadius=$IR ironRadius=$IR*1.2 fieldLength=$L%s'%name + '_%s'%str(elementIndex) + '\\' + '\n')
        fileG4bl.write('ironLength=$L%s'%name + '_%s'%str(elementIndex) + ' ironMaterial=nusaver fringe=0 gradient=$G%s'%name + '_%s'%str(elementIndex)+' ')
        fileG4bl.write('ironColor=0.7,0.3,0.3,0.6 \n')
        fileG4bl.write('place %s'%name + '_%s'%str(elementIndex) + ' z=$Z front=1 \n')
        fileG4bl.write('param Z=$Z+$L%s'%name + '_%s'%str(elementIndex) + '\n\n')
        elementIndex += 1
        
        # Return the half length of this dipole so that the current position can be set:
        return length/2
        
    if kind == 'm':
        ContPrime = Cont
        Cont = re.sub('[ {};]','',Cont)
        ContSplit = string.split(Cont,',')
        name = Cont.split(':')[0]

        # First the drift before this marker:
        if math.fabs(currPos-tempPos)>1:
            fileG4bl.write('# The drift before %s'%name+'\n')
            fileG4bl.write('param LO_to_%s'%name+'_%s'%(str(elementIndex))+'=%s'%(str(math.fabs(currPos-tempPos)))+'\n')
            fileG4bl.write('tubs O_to_%s'%name+'_%s'%(str(elementIndex))+' innerRadius=$IR outerRadius=$IR*1.2 length=$'+'LO_to_%s'%name+'_%s'%(str(elementIndex)))
            fileG4bl.write(' material=nusaver\\'+'\n'+'color=1,1,1,0.5\n')
            fileG4bl.write('place '+'O_to_%s'%name+'_%s'%(str(elementIndex))+' z=$Z front=1\n')
            fileG4bl.write('param Z=$Z+$LO_to_%s'%name+'_%s'%(str(elementIndex))+'\n\n')
            elementIndex += 1

        fileG4bl.write('# %s\n'%ContPrime)
        fileG4bl.write('virtualdetector det_%s'%name + ' length=0.01 radius=$IR material=Vacuum format=ascii file=detect_%s'%name + '_%s'%str(elementIndex) + ' referenceParticle=1\n')
        fileG4bl.write('place det_%s'%name + ' z=$Z+0.02 front=1\n')
        fileG4bl.write('param Z=$Z+0.05 \n\n')
        elementIndex += 1
        return 0.05

def getPDGdata():
    """Get the PDG data from the web server"""
    pdgDataHandle = urllib2.urlopen('http://frankliuao.com/downloads/Codes/pdg_py.dict')
    pdgDataCont = pdgDataHandle.read()
    pdgDataHandle.close()
    exec pdgDataCont
    return pdgData, pdgMap

def Twissdict(TwissCont):
    """Extract the contents in the Twiss file as a dictionary.

    Given a MAD-X Twiss file, extract the Twiss parameters and organize them
    in a dictionary format.
    """
    for ii in range(len(TwissCont)):
        if TwissCont[ii][0:2] == "* ":
            startLine = ii
            paramNames = map(lambda x: x.lower(),
                             TwissCont[ii].lstrip("* ").split()[1:])
            break

    TwissDict = dict.fromkeys(paramNames,[])

    for ii in range(startLine, len(TwissCont)):
        TwissThisLine = TwissCont[ii].split()
        for jj in range(len(paramNames)):
            TwissDict[paramNames[jj]].append(TwissCont[ii][jj])

    # Convert the Twiss table to its corresponding format
    for ii in paramNames:
        TwissDict[ii] = map(lambda x: x.strip('"'),TwissDict[ii]) if ii is \
                                                                     "name" \
                   else map(float,TwissDict[ii])

    return TwissDict

def setparam(inputString):
    """Extract the beam parameter definitions in the Twiss file.

    In the Twiss file of MAD-X, some beam parameters are defined, e.g.
    @ MASS             %le              0.1397
    This code process these parameter definition lines and extract the
    useful information.
    """
    if inputString[0] == "@":
        if inputString.count('"') == 0:
            inputString = inputString.split()
            return "madx_"+inputString[1].lower()+"="+inputString[-1]
        else:
            return "pass"
    else:
        return "pass"

def main(sequence=None, twiss=None):
    """Main function to convert the format.
    """
    global seqCont
    global TwissCont

    currentDir = os.getcwd()    # Current working directory;
    print "Welcome to frankliuao's code!\nwww.frankliuao.com\n"

    # Sequence and Twiss file names with the last extension removed:
    seqName = '.'.join(os.path.basename(sequence).split('.')[:-1])
    g4blName = seqName


    seqHandle = open(sequence,'r')
    seqCont = seqHandle.readlines()
    seqHandle.close()
    seqCont = map(lambda x:x.rstrip('\n'),seqCont)
    TwissHandle = open(Twiss,'r')
    TwissCont = TwissHandle.readlines()
    TwissHandle.close()

    # Find the key parameters in the Twiss file and set them as variables.
    for ii in [setparam(jj) for jj in TwissCont]:
        exec(ii)
    # Convert the units.
    mass, energy, p = madx_mass*1000, madx_energy*1000, madx_pc*1000
    ex, ey = madx_ex*1000, madx_ey*1000
    #
    kE = (energy-mass)
    rigidity = p*3.33564/1000
    sigX = math.sqrt(betx*ex)
    sigXp = math.sqrt((1+alfx*alfx)/betx*ex)
    sigY = math.sqrt(bety*ey)
    sigYp = math.sqrt((1+alfy*alfy)/bety*ey)
    beta = math.sqrt(1-1/(gamma*gamma))

    # Get the dictionary of all the Twiss parameters:
    TwissDict = Twissdict(TwissCont)
    #
    if TwissDict.has_key('betx'): betx_0 = TwissDict['betx'][0]*1000
    if TwissDict.has_key('alfx'): alfx_0 = TwissDict['alfx'][0]
    if TwissDict.has_key('bety'): bety_0 = TwissDict['bety'][0]*1000
    if TwissDict.has_key('alfy'): alfy_0 = TwissDict['alfy'][0]
    if TwissDict.has_key('dx'): dx_0 = TwissDict['dx'][0]*1000
    if TwissDict.has_key('dpx'): dpx_0 = TwissDict['dpx'][0]
    if TwissDict.has_key('dy'): dy_0 = TwissDict['dy'][0]*1000
    if TwissDict.has_key('dpy'): dpy_0 = TwissDict['dpy'][0]
    if TwissDict.has_key('x'): x_0 = TwissDict['x'][0]*1000
    if TwissDict.has_key('y'): y_0 = TwissDict['y'][0]*1000
    if TwissDict.has_key('t'): z_0 = TwissDict['t'][0]*3e8*1000
    if TwissDict.has_key('px'): px_0 = TwissDict['px'][0]
    if TwissDict.has_key('py'): py_0 = TwissDict['py'][0]

    # Temporary use the following:
    betx, bety, alfx, alfy, dx, dy, dpx, dpy = betx_0, bety_0, alfx_0, \
                                               alfy_0, dx_0, dy_0, dpx_0, dpy_0
    initX, initY, initZ = x_0, y_0, z_0
    thetaX, thetaY = px_0, py_0


    # User input: particle name; From the name get the particle charge.
    pdgData, pdgMap = getPDGdata()
    while True:
        particle = raw_input('Please input your particle. See G4bl manual for supported names. e.g. pi+, mu+, e+, ...\n')
        try:
            PDGid = pdgMap[particle]
            break
        except KeyError:
            continue

    particle = ''
    PDGid = {'e+':1, 'e-':-1, 'anti_nu_e':1, 'anti_nu_mu':1, 'anti_nu_tau':1, 'gamma':1, 'gluon':1, 'kaon+':1,'kaon-':-1,
                 'kaon0':1, 'mu+':1, 'mu-':-1, 'neutron':1, 'nu_e':1, 'nu_mu':1, 'nu_tau':1, 'pi+':1, 'pi-':-1, 'pi0':1,
                 'proton':1, 'rho+':1, 'rho-':-1, 'rho0':1, 'tau+':1, 'tau-':-1, 'anti_proton':-1}


    # Set the default aperture radius in mm:
    while True:
        apert = raw_input('Please input the default beam pipe RADIUS: (in mm)\n')
        try:
            apert = float(apert)
            break
        except:
            print 'Sorry, come again.'

    # Set the fringe field flag:
    while True:
        fringeFlag = raw_input('Do you want fringe fields?\n'+'y/n\n')
        if fringeFlag == 'y':
            fringeFlag = True
            break
        elif fringeFlag == 'n':
            fringeFlag = False
            break
        else:
            print 'Sorry, come again. y/n?'

    # Set the beam definition line on line 12:
    while True:
        beamFlag = raw_input('Do you want your beam to be Gaussian distributed over the transverse phase space, based on the TWISS parameters?\n'+'y/n\n')
        if string.lower(beamFlag) == 'y':
            beamFlag = True
            print 'The initial beam parameters will be set based on your TWISS BETA0 configurations. You may change it later on with other beam definitions'
            break
        elif string.lower(beamFlag) == 'n':
            beamFlag = False
            print 'Beam not defined yet, you may do it later on line 12.'
            break
        else:
            print 'Sorry, come again. y/n?'

    # Usually when a genericbend is placed, 5 mm is added to its physical length, but total length should conserve;
    dipolelabel = 0                        # To conserve the length of beamline;
    print "Running the conversion, please wait..."

    # Ready to write the .g4bl file.
    global fileG4bl
    if os.path.exists(seqName+'.g4bl'):
        os.system('''rm %s'''%(seqName+'.g4bl'))
    fileG4bl = open(seqName+'.g4bl','w')

    # Write some pre-definition in it first:
    fileG4bl.write('# G4Beamline version of %s and %s, the rigidity is $rigid T.m.\n\n'%(arg1,arg2)+'''# Ao Liu (ID:frankliuao)'s mad2g4bl.csh; Version 1.0.\n\n''')
    fileG4bl.write('# Use the low energy physics list \n'+'physics QGSP_BERT list=1 disable=Decay\n \n'+'param IR=%s maxStep=1 eventTimeLimit=400\n\n'%(str(apert)))
    fileG4bl.write('# Notice!!!!!!!\n'+'# Define your beam here, on line 12.\n')
    # Defining the beam if the user requested so:
    if beamFlag:
        fileG4bl.write('beam gaussian particle=%s'%particle + ' nEvents=5000 meanP=%s '%str(p)+'\\\n')
        fileG4bl.write('beamX=%s'%str(initX) + ' beamY=%s'%str(initY)+' beamZ=%s'%str(initZ)+' sigmaP=-%s'%str(sigE/beta/beta*p)+' sigmaX=%s'%str(sigX)+' sigmaY=%s'%str(sigY)+' sigmaZ=0 '+'\\'+'\n')
        fileG4bl.write('meanXp=%s'%str(thetaX)+' meanYp=%s'%str(thetaY)+' meanT=0 sigmaXp=%s'%str(sigXp)+' sigmaYp=%s'%str(sigYp)+' sigmaT=0'+'\n\n\n')
    else:
        fileG4bl.write('\n\n\n')
    #
    fileG4bl.write('output %s'%seqName+'.out \n\n')
    fileG4bl.write('# Define a reference particle:\n'+'reference particle=%s'%particle + ' referenceMomentum=%s'%(str(p)) + ' beamZ=0\n\n')
    fileG4bl.write('# Change color for particles:\n'+'particlecolor mu+=1,0,0 pi+=0,1,1 e+=0,1,0\n\n')
    fileG4bl.write('trackcuts keep=%s,nu_mu,anti_nu_mu,nu_e,anti_nu_e\n\n'%particle)
    fileG4bl.write('# Redirect the loss information to a separate file\n'+'beamlossntuple losscheck format=ascii file=%s\n\n'%(arg1Prime))
    fileG4bl.write('# Define a new material for killing particles hitting it except for neutrinos.\n'+'material nusaver Fe,1.0 density=7.874 keep=nu_mu,nu_e,anti_nu_mu,anti_nu_e\n\n')
    fileG4bl.write('param histoFile=%s \n\n'%(seqName))
    fileG4bl.write('# Put a particle filter here to stop particles circulating more than n turns as User defined;\n')
    fileG4bl.write('particlefilter ringfilter radius=$IR length=0.01 kill=mu+,e+,pi+ material=Vacuum color=0,0,1 nWait=2 referenceWait=2 steppingVerbose=1 \n')
    fileG4bl.write('place ringfilter z=0.02\nparam Z=0.05\n\n')

    # Record the element index:
    global elementIndex
    elementIndex = 1

    ########
    # Deal with the sequence file:
    # Find where the sequence definition lines start and end:
    for kk in range(0,len(SeqCont),1):
        lineTemp = string.split(SeqCont[kk])
        # Find where it starts:
        if lineTemp[1]=='sequence,':
            # Set the current longitudinal position to 0:
            currPos = 0
            ll = kk+1
            # From now until the definition ends:
            while SeqCont[ll]!='endsequence;\n':
                # Get the beamline element name and its longitudinal position:
                SeqTemp = SeqCont[ll].replace(';','').replace('\n','')
                SeqTempSplit = string.split(SeqTemp,', at = ')
                SeqTempName = SeqTempSplit[0]
                SeqTempPos = float(SeqTempSplit[1])*1000
                # Then look for this element name in the sequence file:
                for mm in range(0,len(SeqCont),1):
                    NameFind = string.split(SeqCont[mm].replace(';','').replace('\n',''),': ')
                    if NameFind[0]==SeqTempName:
                        # NameFound is the definition line of this element:
                        NameFound = string.split(SeqCont[mm].replace(';','').replace('\n',''),',')
                        NameType = string.split(NameFound[0],': ')
                        NameType = NameType[1].replace(' ','')
                        if NameType == 'quadrupole':
                            halfLength = g4blWrite(SeqCont[mm],'q',currPos,SeqTempPos)
                        elif NameType == 'sbend':
                            halfLength = g4blWrite(SeqCont[mm],'b',currPos,SeqTempPos)
                        elif NameType == 'marker':
                            halfLength = g4blWrite(SeqCont[mm],'m',currPos,SeqTempPos)
                        else:
                            print '*** Error: the element type %s'%NameType+' is not supported. Contact the author please***'
                            print '*** www.frankliuao.com ***'
                        # Move the current position to the tail of this element:
                        currPos = SeqTempPos + halfLength
                        # Break the loop of NameFind:
                        break

                # Go to the next line (while SeqCont[ll]!='endsequence;':)
                ll += 1
            # Once gone through the sequence definitions, the loop can be broken:
            break
    ####
    # Finally, we need to add a virtualdetector whatsoever to the end of this beamline:
    fileG4bl.write('\n\n'+'# Beam tranfer monitoring virtualdetector'+'\n')
    fileG4bl.write('virtualdetector detEnd length=0.01 radius=$IR material=Vacuum format=ascii file=detEnd referenceParticle=1\nplace detEnd z=$Z+0.01 front=1')
    print '\n'+'Conversion complete!'

def mad2g4bl(sequence=None, Twiss=None):
    """Just an alias of the main function.
    """
    if sequence != None and Twiss != None:
        main(sequence,Twiss)
    else:
        print "Usage: mad2g4bl(sequence=FILEPATH, Twiss=FILEPATH2)"+__help__
        print '---Error: MAD-X sequence file and Twiss file must be ' \
              'provided! ---'
        sys.exit(1)