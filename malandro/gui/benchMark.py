'''
Scripts for bench marking the results of the automated assignment.
Normally not used. Just comparing the assignment the algorithm
procuced with the original manually made one.
'''

import csv
from time import localtime


def createBenchmark(results):
    '''
    Output:name of the proh
    Compare the assignment made by malandro to the manual assignment.
    '''

    dirName = '/home/joren/malandroBenchmarks/'
    ts = timeStamp()
    #projectName = project.shortName

    filePath = '%s%s_%s' % (dirName, 'benchMark', ts)

    with open(filePath, 'wb') as csvfile:
        writer = csv.writer(csvfile, dialect='excel')

        writer.writerow([ts])

        spectraNames = [spectrum.getName()
                        for spectrum in results.getSpectra()]

        writer.writerow(spectraNames)
        writer.writerow(['% agree',
                         '% disagree',
                         '% agreeing Joker',
                         '% disagreeing Joker'])

        residues = results.getChain().getResidues()

        for residue in residues:

            writer.writerow(calculatePercentages(residue))


def calculatePercentages(res):
    '''Calculate how often 4 types of spin realSpinSystems
       assignments are made to this residue:
           1) agreeing with solution in project
           2) highest scoring assignment that does
              not agree with the assignment in the
              project.
           3) agreeing joker. A joker is placed on
              a residue that does not have any spin
              system assignment in the project either.
           4) disagreeing joker. A joker is placed
              on a residue that does have a spin
              system assignment in the project.
        args: res (src.cython.malandro.Residue)
        return: list [sequence code,
                      agree percentage,
                      disagree percentage,
                      agreeing joker percentage,
                      disagreeing joker percentage]
    '''

    data = []
    jokers = []
    realSpinSystems = []

    resSeqCode = res.getSeqCode()

    if res.getCcpnResidue().findFirstResonanceGroup():

        manualAssigned = True

    else:

        manualAssigned = False

    solutions = res.getSolutions()
    Nsolutions = len(solutions)
    solutionSet = set(solutions)

    for spinSys in solutionSet:

        if spinSys.getIsJoker():

            jokers.append(spinSys)

        else:

            realSpinSystems.append(spinSys)

    for spinsys in realSpinSystems:

        assignmentPercentage = int(
            float(solutions.count(spinsys)) / Nsolutions * 100.0)

        data.append((assignmentPercentage, spinsys))

    if jokers:

        NumberOfAssignmentsToJoker = 0

        for spinSys in jokers:

            NumberOfAssignmentsToJoker += solutions.count(spinSys)

        assignmentPercentage = int(
            float(NumberOfAssignmentsToJoker) / Nsolutions * 100.0)

        data.append((assignmentPercentage, None))

    data = sorted(data, reverse=True)

    bestWrongPercentage = 0
    agreePercentage = 0
    goodJoker = 0
    badJoker = 0

    foundBad = False

    for percent, spinSys in data:

        if spinSys:
            # print '----'
            # print spinSys.getCcpnResonanceGroup().residue.seqCode
            # print resSeqCode
            if spinSys.getCcpnResonanceGroup().residue and spinSys.getCcpnResonanceGroup().residue.seqCode == resSeqCode:

                agreePercentage = percent

            elif not foundBad:

                bestWrongPercentage = percent
                foundBad = True

        else:  # Joker

            if manualAssigned:

                if not foundBad:

                    badJoker = percent
                    foundBad = True

            else:

                goodJoker = percent

    return [resSeqCode,
            agreePercentage,
            bestWrongPercentage * -1,
            goodJoker, badJoker * -1]


def timeStamp():
    '''Returns string timestamp.'''
    t = localtime()

    return '%s_%s_%s_%s_%s' % (str(t.tm_year), str(t.tm_mon), str(t.tm_mday), str(t.tm_hour), str(t.tm_min))
