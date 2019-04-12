import sys
import os
# cheap hack until I can figure out how to do this properly in
# python 3.6
sys.path.append(os.path.split(os.path.split(os.path.realpath(__file__))[0])[0])

import random
from decimal import Decimal

import backend.Structures.pssmIO as PSSMIO
import backend.Structures.acidMassTableIO as AcidMassTableIO
import backend.Structures.spectrumIO as SpectrumIO
import backend.Structures.pssm as PSSM
import backend.Structures.acidMassTable as AcidMassTable
import backend.Structures.spectrum as Spectrum
import backend.PreProcessing.acidConversion as AcidConversion
import backend.PreProcessing.spectrumConversion as SpectrumConversion
import backend.PostProcessing.processResults as ProcessResults
from backend.constants import *
from backend.userInput import *


def massToleranceMaxDiff(calculatedMass, massTolerance):

  minExp = int(Decimal(calculatedMass) / (Decimal(massTolerance) + Decimal(1)))
  maxExp = int(Decimal(calculatedMass) / (Decimal(1) - Decimal(massTolerance)))

  return (minExp, maxExp)


def getDecoys(mass, peptideDict, massTolerance):
  possibles = []

  minDiff, maxDiff = massToleranceMaxDiff(mass, massTolerance)

  for mass in range(minDiff, maxDiff+1):
    if mass in peptideDict:
      possibles += peptideDict[mass]

  if len(possibles) <= 10:
    return possibles

  else:
    return random.sample(possibles, 10)


def calculatedAminoScore(peptide, aminoScoreList):
  aminoScores = 0
  aminoDict = aminoScoreList[len(peptide) - 9]
  for position, amino in zip(range(len(peptide)), peptide):
    aminoScores += aminoDict[amino][position]
  return (float(aminoScores) / len(peptide))


def convertPeptide(peptide,
                   acidConversion,
                   conversionOrder):
  for conversion in conversionOrder:
    newPep = peptide.replace(conversion,
                             acidConversion[conversion])
  return newPep


def getPepMass(peptide,
               acidMassTable):
  mass = 0
  for acid in peptide:
    mass += AcidMassTable.getMass(acidMassTable, acid)
  return mass


def fileToDict(pepFile,
               precision,
               minP,
               maxP,
               acidConversion,
               conversionOrder,
               acidMassTable,
               reverse):

  pepDict = {}

  for pepLine in pepFile:
    if pepLine == "":
      break
    peptide = pepLine.strip()
    if len(peptide) < minP or len(peptide) > maxP:
      continue
    peptide = convertPeptide(peptide,
                             acidConversion,
                             conversionOrder)
    if reverse == True:
      peptide = peptide[::-1]
    if peptide == '':
      continue
    pepMass = getPepMass(peptide, acidMassTable)
    if pepMass in pepDict:
      pepDict[pepMass].append(peptide)
    else:
      pepDict[pepMass] = [peptide]

  return pepDict


def getPSSMScore(peptide,
                 pssmTitle,
                 allPSSM):

  if len(peptide) in PSSM.getIgnoredLengths(pssmTitle, allPSSM):
    return 0

  score = 0
  matrix = PSSM.getMatrixOfLength(pssmTitle, len(peptide), allPSSM)
  for position, acid in enumerate(peptide):
    score += matrix[acid][position]

  return score


def programUsageOutput():
  print("\nUsage: dbPepToScores.py SPECDIR ACIDMASSFILE PSSMDIR DCYDBDIR REVERSE\n\n")

  padding = '{:<20}{}'
  print("Positional arguments:")
  print(str.format(padding, 'SPECDIR',
        'A directory containing 1 or more spectrum files'))
  print(str.format(padding, 'ACIDMASSFILE',
        'A file containing mass spectrometry data'))
  print(str.format(padding, 'PSSMDIR',
        'A directory containing a positional scoring matrix'))
  print(str.format(padding, 'DCYDBDIR',
        'A directory containing a databases of decoy peptides'))
  print(str.format(padding, 'REVERSE',
        'Should decoy peptides be reversed?'))
  print('{:<20}{}'.format('','"t" or "true" or "f" or "false" (not case sensitive)'))

  printOptionalArguments()
  exit()


def checkInputLength(args):
  if len(args) == 1:
    programUsageOutput()

  elif len(args) < 6:
    print("\nThis script takes at least 6 arguments.\n")
    print("Call this script with no arguments for usage details\n")
    print("Exiting...\n")



if __name__ == '__main__':

  checkInputLength(sys.argv)

  spectrumDirectory, acidMassFile, \
  pssmDirectory, decoyPeptideDirectory, \
  reverse = sys.argv[1:6]
  defaultParameters = parseParameterInput(sys.argv[6:])

  MMTString = str(defaultParameters['MMT'])
  massTolerance = '0.' + ('0' * (6 - len(MMTString))) + MMTString

  reverse = reverse.lower()
  if reverse == 'true' or reverse == 't':
    reverse = True
  elif reverse == 'false' or reverse == 'f':
    reverse = False
  else:
    print("unknown reverse parameter '{}'".format(reverse))

  decoyNames = os.listdir(decoyPeptideDirectory)
  decoyFiles = \
      [os.path.join(decoyPeptideDirectory, name) for name in decoyNames]
  tryFiles(decoyFiles)

  acidMassTable = \
      AcidMassTable.adjustForPrecision(
          AcidMassTableIO.getAminoMasses(acidMassFile),
          defaultParameters['PREC'])

  allPSSM = \
      PSSM.adjustForPrecision(
          PSSMIO.getAllPSSM(pssmDirectory,
                            defaultParameters['minP'],
                            defaultParameters['maxP']),
          defaultParameters['PREC'])

  conversionTable = \
    AcidConversion.createAcidConversionTable([acid for acid in acidMassTable])
  conversionOrder = []
  for conversion in conversionTable:
    conversionOrder.append(conversion)
  conversionOrder.sort(key=lambda x: len(x), reverse=True)

  if conversionTable == {}:
    print("Currently this program only accepts up to 23 amino acids total\n")
    print("The mass file you have provided "
          "has more than 23 amino acids defined\n")
    print("Other files may also have too "
          "many masses but haven't been checked\n"
          "by the program at this time\n")
    exit()

  # quick and dirty way to make sure all our amino masses match for now
  sanityCheck(allPSSM, acidMassTable)

  acidMassTable, allPSSM = \
      AcidConversion.convertAcidModifications(acidMassTable,
                                              allPSSM,
                                              conversionTable)

  H2OMassAdjusted = int(H2OMASS * (10**defaultParameters['PREC']))
  NH3MassAdjusted = int(NH3MASS * (10**defaultParameters['PREC']))
  protonMassAdjusted = int(PROTONMASS * (10**defaultParameters['PREC']))

  allDecoyPeptides = {}

  for decoyFile in decoyFiles:
    with open(decoyFile, 'r') as f:
      decoys = f.read()
    decoys = decoys.split('\n')

    decoys = fileToDict(decoys,
                        defaultParameters['PREC'],
                        defaultParameters['minP'],
                        defaultParameters['maxP'],
                        conversionTable,
                        conversionOrder,
                        acidMassTable,
                        reverse)

    for spectrumFileName in os.listdir(spectrumDirectory):

      if spectrumFileName not in allDecoyPeptides:
        allDecoyPeptides[spectrumFileName] = {}

      spectrumFile = os.path.join(spectrumDirectory, spectrumFileName)

      for spectrum in SpectrumIO.getSpectrums(spectrumFile):
        spectrumTitle = Spectrum.getTitle(spectrum)

        if spectrumTitle not in allDecoyPeptides[spectrumFileName]:
          allDecoyPeptides[spectrumFileName][spectrumTitle] = []

        spectrumMasses, spectrumMassesDouble, \
        spectrumIntensities, spectrumIntensitiesDouble = \
                                                         \
          SpectrumConversion.adjustSpectrumPrecision(
              *SpectrumConversion.processSpectrum(
                            Spectrum.getMasses(spectrum),
                            Spectrum.getIntensities(spectrum),
                            Spectrum.getPrecursorMass(spectrum),
                            H2OMASS,
                            PROTONMASS,
                            Spectrum.getCharge(spectrum),
                            defaultParameters['COMP']),
              defaultParameters['PREC'])

        peptides = getDecoys(spectrumMasses[-1], decoys, massTolerance)

        for peptide in peptides:

          experimentalSpectrum = \
              spectrumMasses + spectrumMassesDouble
          experimentalIntensities = \
              spectrumIntensities + spectrumIntensitiesDouble

          globalScore = \
              ProcessResults.calculateGlobalScore(acidMassTable,
                                                  experimentalSpectrum,
                                                  experimentalIntensities,
                                                  peptide,
                                                  protonMassAdjusted,
                                                  H2OMassAdjusted,
                                                  NH3MassAdjusted,
                                                  defaultParameters['MMT'])

          allDecoyPeptides[spectrumFileName][spectrumTitle].append((peptide,
                                                                   globalScore))


  for spectrumFileName in allDecoyPeptides:
    globalBest = []

    # Code for global best scores is here
    #
    # for spectrumTitle in allDecoyPeptides[spectrumFileName]:

    #   if allDecoyPeptides[spectrumFileName][spectrumTitle] == []:
    #     continue

    #     list.sort(allDecoyPeptides[spectrumFileName][spectrumTitle],
    #               key=lambda x: x[1],
    #               reverse=True)
    #   globalBest.append(allDecoyPeptides[spectrumFileName][spectrumTitle][0])

    # globalBestOut = open(spectrumFileName + ".decoys.globalBest", 'w')

    # for best in globalBest:
    #   globalBestOut.write(str(best[0]) + ',' + str(best[1]) + '\n')

    # globalBestOut.close()

    for pssmTitle in allPSSM:
      ignoredLengths = PSSM.getIgnoredLengths(pssmTitle, allPSSM)
      bestAminos = []
      bestCombined = []

      for spectrumTitle in allDecoyPeptides[spectrumFileName]:
        addedScores = []
        if allDecoyPeptides[spectrumFileName][spectrumTitle] == []:
          continue

        for peptide in allDecoyPeptides[spectrumFileName][spectrumTitle]:
          positionalScore = \
            getPSSMScore(peptide[0], pssmTitle, allPSSM) / \
            (10 ** defaultParameters['PREC'] * len(peptide[0]))

          addedScores.append((peptide[0],
                              positionalScore,
                              positionalScore * peptide[1]))

        addedScores.sort(key=lambda x: x[1], reverse=True)
        bestAminos.append(addedScores[0])
        addedScores.sort(key=lambda x: x[2], reverse=True)
        bestCombined.append(addedScores[0])

      # aminoBestFileName = \
      #     str.format('{}.{}.decoys.aminoBest',
      #                spectrumFileName,
      #                pssmTitle)
      combinedBestFileName = \
          str.format('{}.{}.decoys',
                     spectrumFileName,
                     pssmTitle)

      #aminoBestOut = open(aminoBestFileName, 'w')
      combinedBestOut = open(combinedBestFileName, 'w')

      # for best in bestAminos:
      #   aminoBestOut.write(str(best[0]) + ',' + str(best[1]) + '\n')
      for best in bestCombined:
        combinedBestOut.write(str(best[0]) + ',' + str(best[2]) + '\n')

      #aminoBestOut.close()
      combinedBestOut.close()