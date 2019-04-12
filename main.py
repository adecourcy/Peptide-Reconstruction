#!/usr/bin/env python3
import sys, os

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

###############################################################################
#
# Functions for the handling of temporary files
#
###############################################################################

def writeToFile(acidTable,
                pssm,
                spectrumMasses,
                spectrumMassesDouble,
                spectrumScores,
                spectrumScoresDouble,
                aminoAcids,
                specFileName,
                specDoubleFileName,
                acidFileName,
                massFileName,
                scoreFileName):

  specFile = open(specFileName, "w")
  specDoubleFile = open(specDoubleFileName, "w")
  acidFile = open(acidFileName, "w")
  massFile = open(massFileName, "w")
  scoreFile = open(scoreFileName, "w")

  writeLine = ""
  specFile.write(str(len(spectrumMasses)) + "\n")
  for mass, score in zip(spectrumMasses, spectrumScores):
    writeLine += str(mass) + " " + str(score) + "\n"
  specFile.write(writeLine.strip())

  writeLine = ""
  specDoubleFile.write(str(len(spectrumMassesDouble)) + "\n")
  for mass, score in zip(spectrumMassesDouble, spectrumScoresDouble):
    writeLine += str(mass) + " " + str(score) + "\n"
  specDoubleFile.write(writeLine.strip())

  writeLine = ""
  acidFile.write(str(len(acidTable)) + "\n")
  for acid in aminoAcids:
    if acid in acidTable:
      writeLine += (acid + "\n")
  acidFile.write(writeLine.strip())

  writeLine = ""
  for acid in aminoAcids:
    if acid in acidTable:
      writeLine += str(acidTable[acid]) + "\n"
  massFile.write(writeLine.strip())

  writeLine = ""
  for length in pssm:
    for acid in aminoAcids:
      if acid in pssm[length]:
        for prob in pssm[length][acid]:
          writeLine += str(prob) + " "
        writeLine = writeLine.strip() + "\n"
  scoreFile.write(writeLine.strip())

  specFile.close()
  specDoubleFile.close()
  acidFile.close()
  massFile.close()
  scoreFile.close()


def clean():
  os.remove(specFile)
  os.remove(specDoubleFile)
  os.remove(acidFile)
  os.remove(massFile)
  os.remove(scoreFile)
  os.remove(resultsFile)

###############################################################################
#
# Functions for handling user input
#
###############################################################################

def programUsageOutput():
  print("\nUsage: main.py SPECDIR ACIDMASSFILE PSSMDIR [OPT ARGS]\n\n")

  padding = '{:<20}{}'
  print("Positional arguments:")
  print(str.format(padding, 'SPECDIR',
        'A directory containing 1 or more spectrum files'))
  print(str.format(padding, 'ACIDMASSFILE',
        'A file containing mass spectrometry data'))
  print(str.format(padding, 'PSSMDIR',
        'A directory containing a positional scoring matrix'))

  printOptionalArguments()
  exit()


def improperNumArgumentsOutput():
  print("\nThis script takes a minimum of 3 arguments,")
  print("with up to 13 additional, optional arguments")
  print("Call this script with no arguments for usage details")
  print("Exiting...\n")
  exit()


def checkInputLength(args):
  if len(args) == 1:
    programUsageOutput()

  if len(args) < 4 or len(args) > 17:
    improperNumArgumentsOutput()

###############################################################################

if __name__ == '__main__':

  try:
    clean()
  except:
    pass

  # Take our user input and do some pre-processing to set up our data
  checkInputLength(sys.argv)
  spectrumDirectory, acidMassFile, pssmDirectory = sys.argv[1:4]
  defaultParameters = parseParameterInput(sys.argv[4:])

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

  if conversionTable == {}:
    print("Currently this program only accepts up to 23 amino acids total\n")
    print("The mass file you have provided has more than 23 amino acids defined\n")
    print("Other files may also have too many masses but haven't been checked\n"
          "by the program at this time\n")
    clean()
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

  # Start processing our spectrum files
  for spectrumFileName in os.listdir(spectrumDirectory):
    
    # Set up our output file
    spectrumFile = os.path.join(spectrumDirectory, spectrumFileName)
    # Make sure we can actually write to our output
    try:
      outputFile = open(spectrumFileName + defaultParameters['OUT'], 'w')
    except:
      print(str.format("Error opening {} to write output, skpping {} spectrum file",
                        spectrumFileName + defaultParameters['OUT'],
                        spectrumFileName))
      continue
    outputFile.write("Filename,Spectrum Title,PSSM Name,Peptide,Global Score,"
                     "Positional Score,Combined Score,Precursor Mass,Precursor Error\n")

    for spectrum in SpectrumIO.getSpectrums(spectrumFile):
      # Obtain and convert our spectrum information to a usable format
      spectrumMasses, spectrumMassesDouble, \
      spectrumIntensities, spectrumIntensitiesDouble = \
          SpectrumConversion.adjustSpectrumPrecision(
              *SpectrumConversion.processSpectrum(Spectrum.getMasses(spectrum),
                                                  Spectrum.getIntensities(spectrum),
                                                  Spectrum.getPrecursorMass(spectrum),
                                                  H2OMASS,
                                                  PROTONMASS,
                                                  Spectrum.getCharge(spectrum),
                                                  defaultParameters['COMP']),
              defaultParameters['PREC'])

      if defaultParameters['END'] == 1:
        spectrumMasses[-1] -= H2OMassAdjusted
        spectrumMassesDouble[-1] -= H2OMassAdjusted


      for pssmTitle in allPSSM:

        cPath = os.path.join("./scoring")
        specFile = os.path.join(cPath, "specFile")
        specDoubleFile = os.path.join(cPath, "specDoubleFile")
        acidFile = os.path.join(cPath, "acidFile")
        massFile = os.path.join(cPath, "massFile")
        scoreFile = os.path.join(cPath, "scoreFile")
        resultsFile = os.path.join(cPath, "results")
        aminoAcids = ["G", "A", "S", "P", "V",
                      "T", "L", "I", "N", "D",
                      "Q", "K", "E", "M", "H",
                      "F", "R", "C", "Y", "W",
                      "B", "J", "O"]

        for acid in acidMassTable:
          if acid not in aminoAcids:
            aminoAcids.append(acid)

        writeToFile(acidMassTable,
                    allPSSM[pssmTitle][1],
                    spectrumMasses,
                    spectrumMassesDouble,
                    spectrumIntensities,
                    spectrumIntensitiesDouble,
                    aminoAcids,
                    specFile,
                    specDoubleFile,
                    acidFile,
                    massFile,
                    scoreFile)

        # This is an awful way to handle this. Figure out issues with Boost
        cOutput = os.system("./scoring/findPep {} {} {} {} {} \
                                     {} {} {} {} {} \
                                     {} {} {} {} {} \
                                     {} {} {} {} {} \
                                     {} {} {}".format(specFile,
                                                      specDoubleFile,
                                                      acidFile,
                                                      massFile,
                                                      scoreFile,
                                                      resultsFile,
                                                      defaultParameters["IMC"],
                                                      defaultParameters["TMC"],
                                                      Spectrum.getCharge(spectrum) - 1,
                                                      defaultParameters["PREC"],
                                                      float(defaultParameters["BPEN"]),
                                                      defaultParameters["MMT"],
                                                      defaultParameters["AMT"],
                                                      defaultParameters["MMT"],
                                                      defaultParameters["TOLP"],
                                                      spectrumMasses[-1],
                                                      spectrumMasses[0],
                                                      defaultParameters["minP"],
                                                      defaultParameters["maxP"],
                                                      H2OMassAdjusted,
                                                      protonMassAdjusted,
                                                      defaultParameters["END"],
                                                      defaultParameters["BIN"]))

        if cOutput != 0:
          clean()
          exit()

        # Do some post-processing and write our data
        spectrumMassesDouble = [mass - protonMassAdjusted for mass in spectrumMassesDouble]
        experimentalSpectrum = spectrumMasses + spectrumMassesDouble
        experimentalIntensities = spectrumIntensities + spectrumIntensitiesDouble
        results = ProcessResults.processResults(resultsFile,
                                                acidMassTable,
                                                experimentalSpectrum,
                                                experimentalIntensities,
                                                protonMassAdjusted,
                                                H2OMassAdjusted,
                                                int((NH3MASS * (10 ** defaultParameters["PREC"]))),
                                                defaultParameters["MMT"],
                                                defaultParameters["PREC"],
                                                conversionTable,
                                                defaultParameters['END'],
                                                defaultParameters['BPEN'])


        spectrumTitle = Spectrum.getTitle(spectrum)
        spectrumTitle = spectrumTitle.replace(",", '')
        if len(results) == 0:
          outputString = spectrumFileName + "," + spectrumTitle + \
                         ',' + pssmTitle + ",No Peptide Found" + "\n"
          outputFile.write(outputString)

        else:
          ignoredLengths = PSSM.getIgnoredLengths(pssmTitle, allPSSM)

          for node in results:
            if len(node.peptideString) in ignoredLengths:
              continue

            outputString = \
              str.format('{},{},{},{},{},{},{},{},{}',
                         spectrumFileName,
                         spectrumTitle,
                         pssmTitle,
                         node.peptideString,
                         node.globalScore,
                         node.aminoScore,
                         node.combinedScore,
                         Spectrum.getPrecursorMass(spectrum),
                         node.precursorError)

            outputFile.write(outputString)

        clean()

  outputFile.close()