#include "inclusive.h"
#include "utilities.h"
#include "PWMIO.h"
#include "PWM.h"
#include "BackgroundIO.h"
#include "BackgroundModel.h"
#include "Instance.h"
#include "FastaIO.h"
#include "RandomNumber.h"
#include "SequenceComputation.h"
#include "GFFWriter.h"

// c++ includes
#include <list>
#include <iostream>
#include <fstream>

// c includes
#include <getopt.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstring>

#define OPTIONS "hf:m:b:p:l:o:s:v"
#define VERSION "3.0"

using namespace INCLUSIVE;
void instructions();
void version();
void cleanup();
string *fastaFile = NULL,
  *bgFile = NULL,
  *outputFile = NULL,
  *matrixFile = NULL,
  *listFile = NULL;


int
main(int argc, char *argv[])
{
  bool bFasta = false,
    bBG = false,
    bStrand = true,
    bOut = false,
    bList = false,
    bMatrix = false;
  strand_modes strand = both;
  char c;
  double priorValue = 0.5;
  int wLength = 0,
    i = 0,
    ndx = 0;
  int j = 0;
  if (optopt == 0)
  {
    instructions();
    exit(-1);
  }
  while ((c = getopt(argc, argv, OPTIONS)) != -1)
  {
    switch (c)
    {
    case 'f':

      // set fasta file name
      fastaFile = new string(optarg);
      bFasta = true;
      break;
    case 'm':
      matrixFile = new string(optarg);
      bMatrix = true;
      break;
    case 'b':
      bgFile = new string(optarg);
      bBG = true;
      break;
    case 'p':
      priorValue = (double) atof(optarg);
      break;
    case 'l':
      bList = true;
      listFile = new string(optarg);
      break;
    case 'o':

      // define new output stream 
      bOut = true;
      outputFile = new string(optarg);
      break;
    case 's':
      bStrand = atoi(optarg);
      strand = both;            // default both strands
      if (bStrand == 0)
      {
        strand = plus_strand;
      }
      break;
    case 'v':
      version();
      exit(-1);
    case '?':
    case 'h':
      instructions();
      exit(-1);
    default:
      cerr << "MotifSampler: Error in getopt() function" << endl;
      exit(-1);
    }
  }

  // required elements  
  if (!bFasta || !bBG || !bMatrix)
  {
    instructions();
    cleanup();
    exit(-1);
  }

  // open stream for GFF output
  GFFWriter *pGFFStream;
  if (bOut == true)
  {
    pGFFStream = new GFFWriter(outputFile);
  }
  else
  {
    pGFFStream = new GFFWriter();
  }

  // read list of matrix ID's
  list < string > matrixIdList;
  if (bList)
  {
    ifstream _ifs(listFile->c_str(), ios::in);
    string inLine;
    while (!_ifs.eof())
    {
      getline(_ifs, inLine, '\n');
      matrixIdList.push_back(inLine);
    }
    _ifs.close();
  }
  list < string >::iterator matListIter = matrixIdList.begin();
  i = 1;
  while (matListIter != matrixIdList.end())
  {
    cerr << i++ << ": id = " << (*matListIter) << endl;
    matListIter++;
  }
  PWMIO *pwmIO = new PWMIO(matrixFile, READ);
  if (!pwmIO->IsOpen())
  {
    cerr << "MotifScanner: Unable to open matrix model file: " << *matrixFile
      << endl;
    cleanup();
    exit(-1);
  }

  // some variables to store matrix data
  PWM *myMatrix;
  list < PWM * >matrixList;
  list < PWM * >::iterator matIter;
  while (pwmIO->IsOpen() && (myMatrix = pwmIO->ReadMatrix()))
  {
    bool bMatching = false;
    if (bList)
    {
      matListIter = matrixIdList.begin();
      while (matListIter != matrixIdList.end())
      {
        if ((*matListIter) == *(myMatrix->GetID()))
        {
          bMatching = true;
          break;
        }
        matListIter++;
      }
    }
    else
    {
      bMatching = true;
    }
    if (bMatching)
    {
      matrixList.push_back(myMatrix);
    }
  }

  // close matrix reader
  delete pwmIO;
  if (matrixList.size() == 0)
  {
    cerr << "MotifScanner: No usable matrix found." << endl;
    cleanup();
    delete myMatrix;
  }

  // load the background model
  BackgroundIO *bgIO = new BackgroundIO(*bgFile, READ);
  BackgroundModel *pBgModel = bgIO->ReadBackgroundModel();
  delete bgIO;

  // open fasta file 
  FastaIO *fileIO = new FastaIO(fastaFile);
  if (fileIO == NULL)
  {
    cerr << "MotifScanner: Unable to read sequences" << endl;

    // cleanup the variables
    cleanup();
    delete pBgModel;
    delete myMatrix;
    matIter = matrixList.begin();
    while (matIter != matrixList.end())
    {
      cerr << "deleting matrix: " << *((*matIter)->GetID()) << endl;
      delete(*matIter);
      matIter++;
    }
    exit(-1);
  }

  // define some variables
  int counter = 0,
    L = 0,
    nbr = 0;
  SequenceObject *pSeqObj;
  SequenceComputation *pSeqComp;

  // vector to contain alingment vector
  vector < int >aPos;

  // set source in GFF output
  string *pSource = new string("MotifScanner");
  while (fileIO->HasNext())
  {

    // read next sequence from fasta file
    pSeqObj = fileIO->NextSequence();
    if (pSeqObj == NULL)
    {
      continue;
    }

    // sequence length
    L = pSeqObj->Length();

    // create the computation object
    pSeqComp = new SequenceComputation(pSeqObj);

    // update full sequence background score (P0)
    pSeqComp->UpdateSequenceBackgroundScore(pBgModel, strand);

    // augment sequence counter
    counter++;
    wLength = 0;

    // scroll through matrix list
    matIter = matrixList.begin();
    while (matIter != matrixList.end())
    {
      if (L <= 2 * ((*matIter)->Length()))
      {

        // motif length is larger than sequence length move to next motif
        cerr << "Warning: Sequence too short: " << *(pSeqObj->
                                                     GetID()) << endl;
        matIter++;
        continue;
      }

      // set motif length in computation oject
      pSeqComp->SetMotifLength((*matIter)->Length());

      // score background model => only necessary when motif length changes
      if ((*matIter)->Length() != wLength)
      {
        wLength = (*matIter)->Length();
        pSeqComp->UpdateInstanceBackgroundScore(pBgModel, wLength, strand);
      }

      // score sequences with motif model
      pSeqComp->UpdateInstanceMotifScore(*matIter, strand);

      // compute Wx
      pSeqComp->UpdateInstanceExpWx(strand);

      // compute the distribution of the number of instances
      pSeqComp->UpdateCopyProbability(priorValue, strand);

      // estimate number of copies
      if (strand == plus_strand || strand == both)
      {
        nbr = pSeqComp->GetEstimatedNumberInstances(plus_strand);
        cerr << counter << "| " << *(pSeqObj->
                                     GetID()) << " + " << *((*matIter)->
                                                            GetID()) <<
          "   -> instances:  " << nbr << endl;
        if (nbr)
        {

          // set number of instances in alignment vector
          if (nbr != (int) aPos.size())
            aPos.resize(nbr);

          // select the nbr best scoring instances 
          pSeqComp->SelectBestInstanceStart(aPos, nbr, plus_strand);
          for (j = 0; j < nbr; j++)
          {

            // select position
            ndx = aPos[j];

            // create new instance
            Instance *pInstance =
              new Instance(pSeqObj, plus_strand, ndx, wLength);
            if (pInstance != NULL)
            {
              pInstance->SetScore(pSeqComp->GetWxAt(ndx, plus_strand));

              // set instance ID
              pInstance->SetID((*matIter)->GetID());

              // write instance
              pGFFStream->WriteInstance(pInstance, pSource);
              delete pInstance;
            }
            pInstance = NULL;
          }
        }
      }
      if (strand == both || strand == minus_strand)
      {

        // number of instances found on minus strand
        nbr = pSeqComp->GetEstimatedNumberInstances(minus_strand);
        cerr << counter << "| " << *(pSeqObj->
                                     GetID()) << " - " << *((*matIter)->
                                                            GetID()) <<
          "   -> instances:  " << nbr << endl;
        if (nbr)
        {

          // update alignment vector
          if (nbr != (int) aPos.size())
            aPos.resize(nbr);
          pSeqComp->SelectBestInstanceStart(aPos, nbr, minus_strand);
          for (j = 0; j < nbr; j++)
          {

            // select position
            ndx = aPos[j];

            // create new instance
            Instance *pInstance =
              new Instance(pSeqObj, minus_strand, ndx, wLength);
            if (pInstance != NULL)
            {
              pInstance->SetScore(pSeqComp->GetWxAt(ndx, minus_strand));
              pInstance->SetID((*matIter)->GetID());
              pGFFStream->WriteInstance(pInstance, pSource);
              delete pInstance;
            }
            pInstance = NULL;
          }
        }
      }
      matIter++;
    }
    delete pSeqObj;
    pSeqObj = NULL;
    delete pSeqComp;
    pSeqComp = NULL;
  }

  // close fasta file
  delete fileIO;

  // delete source name
  delete pSource;

  // close outpfile
  delete pGFFStream;

  // clear matrixList
  matIter = matrixList.begin();
  while (matIter != matrixList.end())
  {

    // cerr << "deleting matrix: " << *((*matIter)->GetID()) << endl;
    delete(*matIter);
    matIter++;
  }
  delete pBgModel;
  delete myMatrix;
  cleanup();
  
  // program ended succesful
  exit(0);
}



/*********************
 *  LOCAL FUNCTIONS  *
 *********************/
void
instructions()
{
  cout << endl;
  cout << "Usage:" << " MotifScanner <ARGS>" << endl;
  cout << endl;
  cout << " Required Arguments" << endl;
  cout << "  -f <fastaFile>      Sequences in FASTA format" << endl;
  cout <<
    "  -b <bgFile>         File containing the background model description"
    << endl;
  cout <<
    "  -m <matrixFile>     File containing the matrix model descriptions" <<
    endl;
  cout << endl << endl;
  cout << " Optional Arguments" << endl;
  cout <<
    "  -p <value>          Sets prior probability of 1 motif copy. (default 0.2) "
    << endl;
  cout <<
    "  -o <outFile>        Output file to write results in GFF (default stdout)"
    << endl;
  cout <<
    "  -l <listFile>       File with a list of identifiers to select individual"
    << endl;
  cout <<
    "                      matrices from the matrix file. IDs that are not found are omitted."
    << endl;
  cout << "  -s <0|1>            Select strand. (default both)" << endl;
  cout <<
    "                      0 is only input sequences, 1 include reverse complement."
    << endl;
  cout << endl;
  cout << "  -v                  Version of MotifScanner" << endl;
  cout << endl;
  cout << "Version " << VERSION << endl;
  cout << "Questions and Remarks: <gert.thijs@esat.kuleuven.ac.be>" << endl
    << endl;
} 

void
version()
{
  cout << endl;
  cout << "INCLUSive -- MotifScanner (stand alone C++ version)" << endl;
  cout << "  Version " << VERSION << endl;
} 


void
cleanup()
{
  if (fastaFile != NULL)
    delete fastaFile;
  fastaFile = 0;
  if (bgFile != NULL)
    delete bgFile;
  bgFile = NULL;
  if (matrixFile != NULL)
    delete matrixFile;
  matrixFile = NULL;
  if (outputFile != NULL)
    delete outputFile;
  outputFile = NULL;
  if (listFile != NULL)
    delete listFile;
  listFile = NULL;
  return;
}
