#include "utilities.h"
#include "PWMIO.h"
#include "PWM.h"
#include "BackgroundIO.h"
#include "BackgroundModel.h"
#include "Instance.h"
#include "FastaIO.h"
#include "RandomNumber.h"
#include "GFFWriter.h"
#include "SequenceComputation.h"

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

#define OPTIONS "hf:m:b:t:l:o:s:va"
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
    bMatrix = false,
    bAbsolute = false;
  strand_modes strand = both;
  char c;
  double threshold = 0.85;
  double minX = 0;
  double maxX = 0;
  int wLength = 0,
    i = 0;
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
    case 't':
      threshold = (double) atof(optarg);
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
    case 'a':
      bAbsolute = true;
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
  if (!bFasta || !bBG || !bMatrix)
  {
    instructions();
    cleanup();
    exit(-1);
  }

  // load the background model
  BackgroundIO *bgIO = new BackgroundIO(*bgFile, READ);
  if (bgIO == NULL)
    exit(-1);
  BackgroundModel *pBgModel = bgIO->ReadBackgroundModel();
  delete bgIO;

  // extract pseudo counts
  double value = 0;
  double minSnf = 1;
  double maxSnf = 0;
  for (j = 0; j < 4; j++)
  {
    value = pBgModel->GetSnfValueAt(j);
    if (minSnf > value)
      minSnf = value;
    if (maxSnf < value)
      maxSnf = value;
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

  // list iterator
  list < string >::iterator matListIter = matrixIdList.begin();
  PWMIO *pwmIO = new PWMIO(matrixFile, READ);
  if (!pwmIO->IsOpen())
  {
    cerr << "MotifLocator: Unable to open matrix model file: " << *matrixFile
      << endl;
    cleanup();
    exit(-1);
  }

  // some variables to store matrix data
  PWM *myMatrix;
  list < PWM * >matrixList;
  list < double >minXList;
  list < double >maxXList;
  list < PWM * >::iterator matIter;
  list < double >::iterator minXListIter;
  list < double >::iterator maxXListIter;
  double minElem = 1;
  double maxElem = 0;
  double wx = 0;
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

      // compute minimal and maximal score of matrix
      minX = 0;
      maxX = 0;
      for (i = 0; i < myMatrix->Length(); i++)
      {
        minSnf = 1;
        minElem = 1;
        maxSnf = 0;
        maxElem = 0;
        for (j = 0; j < 4; j++)
        {
          value = myMatrix->GetValueAt(i, j);
          if (value < minElem){
            minElem = value;
            maxSnf = pBgModel->GetSnfValueAt(j);
          }
          if (maxElem < value){
            maxElem = value;
            minSnf = pBgModel->GetSnfValueAt(j);
          }
        }
        minX += log(minElem) - log(maxSnf);
        maxX += log(maxElem) - log(minSnf);
      }
      // minX -= myMatrix->Length() * log(maxSnf);
      // maxX -= myMatrix->Length() * log(minSnf);
      cerr << "+ " << *(myMatrix->
                        GetID()) << ": min Wx = " << minX << " |max Wx = "
        << maxX << endl;
      matrixList.push_back(myMatrix);
      minXList.push_back(minX);
      maxXList.push_back(maxX);

      // myMatrix->StderrPrintMatrix();
    }
  }

  // close matrix reader
  delete pwmIO;
  if (matrixList.size() == 0)
  {
    cerr << "MotifLocator: No usable matrix found." << endl;
    cleanup();
    delete myMatrix;
  }

  // open fasta file
  FastaIO *fileIO = new FastaIO(fastaFile);
  if (fileIO == NULL)
  {
    cerr << "MotifLocator: Unable to read sequences" << endl;

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

  // open file for output writing if necessarry
  GFFWriter *pGFFStream;
  if (bOut)
  {
    pGFFStream = new GFFWriter(outputFile);
  }

  else

  {
    pGFFStream = new GFFWriter();
  }

  // start processing the sequences
  int counter = 0,
    L = 0,
    nbr = 0;
  SequenceObject *pSeqObj = NULL;
  SequenceComputation *pSeqComp = NULL;
  string *Source = new string("MotifLocator");
  while (fileIO->HasNext())
  {
    pSeqObj = fileIO->NextSequence();
    if (pSeqObj == NULL)
    {
      continue;
    }

    // sequence length
    L = pSeqObj->Length();

    // create computation object
    pSeqComp = new SequenceComputation(pSeqObj);

    // augment sequence counter
    counter++;
    wLength = 0;
    Instance *pSite = NULL;

    // scroll through matrix list
    matIter = matrixList.begin();
    maxXListIter = maxXList.begin();
    minXListIter = minXList.begin();
    while (matIter != matrixList.end())
    {

      // (*matIter)->StderrPrintMatrix();
      nbr = 0;
      minX = *minXListIter;
      maxX = *maxXListIter;
      if (L <= 2 * (wLength))
      {

        // motif length is larger than sequence length move to next motif
        cerr << "Warning: Sequence too short: " << *(pSeqObj->
                                                     GetID()) << endl;
        matIter++;
        continue;
      }

      // define motif length in computation object
      pSeqComp->SetMotifLength(wLength);

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

      // select best scoring instances on plus strand
      for (j = 0; j < L - wLength + 1; j++)
      {
        wx =
          ((log(pSeqComp->GetWxAt(j, plus_strand))) - minX) / (maxX - minX);

        // and select motifs with score higher than threshold
        if (wx >= threshold)
        {
          nbr++;
          pSite = new Instance(pSeqObj, plus_strand, j, wLength);

          // go to next instance
          if (pSite == NULL)
            continue;

          // set score 
          if (bAbsolute)
          {
            pSite->SetScore(pSeqComp->GetWxAt(j, plus_strand));
          }
          else
          {
            pSite->SetScore(wx);
          }

          // set instance ID
          pSite->SetID((*matIter)->GetID());

          // output site
          pGFFStream->WriteInstance(pSite, Source);
          delete pSite;
          pSite = NULL;
        }
      }

      // do the reverse complement
      if (bStrand)
      {
        for (j = 0; j < L - wLength + 1; j++)
        {

          // re-normalize wx
          wx =
            (log(pSeqComp->GetWxAt(j,minus_strand)) - minX) / (maxX - minX);

          // and select motifs with score higher than threshold
          if (wx >= threshold)
          {
            nbr++;

            // create a new instance
            pSite = new Instance(pSeqObj, minus_strand, j, wLength);
            if (pSite == NULL)
              continue;

            // set score
            if (bAbsolute)
            {
              pSite->SetScore(pSeqComp->GetWxAt(j, minus_strand));
            }
            else
            {
              pSite->SetScore(wx);
            }

            // set instance ID
            pSite->SetID((*matIter)->GetID());

            // output site
            pGFFStream->WriteInstance(pSite, Source);
            delete pSite;
            pSite = NULL;
          }
        }
      }
      cerr << counter << " |" << *(pSeqObj->
                                   GetID()) << " + " << *((*matIter)->
                                                          GetID()) <<
        "   -> instances:  " << nbr << endl;

      //move on to the next element
      matIter++;
      maxXListIter++;
      minXListIter++;
    }                           // end matrix while loop
    delete pSeqObj;
    pSeqObj = NULL;
    delete pSeqComp;
    pSeqComp = NULL;
  }                             // end next sequence loop

  // close fasta file
  delete fileIO;

  // delete source name
  delete Source;

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
  exit(1);
}



/*********************
 *  LOCAL FUNCTIONS  *
 *********************/
void
instructions()
{
  cout << endl;
  cout << "Usage:" << " MotifLocator <ARGS>" << endl;
  cout << endl;
  cout << " Required Arguments" << endl;
  cout << "  -f <fastaFile>      Sequences in FASTA format" << endl;
  cout <<
    "  -b <bgFile>         File containing the background model description"
    << endl;
  cout <<
    "                      -> Use WriteBackgroundModel to generate such a file."
    << endl;
  cout <<
    "  -m <matrixFile>     File containing the matrix model descriptions" <<
    endl;
  cout <<
    "                      -> Use WriteMatrixModel to generate such a file."
    << endl;
  cout << endl;
  cout << " Optional Arguments" << endl;
  cout <<
    "  -t <value>          Sets threshold above which a motif is selected (default 0.85)."
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
  cout <<
    "  -a                  If selected, the scores are the normalized values "
    << endl;
  cerr << "                      instead of the absolute values." << endl;
  cout << endl;
  cout << "  -v                  Version of MotifLocator" << endl;
  cout << endl;
  cout << "Version " << VERSION << endl;
  cout << "Questions and Remarks: <gert.thijs@esat.kuleuven.ac.be>" << endl
    << endl;
} void

version()
{
  cout << endl;
  cout << "INCLUSive -- MotifLocator " << endl;
  cout << "  version " << VERSION << endl;
} void

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
