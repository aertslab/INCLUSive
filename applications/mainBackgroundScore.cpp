#include "inclusive.h"
#include "BackgroundIO.h"
#include "BackgroundModel.h"
#include "FastaIO.h"
#include "utilities.h"

// c++ includes
#include <iostream>
#include <fstream>

// c includes
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstring>

#define OPTIONS "hf:b:o:s:w:t:v"

// forward declaration function
void instructions();
void version();
void _cleanup();

string *fastaFile = NULL,
  *bgFile = NULL,
  *outputFile = NULL;


int
main(int argc, char *argv[])
{
cerr << endl;
cerr << "At date, this application is not actively used neither maintained. "
<< "Contact <mclaeys@qatar.net.qa> if you have questions "
		<< "or interest to (re)use it. "
<< endl;
exit(0);
	
	
	
  bool bFasta = false,
    bBG = false,
    bOut = false,
    bStepSize = false,
    bWindow = false;
  char c;
	strand_modes bStrand = both;
  int window = 0,
    stepSize = 1;
  
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
    case 'b':
      bgFile = new string(optarg);
      bBG = true;
      break;
    case 'o':
      // define new output stream 
      bOut = true;
      outputFile = new string(optarg);
      break;
    case 's':
      bStrand = both;
      if ( atoi(optarg) == 0)
        bStrand = plus_strand;
      break;
    case 'w':
      window = atoi(optarg);
      bWindow = true;
      break;
    case 't':
      stepSize = atoi(optarg);
      bStepSize = true;
      break;
    case 'v':
      version();
      exit(1);
    case '?':
    case 'h':
      instructions();
      exit(1);
    default:
      cerr << "BackgroundScore: Error in getopt() function" << endl;
      exit(-1);
    }
  }


  if (!bFasta || !bBG)
  {
    instructions();
    _cleanup();
    exit(-1);
  }

  // load the background model
  BackgroundIO *bgIO = new BackgroundIO(*bgFile, READ);
  if (bgIO == NULL)
  {
    cerr << "BackgroundScore: Problems reading background file." << endl;
    _cleanup();
    exit(-1);
  }
  BackgroundModel *pBgModel = bgIO->ReadBackgroundModel();
  delete bgIO;
  if (pBgModel == NULL)
  {
    cerr << "BackgroundScore: Problems reading background file." << endl;
    _cleanup();
    exit(-1);
  }

  // open fasta file 
  FastaIO *fileIO = new FastaIO(fastaFile);
  if (fileIO == NULL)

  {
    cerr << "BackgroundScore: Unable to read sequences" << endl;
    _cleanup();
    if (pBgModel != NULL)
      delete pBgModel;
    pBgModel = NULL;
    exit(-1);
  }

  // open file for output writing if necessarry
  ofstream OFS;
  if (bOut)
  {
    OFS.open(outputFile->c_str(), ios::out);
  }

  // sequence variables
  int L = 0,
    wStart = 0,
    wEnd = 0;
  double p0 = 0;
  SequenceObject *pSeqObj;

  while (fileIO->HasNext())
  {
    pSeqObj = fileIO->NextSequence();
    if (pSeqObj == NULL)
    {
      continue;
    }

    // get sequence length
    L = pSeqObj->Length();
    wStart = 0;
    wEnd = wStart + L - 1;
    
    if ( bWindow )
      wEnd = wStart + window - 1;      
    
    if ( wEnd >= L )
      wEnd = L - 1;
    
    while ( wEnd < L )
    {

      if ( !bWindow )
      {
        p0 = INCLUSIVE::SequenceLogBackgroundScore(pSeqObj, plus_strand, pBgModel);
      }else{
        p0 = INCLUSIVE::SequenceLogBackgroundScore(pSeqObj, wStart, window, plus_strand, pBgModel);
      }
      
      if (bOut == false)
      {
        cout << *(pSeqObj->GetID()) << "\t" 
          << wStart + 1<< "\t" 
          << wEnd + 1<< "\t" 
          << "+" << "\t" 
          << p0 << endl;
      }
      else
      {
        OFS << *(pSeqObj->GetID()) << "\t" 
          << wStart + 1 << "\t" 
          << wEnd + 1 << "\t" 
          << "+" << "\t" 
          << p0 << endl;
      }

      // do the reverse complement
      if ( bStrand == both)
      {
        if ( !bWindow )
        {
          p0 = INCLUSIVE::SequenceLogBackgroundScore(pSeqObj, minus_strand, pBgModel);
        }
        else
        {
          p0 = INCLUSIVE::SequenceLogBackgroundScore(pSeqObj, wStart, window, minus_strand, pBgModel);
        }
          
        if (bOut == false)
        {
          cout << *(pSeqObj->GetID()) << "\t" 
            << wStart + 1<< "\t" 
            << wEnd + 1 << "\t" 
            << "-" << "\t" 
            << p0 << endl;
        }
        else
        {
          OFS << *(pSeqObj->GetID()) << "\t" 
            << wStart + 1 << "\t" 
            << wEnd + 1 << "\t" 
            << "-" << "\t" 
            << p0 << endl;
        }
      }

      // augment windows
      wStart += stepSize;
      wEnd += stepSize;
      
    }
    
    // cleanup variables
    delete pSeqObj;
    pSeqObj = NULL;
  }

  // close fasta file
  delete fileIO;

  // close outpfile
  if (bOut)
    OFS.close();
  _cleanup();  
  if (pBgModel != NULL)
    delete pBgModel;
  pBgModel = NULL;
  
  exit(0);

}



/*********************
 *  LOCAL FUNCTIONS  *
 *********************/
void
instructions()
{
  cout << endl;
  cout << "Usage:" << " BackgroundScore <ARGS>" << endl;
  cout << endl;
  cout << " Required Arguments" << endl;
  cout << "  -f <fastaFile>      Sequences in FASTA format" << endl;
  cout <<
    "  -b <bgFile>         File containing the background model description"
    << endl;
  cout << endl;
  cout << " Optional Arguments" << endl;
  cout <<
    "  -o <outFile>        Output file to write results in GFF (default stdout)"
    << endl;
  cout << "  -s <0|1>            Select strand. (default both)" << endl;
  cout <<
    "                      0 is only input sequences, 1 include reverse complement."
    << endl;
  cout << "  -v                  Version of BackgroundScore" << endl;
  cout << endl;
  cout << "Version " << VERSION << endl;
  cout << "Questions and Remarks: <gert.thijs@esat.kuleuven.ac.be>" << endl
    << endl;
} 

void
version()
{
  cout << endl;
  cout << "INCLUSive -- BackgroundScore (stand alone C++ version)" << endl;
  cout << "   Version " << VERSION << endl;
}

void 
_cleanup()
{
  // cleanup variables
  if (fastaFile != NULL)
    delete fastaFile;
  fastaFile = NULL;
  if (bgFile != NULL)
    delete bgFile;
  bgFile = NULL;
}
