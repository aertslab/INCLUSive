#include "utilities.h"
#include "PWMIO.h"
#include "PWM.h"
#include "BackgroundIO.h"
#include "BackgroundModel.h"
#include "Instance.h"
#include "FastaIO.h"
#include "RandomNumber.h"
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

#define OPTIONS "hf:m:b:s:w:x:t:l:o:s:v"
// #define VERSION "3.0"


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
	
	
cerr << endl;
cerr << "At date, this application is not actively used neither maintained. "
<< "Contact <mclaeys@qatar.net.qa> if you have questions "
		<< "or interest to (re)use it. "
<< endl;
exit(0);
	
	
	
	
  bool bFasta = false,
    bBG = false,
    bStrand = false,
    bOut = false,
    bList = false,
    bMatrix = false,
    bAbsolute = false,
    bWindow = false,
    bStepSize = false,
    bMatching = false;
  strand_modes strand = plus_strand;
  char c;
  double threshold = 0,
    wScore = 0,
    oldScore = 0,
    newScore = 0;
  int wLength = 0,
    window = 0,
    stepSize = 1,
    wStart = 0,
    wEnd = 0,
    index = 0,
    nbr = 0;
  double **pWindowScore;
    
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
      if ( atoi(optarg) == 0 )
			{
				bStrand = false;
        strand = plus_strand;
			}
			else
			{
				bStrand = true;
	      strand = both;
			}
      break;
    case 'w':
      window = atoi(optarg);
      bWindow = true;
      break;
    case 'x':
      stepSize = atoi(optarg);
      bStepSize = true;
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
  {
    cerr << "CireScorer: Problems reading background file." << endl;
    cleanup();
    exit(-1);
  }
  BackgroundModel *pBgModel = bgIO->ReadBackgroundModel();
  if (pBgModel == NULL)
  {
    cerr << "CireScorer: Problems reading background file." << endl;
    cleanup();
    exit(-1);
  }
  delete bgIO;

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
    cerr << "CireScorer: Unable to open matrix model file: " << *matrixFile
      << endl;
    cleanup();
    exit(-1);
  }

  // some variables to store matrix data
  PWM *myMatrix;
  list < PWM * >matrixList;
  list < PWM * >::iterator matIter;
  double wx = 0;
  while ( pwmIO->IsOpen() )
  {
    myMatrix = pwmIO->ReadMatrix(); 
    if ( myMatrix )
    {
	    cerr << " " << *(myMatrix->GetID()) << " " << *(myMatrix->GetConsensus()) << endl;	
      bMatching = true;
      if ( bList )
      {
        bMatching = false;
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
      
      if (bMatching)
        matrixList.push_back(myMatrix);
    }
  }
	cerr << "done reading matrices." << endl;
  // close matrix reader
  delete pwmIO;
	
  if (matrixList.size() == 0)
  {
    cerr << "CireScorer: No usable matrix found." << endl;
    cleanup();
    delete myMatrix;
  }

  uint nbrMotifs = matrixList.size();
  cerr << "Found " << nbrMotifs << " matrices" << endl;
	
  // open fasta file
	cerr << "Open fasta file for reading." << endl;
  FastaIO *fileIO = new FastaIO(fastaFile);
  if (fileIO == NULL)
  {
    cerr << "CireScorer: Unable to read sequences" << endl;

    // cleanup the variables
    cleanup();
    if ( pBgModel != NULL )
      delete pBgModel;
    if ( myMatrix != NULL )
      delete myMatrix;
    matIter = matrixList.begin();
    while (matIter != matrixList.end())
    {
      cerr << "deleting matrix: " << *((*matIter)->GetID()) << endl;
      if ( *matIter != NULL )
        delete *matIter;
      matIter++;
    }
    exit(-1);
  }

  // open file for output writing if necessarry
  ofstream OFS;
  if (bOut)
  {
    OFS.open(outputFile->c_str(), ios::out);
  }

  // start processing the sequences
  int counter = 0,
    L = 0;
  SequenceObject *pSeqObj = NULL;
  SequenceComputation *pSeqComp = NULL;
  while ( fileIO->HasNext() )
  {
    pSeqObj = fileIO->NextSequence();
		cerr << "Read sequence: " << *(pSeqObj->GetID()) << endl;
    if (pSeqObj == NULL)
      continue;

    // sequence length
    L = pSeqObj->Length();

    // define first window size
    wStart = 0;
    wEnd = L - 1;
    if ( bWindow )
    {
      wEnd = wStart + window - 1;      
    }
    else
    {
      window = L;
    }
    
    if ( wEnd >= L )
      continue;
    
    // create computation object
		cerr << "Create computation object..." << endl;
    pSeqComp = new SequenceComputation(pSeqObj);

    // augment sequence counter
    counter++;
    wLength = 0;

    // window score
    pWindowScore = new double*[nbrMotifs];
    for (uint i=0; i<nbrMotifs; i++)
      pWindowScore[i] = new double[L-window+1];

    // scroll through matrix list
    matIter = matrixList.begin();
    nbr = 0;
    while (matIter != matrixList.end())
    {
			cerr << "scoring with matrix: " << *((*matIter)->GetID()) << endl;
      if (L <= wLength )
      {
        // motif length is larger than sequence length move to next motif
        cerr << "Warning: Sequence too short: " << *(pSeqObj->GetID()) << endl;
        matIter++;
        continue;
      }

      // initialize window
      wStart = 0;
      wEnd = wStart + window - 1;
      
      // define motif length in computation object
      // and score background model 
			// => only necessary when motif length changes
      if ( (*matIter)->Length() !=  pSeqComp->GetMotifLength() )
      {
        wLength = (*matIter)->Length();
				pSeqComp->SetMotifLength(wLength);
        pSeqComp->UpdateInstanceBackgroundScore(pBgModel, wLength, strand);
      }

      // score sequences with motif model
      pSeqComp->UpdateInstanceMotifScore(*matIter, strand);

      // cerr << "sequence scored." << endl; 

      // compute score of first window
      wScore = 0;
      for (int j = 0; j < window; j++)
      {
        wx = pSeqComp->GetMotifScoreAt(j, plus_strand) -  pSeqComp->GetBackgroundScoreAt(j, plus_strand);

        // select motifs with score higher than threshold
        if (wx >= threshold)
          wScore += wx;
      }
      if ( bStrand )
      {
        for (int j = 0; j < window; j++)
        {
          wx = pSeqComp->GetMotifScoreAt(L-j-1, minus_strand) -  pSeqComp->GetBackgroundScoreAt(L-j-1, minus_strand);

          // select motifs with score higher than threshold
          if (wx >= threshold)
            wScore += wx;
        }
      }
      pWindowScore[nbr][0] = wScore;

      cerr << wScore << " ";
      
      // proceed with sliding windows
      wStart += stepSize;
      wEnd += stepSize;
      index = 1;
      while ( wEnd < L ) 
      {
        // previous scores
        oldScore = 0;
        for ( int i=0; i<stepSize; i++)
        {
          wx = pSeqComp->GetMotifScoreAt(wStart - i - 1, plus_strand) -  pSeqComp->GetBackgroundScoreAt(wStart - i - 1, plus_strand);
          if ( wx > threshold )
            oldScore += wx;
          cerr << "|" << wx;

          if ( bStrand )
          {
            wx = pSeqComp->GetMotifScoreAt(L - wStart + i, minus_strand) -  pSeqComp->GetBackgroundScoreAt(L - wStart + i, minus_strand);
            if ( wx > threshold )
              oldScore += wx;              
          }          
        }
				cerr << " -" << oldScore;
				
        // new scores
        newScore = 0;
        for ( int i=0; i<stepSize; i++)
        {
          wx = pSeqComp->GetMotifScoreAt(wEnd - i, plus_strand) -  pSeqComp->GetBackgroundScoreAt(wEnd - i, plus_strand);
          if ( wx > threshold )
            newScore += wx;
          cerr << "|" << wx;
					
          if ( bStrand )
          {
            wx = pSeqComp->GetMotifScoreAt(L - wEnd + i,minus_strand) -  pSeqComp->GetBackgroundScoreAt( L - wEnd + i, minus_strand);
            if ( wx > threshold )
              newScore += wx;
          }
        }
        cerr << " +" << newScore << "= ";
				
        // update window score 
        pWindowScore[nbr][index] = pWindowScore[nbr][index - 1] + newScore - oldScore;   

        cerr << " " << pWindowScore[nbr][index] << " ";
        // cerr << nbr << " " << index << "\r";
        // augment counters
        index++;
        wStart += stepSize;
        wEnd += stepSize;
        
      }
      
      cerr << " till " << wEnd << " ?= " << L << endl;
      
      //move on to the next element
      matIter++;
      nbr++;
    }                           // end matrix while loop
    
    // print out results stored in pWindowScore
    wStart = 0;
    wEnd = window - 1;
    index = 0;
    while ( wEnd < L )
    {
      wScore = 0;
      for ( uint i=0; i<nbrMotifs; i++)
        wScore += pWindowScore[i][index];
      
      if ( bOut )
      {
        OFS << *(pSeqObj->GetID()) << "\t" 
          << wStart + 1<< "\t" 
          << wEnd + 1<< "\t" 
          << wScore << endl;
      }
      else
      {
        cout << *(pSeqObj->GetID()) << "\t" 
          << wStart + 1<< "\t" 
          << wEnd + 1<< "\t" 
          << wScore << endl;
      }
      
      wStart += stepSize;
      wEnd += stepSize;
      index++;
    }
    
    for (uint i=0; i<nbrMotifs; i++ )
      delete[] pWindowScore[i]; 
    delete pWindowScore;
    
    delete pSeqObj;
    pSeqObj = NULL;
    delete pSeqComp;
    pSeqComp = NULL;
  }                             // end next sequence loop

  // close fasta file
  delete fileIO;

  if ( bOut )
    OFS.close();
  
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
  cout << "Usage:" << " CireScorer <ARGS>" << endl;
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
    "  -o <outFile>        Output file to write results (default stdout)"
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
