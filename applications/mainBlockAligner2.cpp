#include "inclusive.h"
#include "utilities.h"
#include "PWMIO.h"
#include "PWM.h"
#include "RandomNumber.h"
#include "BlockAlignment.h"

#include <iostream>
#include <fstream>
#include <getopt.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cstring>

#define OPTIONS "d:m:o:M:t:g:w:s:k:r:vh"
#define PI 3.14159265358979

// function prototypes
void instructions();
void version();
void cleanup();

string * outputFile = NULL,
  * dbFile = NULL,
  * matrixFile = NULL,
  * matOutFile = NULL,
  * shuffleFile = NULL;

int
main(int argc, char *argv[])
{
	cerr << endl;
cerr << "At date, this application is not actively used neither maintained. "
<< "Contact <mclaeys@qatar.net.qa> if you have questions "
		<< "or interest to (re)use it. "
<< endl;
exit(0);
	
	

  // initialize random number generator
  RandomNumber::InitRandomNumber();
  srandom((long)time(NULL));
  cerr << "Seed = " << random() << endl;
  
  char c;
  double maxScore = 0,
    maxScoreRev = 0,
    threshold = 0.3,
    gapScore = 0.2;
  int xStart = 0,
    yStart = 0,
    w1 = 0,
    w2 = 0,
    wCommon = 0,
    minLength = 4,
    nbrShuffles = 0,
    nbrSelection = 500;
  bool bDB = false,
    bMatrix = false,
    bOut = false,
    bMatOut = false,
    bShuffle = false,
		bReverse = true,
		bRevMax = false;
  PWMIO* matIO = NULL;
  PWM * pwm1 = NULL,
    * pwm2 = NULL,
    * myMatrix = NULL,
    * pPwmRev = NULL;
  char cid[1024];

  list<PWM *>::iterator matIter;
  
  if (optopt == 0)
  {
    instructions();
    exit(-1);
  }
  while ((c = getopt(argc, argv, OPTIONS)) != -1)
  {
    switch (c)
    {
    case 'd':
      dbFile = new string(optarg);
      bDB = true;
      break;
    case 'm':
      matrixFile = new string(optarg);
      bMatrix = true;
      break;
    case 'M':
      matOutFile = new string(optarg);
      bMatOut = true;
      break;
    case 'o':
      // define new output stream 
      bOut = true;
      outputFile = new string(optarg);
      break;
    case 't':
      threshold = (double) atof(optarg);
      if ( threshold < 0 )
      {
        cerr << "Warning threshold is smaller than 0. Reset to default 0.3" << endl;
        threshold = 0.3;
      }
      break;
    case 'g':
      gapScore = (double) atof(optarg);
      if ( gapScore < 0 )
      {
        cerr << "Warning: gap penalty is smaller than 0. Reset to default 0.2" << endl;
        gapScore = 0.2;
      }
      break;
    case 'w':
      minLength = atoi(optarg);
      if ( minLength < 1 ) 
      {
        cerr << "Warning: minimal length is smaller than 0. Reset to default 4." << endl;
        minLength = 4;
      }
      break;
    case 's':
      nbrShuffles = atoi(optarg);
      if ( nbrShuffles > 0 )
        bShuffle = true;
      break;
    case 'k':
      nbrSelection = atoi(optarg);
      if ( nbrSelection < 0 )
        nbrSelection = 500;
      break;
		case 'r':
			w1 = atoi(optarg);
			if ( w1 == 0 )
			{
				bReverse = false;
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

  if (!bDB || !bMatrix)
  {
    instructions();
    cleanup();
    exit(-1);
  }

  
  // open matrix file for reading
  PWMIO *pwmIO = new PWMIO(dbFile, READ);
  list<PWM *> dbMatrix;
  while (  pwmIO->IsOpen() )
  {
    // cerr << "next" << endl;
    myMatrix = pwmIO->ReadMatrix();
    if ( myMatrix == NULL )
      break;
    
    dbMatrix.push_back(myMatrix);
    cerr << *(myMatrix->GetID()) << endl;
  }
  delete pwmIO;  // close input stream
  delete dbFile;

  // check list of matrices
  if (dbMatrix.size() == 0)
  {
    cerr << "MotifComparison: No usable matrix found." << endl;
    cleanup();
    delete myMatrix;
  }

  // open outputstream, if available
  ofstream _ofs;
  if ( bOut )
    _ofs.open(outputFile->c_str(), ios::out);
    
  // open results file
  if ( bMatOut )
    matIO = new PWMIO(matOutFile, WRITE);
  
  // define alignment objects
  BlockAlignment * pBlock1 = NULL,
    * pBlock2 = NULL;
  
  // open matrix file
  pwmIO = new PWMIO(matrixFile, READ);
  while ( pwmIO->IsOpen() )
  {
    // get the next matrix 
    myMatrix = pwmIO->ReadMatrix();
    
    if ( myMatrix == NULL )
      break;
    
    // length of current motif
    w1 = myMatrix->Length();

    cerr << endl << endl << "++++++ matrix: " << *(myMatrix->GetID()) << endl;
    
    // iterate over all db matrices
    matIter = dbMatrix.begin();
    while (matIter != dbMatrix.end())
    {
      // (*matIter)->StderrPrintMatrix();
      wCommon = 0;
      
      // check if both motifs are the same
      if ( *(myMatrix->GetID()) == *((*matIter)->GetID()) )
      {
        cerr << "--Warning-- Skip same matrix " << *(myMatrix->GetID()) << " == " << *((*matIter)->GetID()) << endl;
        matIter++;
        continue;
      }
      
      cerr << "++++++ aligning with matrix: " << *((*matIter)->GetID()) << endl;
      
      w2 = (*matIter)->Length(); // motif length of block 2

			if ( bReverse )
			{
				// get the reverse complement of the block
				pPwmRev = (*matIter)->ReverseSubMatrix(w2-1,w2);
				pBlock2 = new BlockAlignment(w1,w2);
				pBlock2->UpdateMatrix(myMatrix, pPwmRev, threshold, gapScore);
				maxScoreRev = pBlock2->GetMaxValue();
			} 
			else
			{
				maxScoreRev = -999;
			}
			
      // create a new block alignment object
      pBlock1 = new BlockAlignment(w1,w2);

      // compute alignment scores
      pBlock1->UpdateMatrix(myMatrix, *matIter, threshold, gapScore);
      
      // get maximal scores
      maxScore = pBlock1->GetMaxValue();
      
      // get motif length
      if ( maxScore > maxScoreRev )
      {
        wCommon = pBlock1->GetPathLength();
        pBlock1->GetPathStart(&xStart, &yStart);
				bRevMax = false;
      }
      else
      {
        wCommon = pBlock2->GetPathLength();
        pBlock2->GetPathStart(&xStart, &yStart);
				bRevMax = true;
      }

      // report the corresponding matrices        
      if ( wCommon >= minLength )
      {
        pwm1 = myMatrix->SubMatrix(xStart, wCommon);
        if ( maxScore > maxScoreRev )
        {
          pwm2 = (*matIter)->SubMatrix(yStart, wCommon);
        }
        else
        {
          pwm2 = pPwmRev->SubMatrix(yStart, wCommon);
          maxScore = maxScoreRev;
        }

        
        // define ids of the two common matrices
        sprintf(cid, "motif1-%s+%s", (myMatrix->GetID())->c_str(),((*matIter)->GetID())->c_str());
        string *mID1 = new string(cid);
        pwm1->SetID(mID1);
        pwm1->SetScore(maxScore);
        
        sprintf(cid, "motif2-%s+%s", (myMatrix->GetID())->c_str(),((*matIter)->GetID())->c_str());
        string *mID2 = new string(cid);
        pwm2->SetID(mID2);
        pwm2->SetScore(maxScore);

        cerr << "Motif ids: " << *mID1 << "  " << *mID2 << endl;
        cerr << "Common: " << *(pwm1->GetConsensus()) << endl
             << "        " << *(pwm2->GetConsensus()) << endl;
        
        if ( bMatOut )
        {
          matIO->WriteMatrix(pwm1);
          matIO->WriteMatrix(pwm2);
        }
        
        if ( bOut )
        {
          _ofs << (myMatrix->GetID())->c_str() << "\t" 
            << w1 << "\t" 
            << ((*matIter)->GetID())->c_str() << "\t" 
            << w2 << "\t" 
            << wCommon << "\t" 
            << maxScore << "\t" 
            << *(pwm1->GetConsensus()) << "\t" 
            << *(pwm2->GetConsensus());
        }
        else
        {
          cout << (myMatrix->GetID())->c_str() << "\t" 
            << w1 << "\t" 
            << ((*matIter)->GetID())->c_str() << "\t" 
            << w2 << "\t" 
            << wCommon << "\t" 
            << maxScore << "\t"  
            << *(pwm1->GetConsensus()) << "\t" 
            << *(pwm2->GetConsensus());
        }
        
        // shuffle blocks if asked for
        if ( bShuffle )
        {
          // create vector to store results
          vector<double> *pScores = new vector<double>(nbrShuffles);
          double m = 0,
            temp = 0,
            std = 0,
            b_hat = 0,
            m_hat = 0,
            pvalue = 1,
            score1 = 0,
            score2 = 0;
          int value = 0,
            dummy = 0,
            l1 = 0,
            l2 = 0;
          int *indices1 = new int[w1];
          int *indices2 = new int[w2];
          
          for (int r=0; r<nbrShuffles; r++)
          {

            // create initial indices
            for ( int i=0; i<w1; i++)
              indices1[i] = i;
            // create initial indices
            for ( int i=0; i<w2; i++)
              indices2[i] = i;
          
            score1 = 0;
            for ( int l=0; l<nbrSelection; l++){

              // shuffle indices of matrix 1
              for ( int i=0; i<wCommon; i++ )
              {
                value = (int)(w1 * (random()/(RAND_MAX + 1.0)));
                dummy = indices1[value];
                indices1[value] = indices1[i];
                indices1[i] = dummy;

                value = (int)(w2 * (random()/(RAND_MAX + 1.0)));
                dummy = indices2[value];
                indices2[value] = indices2[i];
                indices2[i] = dummy;
              }
              
              
              // compute score of alignement
              score2 = 0;
              for ( int i=0; i<wCommon; i++)
              {
                l1 = indices1[i];
                l2 = indices2[i];
                temp = threshold;
                for (int k=0; k<4; k++)
                  temp -= ((myMatrix->GetValueAt(l1,k) -  (*matIter)->GetValueAt(l2,k)) * (log(myMatrix->GetValueAt(l1,k)) - log((*matIter)->GetValueAt(l2,k))));
      
                if ( temp > - gapScore ) 
                {
                  score2 += temp;
                }
                else
                {
                  score2 -= gapScore;
                }
              }
              if ( score1 < score2 )
                score1 = score2;              
            }

            // add score to vector
            cerr << "Shuffle\t" << r << "\t" << score1 << "\r";
            (*pScores)[r] = score1;
            m += score1;  // compute mean score
          }

          // delete indices
          delete[] indices1;
          delete[] indices2;
          
          // compute sample mean and standard deviation
          m = m/nbrShuffles;
          std = 0;
          for ( int r=0; r<nbrShuffles; r++ )
            std += pow(m - (*pScores)[r],2);
          std = sqrt(std/nbrShuffles);
                    
          // define EVD parameters
          b_hat = (std * sqrt(6.0)) / PI;
          m_hat = m - (0.5772 * b_hat);
          
          cerr << endl << "parameters: " << m << " " << std << " => " << m_hat << " " << b_hat << endl;

          // pvalue 
          pvalue = 1 - exp(-exp((m_hat-maxScore)/b_hat) );
          cerr << "Score " << maxScore << "  => p-value = " << pvalue << endl << endl;

          if ( bOut )
          {
            _ofs << "\t" << pvalue;
          }
          else
          {
            cout << "\t" << pvalue;
          }
          // delete score vector
          delete pScores;
        }

        
        if ( bOut )
        {
          _ofs << endl;
        }
        else
        {
          cout << endl;
        }
        
        
        // delete 
        if ( mID1 != NULL )
          delete mID1;
        mID1 = NULL;
        if ( mID2 != NULL )
          delete mID2;
        mID2 = NULL;
        
        // clean up internal variables 
        if ( pwm1 != NULL )
          delete pwm1;
        pwm1 = NULL;
        if ( pwm2 != NULL )
          delete pwm2;
        pwm2 = NULL;

      }
        

      // delete reverse complement of block2
      if ( pPwmRev != NULL )
          delete pPwmRev;
      pPwmRev = NULL;
              
              
      // delete local variables
      delete pBlock1;
      pBlock1 = NULL;
      delete pBlock2;
      pBlock2 = NULL;
          
      // go to the next matrix
      // cerr << "move on to the next one from db" << endl;
      matIter++;
    }

    // cerr << endl << " Clean local variables." << endl;
    // clean current matrix
    if ( myMatrix != NULL )
      delete myMatrix;
    myMatrix = NULL;
    
    // cerr << "move on to the next one" << endl;
  }

  // cerr << endl << " Clean local variables." << endl;
  // clean up all variables
  matIter = dbMatrix.begin();
  while (matIter != dbMatrix.end())
  {
    delete *matIter;
    *matIter = NULL;
    matIter++;
  }
  
  // close output stream if open
  if ( bOut )
    _ofs.close();
  delete outputFile;
  
  // motif model output stream
  if ( bMatOut )
  {
    delete matIO;
    delete matOutFile;
  }
  
  // close motif model input stream
  delete pwmIO;
  
  cleanup();
  
  exit(0);
}


void
cleanup()
{
  if (matrixFile != NULL)
    delete matrixFile;
  matrixFile = NULL;
  if (outputFile != NULL)
    delete outputFile;
  outputFile = NULL;
  if ( dbFile != NULL )
    delete dbFile;
  dbFile = NULL;
  if ( matOutFile != NULL )
    delete matOutFile;
  matOutFile = NULL;
  if ( shuffleFile != NULL )
    delete shuffleFile;
  shuffleFile = NULL;

  return;
}


void
version()
{
  cout << endl;
  cout << "INCLUSive -- Block2MotifAligner (stand alone C++ version)" << endl;
  cout << "  Version " << VERSION << endl;
} 


void
instructions()
{
  cout << endl;
  cout << "Usage:" << " Block2MotifAligner <ARGS>" << endl;
  cout << endl;
  cout << " Required Arguments" << endl;
  cout << "  -m <matrixFile>     File containing the query motif models." 
    << endl;
  cout << "  -d <matrixFile>     File containing database of models with which "
    << "                      all query motifs will be compared." 
    << endl;
  cout << endl;
  cout << " Optional Arguments" << endl;
  cout << "  -t <value>          Maximal distance between two motifs to be considered"
    << endl;
  cout << "                      as the same motif (default 0.4)" << endl;
  cout << "  -g <value>          Gap score (default 0.4)" << endl;
  cout << "  -w <value>          Minimal length of reported common motif (default 4)" << endl;
  cout << "  -s <value>          Number of shuffles of blocks to assess significance (default = 0)" << endl;
  cout << "    " << endl;
  cout << "  -o <outFile>        Output file to write results to." << endl;
  cout << "  -M <filename>       File to write common matrices." << endl;
  cout << endl;
  cout << "  -v                  Version of MotifComparison" << endl;
  cout << endl;
  cout << "Version " << VERSION << endl;
  cout << "Questions and Remarks: <gert.thijs@esat.kuleuven.ac.be>" << endl
    << endl;
}
