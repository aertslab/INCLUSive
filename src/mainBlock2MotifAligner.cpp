#include "inclusive.h"
#include "utilities.h"
#include "PWMIO.h"
#include "PWM.h"
#include "RandomNumber.h"


#include <iostream>
#include <fstream>
#include <getopt.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cstring>

#define OPTIONS "d:m:o:M:t:g:w:vh"
#define VERSION "3.0"

// function prototypes
void instructions();
void version();
void cleanup();

string *outputFile = NULL,
  *dbFile = NULL,
  *matrixFile = NULL,
  *matOutFile = NULL;

int
main(int argc, char *argv[])
{

  // initialize random number generator
  RandomNumber::InitRandomNumber();

  char c;
  double score = 0,
    maxScore = 0,
    threshold = 0.3,
    gapScore = 0.2,
    L1 = 0, L2 = 0;
  int maxIndexX = 0,
    maxIndexY = 0,
    w1 = 0,
    w2 = 0,
    wCommon = 0,
    minLength = 4,
    i = 0,
    j = 0,
    k = 0,
    ii = 0,
    jj = 0;
  bool bDB = false,
    bMatrix = false,
    bOut = false,
    bMatOut = false;
  double **pAlingMatrix;
  PWMIO* matIO = NULL;
  PWM * pwm1 = NULL;
  PWM * pwm2 = NULL;
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
  PWM *myMatrix;
  list<PWM *> dbMatrix;
  while ( (myMatrix = pwmIO->ReadMatrix()) )
  {
    dbMatrix.push_back(myMatrix);
    cerr << *(myMatrix->GetID()) << endl;
  }
  delete pwmIO;  // close input stream

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
  {
    _ofs.open(outputFile->c_str(), ios::out);
    // _ofs << "#INCLUSive Block2MotifAligner Output" << endl;
  }
  else
  {
    // cout << "#INCLUSive Block2MotifAligner Output" << endl;
  }
    
  // open results file
  if ( bMatOut )
  {
    matIO = new PWMIO(matOutFile, WRITE);
  }
  
  // open matrix file
  pwmIO = new PWMIO(matrixFile, READ);
  if (!pwmIO->IsOpen())
  {
    cerr << "Block2MotifAligner: Unable to open file with matrix models: " << *matrixFile << endl;
    cleanup();
    exit(-1);
  }

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

      // check if both motifs are the same
      if ( *(myMatrix->GetID()) == *((*matIter)->GetID()) )
      {
        cerr << "--Warning: samer matrix found" << *(myMatrix->GetID()) << " == " << *((*matIter)->GetID()) << endl;
      }
      else
      {
        // create 2dim array to store scores
        w2 = (*matIter)->Length();
        pAlingMatrix = new double*[w1];
        for (i=0; i<w1; i++)
          pAlingMatrix[i] = new double[w2];
        
        maxScore = 0;
        maxIndexX = 0;
        maxIndexY = 0;
        for (i=0; i<w1; i++)
        {
          for (j=0; j<w2; j++ )
          {
            score = threshold;
            for (k=0; k<4; k++)
            {
              score -= ((myMatrix->GetValueAt(i,k) -  (*matIter)->GetValueAt(j,k)) * 
                (log(myMatrix->GetValueAt(i,k)) - log((*matIter)->GetValueAt(j,k))));
            }
  
            if ( i && j )
            {
              L1 = pAlingMatrix[i-1][j-1] + score;
              L2 = pAlingMatrix[i-1][j-1] - gapScore;
            }
            else
            {             
              L1 = score;
              L2 = - gapScore;
            }
            
            if ( L1 < 0 && L2 < 0 )
            {
              pAlingMatrix[i][j] = 0;
            }
            else if ( L1 > L2 )
            {
              pAlingMatrix[i][j] = L1;
            }
            else
            {
              pAlingMatrix[i][j] = L2;
            }
  
            if ( pAlingMatrix[i][j] > maxScore)
            { 
              maxScore = pAlingMatrix[i][j];
              maxIndexX = i;
              maxIndexY = j;
            }
            
            // cerr << pAlingMatrix[i][j] << " ";
          }
          // cerr << endl;
        }  
        
        // cerr << maxScore << " at (" << maxIndexX << "," << maxIndexY << ")" << endl;
        
        // backtrack highest scoring path
        ii = maxIndexX;
        jj = maxIndexY;
        while ( ii && jj && (pAlingMatrix[ii-1][jj-1] > 0) )
        {
          ii--; 
          jj--;
        }        
        
        // report the corresponding matrices
        wCommon = maxIndexX - ii + 1;
        // cerr << "DEBUG: common length found: " << wCommon << " " << ii << " " << jj << endl;
        
        if ( wCommon >= minLength )
        {
          
          pwm1 = myMatrix->SubMatrix(ii, wCommon);
          pwm2 = (*matIter)->SubMatrix(jj, wCommon);
    
          // define ids of the two common matrices
          sprintf(cid, "motif1-%s+%s", (myMatrix->GetID())->c_str(),((*matIter)->GetID())->c_str());
          string *mID1 = new string(cid);
          pwm1->SetID(mID1);
          
          sprintf(cid, "motif2-%s+%s", (myMatrix->GetID())->c_str(),((*matIter)->GetID())->c_str());
          string *mID2 = new string(cid);
          pwm2->SetID(mID2);
  
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
            _ofs << *(myMatrix->GetID()) << "\t" 
              << (myMatrix->Length()) << "\t" 
              << *((*matIter)->GetID()) << "\t" 
              << ((*matIter)->Length()) << "\t" 
              << wCommon << "\t" 
              << pAlingMatrix[maxIndexX][maxIndexY] << "\t" 
              << *(pwm1->GetConsensus()) << "\t" 
              << *(pwm2->GetConsensus()) << endl;
            
            //~ for (i=0; i<wCommon; i++)
            //~ {
              //~ cerr << myMatrix->GetConsensusSymbolAt(ii+i) << " | ";
              //~ cerr << (*matIter)->GetConsensusSymbolAt(jj+i) << endl;  
            //~ }
          }
          else
          {
            cout << *(myMatrix->GetID()) << "\t" 
              << (myMatrix->Length()) << "\t" 
              << *((*matIter)->GetID()) << "\t" 
              << ((*matIter)->Length()) << "\t" 
              << wCommon << "\t" 
              << pAlingMatrix[maxIndexX][maxIndexY] << "\t"  
              << *(pwm1->GetConsensus()) << "\t" 
              << *(pwm2->GetConsensus()) << endl;
            
            //~ for (i=0; i<wCommon; i++)
            //~ {
              //~ cerr << myMatrix->GetConsensusSymbolAt(ii+i) << " | ";
              //~ cerr << (*matIter)->GetConsensusSymbolAt(jj+i) << endl;  
            //~ }
          }
          
          // clean up internal variables
          if ( pwm1 != NULL )
            delete pwm1;
          pwm1 = NULL;
          if ( pwm2 != NULL )
            delete pwm2;
          pwm2 = NULL;
          
          if ( mID1 != NULL )
            delete mID1;
          if ( mID2 != NULL )
            delete mID2;
          mID1 = NULL;
          mID2 = NULL;
          
        }
        
        // delete the score matrixFile
        for (i=0; i<w1; i++ )
          delete[] pAlingMatrix[i];
        delete[] pAlingMatrix;
       
      }        
        
      // go to the next matrix
      matIter++;
    }
    // clean current matrix
    if ( myMatrix != NULL )
      delete myMatrix;
    myMatrix = NULL;
    
  }
  
  _ofs.close();
  
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
  cout << "  -o <outFile>        Output file to write results to." << endl;
  cout << "  -M <filename>       File to write common matrices." << endl;
  cout << endl;
  cout << "  -v                  Version of MotifComparison" << endl;
  cout << endl;
  cout << "Version " << VERSION << endl;
  cout << "Questions and Remarks: <gert.thijs@esat.kuleuven.ac.be>" << endl
    << endl;
}
