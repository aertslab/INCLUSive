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
#include <algorithm>
#include <iterator>

#define OPTIONS "m:o:r:hv"
#define VERSION "3.0"

// using std::random_shuffle;

// function prototypes
void instructions();
void version();
void cleanup();

string *outputFile = NULL,
  *matrixFile = NULL;

int 
main(int argc, char *argv[])
{
  // initialize random number generator
  RandomNumber::InitRandomNumber();
  srandom( 1000 * (long)time(NULL) );
  
  bool bOut = false,
    bMatFile = false;
  int i = 0,
    j = 0, 
    k = 0,
    dummy = 0,
    value = 0,
    w = 0,
    shuffles = 1;
  char c;

  if (optopt == 0)
  {
    instructions();
    exit(-1);
  }
  while ((c = getopt(argc, argv, OPTIONS)) != -1)
  {
    switch (c)
    {
    case 'm':
      // set fasta file name
      matrixFile = new string(optarg);
      bMatFile = true;
      break;
    case 'o':
      // set fasta file name
      outputFile = new string(optarg);
      bOut = true;
      break;
    case 'r':
      shuffles = atoi(optarg);
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
  
  if (!bMatFile || !bOut)
  {
    instructions();
    cleanup();
    exit(-1);
  }
  
  // open file for output writing if necessarry
  ofstream ofs;
  ofs.open(outputFile->c_str(), ios::out);
  
  if ( !ofs.good() )
  {
    cerr << "--ERROR-- Unable to open file for writing" << endl;
    cleanup();
    exit(-1);
  }
  else
  {
    ofs << "#INCLUSive Motif Model v1.0" << endl << "#" << endl << endl;
  }
  
  // process matrix file
  PWMIO *pwmIO = new PWMIO(matrixFile, READ);
  PWM *myMatrix;

  while ( pwmIO->IsOpen() )
  {
    // read matrix
    myMatrix = pwmIO->ReadMatrix();
    
    if ( myMatrix == NULL )
    {
      pwmIO->Close();
      break;
    }
    
    // create index vector
    w = myMatrix->Length();
    vector<int> vIndex;
    for (i=0; i<w; i++)
      vIndex.push_back(i);
    
    // shuffle the vector
    for (int r=0; r<shuffles; r++)
    {

      for ( i=0; i<w; i++ )
      {
        value = (int)(w * (random()/(RAND_MAX + 1.0)));
        dummy = vIndex[value];
        vIndex[value] = vIndex[i];
        vIndex[i] = dummy;
      }
      
      //~for (i=0; i<(int)vIndex.size(); i++)
        //~cerr << vIndex[i] << " ";
      //~cerr << endl;

      ofs << "#ID = shuffled_" << *(myMatrix->GetID()) << "_" << r << endl;
      ofs << "#W = " << w << endl;
      ofs << "#Score = " << myMatrix->Score() << endl;
      
      for ( i=0; i<w; i++ )
      {
        value = vIndex[i];
        ofs << myMatrix->GetValueAt(value, 0) << "\t";
        ofs << myMatrix->GetValueAt(value, 1) << "\t";
        ofs << myMatrix->GetValueAt(value, 2) << "\t";
        ofs << myMatrix->GetValueAt(value, 3) << endl;
      }
      ofs << endl; 
    }
    
    if ( myMatrix != NULL )
      delete myMatrix;
    myMatrix = NULL;
  }
  
  delete pwmIO;
  ofs.close();
  
  cleanup();
  
  exit(0);
}


void
instructions()
{
  cout << endl;
  cout << "Usage:" << " MotifModelShuffler -m <inFile> -o <outFile>" << endl;
  cout << endl;
  cout << endl;
  cout << "Version " << VERSION << endl;
  cout << "Questions and Remarks: <gert.thijs@esat.kuleuven.ac.be>" << endl
    << endl;
} void

version()
{
  cout << endl;
  cout << "INCLUSive -- MotifModelShuffler (stand alone C++ version)" << endl;
  cout << "  Version " << VERSION << endl;
} void

cleanup()
{
  if (matrixFile != NULL)
    delete matrixFile;
  matrixFile = NULL;
  if (outputFile != NULL)
    delete outputFile;
  outputFile = NULL;
}
