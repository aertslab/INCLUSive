#include "utilities.h"
#include "PWMIO.h"
#include "PWM.h"

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

#define OPTIONS "hd:i:r:t:o:s:v"
#define VERSION "3.0"


void instructions();
void version();
void cleanup();
string *matrixFile = NULL,
  *outputFile = NULL;

int
main(int argc, char *argv[])
{
  char c;
  bool bOut = false,
    bMatrix = false;
  double threshold = 0.4,
    score = 0,
    maxScore = 0;
  int shift = 1;
  int rank = 5;
  int mode = 0;
  if (optopt == 0)
  {
    instructions();
    exit(-1);
  }
  while ((c = getopt(argc, argv, OPTIONS)) != -1)
  {
    switch (c)
    {
    case 'i':
      matrixFile = new string(optarg);
      bMatrix = true;
      break;
    case 'r':
      rank = atoi(optarg);
      break;
    case 't':
      threshold = (double) atof(optarg);
      break;
    case 'o':
      // define new output stream 
      bOut = true;
      outputFile = new string(optarg);
      break;
    case 's':
      shift = atoi(optarg);
      break;
    case 'v':
      version();
      exit(-1);
    case '?':
    case 'h':
      instructions();
      exit(-1);
    default:
      cerr << "BlockRanking: Error in getopt() function" << endl;
      exit(-1);
    }
  }

  if (!bMatrix || !bOut)
  {
    instructions();
    cleanup();
    exit(-1);
  }

  // open file for output writing if necessarry
  PWMIO *pwmFile = new PWMIO(outputFile, WRITE);
  
  // process matrix file
  PWMIO *pwmIO = new PWMIO(matrixFile, READ);
  PWM *myMatrix;
  list <PWM*>matrixList;
  list <PWM*>::iterator matIter;

  // open matrix stream reader
  if (!pwmIO->IsOpen())
  {
    cerr << "BlockRanking: Unable to open matrix model file: " << *matrixFile
      << endl;
    cleanup();
    exit(-1);
  }

  // add matrices to list
  while ( pwmIO->IsOpen() )
  {
    myMatrix = pwmIO->ReadMatrix();
    if ( myMatrix )
      matrixList.push_back(myMatrix);
  }
  // close matrix reader
  delete pwmIO;

  if ((int) matrixList.size() < rank)
  {
    cerr <<
      "BlockRanking: Number of matrices in a list is samller than the requested number of motifs."
      << endl;
    cleanup();
    delete myMatrix;
  }

  int count = 0,
    shift1 = 0,
    shift2 = 0;
  for (int i = 0; i < rank; i++)
  {

    // look for best scoring motif 
    matIter = matrixList.begin();
    maxScore = 0;
    while (matIter != matrixList.end())
    {
      score = (*matIter)->Score();
      if (maxScore < score)
      {
        myMatrix = (*matIter);
        maxScore = myMatrix->Score();
      }
      matIter++;
    }

    // write best scoring motif to file
    if (maxScore != 0)
    {
      pwmFile->WriteMatrix(myMatrix);

      // reset score 
      myMatrix->SetScore(0);
      cerr << endl << "= " << *(myMatrix->GetID()) << endl;
      count = 0;
      double mi = 0;
      matIter = matrixList.begin();
      while (matIter != matrixList.end())
      {
        shift1 = (int)(0.4 * myMatrix->Length());
        shift2 = (int)(0.4 * (*matIter)->Length());
        shift = shift1 < shift2 ? shift1 : shift2;

        mi =
          (myMatrix->MutualInformation((*matIter), shift) +
           (*matIter)->MutualInformation(myMatrix, shift)) / 2;
  
        // find similar motif
        if (mi < threshold)
        {
          count++;
          (*matIter)->SetScore(0);
          cerr << "   " << *((*matIter)->GetID()) << "\t" << *((*matIter)->GetConsensus()) << endl;
        }
        matIter++;
      }
      cerr << "==> Total: " << count << endl;
    }
  }

  // let myMatrix point to NULL
  myMatrix = NULL;

  // clear matrix list
  matIter = matrixList.begin();
  while (matIter != matrixList.end())
  {
    delete(*matIter);
    matIter++;
  }
  if (myMatrix != NULL)
  {
    delete myMatrix;
  }
  cleanup();
  
  // program was succesful
  exit(0);
}



/*********************
 *  LOCAL FUNCTIONS  *
 *********************/
void
instructions()
{
  cout << endl;
  cout << "Usage:" << " BlockRanking <ARGS>" << endl;
  cout << endl;
  cout << " Required Arguments" << endl;
  cout <<
    "  -i <matrixFile>     File containing the matrix model descriptions" <<
    endl;
  cout << "  -o <outFile>        Output file to write results to." << endl;
  cout << endl;
  cout << " Optional Arguments" << endl;
  cout << "  -m <0|1|2>          Select the type of score: " << endl;
  cout << "                       0 = Score defined in matrix file (default)"
    << endl;
  cout << "                       1 = Consensus score [2+plog(p)]" << endl;
  cout << "                       2 = Information Content [plog(p/p0)] (requires -b)" <<
    endl;
  cout <<
    "  -b <file>           Load Background model from file. Sets single nucleotide frequency"
    << endl;
  cout <<
    "  -r <value>          Number of topscoring motifs to display (default 5)."
    << endl;
  cout <<
    "  -t <value>          Sets threshold below two motifs are seen as being the same."
    << endl;
  cout << "                      Default value is set to 0.4." << endl;
  cout << "  -s <value>          Maximal allowed shift between motif models."
    << endl;
  cout << endl;
  cout << "  -v                  Version of BlockRanking" << endl;
  cout << endl;
  cout << "Version " << VERSION << endl;
  cout << "Questions and Remarks: <gert.thijs@esat.kuleuven.ac.be>" << endl
    << endl;
} void

version()
{
  cout << endl;
  cout << "INCLUSive -- BlockRanking (stand alone C++ version)" << endl;
  cout << "  version " << VERSION << endl;
} void

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
