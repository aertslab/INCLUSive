#include "inclusive.h"
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

#define OPTIONS "hd:m:t:o:s:v"
#define VERSION "3.0"

void instructions();
void version();
void cleanup();
string *matrixFile = NULL,
  *dbFile = NULL,
  *outputFile = NULL;


int
main(int argc, char *argv[])
{
  char c;
  bool bOut = false,
    bDB = false,
    bMatrix = false;
  double threshold = 0.4;
  int shift = 0;
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

      // set fasta file name
      dbFile = new string(optarg);
      bDB = true;
      break;
    case 'm':
      matrixFile = new string(optarg);
      bMatrix = true;
      break;
    case 't':
      threshold = (double)atof(optarg);
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

  // Read database matrices
  PWMIO *pwmIO = new PWMIO(dbFile, READ);
  if (!pwmIO->IsOpen())
  {
    cerr << "MotifLocator: Unable to open dbfile with matrix models: " <<
      *dbFile << endl;
    cleanup();
    exit(-1);
  }
  PWM *myMatrix;
  list < PWM * >dbMatrix;
  list < PWM * >::iterator matIter;
  while (pwmIO->IsOpen() && (myMatrix = pwmIO->ReadMatrix()))
  {
    dbMatrix.push_back(myMatrix);
  }

  // close matrix reader
  delete pwmIO;
  if (dbMatrix.size() == 0)
  {
    cerr << "MotifComparison: No usable matrix found." << endl;
    cleanup();
    delete myMatrix;
  }

  // open matrix file for reading
  pwmIO = new PWMIO(matrixFile, READ);
  if (!pwmIO->IsOpen())
  {
    cerr << "MotifComparison: Unable to open dbfile with matrix models: " <<
      *matrixFile << endl;
    cleanup();
    exit(-1);
  }

  // open file for output writing if necessarry
  ofstream OFS;
  if (bOut)
  {
    OFS.open(outputFile->c_str(), ios::out);
    OFS << "#query\tquery_consensus\tmatch\tmatch_consensus\tscore" << endl;
  }
  else
  {
    cout << "#query\tquery_consensus\tmatch\tmatch_consensus\tscore" << endl;
  }
  int counter = 0;
  double mi;

  // bool found = false;
  // iterate over all matrices
  while (pwmIO->IsOpen() && (myMatrix = pwmIO->ReadMatrix()))
  {
    matIter = dbMatrix.begin();
    while (matIter != dbMatrix.end())
    {
      cerr << ++counter << "\t" << *(myMatrix->
                                     GetID()) << " <-> " << *((*matIter)->
                                                              GetID()) <<
        "\r";

      // compare matrices
      mi =
        (myMatrix->MutualInformation((*matIter), shift) +
         (*matIter)->MutualInformation(myMatrix, shift)) / 2;

      // check mi
      if (mi < threshold)
      {
        if (bOut == false)
        {
          cout << *(myMatrix->GetID()) << "\t" 
            << *(myMatrix->GetConsensus()) << "\t"
            << *((*matIter)->GetID()) << "\t" 
            << *((*matIter)->GetConsensus()) << "\t" 
            << mi << endl;
        }
        else
        {
          OFS << *(myMatrix->GetID()) << "\t" << *(myMatrix->
                                                   GetConsensus()) << "\t" <<
            *((*matIter)->GetID()) << "\t" << *((*matIter)->
                                                GetConsensus()) << "\t" << mi
            << endl;
        }
      }

      // next matrix
      matIter++;
    }
  }
  if (bOut)
    OFS.close();
  cerr << endl;

  // close matrix reader
  delete pwmIO;

  // clear dbMatrix
  matIter = dbMatrix.begin();
  while (matIter != dbMatrix.end())
  {

    // cerr << "deleting matrix: " << *((*matIter)->GetID()) << endl;
    delete(*matIter);
    matIter++;
  }
  delete myMatrix;
  cleanup();
  
  // program is succesful
  exit(0);
}



/*********************
 *  LOCAL FUNCTIONS  *
 *********************/
void
instructions()
{
  cout << endl;
  cout << "Usage:" << " MotifComparison <ARGS>" << endl;
  cout << endl;
  cout << " Required Arguments" << endl;
  cout << "  -m <matrixFile>     File containing the query motif models." 
    << endl;
  cout << "  -d <matrixFile>     File containing database of models with which "
    << "                      all query motifs will be compared." 
    << endl;
  cout << endl;
  cout << " Optional Arguments" << endl;
  cout <<
    "  -t <value>          Sets threshold below two motifs are seen as being the same."
    << endl;
  cout << "                      default value is set to 0.4." << endl;
  cout << "  -s <value>          Maximal allowed shift between motif models."
    << endl;
  cout << "  -o <outFile>        Output file to write results to." << endl;
  cout << endl;
  cout << "  -v                  Version of MotifComparison" << endl;
  cout << endl;
  cout << "Version " << VERSION << endl;
  cout << "Questions and Remarks: <gert.thijs@esat.kuleuven.ac.be>" << endl
    << endl;
} 


void
version()
{
  cout << endl;
  cout << "INCLUSive -- MotifComparison (stand alone C++ version)" << endl;
  cout << "  Version " << VERSION << endl;
} 


void
cleanup()
{
  if (dbFile != NULL)
    delete dbFile;
  dbFile = 0;
  if (matrixFile != NULL)
    delete matrixFile;
  matrixFile = NULL;
  if (outputFile != NULL)
    delete outputFile;
  outputFile = NULL;
  return;
}
