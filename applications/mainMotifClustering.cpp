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

#define OPTIONS "hf:m:c:t:p:v"

void instructions();
void version();
void cleanup();

string *matrixFile = NULL,
  *outputFile = NULL,
  *resultsFile = NULL;


int
main(int argc, char *argv[])
{
	
		
cerr << endl;
cerr << "At date, this application is not actively used neither maintained. "
<< "Contact <mclaeys@qatar.net.qa> if you have questions "
		<< "or interest to (re)use it. "
<< endl;
exit(0);
	
	
	
  char c;
  bool bOut = false,
    bMatrix = false,
    bIn = true;
  double threshold = 0.4,
    //score = 0,
    //maxScore = 0,
    percentage = 0.8;
  int shift = 0,
    clusters = 3;
  
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
      resultsFile = new string(optarg);
      bIn = true;
      break;
    case 'm':
      matrixFile = new string(optarg);
      bMatrix = true;
      break;
    case 'c':
      clusters = atoi(optarg);
      if ( clusters < 1 )
      {
        cerr << "--Error-- Number of clusters should be a positive number" << endl;
        exit(-1);
      }
      break;
    case 't':
      threshold = (double) atof(optarg);
      if ( threshold < 0 )
      {
        cerr << "--Error-- Threshold should be a positive number" << endl;
        exit(-1);
      }
      break;
    case 'o':
      // define new output stream 
      bOut = true;
      outputFile = new string(optarg);
      break;
    case 'p':
      percentage = (double) atof(optarg);
      if ( percentage < 0 || percentage > 1 )
      {
        cerr << "--Error-- percentage of overlap should between 0 and 1.0" << endl;
        exit(-1);
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
      cerr << "BlockRanking: Error in getopt() function" << endl;
      exit(-1);
    }
  }

  if ( !bMatrix || !bIn )
  {
    instructions();
    cleanup();
    exit(-1);
  }

  // process the output file of BlockAligner
  ifstream _ifs(resultsFile->c_str(), ios::in);
  string inLine;
  while (!_ifs.eof())
  {
    // read next line
    getline(_ifs, inLine, '\n');
    
    // check structure of the line
    
  }
  
  // process matrix file
  PWMIO *pwmIO = new PWMIO(matrixFile, READ);
  PWM *myMatrix;
  list <PWM*>matrixList;
  list <PWM*>::iterator matIter;

  // open matrix stream reader
  if (!pwmIO->IsOpen())
  {
    cerr << "BlockRanking: Unable to open matrix model file: " << *matrixFile << endl;
    cleanup();
    exit(-1);
  }

  // add matrices to list
  while ( pwmIO->IsOpen() )
  {
    myMatrix = pwmIO->ReadMatrix();
    if ( myMatrix != NULL )
      matrixList.push_back(myMatrix);
  }
  // close matrix reader
  delete pwmIO;

  if ( matrixList.size() < 1 )
  {
    cerr <<
      "BlockRanking: Number of matrices in a list is samller than the requested number of motifs."
      << endl;
    cleanup();
    delete myMatrix;
  }

  // compute distance between motifs
  //uint L = matrixList.size();
  double mi = 0;
  
  list<PWM*>::iterator iter2;
  for(matIter = matrixList.begin(); matIter != matrixList.end(); matIter++)
  {
    // cerr << *((*matIter)->GetID()) << endl;
    iter2 = matIter;
    iter2++;
    for (; iter2 != matrixList.end(); iter2++)
    {
      // compute shift
      if ( (*iter2)->Length() < (*matIter)->Length() )
      {
        shift = (int)((1-percentage) * (*iter2)->Length());
      }
      else
      {  
        shift = (int)((1-percentage) * (*matIter)->Length());
      }
        
      mi = ((*iter2)->MutualInformation((*matIter), shift) +
         (*matIter)->MutualInformation((*iter2), shift)) / 2;
      
      cout << mi << endl;
    }
    // cerr << endl;
  }
  
  
 // clear matrix list
  matIter = matrixList.begin();
  while (matIter != matrixList.end())
  {
    delete(*matIter);
    matIter++;
  }
  if (myMatrix != NULL)
    delete myMatrix;
  myMatrix = NULL;
  
  cleanup();
  
  // program was succesful
  exit(0);
}



// local subroutines
void
version()
{
  cout << endl;
  cout << "INCLUSive -- BlockRanking (stand alone C++ version)" << endl;
  cout << "  version " << VERSION << endl;
  return;
} 

void
cleanup()
{
  if (resultsFile != NULL)
    delete resultsFile;
  resultsFile = NULL;
  if (matrixFile != NULL)
    delete matrixFile;
  matrixFile = NULL;
  if (outputFile != NULL)
    delete outputFile;
  outputFile = NULL;
  return;
}

void 
instructions()
{
  cout << endl;
  cout << "Usage:" << " MotifClustering -m <MatrixFile> <OPT ARGS>" << endl;
  cout << endl;
  cout << " Required Arguments" << endl;
  cout << "  -f <results file>   File containing the results of BlockAligner" << endl;
  cout << "  -m <matrixFile>     File containing the matrix model descriptions" << endl;
  cout << endl;
  cout << " Optional Arguments" << endl;
  cout << "  -t <value>          threshold on distance between motifs " << endl;
  cout << "  -p <value>          percentage overlap when comparing to motifs" << endl;
  cout << "  - " << endl;
  cout << endl;
  cout << "  -v                  Version of MotifClustering" << endl;
  cout << endl;
  cout << "Version " << VERSION << endl;
  cout << "Questions and Remarks: <gert.thijs@esat.kuleuven.ac.be>" << endl
    << endl;

  return;
}
