#include "inclusive.h"
#include "utilities.h"
#include "RandomNumber.h"
#include "BackgroundIO.h"
#include "BackgroundModel.h"
#include <getopt.h>

#define OPTIONS "f:b:l:n:L:v"
#define VERSION "3.0"

using namespace INCLUSIVE;

// function prototypes
void instructions();
void version();
void cleanup();

// global file name variables
string *fastaFile = NULL,
  *bgFile = NULL,
  *lengthFile= NULL;
  
int
main(int argc, char *argv[])
{
  // initialize random number generator
  RandomNumber::InitRandomNumber();

  // local variables  
  char c;
  bool bLength = false,
    bBG = false,
    bLengthFile = false,
    bFasta = false;
  int nbr = 1,
    length = 100;
  
  // get options
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
    case 'l':
      length = atoi(optarg);
      bLength = true;
      break;
    case 'L':
      lengthFile = new string(optarg);
      bLengthFile = true;
      break;
    case 'n':
      nbr = atoi(optarg);
      if (nbr < 1)
      {
        cerr <<
          "Warning\n SequenceCreator: number of sequences should be greater than 1. Default = 100."
          << endl;
        instructions();
        exit(-1);
      }
      break;
    case '?':
    case 'h':
      instructions();
      exit(-1);
    case 'v':
      version();
      exit(-1);
    default:
      cerr << "-- SequenceCreator: Error in getopt() function" << endl;
      exit(-1);
    }
  }

  if ( bLength && bLengthFile ){
    cerr << "-- SequenceCreator: Use of -l and -L simultaneously.";
    instructions();
    exit(-1);
  }

  // open background file for reading
  BackgroundIO *bgIO = new BackgroundIO(fileName, READ);
  BackgroundModel *pBgModel;
  if (bgIO != NULL)
  {
    pBgModel = bgIO->ReadBackgroundModel();
    delete bgIO;
  }
  else
  {
    cerr << "Error: Problems opening background model file for reading." <<
      endl;
    cleanup();
    exit(-1);
  }
  
  // loop over all sequences
  for (int s=0; s<nbr; s++)
  {
    
  }
  
  exit(0);
}



void 
instructions()
{
  cout << endl;
  cout << "Usage:" << " SequenceCreator <ARGS>" << endl;
  cout << endl;
  cout << "   " << endl;
  cout << "  -b <bgFile>         File containing the background model description"
    << endl;
  cout << "  -f <fastaFile>      Fasta file with output sequences. When ommitted " 
    << endl;
  cout << "                      sequences are written on stdout." << endl;
  cout << endl;
  cout << "If you like all sequences of the same length, choose  " << endl;
  cout << "  -l <length>         Length of the sequences."
    << endl;
  cout << "  -n <number>         Number of different sequences."
    << endl;
  cout << endl;
  cout << "If you like to create sequences with different lengths, choose  " << endl;
  cout << "  -L <file>           File with sequence lengths on different lines."
    << endl;   
  cout << endl;
  cout << "Version " << VERSION << endl;
  cout << "Questions and Remarks: <gert.thijs@esat.kuleuven.ac.be>" << endl
    << endl;

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
}
