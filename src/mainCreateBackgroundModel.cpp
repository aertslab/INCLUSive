#include "inclusive.h"
#include "BackgroundIO.h"
#include "BackgroundModel.h"
#include "FastaIO.h"
#include "SequenceObject.h"

#include <iostream>
#include <fstream>
#include <getopt.h>
#include <time.h>
#include <cstring>
#include <stdlib.h>
#include <math.h>


#define OPTIONS "hf:b:o:n:v"
#define VERSION "3.0"

// function prototypes
void instructions();
void version();
void cleanup();

string *fastaFile = NULL,
  *bgFile = NULL;


int
main(int argc, char *argv[])
{
  bool bFasta = false,
    bBG = false;
  int order = 1;
  char c;
  string *organism = NULL;

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
      order = atoi(optarg);
      break;
    case 'n':
      organism = new string(optarg);
      break;
    case '?':
    case 'h':
      instructions();
      exit(-1);
    case 'v':
      version();
      exit(-1);
    default:
      cerr << "-- MotifSampler: Error in getopt() function" << endl;
      exit(-1);
    }
  }

  if (!bFasta || !bBG || order < 0)
  {
    instructions();
    cleanup();
    exit(-1);
  }

  // open fasta file 
  FastaIO *fileIO = new FastaIO(fastaFile);
  SequenceObject *pSeq = NULL;

  if (fileIO == NULL)
  {
    cerr << "CreateBackgroundModel: Unable to read sequences from fastaFile:"
      << *fastaFile << endl;
    cleanup();
    exit(-1);
  }

  // define matrices
  int l = (int) pow(4.0, order);

  // create matrices to store data
  double *pSnf = new double[4];
  for (int i = 0; i < 4; i++)
    pSnf[i] = 0;

  int ll = l;
  if (order == 0)
    ll = 4;

  double *pFreq = new double[ll];
  for (int i = 0; i < l; i++)
    pFreq[i] = 0;

  double (*pTrans)[4] = new double[l][4];
  for (int i = 0; i < l; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      pTrans[i][j] = 0;
    }
  }

  int nt,
    index,
    nt1;

  while (fileIO->HasNext())
  {
    pSeq = fileIO->NextSequence();
    cerr << *(pSeq->GetID()) << " ";
    for (int i = 0; i < pSeq->Length(); i++)
    {
      nt = pSeq->GetNucleotideAt(plus_strand, i);
      if (nt != -1)
      {
        pSnf[nt] += 1;

        if (i >= order)
        {
          index = 0;
          for (int j = 0; j < order; j++)
          {
            nt1 = pSeq->GetNucleotideAt(plus_strand, i - order + j);
            if (nt1 == -1)
            {
              // this is a non ACGT symbol
              index = -1;
              break;
            }
            index += (((int) pow(4.0, order - j - 1)) * nt1);
          }

          if (index != -1)
          {
            // add one to frequency matrix
            pFreq[index] += 1;

            // add 1 to right entry in transition matrix
            pTrans[index][nt] += 1;
          }
        }
      }
    }

    // get last oligo of length order
    index = 0;
    for (int j = 0; j < order; j++)
    {
      nt1 = pSeq->GetNucleotideAt(plus_strand, pSeq->Length() - order + j);
      if (nt1 == -1)
      {
        // this is a non ACGT symbol
        index = -1;
        break;
      }
      index += (((int) pow(4.0, order - j - 1)) * nt1);
    }

    if (index != -1)
      pFreq[index]++;

    cerr << "done" << endl;
    // delete sequence
    delete pSeq;
    pSeq = NULL;
  }

  // close fasta file
  if (fileIO != NULL)
    delete fileIO;

  // total number of nucleotides
  double ntCount = 0;
  for (int i = 0; i < 4; i++)
    ntCount += pSnf[i];

  // pseudo count factors
  double pseudoFreq = l / ntCount;
  double pseudoTrans = (4 * l) / ntCount;

  // normalize matrices
  // SNF
  cerr << "SNF: ";
  for (int i = 0; i < 4; i++)
  {
    pSnf[i] = pSnf[i] / ntCount;
    cerr << pSnf[i] << " ";
  }
  cerr << endl;

  // oligo frequency
  double norm = 0;
  cerr << "Pfreq: ";
  for (int i = 0; i < l; i++)
    norm += (pFreq[i] + pseudoFreq);
  for (int i = 0; i < l; i++)
  {
    pFreq[i] = (pFreq[i] + pseudoFreq) / norm;
    cerr << pFreq[i] << " ";
  }
  cerr << endl;

  // transition matrix
  cerr << "Transition Matrix: " << endl;
  for (int i = 0; i < l; i++)
  {
    norm = 0;
    for (int j = 0; j < 4; j++)
      norm += (pTrans[i][j] + pseudoTrans);


    for (int j = 0; j < 4; j++)
    {
      pTrans[i][j] = (pTrans[i][j] + pseudoTrans) / norm;
      cerr << pTrans[i][j] << "  ";
    }
    cerr << endl;
  }
  // create new background model
  BackgroundModel *bg = new BackgroundModel(order, pTrans, pFreq, pSnf);

  if (bg != NULL)
  {
    // define organism name
    if (organism == NULL)
      organism = new string("unknown");

    bg->SetOrganism(organism);
    bg->SetSequences(fastaFile);
    // write background model
    BackgroundIO *bgIO = new BackgroundIO(*bgFile, WRITE);
    bgIO->WriteBackgroundModel(bg);
    delete bgIO;

  }
  else
  {
    cerr <<
      "Error: CreateBackgroundModel -> unable to create background model." <<
      endl;
  }

  // cleanup variables
  cleanup();

  // delete matrices
  delete[]pSnf;
  delete[]pFreq;
  delete[]pTrans;

  if (organism != NULL)
    delete organism;

  exit(0);
}


void
instructions()
{

  cout << endl;
  cout << "Usage:" << " CreateBackgroundModel <ARGS>" << endl;
  cout << endl;
  cout << " Required Arguments" << endl;
  cout << "  -f <fastaFile>      Sequences in FASTA format" << endl;
  cout <<
    "  -b <bgFile>         Output file in which the background model will be written."
    << endl;
  cout << "                      " << endl;
  cout << endl;
  cout << " Optional Arguments" << endl;
  cout << "  -o <value>          Order of the background model. Default = 1;"
    << endl;
  cout << "  -n <name>           Name of the organism." << endl;
  cout << endl;
  cout << "Version " << VERSION << endl;
  cout << "Questions and Remarks: <gert.thijs@esat.kuleuven.ac.be>" << endl <<
    endl;

  return;
}


void
version()
{
  cout << endl;
  cout << "INCLUSive -- CreateBackgroundModel (stand alone C++ version)" << endl;
  cout << "  Version " << VERSION << endl;
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

  return;

}
