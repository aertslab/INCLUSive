// 25 march 2012 : adjust cerr reporting => cerrstr

#include "inclusive.h"
#include "BackgroundIO.h"
#include "BackgroundModel.h"
#include "FastaIO.h"
#include "SequenceObject.h"
#include "GFFWriter.h"

#include <iostream>
#include <fstream>
#include <getopt.h>
#include <time.h>
#include <cstring>
#include <stdlib.h>
#include <math.h>


#define OPTIONS "hf:b:o:n:v"

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

  // 25 march 2012
  ostringstream cerrstr; // to output all error messages (instead of cerr)
  bool wronginput = false;

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
      if (!(order == 0 || order == 1 || order == 2 || order == 3 || order == 4 || order == 5))
      {
        cerr << "--ERROR: -o order should be an integer positive number limited to 5." << endl;
        cerrstr << "--ERROR: -o order should be an integer positive number limited to 5." << endl;
        wronginput = true;
        //cleanup(); exit(-1);
      }
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
      cerr << "-- CreateBackgroundModel: Error in input getopt() function" << endl;
      cerrstr << "-- CreateBackgroundModel: Error in input getopt() function" << endl;
      wronginput = true;
      //cleanup(); exit(-1);
    }
  }

  // server will not proceed if no fasta file supplied
  // output cannot be texted if no bBg file is given
  // therefore no extra cerrstr reporting applicable here
  if (!bFasta || !bBG)
  {
    cerr << "--ERROR: required files (-f,-b) are not provided." << endl;
    instructions();
    wronginput = true;
    //cleanup(); exit(-1);
  }
		
  if (wronginput)
  { 
    // report to bgFile in this case (bBG is surely true here)
    GFFWriter * pGffIO = new GFFWriter(bgFile);
    cerrstr << "#Check our MotifSuite webpage for correct input format. Or contact us, we will be happy to help." << endl;
    pGffIO->AddComment(cerrstr.str());
    cerrstr.flush(); // flush  
    delete pGffIO; pGffIO = NULL;
    // exit
    cleanup(); exit(-1);  
  }

  // open fasta file 
  FastaIO *fileIO = new FastaIO(fastaFile);
  SequenceObject *pSeq = NULL;

  if (fileIO == NULL || !fileIO->IsOpen())
  {
    cerr << "--Error: CreateBackgroundModel: Unable to read sequences from fastaFile: "
      << *fastaFile << endl;
    // report error message to bgFile in this case
    GFFWriter * pGffIO = new GFFWriter(bgFile);
    cerrstr << "--Error: CreateBackgroundModel: Unable to read sequences from fastaFile: "
      << *fastaFile << endl;
    cerrstr << "#Check our MotifSuite webpage for correct input format. Or contact us, we will be happy to help." << endl;
    pGffIO->AddComment(cerrstr.str());
    cerrstr.flush(); // flush  
    delete pGffIO; pGffIO = NULL;
    // exit
    cleanup(); exit(-1);  
  }

  // define matrices
  int l = (int) pow(4.0, order);

  // create matrices to store data
  double *pSnf = new double[4];
  for (int i = 0; i < 4; i++)
    pSnf[i] = 0;

  if (order == 0)
    l = 4;

  double *pFreq = new double[l];
  for (int i = 0; i < l; i++)
    pFreq[i] = 0;

  double **pTrans = new double*[l];
  for (int i = 0; i < l; i++)
  {
    pTrans[i] = new double[4];
    for (int j = 0; j < 4; j++)
      pTrans[i][j] = 0;
  }

  int nt,
    index,
    nt1;

	cerr << "CreateBackgroundModel loading sequences: ";
  int count = 0;
  while (fileIO->HasNext())
  {
    pSeq = fileIO->NextSequence();
    count++;
    cerr << count << ")" << *(pSeq->GetID()) << " ";
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
						if ( order != 0 )
						{
							// add one to frequency matrix
							pFreq[index] += 1;

							// add 1 to right entry in transition matrix
							pTrans[index][nt] += 1;
						}
						else
						{
							pFreq[nt] += 1;
							pTrans[0][nt] += 1;
							pTrans[1][nt] += 1;
							pTrans[2][nt] += 1;
							pTrans[3][nt] += 1;
						}
					}	
        }
      }
    }

    // get last oligo of length order
		if ( order != 0 )
		{
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
		}
		
    //cerr << "done" << endl;
    // delete sequence
    delete pSeq;
    pSeq = NULL;
  }
  cerr << "done." << endl;

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
    if (bgIO->GetError() != NULL)
    {
      cerr << *(bgIO->GetError());
      // report error message to bgFile in this case
      GFFWriter * pGffIO = new GFFWriter(bgFile);
      if (pGffIO != NULL) 
      { 
        pGffIO->AddComment(bgIO->GetError());
        cerrstr.flush(); // flush  
        delete pGffIO; pGffIO = NULL;
      }
      delete bgIO; bgIO = NULL;
    }
    if (bgIO != NULL)
    {
      bgIO->WriteBackgroundModel(bg);
      delete bgIO;
    }
  }
  else
  {
    cerr << "--Error: CreateBackgroundModel -> unable to create background model." << endl;
    // report error message to bgFile in this case
    GFFWriter * pGffIO = new GFFWriter(bgFile);
    cerrstr << "--Error: CreateBackgroundModel -> unable to create background model." << endl;
    cerrstr << "#Check our MotifSuite webpage for correct input format. Or contact us, we will be happy to help." << endl;
    pGffIO->AddComment(cerrstr.str());
    cerrstr.flush(); // flush  
    delete pGffIO; pGffIO = NULL;
  }

  // cleanup variables
  cleanup();

  // delete matrices
  delete[] pSnf;
  delete[] pFreq;
  for (int i = 0; i < l; i++)
    delete[] pTrans[i];
  delete[] pTrans;

  if (organism != NULL)
    delete organism;
  cerr << "CreatBackgroundModel: ending procedure successfully." << endl;
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
  cout << "  -b <bgFile>         Output file in which the background model will be written." << endl;
  cout << "                      " << endl;
  cout << endl;
  cout << " Optional Arguments" << endl;
  cout << "  -o <value>          Order of the background model. Default = 1;" << endl;
  cout << "  -n <name>           Name of the organism." << endl;
  cout << endl;
  cout << "Version " << VERSION << endl;
  cout << "Questions and Remarks: <mclaeys@qatar.net.qa>" << endl << endl;
  return;
}


void
version()
{
  cout << endl;
  cout << "INCLUSive -- CreateBackgroundModel (stand alone C++ version)" << endl;
  cout << "  Version " << VERSION << endl;
  cout << endl;
  cout << "Revision history : " << endl;
  cout << "- (3.1.5) 14/03/2012 : output error messages to user file." << endl;
  cout << "- (3.2.1) 04/02/13 : no changes for CreateBackgroundModel." << endl;
  cout << endl; 
  return;
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
