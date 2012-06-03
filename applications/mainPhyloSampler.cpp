// 20 oct 2009 : make conform changes for 3.1.5 (MotifSamplerRun)
#include "inclusive.h"
#include "utilities.h"
#include "PWMIO.h"
#include "PWM.h"
#include "BackgroundIO.h"
#include "BackgroundModel.h"
#include "Instance.h"
#include "RandomNumber.h"
#include "MotifSamplerRun.h"

#include <iostream>
#include <fstream>
#include <getopt.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cstring>

#define OPTIONS "f:b:p:m:o:w:n:s:x:r:M:v"

// function prototypes
void instructions();
void version();
void cleanup();

string *fastaFile = NULL,
  *bgFile = NULL,
  *outputFile = NULL,
  *matrixFile = NULL;


int
main(int argc, char *argv[])
{

  // initialize random number generator
  RandomNumber::InitRandomNumber();
  srandom((long)time(NULL));
  cerr << "Seed = " << random() << endl;
  
  // define a set of local variables
  bool bFasta = false,
    bBG = false,
    bOut = false,
    bMatrix = false,
    bMax = false,
    bBgComputed = false,
    bBgSet = false;
  char c;
  double priorValue = 0.5;
  strand_modes bStrand = both;
  int wLength = 8,
    overlap = 1,
    nbrMotifs = 1,
    runs = 100,
    maxInstances = 10; //
  BackgroundModel *pBgModel = NULL;
  char pTmpChar[128];

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
    case 's':
      if (atoi(optarg) == 0) // error '!='
      {
        // default single strands
        bStrand = plus_strand;
      }
      else
      {
        bStrand = both;
      }
      break;
    case 'p':
      priorValue = (double) atof(optarg);
      if (priorValue <= 0 || priorValue >= 1)
      {
        cerr <<
          "Warning\n PhyloSampler: prior should be between 0 and 1 -> reset to default value 0.5"
          << endl;
        priorValue = 0.5;
      }
      break;
    case 'w':
      wLength = atoi(optarg);
      if (wLength <= 1)
      {
        cerr <<
          "Warning\n PhyloSampler: motif length should be greater than 1 => set to default 8"
          << endl;
        wLength = 8;
      }
      break;
    case 'n':
      nbrMotifs = atoi(optarg);
      if (nbrMotifs < 1)
      {
        cerr <<
          "Warning\n PhyloSampler: number of different motifs should be greater than 1."
          << endl;
        nbrMotifs = 1;
      }
      break;
    case 'o':

      // define new output stream 
      bOut = true;
      outputFile = new string(optarg);
      break;
      break;
    case 'm':
      matrixFile = new string(optarg);
      bMatrix = true;
      break;
    case 'r':
      runs = atoi(optarg);
      if (runs < 1)
      {
        cerr <<
          "Warning\n PhyloSampler: number of runs should be greater than 1. Default = 100."
          << endl;
        runs = 100;
      }
      break;
    case 'x':
      overlap = atoi(optarg);
      if (overlap < 0)
      {
        cerr <<
          "--Warning-- PhyloSampler: overlap should be greater than 0, use default value (=0)."
          << endl;
        overlap = 0;
      }
      break;
    case 'M':
      bMax = true;
      maxInstances = atoi(optarg);
      if (maxInstances < 1)
      {
        cerr <<
          "--Warning-- PhyloSampler: value added with -M switch is smaller than 1. Algorithm will not use this value."
          << endl;
        bMax = false;
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
      cerr << "-- PhyloSampler: Error in getopt() function" << endl;
      exit(-1);
    }
  }
  if (!bFasta || !bBG)
  {
    instructions();
    cleanup();
    exit(-1);
  }
	
//=========================================================
  // process input from -p (conform code changes for MotifSampler 3.1.5)
  // keep it with prior/kappa=0.25, maxInstances is set to 10 

  // local vector to store distributions
  vector<Distribution *> *_pPriorDistrs = new vector<Distribution *>;
  _pPriorDistrs->push_back(NULL); // do not use vector[0]
  ScoreVector * pProbs = new ScoreVector;
  double value, sum(1);
  pProbs->push_back(1 - priorValue); // Pr(0)
  pProbs->push_back(priorValue); // Pr(1)
  for(int i = 2; i <= maxInstances; i++) 
  { value = 0.25 * (*pProbs)[i-1];
    sum += value;
    if ( value/sum <= 0.001) break; // untill max or <0.001
     pProbs->push_back(value);
  }
  maxInstances = (int)pProbs->size() -1; // adjust -M if so
  for (int n = 1; n <= maxInstances; n++)
   _pPriorDistrs->push_back(new Distribution(pProbs, 0, n+1));
  // cleanup
  delete pProbs; pProbs = NULL;
//=========================================================
	

  // initialize file to write motif models
  PWMIO *pwmFile = NULL;
  if (bMatrix)
  {
    pwmFile = new PWMIO(matrixFile, WRITE);
  }

  // define output stream to write GFF results
  GFFWriter *pGffIO;
  if (bOut)
  {
    pGffIO = new GFFWriter(outputFile);
  }
  else
  {
    pGffIO = new GFFWriter();
  }

  // create a new MotifSampler object
  cerr << "Create PhyloSampler run." << endl;
  MotifSamplerRun *pMainRunner = new MotifSamplerRun(fastaFile, bStrand);

  // set primary parameters of the algorithm
  pMainRunner->SetMotifLength(wLength);
  //pMainRunner->SetOneInstancePrior(priorValue);
  pMainRunner->LinkNbrInstInfo(maxInstances, _pPriorDistrs, 0); //
  pMainRunner->SetOverlap(overlap);

  // check the number of sequences
  if (pMainRunner->NumberOfSequences() < 2)
  {
    cerr <<
      "Error\n Too few sequences in fasta file, there should be at least 2 sequences."
      << endl;
    delete pMainRunner;
    cleanup();
    exit(-1);
  }

  // create map to store background models
  map<string,BackgroundModel*> pBgList;
  map<string,BackgroundModel*>::iterator iBgList;
  
  // read the list of background models
  cerr << "Load the background models." << endl;
  ifstream _ifs(bgFile->c_str(), ios::in);
  while (!_ifs.eof())
  {
    // variables to store the names found in the bgFile
    string inLine, seqID, fileName;
    string::size_type pos = 0;
    string spaces(" \t");

    // read input line
    getline(_ifs, inLine, '\n');

    // split into id and bgFile
    pos = inLine.find_first_of(spaces);
    if ( pos != string::npos )
    {
      seqID = inLine.substr(0,pos);
      fileName = inLine.substr(pos+1,inLine.size()-pos-1);
      cerr << "DEBUG: Read pair: " << seqID << "|" << fileName << "|" << endl;
          
      // open background file for reading
      BackgroundIO *bgIO = new BackgroundIO(fileName, READ);
      if (bgIO != NULL)
      {
        pBgModel = bgIO->ReadBackgroundModel();
        if ( pBgModel != NULL ){
          // cerr << "map :" << seqID << endl;
          pBgList.insert(map<string,BackgroundModel*>::value_type(seqID,pBgModel));
          
          // update background score
          bBgComputed = pMainRunner->UpdateBackgroundScore(seqID, pBgModel);
          if ( !bBgComputed )
          {
            cerr << "--Error-- Problems scoring sequence "
              << seqID << " with bgModel: "
              << fileName << endl;
            delete pMainRunner;
            cleanup();
            exit(-1);
          }
          if ( !bBgSet )
          {
            pMainRunner->SetBackgroundModel(pBgModel);
            bBgSet = true;
          }
        }
        delete bgIO;
      }
      else
      {
        cerr << "--Error-- Problems opening background model file " << fileName
          << " for reading." << endl;
        delete pMainRunner;
        cleanup();
        exit(-1);
      }
    }
    else
    { 
      cerr << "next line" << endl;
    }
  }
  _ifs.close();


  // define string to keep id of motif
  char cid[64];
  string *pSource = new string("PhyloSampler");
  for (int r = 1; r <= runs; r++)
  {

    // reset mask to all 1's
    cerr << "Reset Masks." << endl;
    if (r != 1)
      pMainRunner->ResetMasks();
    for (int motifCounter = 1; motifCounter <= nbrMotifs; motifCounter++)
    {
      cerr << endl << "---- run: " << r << endl;
      cerr << "Searching for motif number: " << r << "." << motifCounter <<
        endl << endl;

      // initialize motif instance from random distribution
      cerr << "Initialization step" << endl;
      pMainRunner->InitializationStep(20);

      // core motif sampler iterations
      cerr << endl;
      cerr << "Starting core sampling steps" << endl;
      //if (!bMax)
      //{
        pMainRunner->CoreSamplingStep(70, 15, 3);
      //}
      //else
      //{
        //pMainRunner->CoreMaxSizeSamplingStep(maxInstances, 70, 15, 3);
      //}
      if (pMainRunner->NumberOfInstances() <= 1
          && pMainRunner->NumberOfSequencesWithInstances() <= 1)
      {
        cerr << "Warning: 0 or 1 instance found => stopping iteration" <<
          endl;
        break;
      }

      // 10 iterations in the convergence step
      cerr << endl;
      cerr << "Starting convergence step" << endl;
      //if (!bMax)
      //{
        pMainRunner->ConvergenceStep(10);
      //}
      //else
      //{
        //pMainRunner->MaxSizeConvergenceStep(maxInstances, 10);
      //}
      if (pMainRunner->NumberOfInstances() > 1
          && pMainRunner->NumberOfSequencesWithInstances() > 1)
      {
        PWM *myMatrix = pMainRunner->GetMotifModel();

        // define name of the motif model
        sprintf(cid, "box_%d_%d_%s", r, motifCounter,
                myMatrix->GetConsensus()->c_str());
        string *motifID = new string(cid);
        pMainRunner->SetInstanceID(motifID);

        // add comment line to GFF
        string *pComment = new string("#id: ");
        pComment->append(*motifID);
        pComment->append("\t");
        pComment->append("consensus: ");
        pComment->append(*(myMatrix->GetConsensus()));
        pComment->append("\t");
        pComment->append("sequences: ");

        // convert number
        sprintf(pTmpChar, "%d",
                pMainRunner->NumberOfSequencesWithInstances());
        pComment->append(pTmpChar);
        pComment->append("\t", 1);
        pComment->append("instances: ", 11);

        // convert number
        sprintf(pTmpChar, "%d", pMainRunner->NumberOfInstances());
        pComment->append(pTmpChar);
        pComment->append("\t", 1);
        pComment->append("cs: ", 4);

        // convert number
        sprintf(pTmpChar, "%1.2f", myMatrix->ConsensusScore());
        pComment->append(pTmpChar);
        pComment->append("\t", 1);
        pComment->append("ic: ", 4);

        // convert number
        sprintf(pTmpChar, "%1.2f",
                myMatrix->InformationContent(pBgModel->GetSNF()));
        pComment->append(pTmpChar);
        pComment->append("\t", 1);
        pComment->append("ll: ", 4);

        // convert number
        sprintf(pTmpChar, "%5.2f", pMainRunner->LogLikelihoodScore());
        pComment->append(pTmpChar);

        // write comment string on GFF output
        pGffIO->AddComment(pComment);

        // write the instance map
        pMainRunner->PrintInstanceMap(pGffIO, pSource);

        // clear the comment string
        delete pComment;
        pComment = NULL;

        // write motif to matrix file
        if (bMatrix)
        {

          // add additional values to matrix model
          myMatrix->SetScore(pMainRunner->LogLikelihoodScore());
          myMatrix->SetID(motifID);
          pwmFile->WriteMatrix(myMatrix);
        }
        myMatrix = NULL;
        delete motifID;
      }

      // apply mask to current instance map
      // not necessarry if last motif found
      if (motifCounter != nbrMotifs)
        pMainRunner->UpdateMasksFromInstanceMap();
    }
  }

  // clear background list
  iBgList = pBgList.begin();
  while ( iBgList != pBgList.end() )
  {
    //cerr << "DEBUG: deleting bgModel: " << (*iBgList).first << endl;
    //delete &((*iBgList).first);
    delete (*iBgList).second;
    iBgList++;
  }
  
  // cleanup local variables
  cleanup();
  if (pwmFile != NULL)
    delete pwmFile;
  if (pGffIO != NULL)
    delete pGffIO;
  pGffIO = NULL;
	
  // cleanup _pPriorDistrs vector
  if (_pPriorDistrs != NULL)
  { // cleanup distributions
    for(int i = 0; i < (int)_pPriorDistrs->size(); i++)
    { if ( (*_pPriorDistrs)[i] != NULL)
        delete (*_pPriorDistrs)[i];
	}
    delete _pPriorDistrs;
  }
  _pPriorDistrs = NULL;

  // delete motifID;
  delete pSource;
  delete pMainRunner;
  //delete pBgModel;

  // end of the program
  exit(0);
}

void
instructions()
{
  cout << endl;
  cout << "Usage:" << " PhyloSampler <ARGS>" << endl;
  cout << endl;
  cout << " Required Arguments" << endl;
  cout << "  -f <fastaFile>      Sequences in FASTA format" << endl;
  cout <<
    "  -b <bgFile>         File containing a list of sequence ids and background model file names."
    << endl;
  cout << "                      " << endl;
  cout << endl;
  cout << " Optional Arguments" << endl;
  cout << "  -s <0|1>            Select strand. (default both)" << endl;
  cout <<
    "                      0 is only input sequences, 1 include reverse complement."
    << endl;
  cout <<
    "  -p <value>          Sets prior probability of 1 motif copy. (default 0.5)."
    << endl;
  cout << "  -M <value>          Maximal number of motif instances per sequence. (default 10)"
    << endl;
  cout <<
    "  -n <value>          Sets number of different motifs to search for (default 1)."
    << endl;
  cout << "     " << endl;
  cout << "  -w <value>          Sets length of the motif (default 8)." <<
    endl;
  cout <<
    "  -x <value>          Sets allowed overlap between different motifs. (default 1)"
    << endl;
  cout <<
    "  -r <runs>           Set number of times the PhyloSampler should be repeated"
    << endl;
  cout << "                      (default = 100)." << endl;
  cout << " Output formatting Arguments" << endl;
  cout <<
    "  -o <outFile>        Output file to write results (default stdout)." <<
    endl;
  cout <<
    "  -m <matrixFile>     Output file to write retrieved motif models." <<
    endl;
  cout << endl;
  cout << "Version " << VERSION << endl;
  cout << "Questions and Remarks: <gert.thijs@esat.kuleuven.ac.be>" << endl
    << endl;
} void

version()
{
  cout << endl;
  cout << "INCLUSive -- PhyloSampler (stand alone C++ version)" << endl;
  cout << "  Version " << VERSION << endl;
  cout << "Revision history : " << endl; 
  cout << "- 06/12/05 : 4 bugs fixed - see ./MotifSampler -v" << endl;
  cout << "- 18/10/06 : 1 bug fixed - see ./MotifSampler -v" << endl;
  cout << "- 02/01/08 : 1 bug fixed - see ./MotifSampler -v" << endl;
  cout << "- 23/10/2009 : no functional change, only code - see ./MotifSampler -v" << endl;
  cout << "Questions and Remarks: <gert.thijs@esat.kuleuven.ac.be>" << endl
       << "                       <mclaeys@qatar.net.qa>" << endl
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
  if (matrixFile != NULL)
    delete matrixFile;
  matrixFile = NULL;
  if (outputFile != NULL)
    delete outputFile;
  outputFile = NULL;
}
