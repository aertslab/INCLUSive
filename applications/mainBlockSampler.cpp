// 20 oct 2009 : make conform changes for 3.1.5 (MotifSamplerRun)
// 17 nov 2009 : check on reading all bg-models

#include "inclusive.h"
#include "utilities.h"
#include "PWMIO.h"
#include "PWM.h"
#include "BackgroundIO.h"
#include "BackgroundModel.h"
#include "Instance.h"
#include "RandomNumber.h"
#include "BlockSamplerRun.h"

#include <iostream>
#include <fstream>
#include <getopt.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cstring>

#define OPTIONS "f:b:p:m:o:w:s:i:r:t:v"

// function prototypes
void instructions();
void version();
void cleanup();

string *fastaFile = NULL,
  *bgFile = NULL,
  *outputFile = NULL,
  *matrixFile = NULL,
  *rootID = NULL;


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
    bBgComputed = false,
    bBgSet = false,
    bRoot = false;
  char c;
  double priorValue = 0.5;
  double T = 1.0,
    ll = 0;
  strand_modes bStrand = plus_strand;
  int wLength = 8,
    wExtended = 8,
    runs = 100,
    pos = 0;
  BackgroundModel *pBgModel = NULL;
  char pTmpChar[128];

  // initialize some pseudocounts
  double pPseudo[4] = { 0.05, 0.05, 0.05, 0.05};
  
  
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
    case 'i':
      rootID = new string(optarg);
      bRoot = true;
      break;
    case 's':
      if (atoi(optarg) != 0)
      {
        // default single strands
        bStrand = both;
      }
      else
      {
        bStrand = plus_strand;
      }
      break;
    case 'p':
      priorValue = (double) atof(optarg);
      if (priorValue <= 0 || priorValue >= 1)
      {
        cerr <<
          "Warning\n BlockSampler: prior should be between 0 and 1 -> reset to default value 0.5"
          << endl;
        priorValue = 0.5;
      }
      break;
    case 'w':
      wLength = atoi(optarg);
      if (wLength <= 1)
      {
        cerr <<
          "Warning\n BlockSampler: motif length should be greater than 1 => set to default 8"
          << endl;
        wLength = 8;
      }
      break;
    case 'o':
      // define new output stream 
      bOut = true;
      outputFile = new string(optarg);
      break;
    case 'm':
      matrixFile = new string(optarg);
      bMatrix = true;
      break;
    case 't':
      T = (double) atof(optarg);
		if ( T < 0 || T >= 2 )
      {
        cerr << "Warning threshold <T> out of range, using default value 1.0 instead." <<endl;
        T = 1.0;
      }
      break;
    case 'r':
      runs = atoi(optarg);
      if (runs < 1)
      {
        cerr <<
          "Warning\n BlockSampler: number of runs should be greater than 1. Default = 100."
          << endl;
        runs = 100;
      }
      break;
    case '?':
    case 'h':
      instructions();
      exit(0);
    case 'v':
      version();
      exit(0);
    default:
      cerr << "-- BlockSampler: Error in getopt() function" << endl;
      exit(-1);
    }
  }

  if (!bFasta || !bBG || !bRoot )
  {
    instructions();
    cleanup();
    exit(-1);
  }
	
//=========================================================
  // process input from -p (conform code changes for MotifSampler 3.1.5)
  // keep it with prior/kappa=0.25
  // maxInstances is 1 (BlockSamplerRun is used like this)

  // local vector to store distributions, use only vector[1]
  vector<Distribution *> *_pPriorDistrs = new vector<Distribution *>; 
  ScoreVector * pProbs = new ScoreVector;
  pProbs->push_back(1 - priorValue); // Pr(0)
  pProbs->push_back(priorValue); // Pr(1)
  _pPriorDistrs->push_back(NULL); // do not use vector[0]
  _pPriorDistrs->push_back(new Distribution(pProbs, 0, 2));
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
  cerr << "Create BlockSampler run." << endl;
  BlockSamplerRun *pMainRunner = new BlockSamplerRun(fastaFile, rootID, bStrand);

  // set primary parameters of the algorithm
  pMainRunner->SetMotifLength(wLength);
  //pMainRunner->SetOneInstancePrior(priorValue);
	
  pMainRunner->LinkNbrInstInfo(1, _pPriorDistrs, 0); // max=1, estimation
  // overwrite to case fixed for root-sequence (always select one instance)
  pMainRunner->FixNbrInstForRoot();//

  pMainRunner->SetOverlap(0);
  pMainRunner->UpdatePseudoCounts(pPseudo);
  
  // check the number of sequences
  if (pMainRunner->NumberOfSequences() < 2)
  {
    cerr << "Error\n Too few sequences in fasta file, there should be at least 2 sequences." << endl;
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
      // cerr << "bgmodel opened " << endl;
      if (bgIO != NULL)
      {
        pBgModel = bgIO->ReadBackgroundModel();
        // cerr << "bgmodel read " << endl;
        if ( pBgModel != NULL ){
          // add model to map
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

  //~ iBgList = pBgList.begin();
  //~ cerr << "DEBUG: background model list" << endl;
  //~ while ( iBgList != pBgList.end() )
  //~ {
    //~ cerr << "    " << (*iBgList).first << "  " << *(((*iBgList).second)->GetOrganism()) << endl;
    //~ iBgList++;
  //~ }
  //~ cerr << endl;
  
  
  // define string to keep id of motif
  char cid[256];
  string *pSource = new string("BlockSampler");

  // start the runs 
  for (int r = 1; r <= runs; r++)
  {
    
    // reset mask to all 1's
    cerr << "Reset Masks." << endl;
    if (r != 1){
      // reset motif length
      pMainRunner->SetMotifLength(wLength);
      pMainRunner->ResetMasks();
    }
    
    cerr << endl << "---- run: " << r << endl << endl;

    // reset motif length
    pMainRunner->SetMotifLength(wLength);
    iBgList = pBgList.begin();
    while ( iBgList != pBgList.end() )
    {
      // cerr << "DEBUG " << *(pBgModel->GetOrganism()) << endl;
      pMainRunner->UpdateBackgroundScore(((*iBgList).first),(*iBgList).second);
      iBgList++;
    }
    
    // initialize motif instance from random distribution
    cerr << "Initialization step" << endl;
    pMainRunner->InitializationStep(15);

    pMainRunner->StderrPrintInstanceMap();
      
    // let us extend the motif model and perform two consecutive core sampling steps
    wExtended = wLength;
    for ( int iii=0; iii<3; iii++ )
    {
      if ( pMainRunner->NumberOfInstances() < 1 || pMainRunner->NumberOfSequencesWithInstances() < 2 )
        break;
      
      // extend sites with extra conserved positions
      cerr << "Check if motif can be extended to left and right hand site." << endl;
      pos = pMainRunner->CheckLeftExtension(T,5);
      if ( pos ){
        if ( pos > 20 )
          pos = 20;
        pMainRunner->ExtendLeftInstanceMap(pos);
        // pMainRunner->StderrPrintInstanceMap();
      }
      wExtended += pos; // add number of added positions to motif length
      pos = pMainRunner->CheckRightExtension(T,5);
      if ( pos ){
        if ( pos > 20 )
          pos = 20;
        pMainRunner->ExtendRightInstanceMap(pos);
        // pMainRunner->StderrPrintInstanceMap();
      }
      wExtended += pos; // add number of added positions to motif length
  
      pMainRunner->StderrPrintInstanceMap();
      
      // change motif length
      pMainRunner->SetMotifLength(wExtended);
      iBgList = pBgList.begin();
      while ( iBgList != pBgList.end() )
      {
        pMainRunner->UpdateBackgroundScore(((*iBgList).first),(*iBgList).second);
        iBgList++;
      }
      pMainRunner->BuildMotifFromInstanceMap();
      pMainRunner->UpdateMotifScores();
  
      // core motif sampler iterations
      cerr << endl;
      cerr << "Starting core sampling steps" << endl;
      pMainRunner->CoreSamplingStep(10);
      // pMainRunner->StderrPrintInstanceMap();
    }
  
    if (pMainRunner->NumberOfInstances() > 1
        && pMainRunner->NumberOfSequencesWithInstances() > 1)
    {

      // 10 iterations in the convergence step
      cerr << endl;
      cerr << "Starting convergence step" << endl;
      pMainRunner->ConvergenceStep(5);
    
      // report results
      if (pMainRunner->NumberOfInstances() > 1
          && pMainRunner->NumberOfSequencesWithInstances() > 1)
      {
  
        cerr << "Preparing results for writing" << endl;        
        PWM *myMatrix = pMainRunner->GetMotifModel();
  
        // define name of the motif model
        sprintf(cid, "block_%s_%d", rootID->c_str(), r);
        string *motifID = new string(cid);
        pMainRunner->SetInstanceID(motifID);
  
        //cerr << "Motif model: " << *motifID << endl;
        
        // add comment line to GFF
        string *pComment = new string("#id: ");
        pComment->append(*motifID);
        pComment->append("\t");
        pComment->append("consensus: ");
        pComment->append(*(myMatrix->GetConsensus()));
        pComment->append("\t");
        pComment->append("sequences: ");
        sprintf(pTmpChar, "%d", pMainRunner->NumberOfSequencesWithInstances());
        pComment->append(pTmpChar);
        pComment->append("\t", 1);
        pComment->append("instances: ", 11);
        sprintf(pTmpChar, "%d", pMainRunner->NumberOfInstances());
        pComment->append(pTmpChar);
        pComment->append("\t", 1);
        pComment->append("cs: ", 4);
        sprintf(pTmpChar, "%1.2f", myMatrix->ConsensusScore());
        pComment->append(pTmpChar);
        pComment->append("\t", 1);
        pComment->append("ic: ", 4);
        sprintf(pTmpChar, "%1.2f", myMatrix->InformationContent(pBgModel->GetSNF()));
        pComment->append(pTmpChar);
        pComment->append("\t", 1);
        pComment->append("ll: ", 4);
        ll = pMainRunner->LogLikelihoodScore();
        sprintf(pTmpChar, "%5.2f", ll);
        pComment->append(pTmpChar);     
  
        // write comment string on GFF output
        pGffIO->AddComment(pComment);
		  
		//cerr << "Motif modelresults :  " << *pComment << endl;
  
        // cerr << "DEBUG: " << *pComment << endl;
        
        // write the instance map
        pMainRunner->PrintInstanceMap(pGffIO, pSource);

        // write motif to matrix file
        if (bMatrix)
        {  
          // add additional values to matrix model
          myMatrix->SetScore(pMainRunner->LogLikelihoodScore());
          myMatrix->SetID(motifID);
          pwmFile->WriteMatrix(myMatrix);  
        }

        // clear the comment string
        if ( pComment !=  NULL )
          delete pComment;
        pComment = NULL;
  
        if ( myMatrix != NULL )
          myMatrix = NULL;
        delete motifID;
      }
      else
      {
        cerr << "Warning: 0 or 1 instance found => stopping iteration" << endl;      
      }
    }
    else
    {
      cerr << "Warning: 0 or 1 instance found => stopping iteration" << endl;      
    }
  }

  cleanup();

  // cleanup local variables
  if (pwmFile != NULL)
    delete pwmFile;
  if (pGffIO != NULL)
    delete pGffIO;
  pGffIO = NULL;
	
  iBgList = pBgList.begin();
  while ( iBgList != pBgList.end() )
  {
    delete (*iBgList).second;
    // cerr << "DEBUG: deleting bgModel: " << (*iBgList).first << endl;
	iBgList++;
  }
  // cerr << "DEBUG: deleting other variables" << endl;
  // delete motifID;
	if ( pSource != NULL )
	  delete pSource;
	pSource = NULL;
	if ( pMainRunner != NULL ) 
	  delete pMainRunner;
	pMainRunner = NULL;

	// debug : already deleted in list.
	//if ( pBgModel != NULL )
	  //delete pBgModel;
	//pBgModel = NULL;

  // cleanup _pPriorDistrs vector
  if (_pPriorDistrs != NULL)
  { if ( (*_pPriorDistrs)[1] != NULL)
      delete (*_pPriorDistrs)[1];
    delete _pPriorDistrs;
  }
  _pPriorDistrs = NULL;

  // end of the program
  exit(0);
}

void
instructions()
{
  cout << endl;
  cout << "Usage:" << " BlockSampler <ARGS>" << endl;
  cout << endl;
  cout << " Required Arguments" << endl;
  cout << "  -f <fastaFile>      Sequences in FASTA format" << endl;
  cout << "  -b <bgFile>         File containing a list of sequence ids and background model file names."
    << endl;
  cout << "  -i <value>          Defines the root sequence of the data set. <value> should be"
    << "                      similar to the identifier of the sequence in the fasta file."
    << endl;
  cout << "                      " << endl;
  cout << endl;
  cout << " Optional Arguments" << endl;
  cout << "  -s <0|1>            Select strand. (default plus strand)" << endl;
  cout <<
    "                      0 is only input sequences, 1 include reverse complement."
    << endl;
  cout <<
    "  -p <value>          Sets prior probability of 1 motif copy. (default 0.5)."
    << endl;
  cout << "     " << endl;
  cout << "  -w <value>          Sets length of the motif (default 8)." <<
    endl;
  cout <<
    "  -r <runs>           Set number of times the BlockSampler should be repeated"
    << endl;
  cout << "                      (default = 100)." << endl;
  cout << "  -t <value>          Sets consensus score threshold to extend motif length (default = 1.0)." << endl;
  cout << " Output formatting Arguments" << endl;
  cout <<
    "  -o <outFile>        Output file to write results (default stdout)." <<
    endl;
  cout <<
    "  -m <matrixFile>     Output file to write retrieved motif models." <<
    endl;
  cout << endl;
  cout << "Version " << VERSION << endl;
  cout << "Questions and Remarks: <mclaeys@qatar.net.qa>" << endl
    << endl;
} void

version()
{
  cout << endl;
  cout << "INCLUSive -- BlockSampler (stand alone C++ version)" << endl;
  cout << "  Version " << VERSION << endl;
  cout << "Revision history : " << endl; 
  cout << "- 06/12/05 : 4 bugs fixed - see ./MotifSampler -v" << endl;
  cout << "- 18/10/06 : 1 bug fixed - see ./MotifSampler -v" << endl;
  cout << "- 02/01/08 : 1 bug fixed - see ./MotifSampler -v" << endl;
  cout << "- 23/10/2009 : no functional change, only code - see ./MotifSampler -v" << endl;
  cout << "- 17/11/2009 : no functional change, some extra checkpoints" << endl;
  cout << "Questions and Remarks: <mclaeys@qatar.net.qa>" << endl
       << endl;
} void

cleanup()
{
  if (fastaFile != NULL)
    delete fastaFile;
  fastaFile = NULL;
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
