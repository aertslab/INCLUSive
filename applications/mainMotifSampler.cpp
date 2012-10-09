// 23 october 2009 - 3.1.5
// include extensions on # instances/sequence (see 3.1.3)
// report extra output to tracking file (-t)
// 20 march 2012 : adjust cerr reporting => cerrstr
// 19 Sept 2012 : implement Position Specific Prior (3.2.0), option -q/-Q 

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
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cstring>
#include <sstream> //
#include <climits>//
#include <ctime>//

#define OPTIONS "f:b:p:q:Q:m:o:t:w:n:s:x:r:M:S:vh?"


// function prototypes
void instructions();
void version();
void cleanup();

// define some global variables to store strings
string *fastaFile = NULL,
  *bgFile = NULL,
  *pspFile = NULL,
  *outputFile = NULL,
  *matrixFile = NULL,
  *trackFile = NULL,
  *priorInput = NULL; // prior input as string


int
main(int argc, char *argv[])
{
  // initialize random number generator
  RandomNumber::InitRandomNumber();
  bool bFasta = false,
    bBG = false,
	bPSP = false,
    bOut = false,
    bMatrix = false;
    //bMax = false;
  bool bSampling = 1; // sampling versus estimation
  char c;
  // double priorValue = 0.5;
  string _pDefault = "0.9_0.25"; // default prior input
  string _pPSPimpact = "n"; // apply prior in update step
  strand_modes bStrand = both;
  int wLength = 8,
    overlap = 1, // 
    nbrMotifs = 1,
    runs = 100,//
    maxInstances = 2; //
  BackgroundModel *pBgModel = NULL;
  char pTmpChar[128];
  // define string to keep id of motif
  char cid[64];

  // 22 march 2012
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
    case 's':
      // default both strands
      bStrand = both;
      if ( atoi(optarg) == 0 )
        bStrand = plus_strand;
      break;
    /*case 'p':
      priorValue = (double) atof(optarg);
      if (priorValue <= 0 || priorValue >= 1)
      {
        cerr << "ERROR: prior should be between 0 and 1" << endl;
        exit(-1);
      }
      break;
	  */
    case 'p':
      priorInput = new string(optarg);
      break;
    case 'q':
      pspFile = new string(optarg);
      bPSP = true;
      break;
    case 'Q': // set PSP correction impact
      _pPSPimpact = string(optarg);
      if ( _pPSPimpact != "s" && _pPSPimpact != "b" && _pPSPimpact != "u" && _pPSPimpact != "n")
      {
        cerr << "ERROR: -Q input should be b|u|s|n." << endl;
        cerrstr << "ERROR: -Q input should be b|u|s|n." << endl;
        wronginput = true;
        //cleanup(); exit(-1);
      }
      break;
    case 'w':
      wLength = atoi(optarg);
      if (wLength <= 1)
      {
        cerr << "ERROR: -w motif length should be higher than 0." << endl;
        cerrstr << "ERROR: -w motif length should be higher than 0." << endl;
        wronginput = true;
        //cleanup(); exit(-1);
      }
      break;
    case 'n':
      nbrMotifs = atoi(optarg);
      if (nbrMotifs < 1)
      {
        cerr << "ERROR: -n number of different motifs should be higher than 0." << endl;
        cerrstr << "ERROR: -n number of different motifs should be higher than 0." << endl;
        wronginput = true;
        //cleanup(); exit(-1);
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
      trackFile = new string(optarg);
      break;
    case 'r':
      runs = atoi(optarg);
      if (runs < 1)
      {
        cerr << "ERROR: -r number of runs should be higher than 0." << endl;
        cerrstr << "ERROR: -r number of runs should be higher than 0." << endl;
        wronginput = true;
        //cleanup(); exit(-1);
      }
      break;
    case 'x':
      overlap = atoi(optarg);
      if (overlap < 0)
      {
        cerr << "ERROR: -x overlap should be greater than 0." << endl;
        cerrstr << "ERROR: -x overlap should be greater than 0." << endl;
        wronginput = true;
        //cleanup(); exit(-1);
      }
      break;
    case 'M':
      //bMax = true;
      maxInstances = atoi(optarg);
      if ( maxInstances < 1 )
      {
        cerr << "ERROR: -M maximal number of instances per sequence should be higher than 0." << endl;
        cerrstr << "ERROR: -M maximal number of instances per sequence should be higher than 0." << endl;
        wronginput = true;
        //cleanup(); exit(-1);
        //bMax = false;
      }
      break;
    case 'S': // 
      bSampling = (bool)atoi(optarg);
      break;
    case '?':
    case 'h':
      instructions();
      exit(-1);
    case 'v':
      version();
      exit(-1);
    default:
      cerr << "-- MotifSampler: Error in input getopt() function" << endl;
      cerrstr << "--MotifSampler: Error in input getopt() function" << endl; 
      wronginput = true;
      //cleanup(); exit(-1);
    }
  }
  // check wrong input if so report/exit
  GFFWriter *pGffIO;
  if (!bFasta || !bBG)
  {
    cerr << "--ERROR: required inputfiles (-f, -b) are not provided." << endl;
    instructions();
    cerrstr << "--ERROR: required inputfiles (-f, -b) are not provided." << endl;
    // report to file if file is available
    if (bOut) 
    { cerrstr << "#Check our MotifSuite webpage for correct input format. Or contact us, we will be happy to help." << endl;
      pGffIO = new GFFWriter(outputFile); 
      pGffIO->AddComment(cerrstr.str());
      cerrstr.flush(); // flush  
      if (pGffIO != NULL) delete pGffIO; pGffIO = NULL; 
    }
    // exit
    cleanup(); exit(-1); 
  }
//=========================================================
  // process input from -p (mainly based on 5.0.4 / partially 3.1.3)
  if (priorInput == NULL){ priorInput = new string(_pDefault);}
  // local vector to store distributions 
  vector<Distribution *> *_pPriorDistrs = new vector<Distribution *>; 
  // (case f) make the fixed distribution---------------->(vector[0])
  // else make prior distributions for n = 1->M(max)----->(vector[1->M])
  ScoreVector * pProbs = NULL;
  string::size_type pos;
  if ((*priorInput)[0] == 'f')
  {
    priorInput->erase(0,1); // erase 'f'
    pProbs = new ScoreVector;
    double p;
    while (!priorInput->empty() && (int)pProbs->size() < maxInstances+1)
    { 
      pos = priorInput->find_first_not_of("0.123456789");
      istringstream istr(priorInput->substr(0,pos)); // convert to double
      istr >> p;
      if (p < 0 || p > 1) // invalid prob
        break;
      priorInput->erase(0,pos);
      while ( priorInput->find_first_of("_ \t") == 0) priorInput->erase(0,1); 
      pProbs->push_back(p);
    }
    // cut low values at the end
    for(int i = (int)pProbs->size()-1; i >= 0; i--) 
    { if ( (*pProbs)[i] <= 0.001) pProbs->pop_back();
      else break;
    }
    // create fixed distribution
    maxInstances = (int)pProbs->size() - 1; // adjust -M if so
    Distribution * dist = new Distribution(pProbs, 0, (int)pProbs->size());
    _pPriorDistrs->push_back(dist); // may be NULL - see end check
  }
  else 
  { // store NULL on vector[0]
    _pPriorDistrs->push_back(NULL);
    // create next distributions
    c = (*priorInput)[0]; // 0,u,b,e
    switch (c)
    {
      case '0' :
      { // read prior p
        double p = 0.5; // prior SEE DEFAULT IN _pDefault
        pos = priorInput->find_first_not_of("0.123456789");
        istringstream istr(priorInput->substr(0,pos)); // convert to double
        istr >> p;
        if (p <= 0 || p >= 1) // invalid prior
          break;
        priorInput->erase(0,pos);
        while ( priorInput->find_first_of("_ \t") == 0) priorInput->erase(0,1); 
        // read kappa 
        double kappa = 0.25; // kappa SEE DEFAULT IN _pDefault
        if ( !priorInput->empty()) 
        { 
          pos = priorInput->find_first_not_of("0.123456789");
          istringstream istr(priorInput->substr(0,pos)); // convert to double
          istr >> kappa;
          if (kappa < 0 || kappa > 1) // invalid kappa
            break;
          priorInput->erase(0,pos);
          while ( priorInput->find_first_of("_ \t") == 0) priorInput->erase(0,1); 
        }
        if (!priorInput->empty()) 
          break; // invalid values or format < >_< >
        // construct distributions
        pProbs = new ScoreVector;
        double value, sum(1);
        pProbs->push_back(1 - p); // Pr(0)
        pProbs->push_back(p); // Pr(1)
        for(int i = 2; i <= maxInstances; i++) 
        { value = kappa * (*pProbs)[i-1];
          sum += value;
          if ( value/sum <= 0.001) break; // untill max or <0.001
          pProbs->push_back(value);
        }
        maxInstances = (int)pProbs->size() -1; // adjust -M if so
        for (int n = 1; n <= maxInstances; n++)
          _pPriorDistrs->push_back(new Distribution(pProbs, 0, n+1));
        break;
      }
      case 'u' : 
      {
        priorInput->erase(0,1); // erase 'u'
        if ( !priorInput->empty()) break;
        // create and store the distributions
        pProbs = new ScoreVector; pProbs->push_back(1);
        for (int n = 1; n <= maxInstances; n++) 
        { pProbs->push_back(1);
          _pPriorDistrs->push_back(new Distribution(pProbs, 0, n+1));
        }
        break;
      }
      case 'b' : 
      {
        priorInput->erase(0,1); // erase 'b'
        // read input
        double p;
        pos = priorInput->find_first_not_of("0.123456789");
        istringstream istr(priorInput->substr(0,pos)); // convert to double
        istr >> p;
        if (p <= 0 || p >= 1) // invalid prior
          break;
        priorInput->erase(0,pos);
        if ( !priorInput->empty()) break;
        // create and store the distributions
        for (int n = 1; n <= maxInstances; n++)
        { pProbs = new ScoreVector; 
          for(int i = 0; i <= n; i++)
            pProbs->push_back(pow(p,i)*pow((1-p),(n-i))*
                   INCLUSIVE::fac(n)/INCLUSIVE::fac(n-i)/INCLUSIVE::fac(i));
          for(int i = n; i >= 0; i--) // cut low values at the end
          { if ( (*pProbs)[i] <= 0.001) pProbs->pop_back();
            else break;
          }
          if ( (int)pProbs->size() == n + 1 )
          { _pPriorDistrs->push_back(new Distribution(pProbs, 0, n+1));
            delete pProbs; pProbs = NULL; 
          }
          else 
          { maxInstances = (int)pProbs->size()-1; // adjust -M if so
            break;
          }
        }
        break;
      }
      case 'e' : 
      {
        priorInput->erase(0,1); // erase 'e'
        pProbs = new ScoreVector;
        double p;
        while (!priorInput->empty() && (int)pProbs->size() < maxInstances+1)
        { 
          pos = priorInput->find_first_not_of("0.123456789");
          istringstream istr(priorInput->substr(0,pos)); // convert to double
          istr >> p;
          if (p < 0 || p >= 1) // invalid prob
            break;
          priorInput->erase(0,pos);
          while ( priorInput->find_first_of("_ \t") == 0) priorInput->erase(0,1); 
          pProbs->push_back(p);
        }
        // cut low values at the end
        for(int i = (int)pProbs->size()-1; i >= 0; i--) 
        { if ( (*pProbs)[i] <= 0.001) { pProbs->pop_back();}
          else break;
        }
        // create distributions
        maxInstances = (int)pProbs->size() - 1; // adjust -M if so
        bool notNull = false;
        for (int n = 1; n <= maxInstances; n++)
        { 
          Distribution * dist = new Distribution(pProbs, 0, n+1);
          if (dist != NULL) notNull = true; // some may be NULL when all-0-subset
          _pPriorDistrs->push_back(dist); 
        }
        if (!notNull) // all distrs are NULL
        { cerr << "ERROR too much zero values in -p e<values> input" << endl;
          for( int i = (int)_pPriorDistrs->size()-1; i >= 1; i--)
          {
            if ((*_pPriorDistrs)[i] != NULL) delete (*_pPriorDistrs)[i]; 
            (*_pPriorDistrs)[i] = NULL;
            _pPriorDistrs->pop_back();
          }
        }
        break;
      }
    }
  }
  // cleanup
  if (pProbs != NULL) delete pProbs; pProbs = NULL;
  // check on correct input
  if (_pPriorDistrs->size() < 2  && (*_pPriorDistrs)[0] == NULL)
  {
    cerr << "ERROR: incorrect format or invalid values for -p. Check guidelines." << endl;
    delete _pPriorDistrs; _pPriorDistrs = NULL;
    instructions();
    cerrstr << "--ERROR: incorrect format or invalid values for -p. Check guidelines." << endl;
    wronginput = true;
    //cleanup(); exit(-1);
  }
//=========================================================

  if (wronginput)
  { 
    // report to file if file is available
    if (bOut) 
    { cerrstr << "#Check our MotifSuite webpage for correct input format. Or contact us, we will be happy to help." << endl;
      pGffIO = new GFFWriter(outputFile); 
      pGffIO->AddComment(cerrstr.str());
      cerrstr.flush(); // flush  
      if (pGffIO != NULL) delete pGffIO; pGffIO = NULL;
    }
    // exit
    cleanup(); exit(-1);  
  }

  cerr << "MotifSampler: Load the background model." << endl;
  // load the background model
  BackgroundIO *bgIO = new BackgroundIO(*bgFile, READ);
  if (bgIO->GetError() != NULL)
  {
    cerr << *(bgIO->GetError());
    if (bOut) { cerrstr << *(bgIO->GetError());}
    if (bgIO != NULL) delete bgIO; bgIO = NULL;
    wronginput = 1;
  }
  else
  { 
    pBgModel = bgIO->ReadBackgroundModel();
    delete bgIO; bgIO = NULL;
    if (pBgModel == NULL)
    { cerr << "--ERROR: problems with background model." << endl;
      if (bOut) cerrstr << "--ERROR: problems with background model." << endl;
      wronginput = 1;
    }
  }
  if (wronginput)
  { 
    // report to file if file is available
    if (bOut) 
    { cerrstr << "#Check our MotifSuite webpage for correct input format. Or contact us, we will be happy to help." << endl;
      pGffIO = new GFFWriter(outputFile); 
      pGffIO->AddComment(cerrstr.str());
      cerrstr.flush(); // flush  
      if (pGffIO != NULL) delete pGffIO; pGffIO = NULL;
    }
    // exit
    cleanup(); exit(-1);  
  }

  // initialize file to write motif models
  PWMIO *pwmFile = NULL;
  if (bMatrix) // this means it is intended to have PWM reporting
  {
    pwmFile = new PWMIO(matrixFile, WRITE);
    if (pwmFile->GetError() != NULL)
    {
      cerr << *(pwmFile->GetError());
      if (bOut)
      {
        pGffIO = new GFFWriter(outputFile); 
        pGffIO->AddComment(pwmFile->GetError()); 
        delete pGffIO; pGffIO = NULL;  
      }
      delete pwmFile; pwmFile = NULL;
      cleanup(); exit(-1);  
    }
  }
  // define output stream to write GFF results // always ok
  if (bOut)
  {
    pGffIO = new GFFWriter(outputFile);
  }
  else
  {
    pGffIO = new GFFWriter();
  }
  GFFWriter *pGffTrack = NULL;
  if ( trackFile != NULL) pGffTrack = new GFFWriter(trackFile);

  // reset some settings depending on the input
  string warning = "";
  if (_pPSPimpact == "n")
  {  if (bPSP == true)
       // give warning that PSP file will not be used
       warning = "WARNING: the supplied PSP file (-q) will not be used as -Q was set to 'n' (=apply no PSP correction).\n";
     bPSP = false; // reset
  }
  else 
  { if (bPSP == false)
      // give warning that no PSP will be used 
      // (i.e. default neutral PSP, all probs are set to '1')
      warning = "WARNING: no PSP correction (yet -Q differs from 'n') will be applied as no PSP file (-q) is provided.\n";
  }
  if (warning != "")
  {
    cerr << warning;
    if (bOut)
    {
      cerrstr << warning;
      pGffIO->AddComment(cerrstr.str());
      cerrstr.flush(); // flush
    }
  }

  int errorSampling = 0;
  int errorMasking = 0;
  int retrievedMotifs = 0;
  int minSeq(INT_MAX), maxSeq(0), minInst(INT_MAX), maxInst(0),minLL(INT_MAX), maxLL(0);
  double minCS(2.1), maxCS(0), minIC(2.1), maxIC(0), score;

  // create a new MotifSampler object
  cerr << "Create MotifSampler run." << endl;
  MotifSamplerRun *pMainRunner = new MotifSamplerRun(fastaFile, bStrand, pGffIO);
  // check the number of sequences
  if (pMainRunner->NumberOfSequences() < 2)
  {
    cerr << "--ERROR: too few sequences in fasta file, there should be at least 2 sequences."
      << endl;
    if (bOut)
    {  
      cerrstr << "--ERROR: too few sequences in fasta file, there should be at least 2 sequences."
      << endl;
      cerrstr << "#Check our MotifSuite webpage for correct input format. Or contact us, we will be happy to help." << endl;
      pGffIO->AddComment(cerrstr.str());
      cerrstr.flush(); // flush  
  
    }
    // cleanup local variables
    delete pGffIO; pGffIO = NULL;
    if (bMatrix) delete pwmFile; pwmFile = NULL;
    if ( trackFile != NULL) delete trackFile; trackFile = NULL;
    if (_pPriorDistrs != NULL)
    { // cleanup distributions
      for(int i = 0; i < (int)_pPriorDistrs->size(); i++)
      { if ( (*_pPriorDistrs)[i] != NULL)
          delete (*_pPriorDistrs)[i];
      }
      delete _pPriorDistrs;
    }
    _pPriorDistrs = NULL;
    delete pMainRunner;
    cleanup();
    exit(-1);
  }
  // set primary parameters of the algorithm
  pMainRunner->SetMotifLength(wLength);
  // pMainRunner->SetOneInstancePrior(priorValue);
  pMainRunner->LinkNbrInstInfo(maxInstances, _pPriorDistrs, bSampling); // 
  //cerr << "debug-bSampling =  " << bSampling << endl;
  pMainRunner->SetOverlap(overlap);
  // pass trackfile to MotifSamplerRun
  pMainRunner->LinkFiles(pGffTrack);
  // set background model and compute background model scores
  pMainRunner->SetBackgroundModel(pBgModel);
  pMainRunner->UpdateBackgroundScores();
  // set the PSP information
  if (bPSP) 
  { if (!(pMainRunner->UpdatePspScores(pspFile, _pPSPimpact)))
    { // exit as incorrect PSP input
      cerr << "--ERROR: incorrect processing of PSP information."<< endl;
      if (bOut)
      {  
        cerrstr << "--ERROR: incorrect processing of PSP information."<< endl;
        cerrstr << "#Check our MotifSuite webpage for correct input format. Or contact us, we will be happy to help." << endl;
        pGffIO->AddComment(cerrstr.str());
        cerrstr.flush(); // flush  
      }
      // cleanup local variables
      delete pGffIO; pGffIO = NULL;
      if (bMatrix) delete pwmFile; pwmFile = NULL;
      if ( trackFile != NULL) delete trackFile; trackFile = NULL;
      if (_pPriorDistrs != NULL)
      { // cleanup distributions
        for(int i = 0; i < (int)_pPriorDistrs->size(); i++)
        { if ( (*_pPriorDistrs)[i] != NULL)
            delete (*_pPriorDistrs)[i];
        }
        delete _pPriorDistrs;
      }
      _pPriorDistrs = NULL;
      delete pMainRunner;
      cleanup();
      exit(-1);
    }
  }

  // 22 juli 2011
  string * intro = new string("#see end of file for a summary of results.\n#The multiple motifs in this file are the result of stochastic motif detection repeats on a given sequences dataset.\n#Please proceed with extraction of the most likely true motif(s) embedded in this file by using MotifRanking (quick assessment) or FuzzyClustering (more refined assessment). Both post processing tools are available in MotifSuite on our webserver on which you can also consult our guidelines.\n");
  pGffIO->AddComment(intro);
  delete intro; intro = NULL;

  string *pSource = new string("MotifSampler");
    time_t t1; time(&t1);
  cerr << "Searching for " << runs << " x " << nbrMotifs << " motifs..." << endl;

  for (int r = 1; r <= runs; r++)
  {

    // reset mask to all 1's
    //cerr << "Reset Masks." << endl;
    if (r != 1)
      pMainRunner->ResetMasks();
    for (int motifCounter = 1; motifCounter <= nbrMotifs; motifCounter++)
    {
      //cerr << endl << "---- run: " << r << endl;
      //cerr << "Searching for motif number: " << r << "." << motifCounter << endl << endl;

      // initialize motif instance from random distribution
      //cerr << "Initialization step" << endl;
      pMainRunner->InitializationStep(19);

      // core motif sampler iterations
      //cerr << endl;
      //cerr << "Starting core sampling steps" << endl;
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
        cerr << "--Warning: 0 or 1 instance found => aborting iteration." << endl;
        /*if (bOut)
        {
          cerrstr << "--Warning: 0 or 1 instance found => aborting iteration." << endl;
          pGffIO->AddComment(cerrstr.str());
          cerrstr.flush(); // flush
        } */
        break;
      }

      // 10 iterations in the convergence step
      //cerr << endl;
      //cerr << "Starting convergence step" << endl;
      //if (!bMax)
      //{
        pMainRunner->ConvergenceStep(10);
      //}
      //else
      //{
       // pMainRunner->MaxSizeConvergenceStep(maxInstances, 10);
      //}
      if (!(pMainRunner->NumberOfInstances() > 1
          && pMainRunner->NumberOfSequencesWithInstances() > 1))
      {
        if (motifCounter > 1) 
        {errorMasking++; motifCounter = nbrMotifs;} // skip further search
        else {errorSampling++;}
      }
      else
      {  
        retrievedMotifs++;
        PWM *myMatrix = pMainRunner->GetMotifModel();

        // define name of the motif model
        // cerr << "DEBUG: create id of motif" << endl;
        sprintf(cid, "box_%d_%d_%s", r, motifCounter,
                myMatrix->GetConsensus()->c_str());
        string *motifID = new string(cid);

        // cerr << "DEBUG: store id of motif in Run" << endl;
        pMainRunner->SetInstanceID(motifID);

        // add comment line to GFF
        // cerr << "DEBUG: make comment string" << endl;
        string *pComment = new string("#id: ");
        if ( pComment != NULL )
        {
          pComment->append(*motifID);
          // cerr << "DEBUG: " << *pComment << endl;

          pComment->append("\t");
          pComment->append("consensus: ");
          pComment->append(*(myMatrix->GetConsensus()));
          // cerr << "DEBUG: " << *pComment << endl;

          pComment->append("\t");
          pComment->append("sequences: ");
          // convert number
          sprintf(pTmpChar, "%d",
                  pMainRunner->NumberOfSequencesWithInstances());
          pComment->append(pTmpChar);
          // cerr << "DEBUG: " << *pComment << endl;

          pComment->append("\t", 1);
          pComment->append("instances: ", 11);
          // convert number
          sprintf(pTmpChar, "%d", pMainRunner->NumberOfInstances());
          pComment->append(pTmpChar);
          // cerr << "DEBUG: " << *pComment << endl;
			
score = (int)pMainRunner->NumberOfSequencesWithInstances();
if (score <= minSeq) { minSeq = score;}
if (score >= maxSeq) { maxSeq = score;}
score = (int)pMainRunner->NumberOfInstances();
if (score <= minInst) { minInst = score;}
if (score >= maxInst) { maxInst = score;}
		
          
          pComment->append("\t", 1);
          pComment->append("cs: ", 4);
          // convert number
          sprintf(pTmpChar, "%1.2f", myMatrix->ConsensusScore());
          pComment->append(pTmpChar);
          // cerr << "DEBUG: " << *pComment << endl;
          pComment->append("\t", 1);
          pComment->append("ic: ", 4);
          // convert number
          sprintf(pTmpChar, "%1.2f",
                  myMatrix->InformationContent(pBgModel->GetSNF()));
          pComment->append(pTmpChar);
          // cerr << "DEBUG: " << *pComment << endl;

          pComment->append("\t", 1);
          pComment->append("ll: ", 4);
          // convert number
          sprintf(pTmpChar, "%5.2f", pMainRunner->LogLikelihoodScore());
          pComment->append(pTmpChar);
          // cerr << "DEBUG: " << *pComment << endl;

score = myMatrix->ConsensusScore();
if (score <= minCS) { minCS = score;}
if (score >= maxCS) { maxCS = score;}
score = myMatrix->InformationContent(pBgModel->GetSNF());
if (score <= minIC) { minIC = score;}
if (score >= maxIC) { maxIC = score;}
score = (int)pMainRunner->LogLikelihoodScore();
if (score <= minLL) { minLL = score;}
if (score >= maxLL) { maxLL = score;}

  
          // write comment string on GFF output
          // cerr << "DEBUG: add comment string to gff writer" << endl;
          pGffIO->AddComment(pComment);

          // clear the comment string
          if ( pComment != NULL )
            delete pComment;
          pComment = NULL;
        }
        else
        {
          cerr << "--Warning-- mainMotifSampler(): Unable to create comment string. No extra line added to GGF." << endl;
          if (bOut)
          {
            cerrstr << "--Warning-- mainMotifSampler(): Unable to create comment string. No extra line added to GGF." << endl;
            pGffIO->AddComment(cerrstr.str());
            cerrstr.flush(); // flush 
          }
        }
        
        // write the instance map
        // cerr << "DEBUG: write instance map" << endl;
        pMainRunner->PrintInstanceMap(pGffIO, pSource);

        // write motif to matrix file
        if (bMatrix)
        {
          // add additional values to matrix model
          myMatrix->SetScore(pMainRunner->LogLikelihoodScore());
          myMatrix->SetID(motifID);
          pwmFile->WriteMatrix(myMatrix);
        }
        myMatrix = NULL;

        if ( motifID != NULL )
          delete motifID;
        motifID = NULL;
        
        cerr << "End run motif_" << r << "_" << motifCounter << "." << endl;

      }

      // apply mask to current instance map
      // not necessarry if last motif found
      if (motifCounter != nbrMotifs)
        pMainRunner->UpdateMasksFromInstanceMap();
    }
  }

string *pSummary = new string("\n#MotifSampler results summary :\n");
pSummary->append("#number of reported motifs :");
// convert number
sprintf(pTmpChar, "%d",retrievedMotifs);
pSummary->append(pTmpChar);
pSummary->append("\n");
pSummary->append("#convergence rate (%):");
// convert number
sprintf(pTmpChar, "%d",100*retrievedMotifs/runs/nbrMotifs);
pSummary->append(pTmpChar);
pSummary->append("\n");
pSummary->append("#maximal/minimal number of sequences with a motif :");
// convert number
sprintf(pTmpChar, "%d",maxSeq);
pSummary->append(pTmpChar);
pSummary->append("/");
// convert number
sprintf(pTmpChar, "%d",minSeq);
pSummary->append(pTmpChar);
pSummary->append("\n");
pSummary->append("#maximal/minimal number of instances in a motif :");
// convert number
sprintf(pTmpChar, "%d",maxInst);
pSummary->append(pTmpChar);
pSummary->append("/");
// convert number
sprintf(pTmpChar, "%d",minInst);
pSummary->append(pTmpChar);
pSummary->append("\n");
pSummary->append("#maximal/minimal motif Consensus score :");
// convert number
sprintf(pTmpChar, "%1.2f",maxCS);
pSummary->append(pTmpChar);
pSummary->append("/");
// convert number
sprintf(pTmpChar, "%1.2f",minCS);
pSummary->append(pTmpChar);
pSummary->append("\n");
pSummary->append("#maximal/minimal motif Information Content score :");
// convert number
sprintf(pTmpChar, "%1.2f",maxIC);
pSummary->append(pTmpChar);
pSummary->append("/");
// convert number
sprintf(pTmpChar, "%1.2f",minIC);
pSummary->append(pTmpChar);
pSummary->append("\n");
pSummary->append("#maximal/minimal motif LogLikelihood score :");
// convert number
sprintf(pTmpChar, "%d",maxLL);
pSummary->append(pTmpChar);
pSummary->append("/");
// convert number
sprintf(pTmpChar, "%d",minLL);
pSummary->append(pTmpChar);
pSummary->append("\n");
if (errorSampling>= 10 && errorSampling >= (int) 0.2 * (runs-errorMasking)*nbrMotifs)
{ pSummary->append("# ");
  // convert number
  sprintf(pTmpChar, "%d", errorSampling);
  pSummary->append(pTmpChar);
  pSummary->append(" predictions were rejected because of too few retrieved instances. You may want to repeat the program with a higher probability for number of instances per sequence (-p).\n");
}
if (errorMasking>= 10 && errorMasking >= (int)0.2 * runs)
{ pSummary->append("# ");
  // convert number
  sprintf(pTmpChar, "%d", errorMasking);
  pSummary->append(pTmpChar);
  pSummary->append(" runs were aborted for searching a next motif most probably because of too much sequence masking by earlier retrieved motifs in this run. Alternatively, you may want to repeat the program with a lower -n (number of motifs) yet higher -r (number of runs) and explore more of the top ranked motifs (use MotifRanking).\n");
}


    // ending
    time_t t2; time(&t2);
    // add time to outputFile
pSummary->append("#MotifSampler runtime information : ");
    //out_3 << "# Date :"  << ctime(&t1);
    double T = difftime(t2, t1);
  // convert number
  sprintf(pTmpChar, "%d", (int)T);
  pSummary->append(pTmpChar);
  pSummary->append("sec (=~ ");
  sprintf(pTmpChar, "%d", (int)T/3600);
  pSummary->append(pTmpChar);
  pSummary->append("hrs, ");
  sprintf(pTmpChar, "%d", (int)((T - ((int)T/3600)*3600)/60));
  pSummary->append(pTmpChar);
  pSummary->append("min).");

pGffIO->AddComment(pSummary);
delete pSummary; pSummary = NULL;

  // cleanup local variables
  cleanup();
  if (pwmFile != NULL)
    delete pwmFile;
  if (pGffIO != NULL)
    delete pGffIO;
  pGffIO = NULL;
  if (pGffTrack != NULL)
    delete pGffTrack;
  pGffTrack = NULL;
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
  delete pBgModel;
  cerr << "MotifSampler: ending procedure successfully." << endl;
  // end of the program
  exit(0);
}

void
instructions()
{
  cout << endl;
  cout << "Usage:" << " MotifSampler <ARGS>" << endl;
  cout << endl;
  cout << " Required Arguments" << endl;
  cout << "  -f <fastaFile>      Sequences in FASTA format" << endl;
  cout << "  -b <bgFile>         File containing the background model description."<< endl;
  cout << endl;
  cout << " Output Arguments" << endl;
  cout << "  -o <outFile>        Output file to write putative motifs in instances format (default stdout)." << endl;
  cout << "  -m <matrixFile>     Output file to write putative motifs in pwm format (default not used)." << endl;
  //cout << "  -t <trackFile>      Output file to track sampled positions (default not used)." << endl;
  cout << endl;
  cout << " Optional Arguments" << endl;
  cout << "  -s <0|1>            Select strand (default both)." << endl;
  cout << "                      0 is only input sequences, 1 includes reverse complement."<< endl;
  cout << "  -p <(u|b|e|f)values> Specifies the prior distribution on number of instances per sequence (default 0.9_0.25)."<< endl;
  cout << "  -q <pspFile>        File containing Position Specific Prior information."<< endl;
  cout << "  -Q <u|s|b|n>        Apply Position Specific Prior correction in (u)updating, (s)sampling, (b)both or (n)none steps.\n";
  cout << "  -M <value>          Sets the maximal number of instances per sequence (default 2)." << endl;
  cout << "  -n <value>          Sets number of different motifs to search for (default 1)."<< endl;
  cout << "  -w <value>          Sets the motif width (default 8)." << endl;
  cout << "  -x <value>          Sets allowed overlap between different motifs (default 1)." << endl;
  cout << "  -r <runs>           Set number of times the algorithm should be repeated (default 100)." << endl;
  cout << "  -S <0|1>            Sampling (1) versus estimation (0) of the number of instances per sequence to be searched for (default 1)." << endl;
  cout << endl;
  cout << "Version " << VERSION << endl;
  cout << "Questions and Remarks: <mclaeys@qatar.net.qa>" << endl;
  cout << endl;
} 

void
version()
{
  cout << endl;
  cout << "INCLUSive -- MotifSampler (stand alone C++ version)" << endl;
  cout << "  Version " << VERSION << endl;
  cout << endl;
  cout << "Revision history : " << endl; 
  cout << "- (3.1.2) 06/12/2005 : debug : Utilities/LogDirichlet()"
       << "/ComputeNInstancesProbability(), RandomNumber/GetUniform()."
       << endl;
  cout << "- (3.1.2) 18/10/2006 : debug : MotifSamplerRun/"
       << "LogLikelihoodScore()." << endl;
  cout << "- (3.1.2) 02/01/2008 : debug : MotifSamplerRun/"
       << "UpdateMasksFromInstanceMape()." << endl;
  cout << "- (3.1.5) 23/10/2009 : upgrade option -p (prior) and implement -S (sampling)." << endl;
//<< "    -t : extra output to trackfile." << endl;
  cout << "- (3.1.5) 07/04/2011 : minor revisions in output formatting" << endl;
  //cout << " - (3.1.5) 22/07/2011 : again very minor revisions in output formatting" << endl;
  cout << "- (3.1.5) 24/10/2011 : initial iterations = 19 instead of 20." << endl; 
  cout << "- (3.1.5) 14/03/2012 : output error messages to user file." << endl;
  cout << "- (3.2.0) 17/09/2012 : (test-phase) implement Position Specific Prior (-q/-Q)." << endl;
  cout << endl;
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
  if (pspFile != NULL)
    delete pspFile;
  pspFile = NULL;
  if (matrixFile != NULL)
    delete matrixFile;
  matrixFile = NULL;
  if (outputFile != NULL)
    delete outputFile;
  outputFile = NULL;
  if (trackFile != NULL)
    delete trackFile;
  trackFile = NULL;
  if (priorInput != NULL)
    delete priorInput;
  priorInput = NULL;
}
