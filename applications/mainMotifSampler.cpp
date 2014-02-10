
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

#define OPTIONS "f:b:p:q:Q:m:o:w:n:s:x:r:M:S:vh?"


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
  *priorInput = NULL; // prior input as string


int
main(int argc, char *argv[])
{
  // initialize random number generator
  RandomNumber::InitRandomNumber();
  bool bFasta = false,
    bBG = false,
    bOut = false,
    bMatrix = false;
    //bMax = false;
  bool bSampling = 1; // sampling versus estimation
  char c;
  // double priorValue = 0.5;
  string _pDefault = "0.9_0.25"; // default prior input
  bool _bPSP_s = 0; // default only PSP for updating step
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
  // to output all error messages (instead of cerr)
  ostringstream cerrstr; 
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
    case 'p':
      priorInput = new string(optarg);
      break;
    case 'q':
      pspFile = new string(optarg);
      break;
    case 'Q': // for test purposes
      _bPSP_s = (bool)atoi(optarg);
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
    delete pMainRunner;
    cleanup();
    exit(-1);
  }

  // set primary parameters of the algorithm
  pMainRunner->SetMotifLength(wLength);
  pMainRunner->SetOverlap(overlap);

  // set background model and compute background model scores
  pMainRunner->SetBackgroundModel(pBgModel);
  pMainRunner->UpdateBackgroundScores();

  // pMainRunner->SetOneInstancePrior(priorValue);
  string warning = ""; 
  if (priorInput == NULL)  { priorInput = new string(_pDefault);}
  if (!pMainRunner->LoadNbrInstInfo(maxInstances, priorInput, bSampling))
  { warning = "--ERROR: incorrect processing of -p (#inst/seq) information.";}
  else if (pspFile != NULL && !pMainRunner->LoadPspScores(pspFile, _bPSP_s))
  { warning = "--ERROR: incorrect processing of -q (PSP) information.";}
  if (warning != "")
  {  
    cerr << warning << endl;
    if (bOut)
    {  
      cerrstr << warning << endl;
      cerrstr << "#Check our MotifSuite webpage for correct input format. Or contact us, we will be happy to help." << endl;
      pGffIO->AddComment(cerrstr.str());
      cerrstr.flush(); // flush  
    }
    // cleanup local variables
    delete pGffIO; pGffIO = NULL;
    if (bMatrix) delete pwmFile; pwmFile = NULL;
    delete pMainRunner;
    cleanup();  exit(-1);
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
      pMainRunner->CoreSamplingStep(70, 15, 3);
      if (pMainRunner->NumberOfInstances() <= 1
          && pMainRunner->NumberOfSequencesWithInstances() <= 1)
      {
        cerr << "--Warning: 0 or 1 instance found => aborting iteration." << endl;
        break;
      }

      // 10 iterations in the convergence step
      //cerr << endl;
      //cerr << "Starting convergence step" << endl;
      pMainRunner->ConvergenceStep(10);
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
  cout << endl;
  cout << " Optional Arguments" << endl;
  cout << "  -s <0|1>            Select strand (default both)." << endl;
  cout << "                      0 is only input sequences, 1 includes reverse complement."<< endl;
  cout << "  -p <(u|b|e|f)values> Specifies the prior distribution on number of instances per sequence (default 0.9_0.25)."<< endl;
  cout << "  -q <pspFile>        File containing Position Specific Prior information."<< endl;
  cout << "  -Q <0|1>            To apply psp-correction in sampling step (1=yes) (default 0).\n";
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
  cout << "- (3.1.5) 07/04/2011 : minor revisions in output formatting" << endl;
  //cout << " - (3.1.5) 22/07/2011 : again very minor revisions in output formatting" << endl;
  cout << "- (3.1.5) 24/10/2011 : initial iterations = 19 instead of 20." << endl; 
  cout << "- (3.1.5) 14/03/2012 : output error messages to user file." << endl;
  cout << "- (3.2.0) 17/09/2012 : (test-phase) implement Position Specific Prior (-q/-Q)." << endl;
  cout << "- (3.2.1) 04/02/2013 : no changes for MotifSampler." << endl;
  cout << "- (3.2.1) 12/03/2013 : checkup/revision PSP code." << endl;
  cout << "- (3.2.1) 14/05/2013 : debug segm fault for too short sequences." << endl;
  cout << "- (3.2.2) 16/01/2014 : lean PSP implementation." << endl;
  cout << "- (3.2.2) 30/01/2014 : document -Q to be boolean." << endl;
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
  if (priorInput != NULL)
    delete priorInput;
  priorInput = NULL;
}
