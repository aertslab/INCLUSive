#include "inclusive.h"
#include "utilities.h"
#include "PWMIO.h"
#include "PWM.h"
#include "BackgroundIO.h"
#include "BackgroundModel.h"
#include "Instance.h"
#include "RandomNumber.h"
#include "PriorMotifSamplerRun.h"

#include <iostream>
#include <fstream>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cstring>

#define OPTIONS "f:b:l:p:m:o:s:x:r:t:v"
#define VERSION "3.0"

using namespace INCLUSIVE;

// function prototypes
void instructions();
void version();
void cleanup();

string *fastaFile = NULL,
  *bgFile = NULL,
  *outputFile = NULL,
  *matrixFile = NULL,
  *priorFile = NULL;
  

int
main(int argc, char *argv[])
{
  // initialize random number generator
  RandomNumber::InitRandomNumber();

  bool bFasta = false,
    bBG = false,
    bOut = false,
    bMatrix = false,
    bPrior = false,
    bFirstFlag = false;
  char c;
  double weight = 1,
    priorValue = 0.5,
    threshold = 0.2;
  strand_modes bStrand = both;
  int overlap = 0,
    runs = 1,
    wLength = 0;
  BackgroundModel *pBgModel = NULL;
  char pTmpChar[128];
  // define string to keep id of motif
  char cid[64];
  string *pSource = new string("MotifSampler");
  
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
      priorFile = new string(optarg);
      bPrior = true;
      break;
    case 's':
      if (atoi(optarg) != 0)
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
          "--Warning-- MotifSampler: prior should be between 0 and 1 -> reset to default value 0.5"
          << endl;
        priorValue = 0.5;
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
          "Warning\n MotifSampler: number of runs should be greater than 1. Default = 100."
          << endl;
        runs = 100;
      }
      break;
    case 't':
      threshold = (double)atof(optarg);
      break;
    case 'x':
      overlap = atoi(optarg);
      if (overlap < 0)
      {
        cerr <<
          "--Warning-- mainMotifSampler: overlap should be greater than 0, use default value (=0)."
          << endl;
        overlap = 0;
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
      cerr << "-- MotifSampler: Error in getopt() function" << endl;
      exit(-1);
    }
  }
  if (!bFasta || !bBG || !bPrior )
  {
    instructions();
    cleanup();
    exit(-1);
  }

  // initialize file to write motif models
  PWMIO* pwmFile = NULL;
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

  // load the background model
  cerr << "Load the background model." << endl;
  BackgroundIO *bgIO = new BackgroundIO(*bgFile, READ);
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


  // open file with prior matrix information
  PWMIO* motifFile = new PWMIO(priorFile, READ);
  PWM* pPriorMotif = NULL;
  if ( motifFile == NULL )
  {
    cerr << "--Error-- Unable to open file with prior motif model."
      << endl;
    cleanup();
    exit(-1);
  }  
  
  
  // create a list to store all found matrices
  // list<PWM*> *pPWMList = new list<PWM*>;
  
  // create a list to store all found instances
  // list<Instance*> *pInstanceList = new list<Instance*>;
  
  PriorMotifSamplerRun* pMainRunner = NULL; 
  // check matrices
  while (motifFile->IsOpen() && ( pPriorMotif = motifFile->ReadMatrix()) ){

    if ( pPriorMotif == NULL )
    {
      cerr << "--Error-- Unable to read/create prior motif model."
        << endl;
      cleanup();
      exit(-1);
    }  
  
    if ( !bFirstFlag )
    {
      // create a new MotifSampler object
      cerr << "Create new PriorMotifSampler run." << endl;
      pMainRunner = new PriorMotifSamplerRun(fastaFile, bStrand, pPriorMotif);
    
      // set primary parameters of the algorithm
      pMainRunner->SetOneInstancePrior(priorValue);
      pMainRunner->SetOverlap(overlap);
    
      // set the weighting factor of the pseudo counts
      weight = sqrt((double)pMainRunner->NumberOfSequences());
      pMainRunner->UpdateWeightingFactor(weight);
      
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
    
      // set background model and compute background model scores
      cerr << "Initialize background model scores." << endl;
      if (pBgModel != NULL)
      {
        pMainRunner->SetBackgroundModel(pBgModel);
        pMainRunner->UpdateBackgroundScores();
        pMainRunner->UpdatePseudoCounts(pBgModel->GetSNF());
      }
      else
      {
        cerr << "Error: Problems with the background models." << endl;
        delete pMainRunner;
        cleanup();
        exit(-1);
      }
      
      cerr << "\n\nSet new prior motif model in PriorMotifSamplerRun(): " 
        << *(pPriorMotif->GetID()) << endl;
      
      bFirstFlag = true;
    }
    else
    {
      // update motif model
      cerr << "\n\nSet new prior motif model in PriorMotifSamplerRun(): " 
        << *(pPriorMotif->GetID()) << endl;
      pMainRunner->SetNewPriorMotifModel(pPriorMotif);
    }

    for (int r = 1; r <= runs; r++)
    {
      
      // reset mask to all 1's
      cerr << endl 
        << "Start run: " << r << endl << endl 
        << "Reset Masks." << endl;

      wLength = pPriorMotif->Length();
      
      pMainRunner->ResetMasks();
      pMainRunner->ResetMixtureCoefficients();
      pMainRunner->StderrPrintMixtureCoefficients();
  
      // initialize motif instance from random distribution
      cerr << "Initialization step" << endl;
      pMainRunner->InitializationStep(20);
      pMainRunner->StderrPrintMixtureCoefficients();
  
      // core motif sampler iterations
      cerr << endl;
      cerr << "Starting core sampling steps" << endl;
        
      pMainRunner->CoreSamplingStep(30);      
      pMainRunner->StderrPrintMixtureCoefficients();
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
      pMainRunner->ConvergenceStep(10);
      pMainRunner->StderrPrintMixtureCoefficients();
  
      if (pMainRunner->NumberOfInstances() > 1
          && pMainRunner->NumberOfSequencesWithInstances() > 1)
      {
        // POST-PROCESSING
        // get the resulting matrix
        PWM *myMatrix = pMainRunner->GetMotifModel();
        double D1 = 0, D2 = 0, value = 0, value2 = 0;
        vector<double>* pMix = pMainRunner->MixtureCoefficients();
        vector<double>* pParam1 = new vector<double>(4,0);
        vector<double>* pParam2 = new vector<double>(4,0);
        vector<double>* pData = new vector<double>(4,0);
        string* str1 = new string("");
        string* str2 = new string("");
        string* str3 = new string("");
        
        for(int ii=0; ii<wLength; ii++)
        {
          for (uint jj=0; jj<4; jj++)
          {
            (*pParam1)[jj] = weight * pPriorMotif->GetValueAt(ii,jj);
            (*pParam2)[jj] = weight * pBgModel->GetSnfValueAt(jj);
            // (*pParam1)[jj] = (*pMix)[ii] * weight * pPriorMotif->GetValueAt(ii,jj);
            // (*pParam2)[jj] = (1 - (*pMix)[ii]) * weight * pBgModel->GetSnfValueAt(jj);
            (*pData)[jj] = myMatrix->GetValueAt(ii,jj);;
          }
          
          value = LogDirichlet(pParam1, pData);
          D1 += value;
          sprintf(pTmpChar,"%3.2f|",value);
          str1->append(pTmpChar);
          
          value = LogDirichlet(pParam2, pData);
          D2 += value;
          sprintf(pTmpChar,"%3.2f|",value);
          str2->append(pTmpChar);
          
          sprintf(pTmpChar,"%3.2f|",(*pMix)[ii]);
          str3->append(pTmpChar);
        }          
        
        value = myMatrix->WeightedKullbackLeiberDistance(pPriorMotif, pMix);
        cerr << endl << "Weigthed Kullback-Leiber distance motif models (p): " << value << endl;
         
        vector<double>* pMixRev = new vector<double>(wLength);
        for (int ii=0; ii < wLength; ii++ )
          (*pMixRev)[ii] = 1 - (*pMix)[ii];
        
        value2 = myMatrix->WeightedKullbackLeiberDistance(pPriorMotif, pMixRev);
        cerr <<  "Weigthed Kullback-Leiber distance motif models (1-p): " << value2 << endl;
       
        cerr << " prior motif model:    " << *str1 << "= " << D1 << endl;
        cerr << " SNF model:            " << *str2 << "= " << D2 << endl;
        cerr << " mixture coeff:        " << *str3 << endl;
        
        // sum of mixturecoefficients
        value2 = SumArray(pMix,0, wLength);
        
        // delete these local variables
        delete pParam1;
        delete pParam2;
        delete pData;
        delete str1;
        delete str2;
        delete str3;
        
        // if ( exp(D1 - D2) >= 1 )
        if ( value2 > wLength/2 )
        {
          // write results to file        
          // define name of the motif model
          sprintf(cid, "box_%s_%d", pPriorMotif->GetID()->c_str(), r);
          string *motifID = new string(cid);
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
    
            pComment->append("\t", 1);
            pComment->append("kl: ", 4);
            // convert number
            sprintf(pTmpChar, "%1.4f", value);
            pComment->append(pTmpChar);
            
            pComment->append("\t", 1);
            pComment->append("sum_weights: ", 4);
            // convert number
            sprintf(pTmpChar, "%2.4f", value2);
            pComment->append(pTmpChar);
            
            // write comment string on GFF output
            pGffIO->AddComment(pComment);
    
            // clear the comment string
            if ( pComment != NULL )
              delete pComment;
            pComment = NULL;
          }
          else
          {
            cerr << "--Warning-- mainMotifSampler(): Unable to create comment string. No extra line added to GGF." << endl;
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
            cerr << "Write motif model " << *(myMatrix->GetID()) << " to file." << endl;
            pwmFile->WriteMatrix(myMatrix);
          }

          if ( motifID != NULL )
            delete motifID;
          motifID = NULL;

        }
        
        // clear local variables
        myMatrix = NULL;
        
      }
    }
    
    if ( pPriorMotif != NULL )
      delete pPriorMotif;
    pPriorMotif = NULL;
    
  }

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
  cout << "  -b <bgFile>         File containing the background model description"
    << endl;
  cout << "  -l <priorMatrix>    File containing the prior motif model"
    << endl;
  cout << "                      " << endl;
  cout << endl;
  cout << " Optional Arguments" << endl;
  cout << "  -s <0|1>            Select strand. (default both)" << endl;
  cout <<
    "                      0 is only input sequences, 1 include reverse complement."
    << endl;
  cout <<
    "  -p <value>          Sets prior probability of 1 motif copy. (default 0.2)."
    << endl;
  cout <<
    "  -x <value>          Sets allowed overlap between different motifs. (default 1)"
    << endl;
  cout <<
    "  -r <runs>           Set number of times the MotifSampler should be repeated (default = 1)"
    << endl;
  cout << "                      (default = 1)." << endl;
  cout << 
    "  -t <threshold>      Set the threshold on the selected motif models (default = 0.2)." << endl; 
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
} 

void
version()
{
  cout << endl;
  cout << "INCLUSive -- PriorMotifSampler (stand alone C++ version)" << endl;
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
  if (matrixFile != NULL)
    delete matrixFile;
  matrixFile = NULL;
  if (outputFile != NULL)
    delete outputFile;
  outputFile = NULL;
}
