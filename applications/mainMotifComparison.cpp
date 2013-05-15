// 21 july 2009 - 3.1.5
// inclusive.h FastaIO.h/cpp PWMIO.cpp PWM.h/cpp
// report extra output to tracking file (-t)
// 20 march 2012 : adjust cerr reporting => cerrstr
// 28 sept 2012 : error (and exit) for incorrect PWM format reading



#include "inclusive.h"
#include "utilities.h"
#include "PWMIO.h"
#include "PWM.h"
#include "BackgroundIO.h" //
#include "BackgroundModel.h" //
#include "FastaIO.h" //
#include "GFFWriter.h" //
#include <iomanip>

// c++ includes
#include <list>
#include <iostream>
#include <fstream>
#include <sstream> //

// c includes
#include <getopt.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstring>

#define OPTIONS "hd:m:t:o:s:x:vl:p:b:n:Y:Z:" //
#define PI 3.14159265358979 //

void instructions();
void version();
void cleanup();
string *matrixFile = NULL,
  *dbFile = NULL,
  *outputFile = NULL,
  *bgFile = NULL,//
  *priorFile = NULL;//

int
main(int argc, char *argv[])
{
  char c;
  bool bOut = false,
    bDB = false,
    bMatrix = false;
  double KLtreshold(0.4), Ptreshold(0.001); //
  double treshold = -1; //
  int shift = 1;
  int overlap = 6; // -x
  int ScoreType = 0; // -l 
  int nbrShuffles = 20; // -n
  bool avgA2A1 = false; // -Z (take avgA2A1 yes or no)
  bool bPrint = 1; // to print ## lines yes or no

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
      treshold = (double)atof(optarg);
      if ( treshold <=0)
      {
        cerr << "ERROR: -t threshold should be higher than 0." << endl;
        cerrstr << "ERROR: -t threshold should be higher than 0." << endl;
        wronginput = true;
        //cleanup(); exit(-1);
      }
      break;
    case 'o':
      // define new output stream 
      bOut = true;
      outputFile = new string(optarg);
      break;
    case 's':
      shift = atoi(optarg);
      if ( shift < 0 )
      {
        cerr << "ERROR: -s shift should be higher than 0." << endl;
        cerrstr << "ERROR: -s shift should be higher than 0." << endl;
        wronginput = true;
        //cleanup(); exit(-1);
      }
      break;
    case 'x': //
      overlap = atoi(optarg);
      if ( overlap <= 0 )
      {
        cerr << "ERROR: -x overlap should be higher than 0." << endl;
        cerrstr <<  "ERROR: -x overlap should be higher than 0." << endl;
        wronginput = true;
        //cleanup(); exit(-1);
      }
      break;
    case 'l': // KL=0, BLiC=1
      ScoreType = atoi(optarg);
      if ( ScoreType != 0 && ScoreType != 1)
      {
        cerr << "ERROR: -l score type should be 0(KL) or 1(BLiC)." << endl;
        cerrstr << "ERROR: -l score type should be 0(KL) or 1(BLiC)." << endl;
        wronginput = true;
        //cleanup(); exit(-1);
      }
      break;
    case 'p': // Dirichlet prior
      priorFile = new string(optarg);
      break;
    case 'b': // BackgroundModel
      bgFile = new string(optarg);
      break;
    case 'n': // nbr Shuffles
      nbrShuffles = atoi(optarg);
      if ( nbrShuffles <= 0 )
      {
        cerr << "ERROR: -n number of shuffles should be higher than 0." << endl;
        cerrstr << "ERROR: -n number of shuffles should be higher than 0." << endl;
        wronginput = true;
        //cleanup(); exit(-1);
      }
      break;				
    case 'Y': // to print ## lines yes or no
      bPrint = (bool)atoi(optarg);
      break;
				
    case 'Z': // TEST (BLiC avg (a2-a1) y/n)
//For KL, the score is averaged over the length of the overlapping region : score = score / (a2-a1);
//This is not defined for BLiC, but if you want to test if this would be better, set -Z = 1
      avgA2A1 = (bool)atoi(optarg);
      break;
    case 'v':
      version();
      exit(-1);
    case '?':
    case 'h':
      instructions();
      exit(-1);
    default:
      cerr << "MotifComparison: Error in input getopt() function" << endl;
      cerrstr << "MotifComparison: Error in input getopt() function" << endl; 
      wronginput = true;
      //cleanup(); exit(-1);
    }
  }

  // no cerrstr reporting, error is covered in website interface
  if (!bDB || !bMatrix)
  {
    instructions();
    cleanup();
    exit(-1);
  }

  GFFWriter *pgffText = NULL;
  // open file for output writing 
  if (bOut) pgffText = new GFFWriter(outputFile);
  else pgffText = new GFFWriter(); // screen

  if (wronginput)
  { 
    // report 
    cerrstr << "#Check our MotifSuite webpage for correct input format. Or contact us, we will be happy to help." << endl; 
    pgffText->AddComment(cerrstr.str());
    cerrstr.flush(); // flush  
    delete pgffText; pgffText = NULL;
    // exit
    cleanup(); exit(-1);  
  }

  // open input matrix file for reading
  PWMIO *pwmIO1 = new PWMIO(matrixFile, READ);
  if (pwmIO1->GetError() != NULL) // comes to the same as (!pwmIO->IsOpen())
  {
    pgffText->AddComment(pwmIO1->GetError()); 
    delete pgffText; pgffText = NULL;  
    delete pwmIO1; pwmIO1 = NULL;
    cleanup(); exit(-1);  
  }

  // Read database matrices
  cerr << "Loading the database matrices..." << endl;
  PWMIO * pwmIO = new PWMIO(dbFile, READ);
  // open matrix stream reader
  if (pwmIO->GetError() != NULL) // comes to the same as (!pwmIO->IsOpen())
  {
    pgffText->AddComment(pwmIO->GetError()); 
    delete pgffText; pgffText = NULL;  
    delete pwmIO; pwmIO = NULL;
    delete pwmIO1; pwmIO1 = NULL;
    cleanup(); exit(-1);  
  }
  PWM * myMatrix = NULL;
  list < PWM * > dbMatrix;
  list < PWM * >::iterator matIter;
  int count = 0;
  int similars = 0; int nonsimilars = 0;

  while (pwmIO->IsOpen())
  { myMatrix = pwmIO->ReadMatrix(); 
    if (myMatrix == NULL && (pwmIO->GetError()) != NULL) 
    { 
      cerr << *(pwmIO->GetError()) << endl;
      pgffText->AddComment(pwmIO->GetError()); 
      cerr << "--Error: MotifComparison : incorrect matrix format found in matrix database file." << endl;
      cerrstr << "--Error: MotifComparison : incorrect matrix format found in matrix database file." << endl;
      pgffText->AddComment(cerrstr.str());
      cerrstr.flush(); // flush
      delete pgffText; pgffText = NULL;  
      delete pwmIO; pwmIO = NULL;
      delete pwmIO1; pwmIO1 = NULL;
      cleanup(); exit(-1);
    }
    if (myMatrix != NULL) // if NULL: processing not-matrix lines
    {
      dbMatrix.push_back(myMatrix); 
      count++;
    }
  }
  delete pwmIO; pwmIO = NULL;  // close matrix reader
  if (count == 0)
  {
    cerrstr << "--Error: MotifComparison : no usable matrix found in database input file." << endl;
    pgffText->AddComment(cerrstr.str());
    cerrstr.flush(); // flush  
    delete pgffText; pgffText = NULL;  
    delete pwmIO; pwmIO = NULL;
    delete pwmIO1; pwmIO1 = NULL;
    delete dbFile; dbFile = NULL;
    // exit
    cleanup(); exit(-1); 
  }
  cerr << dbMatrix.size() << " database matrices loaded." << endl;

  // assign user-input treshold 
  if (treshold != -1) // 
  { if (ScoreType == 0) KLtreshold = treshold;
    else 
    { Ptreshold = treshold; }
  }
  // load the backgroundModel for BLiC
  double * pDirichlet(NULL), *pBgSNF(NULL);
  wronginput = 0;
  if (ScoreType == 1)
  {
    pBgSNF = new double[4];
    for(int i = 0; i < 4; i++) pBgSNF[i] = 0.25; // default -b
    if (bgFile != NULL)
    {
      BackgroundModel * pBgModel = NULL;
      cerr << "Loading the background model... ";
      BackgroundIO *bgIO = new BackgroundIO(*bgFile, READ);
      if (bgIO->GetError() != NULL) 
      { wronginput = 1; cerrstr << "--Warning: " << *(bgIO->GetError()) << " Default values will be used instead." << endl; }
      else
      { pBgModel = bgIO->ReadBackgroundModel();
        if (pBgModel == NULL)
        { wronginput = 1; cerrstr << "--Warning: background model not suitable. Default values will be used instead." << endl; }
        else
        { // store snf
          pBgSNF = new double[4];
          for (int i = 0; i < 4; i++) pBgSNF[i] = pBgModel->GetSnfValueAt(i);
          delete pBgModel; pBgModel = NULL;
        }
      }
      delete bgIO; bgIO = NULL;
      delete bgFile; bgFile = NULL;
    }
    // report to screen
    cerr << "#Background snf : [";
    for(int i = 0; i < 4; i++) { cerr << pBgSNF[i] << " ";} cerr << "]" << endl;

    pDirichlet = new double[4];
    for(int i = 0; i < 4; i++) pDirichlet[i] = 0.01;// default -p
    if (priorFile != NULL)
    {
      // load the Dirichlet prior
      FastaIO * pFasta = new FastaIO(priorFile);
      cerr << "Loading the Dirichlet prior... ";
      if (pFasta->IsOpen() && pFasta->HasNext())
      {
        pDirichlet = pFasta->ReadDirichlet(); // simple case
        if (pDirichlet == NULL)
        { wronginput = 1; cerrstr << "--Warning: Problems reading Dirichlet prior (-p). Default values will be used instead." << endl; }
      }
      delete pFasta; pFasta = NULL;
      delete priorFile; priorFile = NULL;
    }
    cerr << "#Dirichlet prior : [";
    for(int i = 0; i < 4; i++) { cerr << pDirichlet[i] << " ";} cerr << "]" << endl;

    // add warnings to text outputfile
    if (wronginput)
    { pgffText->AddComment(cerrstr.str());
      cerrstr.flush(); // flush
    }  
  }

  string * comment = new string("");
  // 07 april 2011
  if (ScoreType == 0)
  { *comment = "#Please read MotifComparison Guidelines on our webserver to optimally evaluate below results.\n\n#query:query_consensus\tmatch:match_consensus\tscore;(shiftm-d,shiftd-m)\n"; }
  else
  { *comment = "#Please read MotifComparison Guidelines on our webserver to optimally evaluate below results.\n\n#query:query_consensus\tmatch:match_consensus\tscore;(shiftm-d,shiftd-m);shift;overlap;strand\tp-value;common;(n:avg,stdev)\n"; }
  pgffText->AddComment(comment);
  *comment = "";

  cerr << "Start comparing matrices with matrices from database." << endl;    
  double score;
  double * result; // score/shift/overlap/strand
  double * resultID; // BLiC for matrix with itself
  int shiftM; // maximum shift between two matrices, respecting minumum overlap
  double pvalue(1), common(0);

  ostringstream report;

  int L1(0), L2(0); //
  double avg(0), stdev(0);
  // iterate over all matrices
  int counter = 0;
  int maximalS = 0;
  while (pwmIO1->IsOpen())
  {
    myMatrix = pwmIO1->ReadMatrix();
    if (myMatrix == NULL && (pwmIO1->GetError()) != NULL) 
    { 
      cerr << *(pwmIO1->GetError()) << endl;
      pgffText->AddComment(pwmIO1->GetError()); 
      cerr << "--Error: MotifComparison : incorrect matrix format found in matrix input file." << endl;
      cerrstr << "--Error: MotifComparison : incorrect matrix format found in matrix input file." << endl;
      pgffText->AddComment(cerrstr.str());
      cerrstr.flush(); // flush  
      // do not exit here, but break
      break;
      //delete pgffText; pgffText = NULL;  
      //delete pwmIO; pwmIO = NULL;
      //delete pwmIO1; pwmIO1 = NULL;
      //cleanup(); exit(-1);
    }
    if (myMatrix == NULL) {break;} // processing not-matrix lines, end of file
    // else this is a good matrix, continue
    counter++;
    maximalS = 0;
    matIter = dbMatrix.begin();
    while (matIter != dbMatrix.end())
    {
      // take most stringent value from shift/overlap 
      shiftM = ( (L1 = myMatrix->Length()) < (L2 = (*matIter)->Length()) ?
                    (L1-overlap) : (L2-overlap) );
      if (shiftM < 0) {shiftM = 0;} 
      else {if (shiftM > shift) {shiftM = shift;}}
      // intermediate reporting // see below
      if (maximalS < shiftM){maximalS = shiftM;}

      // compare matrices 
      if (ScoreType == 0) // MutualInformation (Kullback-Leiber)
      { score =
         (myMatrix->MutualInformation((*matIter), shiftM) +
         (*matIter)->MutualInformation(myMatrix, shiftM)) / 2;
        // check score against treshold
        if( score >= KLtreshold) {nonsimilars++;} else {similars++;}
        if( score >= KLtreshold && bPrint) report << "##"; // non-similars
        if (score < KLtreshold || (score >= KLtreshold && bPrint))
        { report << *(myMatrix->GetID()) << ":" 
            << *(myMatrix->GetConsensus()) << "\t"
            << *((*matIter)->GetID()) << ":" 
            << *((*matIter)->GetConsensus()) 
            << "\t" << score
		    // get the locally stored applied shift from PWM
		    << ";(" << myMatrix->GetMIshift() << "," 
				     << (*matIter)->GetMIshift() << ")" 
            << endl;
        }
      }
      else // BLiC
      {
        result = myMatrix->BLiC(*matIter, shiftM, pDirichlet, pBgSNF, avgA2A1);
        resultID = myMatrix->BLiC(myMatrix, shiftM, pDirichlet, pBgSNF, avgA2A1);
        score = result[0];
        if (1) // always compute the pvalue
        {
          PWM * pShuffle1 = myMatrix->NewShuffledMatrix();
          PWM * pShuffle2 = (*matIter)->NewShuffledMatrix();
          pShuffle1->SetWeight(myMatrix->GetWeight());
          pShuffle2->SetWeight((*matIter)->GetWeight());
          double * shuffle; double b_hat, m_hat; avg = 0; stdev = 0; 
          vector<double>* pScores = new vector<double>;
          common = 0; // track same shift/overlap/strand
          for(int n = 0; n < nbrShuffles; n++)
          { // compute score for 2 shuffled matrices
            shuffle = pShuffle1->BLiC(pShuffle2, shiftM, pDirichlet, pBgSNF, avgA2A1);
            pScores->push_back(shuffle[0]); avg += shuffle[0];
            if( result[1] == shuffle[1] && result[2] == shuffle[2] &&
                result[3] == shuffle[3]) common++;
            // shuffle again
            pShuffle1->ShuffleMatrix();
            pShuffle2->ShuffleMatrix();
            delete shuffle; shuffle = NULL;
          }
          avg /= nbrShuffles; // average BLiC-score
          common /= nbrShuffles;
          for(int n = 0; n < nbrShuffles; n++) 
            stdev += pow(avg-(*pScores)[n], 2);
          stdev = sqrt(stdev/nbrShuffles); // stdev BLiC-score
          // compute p-value
          b_hat =(stdev * sqrt(6.0))/PI;
          m_hat = avg - (0.5772 * b_hat);
          pvalue = 1-exp(-exp((m_hat-score)/b_hat));
          // cleanup
          delete pShuffle1; pShuffle1 = NULL;
          delete pShuffle2; pShuffle2 = NULL;
          delete pScores; pScores = NULL;
        }
        // output the results 
        if( pvalue >= Ptreshold) {nonsimilars++;} else {similars++;}
        if ( pvalue >= Ptreshold && bPrint) report << "##";
        if ( pvalue < Ptreshold || (pvalue >= Ptreshold && bPrint)) 
        {     
          report << *(myMatrix->GetID()) << ":" 
            << *(myMatrix->GetConsensus()) << "\t"
            << *((*matIter)->GetID()) << ":" 
            << *((*matIter)->GetConsensus())
            << "\t" << score << ";" << result[1] << ";" << result[2];
          if (result[3] == 1) report << ";fwd"; else report << ";rev";
          report << "\t" << setprecision(5) << pvalue << ";" << common 
                 << ";(n=" << nbrShuffles << ":" << avg << "," << stdev << ")";
          report << endl;
        }
        // cleanup
        delete result; result = NULL; //
        delete resultID; resultID = NULL; //
      }
      // next matrix
      matIter++;
    }
    // report to screen for this matrix
    cerr << counter << ") " << *(myMatrix->GetID()) << " <-> " << count 
           << " db-matrices, maximal-shift was " << maximalS << "\n"; // 
    delete myMatrix; myMatrix = NULL; //
  }

  // close matrix reader
  delete pwmIO1; pwmIO1 = NULL;
	
	
	report << endl << "#SUMMARY of " << count << "(DB-matrices)*" 
		<< counter << "(INPUT-matrices) = " << count*counter << " comparisons:" << endl;
	report << "#SIMILARS: " << similars << "." << endl;
	report << "#NON-SIMILARS: " << nonsimilars << "." << endl;
	
	
  // print results
  pgffText->AddComment(report.str());
  report.flush(); // flush

  // clear dbMatrix
  matIter = dbMatrix.begin();
  while (matIter != dbMatrix.end())
  {
    // cerr << "deleting matrix: " << *((*matIter)->GetID()) << endl;
    delete(*matIter);
    matIter++;
  }
  cleanup();
  delete[] pDirichlet; pDirichlet = NULL; //
  delete[] pBgSNF; pBgSNF = NULL; // 
  if (pgffText != NULL) delete pgffText; pgffText = NULL;
  if (pwmIO != NULL) delete pwmIO; pwmIO = NULL;
  if (pwmIO1 != NULL) delete pwmIO1; pwmIO1 = NULL;

  // program is succesful
  cerr << "MotifComparison: ending procedure successfully." << endl;
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
    << "all query motifs will be compared." 
    << endl;
  cout << "  -o <outFile>        Output file to write results to." << endl;
  cout << endl;
  cout << " Optional Arguments" << endl;
  cout << "  -t <value>          Treshold where two motifs are seen as being "
       << "the same. Default KL<0.4> pvalue<0.001>" << endl;
  cout << "  -s <value>          Maximal allowed shift between motif models. "
       << "Default <1>." << endl;
  cout << "  -x <value>          Minimum required overlap between motif models. "
	   << "Default <6>. " << endl;
  cout << "  -l <0|1>          Type of score, KL(0) or BLiC(1). Default <0>" 
	   << endl;
  cout << "  -p <priorFile>      File containing the Dirichlet prior counts. "
	   << "Default <1_1_1_1>." << endl;
  cout << "  -b <bgFile>         File containing the BackgroundModel description. "
       << "Default <order0 : 0.25_0.25_0.25_0.25>."<< endl;	
  cout << "  -n <value>          Number of shuffles for p-value computation. "
       << "Default<20>" << endl;
  cout << endl;
  cout << "  -v                  Version of MotifComparison" << endl;
  cout << "Version " << VERSION << endl;
  cout << "Questions and Remarks: <mclaeys@qatar.net.qa>" << endl;
  cout << endl;
} 

void
version()
{
  cout << endl;
  cout << "INCLUSive -- MotifComparison (stand alone C++ version)" << endl;
  cout << "  Version " << VERSION << endl;
  cout << endl; 
  cout << "  Revision history : " << endl; 
  cout << "  - (3.1.2) 14/08/07 : explicit x-parameter for minimum overlap." << endl;
  cout << "  - (3.1.4) 14/01/08 : debug error input x-parameter." << endl;
  cout << "  - (3.1.4) 08/06/09 : extended output (personal use)." << endl;
  cout << "  - (3.1.5) 21/07/09 : implementation BLiC-score." << endl;
  cout << "  - (3.1.5) 08/09/09 : some updates on BLiC." << endl;
  cout << "  - (3.1.5) 07/04/11 : minor revisions in output formatting" << endl;
  cout << "  - (3.1.5) 14/03/12 : output error messages to user file." << endl;
  cout << "  - (3.2.0) 28/09/12 : exit on PWM-reading error." << endl;
  cout << "  - (3.2.1) 04/02/13 : no changes for MotifComparison." << endl;
  cout << "  - end." << endl;
  cout << endl; 

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
  if (bgFile != NULL)
    delete bgFile;
  bgFile = NULL;
  if (priorFile != NULL)
    delete priorFile;
  priorFile = NULL;
  return;
}
