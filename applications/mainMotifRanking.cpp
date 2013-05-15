// 28 october 2009 - 3.1.5
// GFFWriter.h/cpp
// report extra output to tracking file (-t)
// 20 march 2012 : adjust cerr reporting => cerrstr
// 28 sept 2012 : error (and exit) for incorrect PWM format reading

#include "utilities.h"
#include "PWMIO.h"
#include "PWM.h"
#include "BackgroundIO.h"
#include "BackgroundModel.h"
#include "GFFWriter.h" //

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

#define OPTIONS "hd:i:m:b:r:t:o:s:x:O:v"

void instructions();
void version();
void cleanup();
string *matrixFile = NULL,
  *dbFile = NULL,
  *outputFile = NULL,
  *textFile = NULL,
  *bgFile = NULL;


int
main(int argc, char *argv[])
{
  char c;
  bool bOut = false, bOut2 = false,
    bMatrix = false,
    bBackground = false;
  double threshold = 0.4;
  int shift = 1;
  int overlap = 6; 
  int rank = 5;
  int mode = 0;

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
    case 'i':
      matrixFile = new string(optarg);
      bMatrix = true;
      break;
    case 'b':
      bgFile = new string(optarg);
      bBackground = true;
      break;
    case 'r':
      rank = atoi(optarg);
      if (rank < 1)
      {
        cerr << "ERROR: -r rank should be greater than 0." << endl;
        cerrstr << "ERROR: -r rank should be greater than 0." << endl;
        wronginput = true;
        //cleanup(); exit(-1);
      }
      break;
    case 'm':
      mode = atoi(optarg);
      if (mode < 0 || mode > 2)
      {
        cerr << "ERROR: -m mode should be 0, 1 or 2." << endl;
        cerrstr << "ERROR: -m mode should be 0, 1 or 2." << endl;
        wronginput = true;
        //cleanup(); exit(-1);
      }
      break;
    case 't':
      threshold = (double) atof(optarg);
      if ( threshold <=0)
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
    case 'O':
      // define new output text stream 
      textFile = new string(optarg);
      bOut2 = 1;
      break;
    case 's':
      shift = atoi(optarg);
      if ( shift < 0 )
      {
        cerr << "ERROR: -s shift should be greater than 0." << endl;
        cerrstr << "ERROR: -s shift should be greater than 0." << endl;
        wronginput = true;
        //cleanup(); exit(-1);
      }
      break;
    case 'x': //
      overlap = atoi(optarg);
      if ( overlap <= 0 )
      {
        cerr << "ERROR: -x overlap should be greater than 0." << endl;
        cerrstr << "ERROR: -x overlap should be greater than 0." << endl;
        wronginput = true;
        //cleanup(); exit(-1);
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
      cerr << "-- MotifRanking: Error in input getopt() function" << endl;
      cerrstr << "MotifRanking: Error in input getopt() function" << endl; 
      wronginput = true;
      //cleanup(); exit(-1);
    }
  }

  // no cerrstr reporting, error is covered in website interface
  if (!bMatrix || !bOut) 
  {
    instructions();
    cleanup();
    exit(-1);
  }

  GFFWriter *pgffText = NULL;
  // open file for output writing 
  if (bOut2) pgffText = new GFFWriter(textFile);
  else pgffText = new GFFWriter();


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

  // process input matrix file
  PWMIO *pwmIO = new PWMIO(matrixFile, READ);
  // open matrix stream reader
  if (pwmIO->GetError() != NULL) // (file is not open)
  {
    pgffText->AddComment(pwmIO->GetError()); 
    delete pgffText; pgffText = NULL;  
    delete pwmIO; pwmIO = NULL;
    cleanup(); exit(-1);  
  }

  // open file for output writing if necessarry
  PWMIO *pwmFile = new PWMIO(outputFile, WRITE);
  if (pwmFile->GetError() != NULL)
  {
    pgffText->AddComment(pwmFile->GetError()); 
    delete pgffText; pgffText = NULL;  
    delete pwmFile; pwmFile = NULL;
    cleanup(); exit(-1);  
  }

  // add matrices to list
  PWM *myMatrix = NULL;
  list <PWM*>matrixList;
  list <PWM*>::iterator matIter;
  int inputcount = 0; //and count how many motifs are loaded
  while (pwmIO->IsOpen())
  { myMatrix = pwmIO->ReadMatrix(); 
    if (myMatrix == NULL && (pwmIO->GetError()) != NULL) 
    { 
      cerr << *(pwmIO->GetError()) << endl;
      pgffText->AddComment(pwmIO->GetError()); 
      cerr << "--Error: MotifRanking : incorrect matrix format found in matrix input file." << endl;
      cerrstr << "--Error: MotifRanking : incorrect matrix format found in matrix input file." << endl;
      pgffText->AddComment(cerrstr.str());
      cerrstr.flush(); // flush
      delete pgffText; pgffText = NULL;  
      delete pwmFile; pwmFile = NULL;
      cleanup(); exit(-1);
    }
    if (myMatrix == NULL){break;} // processing not-matrix lines end-of-file
    //
    matrixList.push_back(myMatrix); 
    inputcount++;
  }
  // close matrix reader
  delete pwmIO; pwmIO = NULL;
  if (inputcount == 0)
  {
    cerr << "--Error: MotifRanking : no usable matrix found in matrix input file." << endl;
    cerrstr << "--Error: MotifRanking : no usable matrix found in matrix input file." << endl;
    pgffText->AddComment(cerrstr.str());
    cerrstr.flush(); // flush  
    delete pgffText; pgffText = NULL;
    delete pwmFile; pwmFile = NULL;
    // exit
    cleanup(); exit(-1); 
  }
  string * comment = new string("");
  // 07 april 2011
 *comment = "#For each reported motif, the identifiers of similar motifs are listed in order as they have been loaded by MotifRanking and NOT by their absolute similarity to the reported motif. The higher the motif score and the motif count ('==>Total') of a top-ranked motif, the more likely it represents a true motif. We advise a minimal count of 10 similars for a top-ranked motif to be a reliable prediction. If you doubt that top-ranked motifs are truly different, we advise to re-evaluate their (non-)similarity with MotifComparison (using p-BLiC metric).\n#Please consult the MotifRanking guidelines available on MotifSuite to optimally interpretate below results.\n";
  pgffText->AddComment(comment);
  *comment = "";

  // open background file to read single nucleotide frequency
  wronginput = false;
  double snf[4];
  for (int i = 0; i < 4; i++) snf[i] = 0.25; // default
  if (mode == 2)
  {  wronginput = 0;
    if ( bBackground )
    {
      BackgroundIO *bgIO = new BackgroundIO(*bgFile, READ);
      if (bgIO->GetError() != NULL) 
      { wronginput = 1; cerrstr << "--Warning: " << *(bgIO->GetError()) << " Default values will be used instead." << endl; }
      else
      { BackgroundModel *pBgModel = bgIO->ReadBackgroundModel();
        if (pBgModel == NULL)
        { wronginput = 1; cerrstr << "--Warning: background model not suitable. Default values will be used instead." << endl; }
        else
        { // store snf
          for (int i = 0; i < 4; i++) snf[i] = pBgModel->GetSnfValueAt(i);
          delete pBgModel; pBgModel = NULL;
        }
      }
      delete bgIO; bgIO = NULL;
    }
    if (wronginput)
    {
      pgffText->AddComment(cerrstr.str());
      cerrstr.flush(); // flush
    }
    // no exit, just go on with default values
  } 
  double score = 0,
    maxScore = 0;
  // set scores if asked for CS or LL
  matIter = matrixList.begin();
  if (mode != 0)
  {
    while (matIter != matrixList.end())
    {

      // cerr << "Updating score to mode: " << mode << endl;
      if (mode == 1)
      {
        score = (*matIter)->ConsensusScore();
      }
      else if (mode == 2)
      {
        score = (*matIter)->InformationContent(snf);
      }
      else
      {
        break;
      }
      (*matIter)->SetScore(score);
      matIter++;
    }
  }
		
  int count = 0;
  int L1(0), L2(0); //
  for (int i = 0; i < rank; i++)
  {

    // look for best scoring motif 
    matIter = matrixList.begin();
    maxScore = 0;
    while (matIter != matrixList.end())
    {
      score = (*matIter)->Score();
      if (maxScore < score)
      {
        myMatrix = (*matIter);
        maxScore = myMatrix->Score();
      }
      matIter++;
    }

    // write best scoring motif to file
    int shiftM; // maximum shift between two matrices, respecting minumum overlap
    if (maxScore != 0)
    {
      pwmFile->WriteMatrix(myMatrix);
      // add to file/screen    i)id.. score
      ostringstream rankstr; rankstr << i+1 << ")";
      *comment = rankstr.str();
      comment->append(*(myMatrix->GetID())); comment->append("\t");
      ostringstream scorestr; scorestr << myMatrix->Score();
      comment->append(scorestr.str());
      pgffText->AddComment(comment);
      // also add to screen (will appear double on screen when no txt-file is given)
      cerr << *comment << endl;

      // reset score 
      myMatrix->SetScore(0);
      //cerr << endl << "= " << *(myMatrix->GetID()) << endl;

      count = 1; // also count the ranked motif itself
      double mi = 0;
      matIter = matrixList.begin();
      while (matIter != matrixList.end())
      {
		  
		  // 2011/04/29 : do not re-evaluate motifs already assigned
		if ((*matIter)->Score() != 0)
		{
		  
        // take most stringent value from shift/overlap 
        shiftM = ( (L1 = myMatrix->Length()) < (L2 = (*matIter)->Length()) ?
                    (L1-overlap) : (L2-overlap) );
        if (shiftM < 0) {shiftM = 0;} 
        else {if (shiftM > shift) {shiftM = shift;}}

        mi =
          (myMatrix->MutualInformation((*matIter), shiftM) +
           (*matIter)->MutualInformation(myMatrix, shiftM)) / 2;
  
        // find similar motif
        if (mi < threshold)
        {
          count++;
          //cerr << "   " << *((*matIter)->GetID()) << "\t" 
          //     << *((*matIter)->GetConsensus()) << endl;
          *comment = "   "; comment->append(*((*matIter)->GetID()));
          comment->append("\t"); 
          ostringstream score2str; score2str << (*matIter)->Score();
          comment->append(score2str.str());
          pgffText->AddComment(comment);
			
          (*matIter)->SetScore(0);
        }
		// end of added restriction 
		}
        matIter++;
      }
      *comment = "==> Total: "; 
      ostringstream totalstr; totalstr << count;
      // and add also RR as statistical significance 22 july 2011
      totalstr << "\t(RR = " << int(100*count/inputcount) << "%)";
      comment->append(totalstr.str());
      pgffText->AddComment(comment);
      // also add to screen (will appear double on screen when no txt-file is given)
      cerr << *comment << endl;
    }
  }

  // let myMatrix point to NULL
  myMatrix = NULL;

  // clear matrix list
  matIter = matrixList.begin();
  while (matIter != matrixList.end())
  {
    delete(*matIter);
    matIter++;
  }
		
  delete pgffText; pgffText = NULL; //
  delete pwmFile; pwmFile = NULL; //
  if (comment != NULL) delete comment; comment = NULL;
		
  cleanup();
  
  // program was succesful
  cerr << "MotifRanking: ending procedure successfully." << endl;
  exit(0);
}



/*********************
 *  LOCAL FUNCTIONS  *
 *********************/
void
instructions()
{
  cout << endl;
  cout << "Usage:" << " MotifRanking <ARGS>" << endl;
  cout << endl;
  cout << " Required Arguments" << endl;
  cout << "  -i <matrixFile>     File containing the matrix model descriptions" << endl;
  cout << endl;
  cout << " Output Arguments" << endl;
  cout << "  -o <outFile>        Output file to write matrix results to." << endl;
  cout << "  -O <outFile>        Output file to write txt results to (default screen)." << endl;
  cout << endl;
  cout << " Optional Arguments" << endl;
  cout << "  -m <0|1|2>          Select the type of score: " << endl;
  cout << "                       0 = Score defined in matrix file (default)" << endl;
  cout << "                       1 = Consensus score [2+plog(p)]" << endl;
  cout << "                       2 = Information Content [plog(p/p0)] (requires -b)" << endl;
  cout << "  -b <file>           Load Background model from file to set single nucleotide frequency" << endl;
  cout << "  -r <value>          Number of topscoring motifs to display (default 5)." << endl;
  cout << "  -t <value>          Sets threshold below two motifs are seen as being the same." << endl;
  cout << "                      Default value is set to 0.4." << endl;
  cout << "  -s <value>          Maximal allowed shift between motif models." << endl;
  cout << "                      Default<1>." << endl;
  cout << "  -x <value>          Minimum required overlap between motif models." << endl;
  cout << "                      Default<6>. " << endl;
  cout << endl;
  cout << "  -v                  Version of MotifRanking" << endl;
  cout << "Version " << VERSION << endl;
  cout << "Questions and Remarks: <mclaeys@qatar.net.qa>" << endl;
  cout << endl;
} 


void
version()
{
  cout << endl;
  cout << "INCLUSive -- MotifRanking (stand alone C++ version)" << endl;
  cout << "  version " << VERSION << endl;
  cout << endl;
  cout << "Revision history : " << endl; 
  cout << "- (3.1.2) 14/08/07 : explicit x-parameter for minimum overlap." << endl;
  cout << "- (3.1.4) 14/01/08 : debug error input x-parameter." << endl;  
  cout << "- (3.1.5) 21/10/09 : always use most stringent of -x and -s." << endl;
  cout << "- (3.1.5) 28/10/09 : optional -O txt.output to file." << endl;
  cout << "- (3.1.5) 07/04/11 : minor revisions in output formatting" << endl;
  cout << "- (3.1.5) 29/04/11 : debug wrongly re-assignment of assigned motifs." << endl;
  cout << "- (3.1.5) 22/07/11 : add RR in output formatting" << endl;
  cout << "- (3.1.5) 14/03/12 : output error messages to user file." << endl;
  cout << "- (3.2.0) 28/09/12 : exit on PWM-reading error." << endl;
  cout << "- (3.2.1) 04/02/13 : no changes for MotifRanking." << endl;
  cout << "- end." << endl;
  cout << endl;
} 



void
cleanup()
{
  if (matrixFile != NULL)
    delete matrixFile;
  matrixFile = NULL;
  if (outputFile != NULL)
    delete outputFile;
  outputFile = NULL;
  if (textFile != NULL)
    delete textFile;
  textFile = NULL;
  return;
}
