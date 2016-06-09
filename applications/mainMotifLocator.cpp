// 27 march 2012 : adjust cerr reporting => cerrstr
// 28 sept 2012 : error (and exit) for incorrect PWM format reading

// 03 feb 2013 : compute minX/maxX against higher order bgmodel instead of SNF !
// --> then rescaled score will never be higher than 1 which may now be the case.


#include "utilities.h"
#include "PWMIO.h"
#include "PWM.h"
#include "BackgroundIO.h"
#include "BackgroundModel.h"
#include "Instance.h"
#include "FastaIO.h"
#include "RandomNumber.h"
#include "GFFWriter.h"
#include "SequenceComputation.h"

// c++ includes
#include <list>
#include <iostream>
#include <fstream>

// c includes
#include <getopt.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstring>

#define OPTIONS "hf:m:b:t:l:o:s:va:"


using namespace INCLUSIVE;
void instructions();
void version();
void cleanup();

string *fastaFile = NULL,
  *bgFile = NULL,
  *outputFile = NULL,
  *matrixFile = NULL,
  *listFile = NULL;


int
main(int argc, char *argv[])
{
  bool bFasta = false,
    bBG = false,
    bStrand = true,
    bOut = false,
    bList = false,
    bMatrix = false,
    bAbsolute = false;
  int iAbsolute = 0;
  strand_modes strand = both;
  char c;
  double threshold = 0.85;
  double minX = 0;
  double maxX = 0;
  int wLength = 0,
    i = 0;
  int j = 0;

  // 27 march 2012
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
    case 'm':
      matrixFile = new string(optarg);
      bMatrix = true;
      break;
    case 'b':
      bgFile = new string(optarg);
      bBG = true;
      break;
    case 't':
      threshold = (double) atof(optarg);
      if (threshold < 0 || threshold > 1 )
      {
        cerr << "--ERROR: -t threshold should be between 0 and 1." << endl;
        cerrstr << "--ERROR: -t threshold should be between 0 and 1." << endl;
        wronginput = true;
        //cleanup(); exit(-1);
      }
      break;
    case 'l':
      bList = true;
      listFile = new string(optarg);
      break;
    case 'o':
      // define new output stream
      bOut = true;
      outputFile = new string(optarg);
      break;
    case 's':
      bStrand = atoi(optarg);
      strand = both;            // default both strands
      if (bStrand == 0){ strand = plus_strand;}
      break;
    case 'a':
      iAbsolute = atoi(optarg);
      if (!(iAbsolute == 0 || iAbsolute == 1))
      {
        cerr << "--ERROR: -a score type should be 0 or 1." << endl;
        cerrstr << "--ERROR: -a score type should be 0 or 1." << endl;
        wronginput = true;
        //cleanup(); exit(-1);
      } 
      else { bAbsolute = iAbsolute;}
      break;
    case 'v':
      version();
      exit(-1);
    case '?':
    case 'h':
      instructions();
      exit(-1);
    default:
      cerr << "MotifLocator: Error in input getopt() function" << endl;
      cerrstr << "MotifLocator: Error in input getopt() function" << endl; 
      wronginput = true;
      //cleanup(); exit(-1);
    }
  }

  // no extra cerrstr reporting (no entry for these inputs is captured by the server)
  if (!bFasta || !bBG || !bMatrix)
  {
    instructions();
    cleanup();
    exit(-1);
  }
	
  GFFWriter *pGFFStream;
  if (wronginput)
  { 
    // report to file if file is available
    if (bOut) 
    { cerrstr << "#Check our MotifSuite webpage for correct input format. Or contact us, we will be happy to help." << endl;
      pGFFStream = new GFFWriter(outputFile); 
      pGFFStream->AddComment(cerrstr.str());
      cerrstr.flush(); // flush  
      delete pGFFStream; pGFFStream = NULL;
    }
    // exit
    cleanup(); exit(-1);  
  }

  cerr << "MotifLocator: Load the background model." << endl;
  // load the background model
  BackgroundModel * pBgModel = NULL;
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


  PWMIO *pwmIO = new PWMIO(matrixFile, READ);
  if (!pwmIO->IsOpen())
  {
    cerr << "--ERROR: MotifLocator: Unable to open matrix model file: " << *matrixFile
      << endl;
    cerrstr << "--ERROR: MotifLocator: Unable to open matrix model file: " << *matrixFile
      << endl;
    wronginput = 1;
  }

  if (wronginput)
  { 
    // report to file if file is available
    if (bOut) 
    { cerrstr << "#Check our MotifSuite webpage for correct input format. Or contact us, we will be happy to help." << endl;
      pGFFStream = new GFFWriter(outputFile); 
      pGFFStream->AddComment(cerrstr.str());
      cerrstr.flush(); // flush  
      delete pGFFStream; pGFFStream = NULL;
    }
    // exit
    cleanup(); exit(-1);  
  }
  // read list of matrix ID's
  list < string > matrixIdList;
  list < string >::iterator matListIter = matrixIdList.begin();
  if (bList)
  {
    ifstream _ifs(listFile->c_str(), ios::in);
    string inLine;
    while (!_ifs.eof())
    {
      getline(_ifs, inLine, '\n');
      matrixIdList.push_back(inLine);
    }
    _ifs.close();
  }

  // open file for output writing if necessarry
  if (bOut)
  {
    pGFFStream = new GFFWriter(outputFile);
  }

  else

  {
    pGFFStream = new GFFWriter();
  }

  // some variables to store matrix data
  PWM *myMatrix;
  list < PWM * >matrixList;
  list < double >minXList;
  list < double >maxXList;
  list < PWM * >::iterator matIter;
  list < double >::iterator minXListIter;
  list < double >::iterator maxXListIter;
  double value = 0;
  double minElem = 1;
  double maxElem = 0;
  double wx = 0;

  cerr << "MotifLocator : Loading matrices... ";
  int count = 0; 
  while (pwmIO->IsOpen())
  {
    myMatrix = pwmIO->ReadMatrix();
    if (myMatrix == NULL && (pwmIO->GetError()) != NULL) 
    { 
      cerr << *(pwmIO->GetError()) << endl;
      pGFFStream->AddComment(pwmIO->GetError()); 
      cerrstr << "--Error: MotifLocator : incorrect matrix format found in matrix input file." << endl;
      pGFFStream->AddComment(cerrstr.str());
      cerrstr.flush(); // flush 
      wronginput = true;
      break;
    }
    if (myMatrix == NULL){break;} // processing not-matrix lines... end of file
    //
    bool bMatching = false;
    if (bList)
    {
      matListIter = matrixIdList.begin();
      while (matListIter != matrixIdList.end())
      {
        if ((*matListIter) == *(myMatrix->GetID()))
        {
          bMatching = true;
          break;
        }
        matListIter++;
      }
    }
    else
    {
      bMatching = true;
    }
    if (bMatching)
    {

      // compute minimal and maximal score of matrix
      minX = 0;
      maxX = 0;
      for (i = 0; i < myMatrix->Length(); i++)
      {
        minElem = 1;
        maxElem = 0;

        for (j = 0; j < 4; j++)
        {
          //value = myMatrix->GetValueAt(i, j) / pBgModel->GetSnfValueAt(j);
          value = myMatrix->GetValueAt(i, j) / pBgModel->GetMaxTransitionMatrixValueAt(j);
          if (value < minElem){ minElem = value;}
          value = myMatrix->GetValueAt(i, j) / pBgModel->GetMinTransitionMatrixValueAt(j);
          if (value > maxElem){ maxElem = value;}
        }
        minX += log(minElem);
        maxX += log(maxElem);
      }
      // minX -= myMatrix->Length() * log(maxSnf);
      // maxX -= myMatrix->Length() * log(minSnf);

      /*cerr << "+ " << *(myMatrix->
                        GetID()) << ": min Wx = " << minX << " |max Wx = "
        << maxX << endl;
      */
      count++;
      matrixList.push_back(myMatrix);
      minXList.push_back(minX);
      maxXList.push_back(maxX);

      // myMatrix->StderrPrintMatrix();
    }
  }
	cerr << endl;
  cerr << count << " matrices loaded." << endl;

  if (!wronginput && count == 0) // count == 0
  { cerr << "--ERROR: MotifLocator: No usable input matrix found." << endl;
    cerrstr << "--ERROR: MotifLocator: No usable input matrix found." << endl;
    pGFFStream->AddComment(cerrstr.str());
    cerrstr.flush(); // flush
    wronginput = true;   
  }

  FastaIO *fileIO = NULL;
  if (!wronginput)
  {
    // open fasta file
    fileIO = new FastaIO(fastaFile);
    cerr << "MotifLocator : Loading sequences... " << endl;
    if (fileIO == NULL || !fileIO->IsOpen() || !fileIO->HasNext())
    {
      cerr << "--ERROR: MotifLocator: Unable to read sequences from fasta file." << endl;
      cerrstr << "--ERROR: MotifLocator: Unable to read sequences from fasta file." << endl;
      pGFFStream->AddComment(cerrstr.str());
      cerrstr.flush(); // flush 
      wronginput = true;
    }
  }
  if (wronginput)
  {
    // cleanup
    delete pGFFStream; pGFFStream = NULL;
    delete pwmIO; pwmIO = NULL;
    delete pBgModel;
    if (myMatrix != NULL) delete myMatrix; myMatrix = NULL;
    // cleanup matrices
    if (count != 0)
    {
      matIter = matrixList.begin();
      while (matIter != matrixList.end())
      { if ((*matIter) != NULL) delete (*matIter); 
        matIter++;
      }
    }
    cleanup(); exit(-1);
  }

  // close matrix reader
  delete pwmIO; pwmIO = NULL;

  // 22 july 2011
  string * intro = new string("#Please read MotifLocator Guidelines on our webserver to optimally evaluate below results.\n");  
  pGFFStream->AddComment(intro);
  delete intro; intro = NULL;

  // start processing the sequences
  int counter = 0,
    L = 0,
    nbr = 0;
  SequenceObject *pSeqObj = NULL;
  SequenceComputation *pSeqComp = NULL;
  string *Source = new string("MotifLocator");
  while (fileIO->HasNext())
  {
    pSeqObj = fileIO->NextSequence();
    if (pSeqObj == NULL)
    {
      continue;
    }

    // sequence length
    L = pSeqObj->Length();

    // create computation object
    pSeqComp = new SequenceComputation(pSeqObj);

    // augment sequence counter
    counter++;
    wLength = 0;
    Instance *pSite = NULL;

    // scroll through matrix list
    matIter = matrixList.begin();
    maxXListIter = maxXList.begin();
    minXListIter = minXList.begin();
    while (matIter != matrixList.end())
    {
      // (*matIter)->StderrPrintMatrix();
      nbr = 0;
      minX = *minXListIter;
      maxX = *maxXListIter;
      if (L <= 2 * ((*matIter)->Length()) )
      {

        // motif length is larger than sequence length move to next motif
        
		  cerr << "Warning: Sequence too short: " << *(pSeqObj->
                                                     GetID()) << endl;
        matIter++;
        minXListIter++;
        maxXListIter++;
        continue;
      }

      // score background model => only necessary when motif length changes
      //if ((*matIter)->Length() != wLength)
      //{
        wLength = (*matIter)->Length();
	      // define motif length in computation object
  	    pSeqComp->SetMotifLength(wLength);
        pSeqComp->UpdateInstanceBackgroundScore(pBgModel, wLength, strand);
      //}

      // score sequences with motif model
      pSeqComp->UpdateInstanceMotifScore(*matIter, strand);

      // compute Wx
      pSeqComp->UpdateInstanceExpWx(strand);

      // select best scoring instances on plus strand
      for (j = 0; j < L - wLength + 1; j++)
      {
        wx =
          ((log(pSeqComp->GetWxAt(j, plus_strand))) - minX) / (maxX - minX);

        // and select motifs with score higher than threshold
        if (wx >= threshold)
        {
          nbr++;
          pSite = new Instance(pSeqObj, plus_strand, j, wLength);

          // go to next instance
          if (pSite == NULL)
            continue;

          // set score 
          if (bAbsolute)
          {
            pSite->SetScore(pSeqComp->GetWxAt(j, plus_strand));
          }
          else
          {
            pSite->SetScore(wx);
          }

          // set instance ID
          pSite->SetID((*matIter)->GetID());

          // output site
          pGFFStream->WriteInstance(pSite, Source);
          delete pSite;
          pSite = NULL;
        }
      }

      // do the reverse complement
      if (bStrand)
      {
        for (j = 0; j < L - wLength + 1; j++)
        {

          // re-normalize wx
          wx =
            (log(pSeqComp->GetWxAt(j,minus_strand)) - minX) / (maxX - minX);
          // and select motifs with score higher than threshold
          if (wx >= threshold)
          {
            nbr++;

            // create a new instance
            pSite = new Instance(pSeqObj, minus_strand, j, wLength);
            if (pSite == NULL)
              continue;

            // set score
            if (bAbsolute)
            {
              pSite->SetScore(pSeqComp->GetWxAt(j, minus_strand));
            }
            else
            {
              pSite->SetScore(wx);
            }

            // set instance ID
            pSite->SetID((*matIter)->GetID());

            // output site
            pGFFStream->WriteInstance(pSite, Source);
            delete pSite;
            pSite = NULL;
          }
        }
      }
      if (nbr > 0)
		{
      cerr << counter << " |" << *(pSeqObj->
                                   GetID()) << " + " << *((*matIter)->
                                                          GetID()) <<
        "   -> instances:  " << nbr << endl;
		}
      
      //move on to the next element
      matIter++;
      maxXListIter++;
      minXListIter++;
    }                           // end matrix while loop
    delete pSeqObj;
    pSeqObj = NULL;
    delete pSeqComp;
    pSeqComp = NULL;
  }                             // end next sequence loop

  // close fasta file
  delete fileIO;

  // delete source name
  delete Source;

  // close outpfile
  delete pGFFStream;

  // clear matrixList
  matIter = matrixList.begin();
  while (matIter != matrixList.end())
  {

    // cerr << "deleting matrix: " << *((*matIter)->GetID()) << endl;
    delete(*matIter);
    matIter++;
  }
  delete pBgModel;
  delete myMatrix;
  cleanup();
  cerr << "MotifLocator: ending procedure successfully." << endl;
  exit(0);
}



/*********************
 *  LOCAL FUNCTIONS  *
 *********************/
void
instructions()
{
  cout << endl;
  cout << "Usage:" << " MotifLocator <ARGS>" << endl;
  cout << endl;
  cout << " Required Arguments" << endl;
  cout << "  -f <fastaFile>      Sequences in FASTA format" << endl;
  cout << "  -b <bgFile>         File containing the background model description" << endl;
  cout << "  -m <matrixFile>     File containing the matrix model descriptions" << endl;
  cout << "  -o <outFile>        Output file to write results in GFF (default stdout)" << endl;
  cout << endl;
  cout << " Optional Arguments" << endl;
  cout << "  -t <value>          Sets threshold above which a motif is selected (default 0.85)." << endl;
  cout << "  -l <listFile>       File with a list of identifiers to select individual" << endl;
  cout << "                      matrices from the matrix file. IDs that are not found are omitted." << endl;
  cout << "  -s <0|1>            Select strand. (default both)" << endl;
  cout << "                      0 is only input sequences, 1 include reverse complement." << endl;
  cout << "  -a <0|1>            1 will report the absolute scores, 0(=default) "
       << " will report the normalized scores as used in treshold comparison." << endl;
  cout << endl;
  cout << "  -v                  Version of MotifLocator" << endl;
  cout << endl;
  cout << "Version " << VERSION << endl;
  cout << "Questions and Remarks: <mclaeys@qatar.net.qa>" << endl
       << endl;
} void

version()
{
  cout << endl;
  cout << "INCLUSive -- MotifLocator " << endl;
  cout << "  version " << VERSION << endl;
  cout << endl; 
  cout << "Revision history : " << endl; 
  cout << "- (3.1.5) 02/03/10 : bug fixed in -a user parameter. " << endl;
  cout << "- (3.1.5) 07/04/11 : minor revisions in output formatting" << endl;
  cout << "- (3.1.5) 14/03/12 : output error messages to user file." << endl;
  cout << "- (3.2.0) 28/09/12 : exit on PWM-reading error." << endl;
  cout << "- (3.2.1) 04/02/13 : scaling factors against higher-order bg freqs." << endl;
  cout << "- (3.2.2) 16/01/14 : fully reset MaskVector at every next PWM assessment." << endl;
  cout << "- end." << endl;
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
  if (matrixFile != NULL)
    delete matrixFile;
  matrixFile = NULL;
  if (outputFile != NULL)
    delete outputFile;
  outputFile = NULL;
  if (listFile != NULL)
    delete listFile;
  listFile = NULL;
  return;
}
