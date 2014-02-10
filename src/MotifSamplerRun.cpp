#include "MotifSamplerRun.h"
#include "FastaIO.h"
#include <math.h>
// release version 3.1.2 : 1 bug fixed in LogLikelihoodScore
// release version 3.1.4 : 1 bug fixed in UpdateMasksFromInstanceMap
// release version 3.1.5 : multiple priors and sampling extensions
// release version 3.2.0 : implement PSP
// release version 3.2.2 : revise PSP

/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
MotifSamplerRun::MotifSamplerRun(string * pFastaFile, strand_modes strand, GFFWriter * pgff)
{
  // cerr << "DEBUG: MotifSamplerRun::MotifSamplerRun" << endl;
  _pgffio = pgff;
  // initialize sequence set
  _strand = strand;
  _pSequenceList = new SequenceSet;
  _pComputationMap = new map < SequenceObject *, SequenceComputation * >;
  _nbrSequences = 0;
  _w = 0; // initialize the motif length to 0
  _pPriorDistributions = NULL; // to be set later
  ostringstream cerrstr;

  // internal sequence pointer
  SequenceObject *pSeq = NULL;
  SequenceComputation *pSeqComp = NULL;

  // open fasta file
  cerr << "Open sequence file for reading." << endl;
  FastaIO *fileIO = new FastaIO(pFastaFile);
  if (fileIO == NULL)
  {
    cerr << "--Error-- MotifSamplerRun::_LoadSequences(): "
      << "Unable to open sequence fasta file" << endl;
    cerrstr << "--Error-- MotifSamplerRun::_LoadSequences(): "
      << "Unable to open sequence fasta file" << endl;
    _pgffio->AddComment(cerrstr.str());
    cerrstr.flush(); // flush     
    return;
  }
  else
  {
    while (fileIO->HasNext())
    {
      // read next sequence
      cerr << "Next sequence: ";
      pSeq = fileIO->NextSequence();
      if (pSeq != NULL)
      {
        // add sequence to sequence list
        _pSequenceList->push_back(pSeq);
        _nbrSequences++;
        // initialize computation object
        cerr << *(pSeq->GetID()) << endl;
        pSeqComp = new SequenceComputation(pSeq);
        // add key/value pair to computation map
        _pComputationMap->insert(MapValue(pSeq, pSeqComp));
      }
      else
      {
        cerr << "--Warning-- MotifSamplerRun::_LoadSequences(): "
          << "Empty sequence." << endl;
        cerrstr << "--Warning-- MotifSamplerRun::_LoadSequences(): "
          << "Empty sequence." << endl;
        _pgffio->AddComment(cerrstr.str());
        cerrstr.flush(); // flush
      }
    }
  }

  // set internal pointers to NULL
  pSeq = NULL;
  pSeqComp = NULL;

  // close fasta file
  delete fileIO;

  // number of instances 
  _nbrInstances = 0;
  _pInstanceMap = new InstanceMap();
  _nbrSequencesWithInstances = 0;

  // cerr << "DEBUG: MotifSamplerRun::MotifSamplerRun" << endl;
  // set pseudo counts to [0.25, 0.25, 0.25, 0.25]
  _pPseudoCounts = new double[4];
  for (int i=0; i<4; i++)
    _pPseudoCounts[i] = 0.25;

  // dividing constant
  _scale = 1.0;
  if ( _nbrSequences > 0 )
    _scale = 1/sqrt((double) _nbrSequences);

  _localMotif = NULL;
  _pLocalCounts = NULL;


}


/****************************************************************************
  Method:       ~MotifSamplerRun ()
  Class:        
  Arguments:    none
  
  Description:  destructor
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
****************************************************************************/
MotifSamplerRun::~MotifSamplerRun()
{
  // delete individual sequences
  // cerr << "DEBUG: MotifSamplerRun::~MotifSamplerRun clear sequence list" << endl;
  if (_pSequenceList != NULL)
  {
    SeqIterator i = _pSequenceList->begin();
    for (; i != _pSequenceList->end(); i++)
    {
	  if ( *i !=NULL )
		delete *i;
      *i = NULL;
    }

    // delete vector itself
    delete _pSequenceList;
  }
  _pSequenceList = NULL;

  // clear the computation map
  // cerr << "DEBUG: MotifSamplerRun::~MotifSamplerRun clear computation list" << endl;
  if (_pComputationMap != NULL)
  {
    MapIterator i = _pComputationMap->begin();
    for (; i != _pComputationMap->end(); i++)
    {
	  if ( (*i).second != NULL )
        delete((*i).second);
      (*i).second = NULL;
    }
    delete _pComputationMap;
  }
  _pComputationMap = NULL;
  
  
  // cerr << "DEBUG: MotifSamplerRun::~MotifSamplerRun clear instance map" << endl;
  if (_pInstanceMap != NULL)
  {
    _pInstanceMap->ClearMap();
    delete _pInstanceMap;
  }
  _pInstanceMap = NULL;
  if (_localMotif != NULL)
    delete _localMotif;
  _localMotif = NULL;
  
  // delete local matrices 
  delete[] _pPseudoCounts;
  for (int i=0; i<_w; i++)
    delete[] _pLocalCounts[i];
  delete[] _pLocalCounts;

  // _pPriorCopyDistribution 
  if (_pPriorDistributions != NULL)
  { // cleanup distributions
    for(int i = 0; i < (int)_pPriorDistributions->size(); i++)
    { if ( (*_pPriorDistributions)[i] != NULL)
        delete (*_pPriorDistributions)[i];
      (*_pPriorDistributions)[i] = NULL;
    }
    delete _pPriorDistributions;
  }
  _pPriorDistributions = NULL;

  // file unlink
  _pgffio = NULL;
  // cerr << "DEBUG: MotifSamplerRun::~MotifSamplerRun" << endl;
}



/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
void
  MotifSamplerRun::_ClearInstanceMap()
{
  // cerr << "DEBUG: MotifSamplerRun::_ClearInstanceMap" << endl;
  if (_pInstanceMap != NULL)
    _pInstanceMap->ClearMap();
  _nbrInstances = 0;
  _nbrSequencesWithInstances = 0;
  // cerr << "DEBUG: MotifSamplerRun::_ClearInstanceMap" << endl;
  return;
}


/*****************************************************************************
  Method:       SetPseudoCounts
  Class:        MotifSamplerRun 
  Arguments:    double* psCounts
  
  Description:  Update pseudo counts

  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
*****************************************************************************/
void
  MotifSamplerRun::UpdatePseudoCounts(double *psCounts)
{
  //cerr << "Pseudo Counts = ";
  for(int i=0; i<4; i++){
    _pPseudoCounts[i] = _scale * psCounts[i];
    //cerr << _pPseudoCounts[i] << "  ";
  }
  //cerr << endl;  
  
  return;
}


/*****************************************************************************
  Method:       UpdatePseudoCounts
  Class:        MotifSamplerRun 
  Arguments:    double* psCounts
                double scale
  Description:  Update pseudo counts

  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
*****************************************************************************/
void
  MotifSamplerRun::UpdatePseudoCounts(double *psCounts, double scale)
{
  //cerr << "Pseudo Counts = ";
  _scale = scale;
  for(int i=0; i<4; i++){
    _pPseudoCounts[i] = _scale * psCounts[i];
    //cerr << _pPseudoCounts[i] << "  ";
  }
  //cerr << endl;  
  
  return;
}


/*****************************************************************************
  Method:       InitFixedSizeInstanceMap
  Class:        MotifSamplerRun 
  Arguments:    none
  
  Description:  Generate an InstanceMap by sampling exactly one instance
                in each sequence.
                Results are stored in _pInstanceMap

  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
*****************************************************************************/
void
  MotifSamplerRun::InitFixedSizeInstanceMap(int nbr)
{

  // initialize local variables
  Instance *
    pSite = NULL;
  SequenceObject *
    pSeq = NULL;
  SequenceComputation *
    pComp = NULL;
	  ostringstream cerrstr;
  bool foundPlus = false;
  // cerr << "DEBUG: MotifSamplerRun::InitFixedSizeInstanceMap" << endl;
  if (nbr < 1)
  {
    cerr <<
      "--Error-- MotifSamplerRun::InitFixedSizeInstanceMap(): negative number of instances asked for."
      << endl;
    cerrstr << "--Error-- MotifSamplerRun::InitFixedSizeInstanceMap(): negative number of instances asked for."
      << endl;
    _pgffio->AddComment(cerrstr.str());
    cerrstr.flush(); // flush
    if (_pInstanceMap != NULL)
      _pInstanceMap->ClearMap();
    _pInstanceMap = NULL;
    return;
  }

  // vector to store start position
  vector < int >
  pVec(nbr);

  // check current state of the instance map
  if (_pInstanceMap != NULL)
  {
    _pInstanceMap->ClearMap();
  }
  else
  {
    _pInstanceMap = new InstanceMap;
  }
  _nbrInstances = 0;
  _nbrSequencesWithInstances = 0;
  MapIterator mi = _pComputationMap->begin();
  int
    i = 0;
  while (mi != _pComputationMap->end())
  {

    // get key and value of current map entry
    pSeq = (*mi).first;         // sequence object 
    pComp = (*mi).second;       // sequence computation object

    // reset boolean to indicate if an instance is selected on the plus strand
    foundPlus = false;
    if (_strand == plus_strand || _strand == both)
    {

      // sample nbr start position
      pComp->SampleUniformInstanceStart(pVec, nbr, plus_strand);
      for (int j = 0; j < nbr; j++)
      {
        if (pVec[j] != -1)
        {

          // create instance from selected start position
          pSite = new Instance(pSeq, plus_strand, pVec[j], _w);
          if (pSite != 0)
          {
            // cerr << "DEBUG " << *(pSite->ParentSequence()->GetID()) << " sampled site: " << *(pSite->PrintSite()) << endl;
            pSite->SetScore(pComp->GetWxAt(pVec[j], plus_strand));
            _nbrInstances++;
            if (!foundPlus)
              _nbrSequencesWithInstances++;
            foundPlus = true;

            // append sites to InstanceMap
            _pInstanceMap->AddInstance(pSite);
          }
        }
      }
    }
    if (_strand == minus_strand || _strand == both)
    {

      // sample nbr start positions
      pComp->SampleUniformInstanceStart(pVec, nbr, minus_strand);
      for (int j = 0; j < nbr; j++)
      {
        if (pVec[j] != -1)
        {

          // create instance
          pSite = new Instance(pSeq, minus_strand, pVec[j], _w);
          if (pSite != 0)
          {
            // cerr << "DEBUG " << *(pSite->ParentSequence()->GetID()) << " sampled site: " << *(pSite->PrintSite()) << endl;
            pSite->SetScore(pComp->GetWxAt(pVec[j], minus_strand));
            _nbrInstances++;
            if (!foundPlus)
              _nbrSequencesWithInstances++;

            // append sites to InstanceMap
            _pInstanceMap->AddInstance(pSite);
          }
        }
      }
    }

    // augment counter and iterator
    i++;
    mi++;
  }
  //cerr << "DEBUG: MotifSamplerRun::InitFixedSizeInstanceMap" << endl;
  return;
}



/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
void
  MotifSamplerRun::SampleFixedSizeInstanceMap(int nbr)
{

  // initialize local variables
  Instance *
    pSite = NULL;
  SequenceObject *
    pSeq = NULL;
  SequenceComputation *
    pComp = NULL;
  bool foundInstance = false;
  // cerr << "DEBUG: MotifSamplerRun::SampleFixedSizeInstanceMap" << endl;
  if (nbr < 1)
  {
    if (_pInstanceMap != NULL)
      _pInstanceMap->ClearMap();
    _pInstanceMap = NULL;
    return;                     // empty instance map
  }

  // vector to store start position
  vector < int >
  pVec(nbr);

  // check current state of the instance map
  if (_pInstanceMap != NULL)
  {

    // clear all sites from the current instance map
    _pInstanceMap->ClearMap();
  }
  else
  {

    // instance map is not defined create a new one
    _pInstanceMap = new InstanceMap;
  }
  _nbrInstances = 0;
  _nbrSequencesWithInstances = 0;
  MapIterator mi = _pComputationMap->begin();
  int
    i = 0;
  while (mi != _pComputationMap->end())

  {

    // get key and value of current map entry
    pSeq = (*mi).first;         // sequence object 
    pComp = (*mi).second;       // sequence computation object

    // reset boolean 
    foundInstance = false;
    if (_strand == plus_strand || _strand == both)
    {

      // sample nbr start positions
      pComp->SampleInstanceStart(pVec, nbr, plus_strand);
      for (int j = 0; j < nbr; j++)
      {
        if (pVec[j] != -1)
        {

          // create instance
          // cerr << "DEBUG: Create new instance :" << pVec[j] << endl;
          pSite = new Instance(pSeq, plus_strand, pVec[j], _w);
          if (pSite != 0)
          {
              // cerr << "DEBUG sampled site: " << *(pSite->PrintSite()) << endl;
            pSite->SetScore(pComp->GetWxAt(pVec[j], plus_strand));
            _nbrInstances++;
            if (!foundInstance)
              _nbrSequencesWithInstances++;
            foundInstance = true;

            // append sites to InstanceMap
            // cerr << "DEBUG: Append site to instance map" << endl;
            _pInstanceMap->AddInstance(pSite);
          }
        }
      }
    }
    if (_strand == minus_strand || _strand == both)
    {

      // sample nbr start positions
      pComp->SampleInstanceStart(pVec, nbr, minus_strand);
      for (int j = 0; j < nbr; j++)
      {
        if (pVec[j] != -1)
        {

          // create instance
          pSite = new Instance(pSeq, minus_strand, pVec[0], _w);
          if (pSite != 0)
          {
              // cerr << "DEBUG sampled site: " << *(pSite->PrintSite()) << endl;
            pSite->SetScore(pComp->GetWxAt(pVec[j], minus_strand));
            _nbrInstances++;
            if (!foundInstance)
              _nbrSequencesWithInstances++;
            foundInstance = true;

            // append sites to InstanceMap
            _pInstanceMap->AddInstance(pSite);
          }
        }
      }
    }

    // augment counter and iterator
    i++;
    mi++;
  }
  // cerr << "DEBUG: MotifSamplerRun::SampleFixedSizeInstanceMap" << endl;
  return;
}


/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
void
  MotifSamplerRun::SampleInstanceMap()
{

  // initialize local variables
  Instance *
    pSite = NULL;
  SequenceObject *
    pSeq = NULL;
  SequenceComputation *
    pComp = NULL;
  bool foundInstance = false;
  int
    nbr = 0;
  // cerr << "DEBUG: MotifSamplerRun::SampleInstanceMap" << endl;

  // vector to store start position
  vector < int >
  pVec(1);

  // check current state of the instance map
  if (_pInstanceMap != NULL)
  {
    // clear all sites from the current instance map
    _pInstanceMap->ClearMap();
  }
  else
  {

    // instance map is not defined create a new one
    _pInstanceMap = new InstanceMap;
  }
  _nbrInstances = 0;
  _nbrSequencesWithInstances = 0;
  MapIterator mi = _pComputationMap->begin();
  while (mi != _pComputationMap->end())
  {
    // get key and value of current map entry
    pSeq = (*mi).first;         // sequence object 
    pComp = (*mi).second;       // sequence computation object

    // reset boolean 
    foundInstance = false;
    if (_strand == plus_strand || _strand == both)
    {

      // get number of estimated motif instances
      //nbr = pComp->GetEstimatedNumberInstances(plus_strand);
      nbr = pComp->GetNumberInstances(plus_strand, 0);//
      // cerr << "DEBUG: Number of instances found (+): " << nbr << endl;

      if (nbr > 0)
      {

        // update vector to store start positions
        if ( nbr > (int)pVec.size() )
          pVec.resize(nbr, -1);

        // sample nbr start positions
        pComp->SampleInstanceStart(pVec, nbr, plus_strand);
        for (int j = 0; j < nbr; j++)
        {
          if (pVec[j] != -1)
          {

            // create instance
            pSite = new Instance(pSeq, plus_strand, pVec[j], _w);
            if (pSite != 0)
            {
              // cerr << "DEBUG sampled site: " << *(pSite->PrintSite()) << endl;
              pSite->SetScore(pComp->GetWxAt(pVec[j], plus_strand));
              _nbrInstances++;
              if (!foundInstance)
                _nbrSequencesWithInstances++;
              foundInstance = true;

              // append sites to InstanceMap
              // cerr << "Append site to instance map" << endl;
              _pInstanceMap->AddInstance(pSite);
            }
          }
        }
      }
    }
    if (_strand == minus_strand || _strand == both)
    {

      // get number of estimated motif instances
      //nbr = pComp->GetEstimatedNumberInstances(minus_strand);
      nbr = pComp->GetNumberInstances(minus_strand, 0);//
      // cerr << "DEBUG: Number of instances found (-): " << nbr << endl;

      if (nbr > 0)
      {

        // update vector to store start positions
        if (nbr > (int) pVec.size())
          pVec.resize(nbr, -1);

        // sample nbr start positions
        pComp->SampleInstanceStart(pVec, nbr, minus_strand);
        for (int j = 0; j < nbr; j++)
        {
          if (pVec[j] != -1)
          {

            // create instance
            pSite = new Instance(pSeq, minus_strand, pVec[0], _w);
            if (pSite != 0)
            {
              // cerr << "DEBUG sampled site: " << *(pSite->PrintSite()) << endl;
              pSite->SetScore(pComp->GetWxAt(pVec[j], minus_strand));
              _nbrInstances++;
              if (!foundInstance)
                _nbrSequencesWithInstances++;
              foundInstance = true;

              // append sites to InstanceMap
              _pInstanceMap->AddInstance(pSite);
            }
          }
        }
      }
    }

    // augment counter and iterator
    mi++;
  }
  // cerr << "DEBUG: MotifSamplerRun::SampleInstanceMap" << endl;
  return;
}



/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
void
  MotifSamplerRun::SelectFixedSizeBestInstanceMap(int nbr)
{

  // initialize local variables
  Instance *
    pSite = NULL;
  SequenceObject *
    pSeq = NULL;
  SequenceComputation *
    pComp = NULL;
  bool foundInstance = false;
  // cerr << "DEBUG: MotifSamplerRun::SelectFixedSizeBestInstanceMap" << endl;

  // check number of instances asked for
  if (nbr < 1)
  {
    if (_pInstanceMap != NULL)
      _pInstanceMap->ClearMap();
    _pInstanceMap = NULL;
    return;
  }

  // create vector to store start positions
  vector < int >
  pVec(nbr);

  // clear instance map
  if (_pInstanceMap != NULL)
  {
    _pInstanceMap->ClearMap();
  }
  else
  {
    _pInstanceMap = new InstanceMap();
  }
  _nbrInstances = 0;
  _nbrSequencesWithInstances = 0;

  // iterate through sequences and sample initial site
  MapIterator mi = _pComputationMap->begin();
  while (mi != _pComputationMap->end())
  {
    pSeq = (*mi).first;
    pComp = (*mi).second;

    // So far, no instance was found 
    foundInstance = false;
    if (_strand == plus_strand || _strand == both)
    {

      // select nbr best scoring start positions
      pComp->SelectBestInstanceStart(pVec, nbr, plus_strand);
      for (int j = 0; j < nbr; j++)
      {
        if (pVec[j] != -1)
        {
          // create instance
          pSite = new Instance(pSeq, plus_strand, pVec[j], _w);
          if (pSite != 0)
          {
              // cerr << "DEBUG sampled site: " << *(pSite->PrintSite()) << endl;
            pSite->SetScore(pComp->GetWxAt(pVec[j], plus_strand));
            _nbrInstances++;
            if (!foundInstance)
              _nbrSequencesWithInstances++;
            foundInstance = true;

            // append sites to InstanceMap
            _pInstanceMap->AddInstance(pSite);
          }
        }
      }
    }
    if (_strand == minus_strand || _strand == both)
    {

      // sample 1 start position
      pComp->SelectBestInstanceStart(pVec, nbr, minus_strand);
      for (int j = 0; j < nbr; j++)
      {
        if (pVec[j] != -1)
        {
          // create instance
          pSite = new Instance(pSeq, minus_strand, pVec[j], _w);
          if (pSite != 0)
          {
              // cerr << "DEBUG sampled site: " << *(pSite->PrintSite()) << endl;
            pSite->SetScore(pComp->GetWxAt(pVec[j], minus_strand));
            _nbrInstances++;
            if (!foundInstance)
              _nbrSequencesWithInstances++;
            foundInstance = true;

            // append sites to InstanceMap
            _pInstanceMap->AddInstance(pSite);
          }
        }
      }
    }

    // go to the next sequence in the set
    ++mi;
  }
  // cerr << "DEBUG: MotifSamplerRun::SelectFixedSizeBestInstanceMap" << endl;
  return;
}



/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
/*void
  MotifSamplerRun::SelectMaxSizeBestInstanceMap(int nbr)
{

  // initialize local variables
  Instance *
    pSite = NULL;
  SequenceObject *
    pSeq = NULL;
  SequenceComputation *
    pComp = NULL;
  bool foundInstance = false;
  int
    n = 0;
  // cerr << "DEBUG: MotifSamplerRun::SelectMaxSizeBestInstanceMap" << endl;

  // check number of instances asked for
  if (nbr < 1)
  {
    if (_pInstanceMap != NULL)
      _pInstanceMap->ClearMap();
    _pInstanceMap = NULL;
    return;
  }

  // create vector to store start positions
  vector < int >
  pVec(nbr);

  // clear instance map
  if (_pInstanceMap != NULL)
  {
    _pInstanceMap->ClearMap();
  }
  else
  {
    _pInstanceMap = new InstanceMap();
  }
  _nbrInstances = 0;
  _nbrSequencesWithInstances = 0;

  // iterate through sequences and sample initial site
  MapIterator mi = _pComputationMap->begin();
  while (mi != _pComputationMap->end())
  {
    pSeq = (*mi).first;
    pComp = (*mi).second;

    // So far, no instance was found 
    foundInstance = false;
    if (_strand == plus_strand || _strand == both)
    {

      // get number of estimated motif instances
      n = pComp->GetEstimatedNumberInstances(plus_strand);
      if (n < 1)
      {

        // go to the next sequence
        mi++;
        continue;
      }

      // check if the estimated number of instances is larger than nbr
      if (n > nbr)
        n = nbr;

      // select min(n,nbr) best scoring start positions
      pComp->SelectBestInstanceStart(pVec, n, plus_strand);
      for (int j = 0; j < n; j++)
      {
        if (pVec[j] != -1)
        {

          // create instance
          pSite = new Instance(pSeq, plus_strand, pVec[j], _w);
          if (pSite != 0)
          {
              // cerr << "DEBUG sampled site: " << *(pSite->PrintSite()) << endl;
            pSite->SetScore(pComp->GetWxAt(pVec[j], plus_strand));
            _nbrInstances++;
            if (!foundInstance)
              _nbrSequencesWithInstances++;
            foundInstance = true;

            // append sites to InstanceMap
            _pInstanceMap->AddInstance(pSite);
          }
        }
      }
    }
    if (_strand == minus_strand || _strand == both)
    {

      // get number of estimated motif instances
      n = pComp->GetEstimatedNumberInstances(minus_strand);
      if (n < 1)
      {

        // go to the next sequence
        mi++;
        continue;
      }

      // check if the estimated number of instances is larger than nbr
      if (n > nbr)
        n = nbr;

      // select nbr best start position
      pComp->SelectBestInstanceStart(pVec, n, minus_strand);
      for (int j = 0; j < n; j++)
      {
        if (pVec[j] != -1)
        {

          // create instance
          pSite = new Instance(pSeq, minus_strand, pVec[j], _w);
          if (pSite != 0)
          {
              // cerr << "DEBUG sampled site: " << *(pSite->PrintSite()) << endl;
            pSite->SetScore(pComp->GetWxAt(pVec[j], minus_strand));
            _nbrInstances++;
            if (!foundInstance)
              _nbrSequencesWithInstances++;
            foundInstance = true;

            // append sites to InstanceMap
            _pInstanceMap->AddInstance(pSite);
          }
        }
      }
    }

    // go to the next sequence in the set
    ++mi;
  }
  // cerr << "DEBUG: MotifSamplerRun::SelectMaxSizeBestInstanceMap" << endl;
  return;
}

*/

/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
void
  MotifSamplerRun::SelectBestInstanceMap()
{
  //cerr << "debug:MotifSamplerRun::SelectBestInstanceMap - begin" << endl;
  // initialize local variables
  Instance *
    pSite = NULL;
  SequenceObject *
    pSeq = NULL;
  SequenceComputation *
    pComp = NULL;
  bool foundInstance = false;
  int
    nbr = 1;
  // cerr << "DEBUG: MotifSamplerRun::SelectBestInstanceMap" << endl;

  // create vector to store start positions
  vector < int >
  pVec(1);

  // clear instance map
  if (_pInstanceMap != NULL)
  {
    _pInstanceMap->ClearMap();
  }
  else
  {
    _pInstanceMap = new InstanceMap();
  }
  _nbrInstances = 0;
  _nbrSequencesWithInstances = 0;

  // iterate through sequences and sample initial site
  MapIterator mi = _pComputationMap->begin();
  while (mi != _pComputationMap->end())
  {
    pSeq = (*mi).first;
    pComp = (*mi).second;

    // So far, no instance was found 
    foundInstance = false;
    if (_strand == plus_strand || _strand == both)
    {

      // get number of estimated motif instances
      //nbr = pComp->GetEstimatedNumberInstances(plus_strand);
      nbr = pComp->GetNumberInstances(plus_strand, 1);//
      if (nbr < 1)
      {

        // go to the next sequence
        mi++;
        continue;
      }

      // check if the estimated number of instances is larger than nbr
      if (nbr > (int) pVec.size())
        pVec.resize(nbr, -1);

      // select nbr best scoring start positions
      pComp->SelectBestInstanceStart(pVec, nbr, plus_strand);
      for (int j = 0; j < nbr; j++)
      {
        if (pVec[j] != -1)
        {

          // create instance
          pSite = new Instance(pSeq, plus_strand, pVec[j], _w);
          if (pSite != 0)
          {
              // cerr << "DEBUG sampled site: " << *(pSite->PrintSite()) << endl;
            pSite->SetScore(pComp->GetWxAt(pVec[j], plus_strand));
            _nbrInstances++;
            if (!foundInstance)
              _nbrSequencesWithInstances++;
            foundInstance = true;

            // append sites to InstanceMap
            _pInstanceMap->AddInstance(pSite);
          }
        }
      }
    }
    if (_strand == minus_strand || _strand == both)
    {

      // get number of estimated motif instances
      //nbr = pComp->GetEstimatedNumberInstances(minus_strand);
      nbr = pComp->GetNumberInstances(minus_strand, 1);//
      if (nbr < 1)
      {

        // go to the next sequence
        mi++;
        continue;
      }

      // check if the estimated number of instances is larger than nbr
      if (nbr > (int) pVec.size())
        pVec.resize(nbr, -1);

      // select nbr best start position
      pComp->SelectBestInstanceStart(pVec, nbr, minus_strand);
      for (int j = 0; j < nbr; j++)
      {
        if (pVec[j] != -1)
        {

          // create instance
          pSite = new Instance(pSeq, minus_strand, pVec[j], _w);
          if (pSite != 0)
          {
              // cerr << "DEBUG sampled site: " << *(pSite->PrintSite()) << endl;
            pSite->SetScore(pComp->GetWxAt(pVec[j], minus_strand));
            _nbrInstances++;
            if (!foundInstance)
              _nbrSequencesWithInstances++;
            foundInstance = true;

            // append sites to InstanceMap
            _pInstanceMap->AddInstance(pSite);
          }
        }
      }
    }

    // go to the next sequence in the set
    ++mi;
  }
  //cerr << "DEBUG: MotifSamplerRun::SelectBestInstanceMap - end" << endl;
  return;
}



/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
void
  MotifSamplerRun::SetBackgroundModel(BackgroundModel * bgM)
{
  // cerr << "DEBUG: MotifSamplerRun::SetBackgroundModel" << endl;

  // set pointer to input background model
  _pBgModel = bgM;

  // set pseudo motif 
  // cerr << "MotifSamplerRun::SetBackgroundModel => set pseudo counts: " << endl;
  // UpdatePseudoCounts(_pBgModel->GetSNF(), constant);
  UpdatePseudoCounts(_pBgModel->GetSNF());
  
  // cerr << "DEBUG: MotifSamplerRun::SetBackgroundModel" << endl;
  return;
}



/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
void
  MotifSamplerRun::UpdateBackgroundScores()
{
  // cerr << "DEBUG: MotifSamplerRun::UpdateBackgroundScores" << endl;

  // iterate through sequences 
  MapIterator seqIter = _pComputationMap->begin();
  while (seqIter != _pComputationMap->end())
  {
    // update background model score
    (*seqIter).second->UpdateInstanceBackgroundScore(_pBgModel, _w, _strand);
    (*seqIter).second->UpdateSequenceBackgroundScore(_pBgModel, _strand);
    seqIter++;
  }
  // cerr << "DEBUG: MotifSamplerRun::UpdateBackgroundScores" << endl;
  return;
}


/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
bool
  MotifSamplerRun::UpdateBackgroundScore(const string &pID, BackgroundModel *pBgM)
{
  // iterate through sequences 
  MapIterator seqIter = _pComputationMap->begin();
  bool found = false;
  while (seqIter != _pComputationMap->end())
  {
    // cerr << "Check sequence : " << *((*seqIter).first->GetID()) << endl;
    if ( ! pID.compare(*((*seqIter).first->GetID()) ) ){
      // cerr << "Found: " << pID << "|" << *((*seqIter).first->GetID()) << "|" << endl;
      // update background model score
      (*seqIter).second->UpdateInstanceBackgroundScore(pBgM, _w, _strand);
      (*seqIter).second->UpdateSequenceBackgroundScore(pBgM, _strand);
      found = true;
      break;
    }
    seqIter++;
  }
  // cerr << "DEBUG: MotifSamplerRun::UpdateBackgroundScores" << endl;
  return found;
}

/******************************************************************************
  Description:  LoadPspScores(string * pPspFile, string bPSPimpact)
                Read the scores from file (L or 2L entries) and linearize
  Date:         2013/12/30
  Author:       Marleen Claeys <mclaeys@qatar.net.qa>
  
******************************************************************************/
bool
  MotifSamplerRun::LoadPspScores(string * pPspFile, bool bPSPs)
{
  // cerr << "DEBUG: MotifSamplerRun::LoadPspScores" << endl;
  bool implemented = true;
  string error = "";
  string *id = NULL;
  SequenceComputation *pSeqComp = NULL;
  ScoreVector *pPSP = NULL;
  ostringstream cerrstr;
  int nbrPSPseq = 0;
  string trackSeq = "";

  // open fasta file
  cerr << "Open psp file for reading." << endl;
  FastaIO *fileIO = new FastaIO(pPspFile);
  if (fileIO == NULL)
  {
    cerr << "--Error-- MotifSamplerRun::LoadPspScores: "
      << "Unable to open psp file (-q)" << endl;
    cerrstr << "--Error-- MotifSamplerRun::LoadPspScores: "
      << "Unable to open psp file (-q)" << endl;
    _pgffio->AddComment(cerrstr.str());
    cerrstr.flush(); // flush     
    return false;
  }
  else
  { 
    while (fileIO->HasNext()) // pointer set at start of > line
    {
      // read the sequence id (must be same as previously loaded in FASTA)
      id = fileIO->ReadSeqID();
      if (id == NULL)
      { // incorrect format
        error = "--Error-- MotifSamplerRun::LoadPspScores: incorrect format for sequence identification (> lines).";
        implemented = false; break;
      }
      pSeqComp = FindSequence(id);
      if (pSeqComp == NULL)
      { // sequence id not found, skip this information
        error = "--Warning--MotifSamplerRun::LoadPspScores: sequence not found in FASTA file (-f). Skip sequence >";
        error.append(*id);
        cerr << error << endl;
        cerrstr << error << endl;
        _pgffio->AddComment(cerrstr.str());
        cerrstr.flush(); // flush     
        // set pointer ready on next >line
        fileIO->LoadPspData(0, 0,1);
        delete id; id = NULL;
        continue;
      }
      pPSP = fileIO->LoadPspData((pSeqComp->ParentSequence())->Length(), _w, 0); 
      if (pPSP == NULL)
      { // incorrect format
        error = "--Error-- MotifSamplerRun::LoadPspScores: incorrect format (check number of entries (L or 2L), white spaces, no kommas) for sequence >";
        error.append(*id);
        implemented = false; break;
      }
      pSeqComp->LoadPspScores(pPSP, _strand, bPSPs);
      nbrPSPseq++;
      trackSeq.append(">");trackSeq.append(*id);trackSeq.append(" ");
      delete id; id = NULL;
      pSeqComp = NULL;
      delete pPSP; pPSP = NULL;
    }
    if (!implemented)
    {
      cerr << error << endl;
      cerrstr << error << endl;
      _pgffio->AddComment(cerrstr.str());
      cerrstr.flush(); // flush     
    }
  }
  // close psp file
  delete fileIO;
  // cleanup
  if (id != NULL) delete id; id = NULL;
	
  // check if all fasta sequences have been covered
  if (implemented && nbrPSPseq != _nbrSequences)
  {
    char pTmp[128];
    string warning = "--Warning--MotifSamplerRun::LoadPspScores: The number of sequences in PSP file (";
    sprintf(pTmp, "%d", nbrPSPseq);
    warning.append(pTmp);
    warning.append(") does not equal the number of sequences in FASTA file (");
    sprintf(pTmp, "%d", _nbrSequences);
    warning.append(pTmp);
    warning.append("). PSP was stored for [");
    warning.append(trackSeq);
	warning.append("].");
    cerr << warning << endl;
    cerrstr << warning << endl;
    _pgffio->AddComment(cerrstr.str());
    cerrstr.flush(); // flush     
  }
  // cerr << "DEBUG: MotifSamplerRun::LoadPspScores" << endl;
  return implemented;
}

/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
void
  MotifSamplerRun::UpdateMotifScoresFromReducedInstanceMap()
{
  // cerr << "DEBUG: MotifSamplerRun::UpdateMotifScoresFromReducedInstanceMap" << endl;

  // loop through the sequence set
  MapIterator ci = _pComputationMap->begin();
  while (ci != _pComputationMap->end())
  {
    // create motif model ( model is stored in _localMotif)
    BuildMotifFromReducedInstanceMap((*ci).first);

    // score sequence
    (*ci).second->UpdateInstanceMotifScore(_localMotif, _strand);
    (*ci).second->UpdateInstanceExpWx(_strand);
    // next sequence
    ci++;
  }
  // cerr << "DEBUG: MotifSamplerRun::UpdateMotifScoresFromReducedInstanceMap" << endl;
  return;
}



/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
void
  MotifSamplerRun::StderrPrintMotifInfo()
{
  if (_nbrInstances == 0)

  {
    cerr << "Warning: No instances in InstanceMap.\n";
    return;
  }

  // build motif model
  BuildMotifFromInstanceMap(); // NO psp applicable
  if (_localMotif == NULL)

  {
    cerr << "Warning: Unable to build motif from current instance map." <<
      endl;
    return;
  }

  cerr << *(_localMotif->GetConsensus())
    << "\t" << _nbrSequencesWithInstances 
    << "\t" << _nbrInstances
    << "\tscore: " << _localMotif->ConsensusScore();
  if ( _pBgModel != NULL )
    cerr << "  " << _localMotif->InformationContent(_pBgModel->GetSNF());
  cerr  << "  " <<    LogLikelihoodScore() 
    << endl;

  return;
}



/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
void
  MotifSamplerRun::ShiftInstanceMap(int maxShift)
{
  // cerr << "DEBUG: MotifSamplerRun::ShiftInstanceMap" << endl;
  // first get score of latest motif from instance map
  BuildMotifFromInstanceMap(); // NO psp applicable
  // initially there is no shift
  double
    ic = _localMotif->InformationContent(_pBgModel->GetSNF());
  double
    ic2 = ic;
  int
    bestShift = 0,
    nt = -1; 
  //cerr << "  Shift Instance Map (init): " << bestShift << ": " << ic << endl;
  for (int shift = -maxShift; shift <= maxShift; shift++)

  {
    if (shift == 0)
      continue;

    // empty local counts
    for (int i = 0; i < _w; i++)
      for (int j = 0; j < 4; j++)
        _pLocalCounts[i][j] = 0;
    
    list < Instance * >::iterator ii = _pInstanceMap->begin();
    int
      start = 0;
    while (ii != _pInstanceMap->end())
    {

      // some local variables
      SequenceObject *
        pSeq = (*ii)->ParentSequence();
      strand_modes str = (*ii)->Strand();
      start = (*ii)->Start() + shift;
      if (start < 0)
      {
        ii++;
        continue;
      }
      if (start + _w > pSeq->Length())
      {
        ii++;
        continue;
      }
      for (int j = 0; j < _w; j++)
      {
        nt = pSeq->GetNucleotideAt(str, start + j);
        if (nt >= 0 && nt < 4)
          _pLocalCounts[j][nt] += 1;
      }
      // next instance
      ++ii;
      pSeq = NULL;
    }
    // create temporary matrix from count matrix
    PWM *
      tmpMatrix = new PWM(_w, _pLocalCounts, _pPseudoCounts);

    // get matrix scores
    ic = tmpMatrix->InformationContent(_pBgModel->GetSNF());

    // delete temporary matrixFile
    delete tmpMatrix;
    if (ic > ic2)
    {

      // reset bestShift and best score
      ic2 = ic;
      bestShift = shift;

      // update local motif
      _localMotif->RebuildMatrix(_w, _pLocalCounts, _pPseudoCounts);
    }

    // give some information on current shift
    //cerr << "  Shift Instance Map: " << shift << ": " << ic << " |  " <<
      //bestShift << ": " << ic2 << endl;
  }
  if (bestShift != 0)
  {

    // create a temporary instance map to store shifted instances
    InstanceMap *
      pTmp = new InstanceMap();

    // update all instances in the instance map
    list < Instance * >::iterator mIter = _pInstanceMap->begin();
    while (mIter != _pInstanceMap->end())
    {
      if ((*mIter)->Start() + bestShift < 0
          || (*mIter)->Start() + bestShift + _w >=
          (*mIter)->ParentSequence()->Length())
      {
        //cerr <<
         // "--Warning-- ShiftInstanceMap() shifted instance out of range: " <<
         // (*mIter)->Start() + bestShift << " - " << (*mIter)->Start() +
         // bestShift + _w << endl;
      }
      else
      {

        // create new instance
        // cerr << "ShiftInstanceMap()  make new instance" << endl;
        Instance *
          pInst = new Instance((*mIter)->ParentSequence(), (*mIter)->Strand(),
                               (*mIter)->Start() + bestShift, _w);

        // add new instance
        pTmp->AddInstance(pInst);
      }
      mIter++;
    }

    // delete old instance map
    delete _pInstanceMap;

    // refer to new instance map
    _pInstanceMap = pTmp;

    // reset pTmp
    pTmp = NULL;

    // create motif model from current instance map and update motif scores
    BuildMotifFromInstanceMap(); // no PSP needed here 
    // iterate through sequences 
    MapIterator seqIter = _pComputationMap->begin();
    while (seqIter != _pComputationMap->end())
    {
      if (_localMotif != NULL)
      {
        // score sequence with motif model
        (*seqIter).second->UpdateInstanceMotifScore(_localMotif, _strand);
        (*seqIter).second->UpdateInstanceExpWx(_strand); // mark : psp will apply in case of 's'
      }
      seqIter++;
    }
  }
  // cerr << "DEBUG: MotifSamplerRun::ShiftInstanceMap" << endl;
  return;
}

/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2012/09/19
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
		  Revised: Marleen Claeys <mclaeys@qatar.net.qa>
  
******************************************************************************/
void
  MotifSamplerRun::BuildMotifFromInstanceMap()
{
  // cerr << "DEBUG: MotifSamplerRun::BuildMotifFromInstanceMap" << endl;
  ostringstream cerrstr;
  // check instance map
  if (_pInstanceMap == NULL)
  {
    cerr <<
      "--Error-- MotifSamplerRun::MotifFromInstanceMap  Trying to make motif model from empty instance map, returning NULL Pointer."
      << endl;
    cerrstr << "--Error-- MotifSamplerRun::MotifFromInstanceMap  Trying to make motif model from empty instance map, returning NULL Pointer."
      << endl;
    _pgffio->AddComment(cerrstr.str());
    cerrstr.flush(); // flush

    // clear local motif
    if (_localMotif != NULL)
      delete _localMotif;
    _localMotif = NULL;
    return;
  }

  // create empty motif count
  // cerr << "DEBUG: initialize local count matrix" << endl;
  for (int i = 0; i < _w; i++)
    for (int j = 0; j < 4; j++)
      _pLocalCounts[i][j] = 0;

  // loop over all instances in the InstanceMap and add them to the count matrix
  // no corrective weighting by PSP here
  list < Instance * >::iterator iter = _pInstanceMap->begin();

  int
    nt;
  while (iter != _pInstanceMap->end())
  {
    // get site
    if ((*iter)->Site() != NULL)
    {
      for (int j = 0; j < _w; j++)
      {
        nt = (*((*iter)->Site()))[j];
        if (nt >= 0 && nt < 4)
          _pLocalCounts[j][nt] += 1;
      }
    }
    iter++;                     // next
  }

  // rebuild motif model from new counts
  if (_localMotif != NULL)
  {
    _localMotif->RebuildMatrix(_w, _pLocalCounts, _pPseudoCounts);
  }
  else
  {
    _localMotif = new PWM(_w, _pLocalCounts, _pPseudoCounts);
  }
  
  // cerr << "DEBUG: motif model = " << endl;
  // _localMotif->StderrPrintMatrix();
  
  // cerr << "DEBUG: MotifSamplerRun::BuildMotifFromInstanceMap" << endl;
  return;
}



/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2012/09/19
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
                Revised: Marleen Claeys <mclaeys@qatar.net.qa>
  
******************************************************************************/
void
  MotifSamplerRun::BuildMotifFromReducedInstanceMap(SequenceObject * pSeq)
{
  // cerr << "DEBUG: MotifSamplerRun::BuildMotifFromReducedInstanceMap" << endl;
  ostringstream cerrstr;
 
  // get id of sequence to be skipped in 'reduced'
  // 'dummmyy' is a work around for full instance map in convergence step, 
  // this assumes a seqid will never be named dummmyy
  string * pSeqId = new string("dummmyy"); 
  if (pSeq != NULL) {*pSeqId = *pSeq->GetID();}
  
  if (_pInstanceMap == NULL)
  {
    cerr <<
      "--Error-- MotifSamplerRun::MotifFromInstanceMap():  Trying to make motif model from empty instance map."
      << endl;
    cerrstr << "--Error-- MotifSamplerRun::MotifFromInstanceMap():  Trying to make motif model from empty instance map."
      << endl;
    _pgffio->AddComment(cerrstr.str());
    cerrstr.flush();
    if (_localMotif != NULL)
      delete _localMotif;
    _localMotif = NULL;
    return;
  }
  
  // create empty motif count
  // cerr << "DEBUG: initialize local count matrix" << endl;
  for (int i = 0; i < _w; i++)
    for (int j = 0; j < 4; j++)
      _pLocalCounts[i][j] = 0;

  // loop over all instances in the InstanceMap and add them to the count matrix
  // cerr << "DEBUG: Loop over all instances." << endl;
  double prob; // refers to psp priorization
	
  list < Instance * >::iterator iter = _pInstanceMap->begin();
  int
    nt;
  while (iter != _pInstanceMap->end())
  {
    if ((*iter) == NULL)
    {
      cerr << "--Warning-- MotifSamplerRun::BuildMotifFromReducedInstanceMap(): found empty instance." << endl;
      cerrstr << "--Warning-- MotifSamplerRun::BuildMotifFromReducedInstanceMap(): found empty instance." << endl;
      _pgffio->AddComment(cerrstr.str());
      cerrstr.flush();
      iter++;
      continue;
    }

    // get ID of ParentSequence
    string *
      pInstanceID = (*iter)->PrintSeqName();
    // cerr << "DEBUG: next instance = " << *pInstanceID << endl;
    // if (pInstanceID->compare(0, pInstanceID->size(), *pSeqId))
    if (pInstanceID->compare(*pSeqId))
    { // compare is 0 (thus skipped) if both strings are equal
      // cerr << "DEBUG selected instance " << *((*iter)->PrintSite())  << " from " << *pInstanceID << " - " << *pSeqId << endl;
      // get the PSP(x) of this instance
      prob = (FindSequence(((*iter)->ParentSequence())->GetID()))->GetPspScoreAt((*iter)->Start(), (*iter)->Strand());
      // get site
      if ((*iter)->Site() != NULL)
      {
        for (int j = 0; j < _w; j++)
        {
          nt = (*((*iter)->Site()))[j];
          if (nt >= 0 && nt < 4)
            _pLocalCounts[j][nt] += prob; // prob instead of '1'
        }
      }
    }
    iter++;                     // next
  }

  // rebuild motif model from new counts
  // cerr << "DEBUG: Rebuild _localMotif from pC." << endl;
  if (_localMotif != NULL)
  {
    _localMotif->RebuildMatrix(_w, _pLocalCounts, _pPseudoCounts);
  }
  else
  {
    _localMotif = new PWM(_w, _pLocalCounts, _pPseudoCounts);
  }

  //cerr << "DEBUG: motif model (reduced) = " << endl;
  //_localMotif->StderrPrintMatrix();
  // cleanup
  delete pSeqId; pSeqId = NULL;
  
  // cerr << "DEBUG: MotifSamplerRun::BuildMotifFromReducedInstanceMap" << endl;
  return;
}



/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>

******************************************************************************/
PWM *
MotifSamplerRun::GetMotifModel()
{
  // cerr << "DEBUG: MotifSamplerRun::GetMotifModel" << endl;

  // first get score of latest motif from instance map
  BuildMotifFromInstanceMap(); // no PSP applicable
  // cerr << "DEBUG: MotifSamplerRun::GetMotifModel" << endl;
  return _localMotif;
}



/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
    				// debug 02/01/2008 - MC
******************************************************************************/
void
  MotifSamplerRun::UpdateMasksFromInstanceMap()
{
  // cerr << "DEBUG: MotifSamplerRun::UpdateMasksFromInstanceMap" << endl;
  if (_pInstanceMap == NULL)
    return;
  SequenceObject *
    pSeq = NULL;
  ostringstream cerrstr;
  int
    s = 0,
    W = 0,
    L = 0;
  strand_modes str = plus_strand;
  list < Instance * >::iterator ii = _pInstanceMap->begin();
  MapIterator mi = _pComputationMap->begin();
  while (ii != _pInstanceMap->end())
  {

    // get sequence
    pSeq = (*ii)->ParentSequence();

    // find the computation object 
    mi = _pComputationMap->find(pSeq);
    if (mi != _pComputationMap->end())
    {
      // cerr << "DEBUG: Found sequence: " << *((*mi).first->GetID()) << endl;
      s = (*ii)->Start();
      W = (*ii)->Length();
      L = pSeq->Length();
      str = (*ii)->Strand();
      if (str == plus_strand)
      {
		// debug '+1' in start-argument
		// debug '-1' in length-argument
        (*mi).second->UpdateMask(s - W + _overlap +1, 2 * (W - _overlap) -1, 0,
                                 plus_strand);
		// debug '2*W' and '+1' in start-argument
		// debug '-1' in length-argument
        (*mi).second->UpdateMask(L - s - 2*W + _overlap +1, 
                                 2 * (W - _overlap) -1, 0, minus_strand);
      }
      else
      {
		// debug idem
        (*mi).second->UpdateMask(s - W + _overlap +1, 2 * (W - _overlap) -1, 0,
                                 minus_strand);
		// debug idem
        (*mi).second->UpdateMask(L - s - 2* W + _overlap +1, 
                                 2 * (W - _overlap) -1, 0, plus_strand);
      }
    }
    else
    {
      cerr <<
        "--Warning-- MotifSamplerRun::UpdateMasksFromInstanceMap(): not found in map, sequence "
        << *(pSeq->GetID()) << endl;
      cerrstr << "--Warning-- MotifSamplerRun::UpdateMasksFromInstanceMap(): not found in map, sequence "
        << *(pSeq->GetID()) << endl;
      _pgffio->AddComment(cerrstr.str());
      cerrstr.flush();
    }
    ii++;
  }
  // cerr << "DEBUG: MotifSamplerRun::UpdateMasksFromInstanceMap" << endl;
  return;
}



/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
void
MotifSamplerRun::ResetMasks()
{
  // cerr << "DEBUG: MotifSamplerRun::ResetMasks" << endl;

  // iterate through sequences 
  MapIterator seqIter = _pComputationMap->begin();
  while (seqIter != _pComputationMap->end())

  {

    // reset mask and set last W-1 positions to zero
    (*seqIter).second->GetMask(plus_strand)->ResetMask();
    (*seqIter).second->GetMask(plus_strand)->UpdateMask((*seqIter).first->
                                                        Length() - _w + 1,
                                                        _w - 1, 0);
    (*seqIter).second->GetMask(minus_strand)->ResetMask();
    (*seqIter).second->GetMask(minus_strand)->UpdateMask((*seqIter).first->
                                                         Length() - _w + 1,
                                                         _w - 1, 0);
    seqIter++;
  }
  // cerr << "DEBUG: MotifSamplerRun::ResetMasks" << endl;
  return;
}


/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
void
MotifSamplerRun::ResetMasks(int wNew)
{
  if ( wNew < 1 )
    wNew = 1;
  
  // iterate through sequences 
  MapIterator seqIter = _pComputationMap->begin();
  while (seqIter != _pComputationMap->end())
  {
    // reset mask and set last W-1 positions to zero
    (*seqIter).second->GetMask(plus_strand)->ResetMask();
    (*seqIter).second->GetMask(plus_strand)->UpdateMask((*seqIter).first->Length() - wNew + 1,wNew - 1, 0);
    (*seqIter).second->GetMask(minus_strand)->ResetMask();
    (*seqIter).second->GetMask(minus_strand)->UpdateMask((*seqIter).first->
                                                         Length() - wNew + 1,
                                                         wNew - 1, 0);
    seqIter++;
  }
  return;
}



/****************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
****************************************************************************/
double
  MotifSamplerRun::LogLikelihoodScore()
{
  //cerr << "DEBUG: MotifSamplerRun::LogLikelihoodScore - begin" << endl;
  double
    ll = 0;
  MapIterator mi = _pComputationMap->begin();
  list < Instance * >::iterator ii = _pInstanceMap->begin();
  vector < int >*
    pVec = NULL,
    *pRVec = NULL;
  string *
    pStr = NULL;
  while (mi != _pComputationMap->end())
  {
    pVec = new vector < int >;
    pRVec = new vector < int >;
    ii = _pInstanceMap->begin();
    while (ii != _pInstanceMap->end())
    {

      // check if instance belongs to sequence
      pStr = (*ii)->ParentSequence()->GetID();
      if (pStr->compare(*((*mi).first->GetID())))
      {
        ii++;
        continue;
      }
      if ((*ii)->Strand() == plus_strand)
      {
        pVec->push_back((*ii)->Start());
      }
      else if ((*ii)->Strand() == minus_strand)
      {
        pRVec->push_back((*ii)->Start());
      }
      ii++;
    }

    // add scores
    // version 3.1.2 : 'else' verwijderd zodat addititief igv _strand == both
    if ( _strand == both || _strand == plus_strand )
    {
      ll += (*mi).second->LogLikelihoodScore(pVec, plus_strand);
    }
    /*else */if ( _strand == both || _strand == minus_strand)
    {
      ll += (*mi).second->LogLikelihoodScore(pRVec, minus_strand);
    }
    
    // clear local vectors
    delete pVec;
    delete pRVec;
    mi++;
  }
//cerr << "ll-all=" << ll << endl;
  //cerr << "DEBUG: MotifSamplerRun::LogLikelihoodScore - end" << endl;
  return ll;
}



/****************************************************************************
  Method:       SetMotifLength
  Class:        MotifSamplerRun
  Arguments:    int wLength
  
  Description:  
  
  Date:         2013/05/14
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
****************************************************************************/
void
MotifSamplerRun::SetMotifLength(int wLength)
{
  //cerr << "DEBUG: MotifSamplerRun::SetMotifLength" << wLength << endl;
  ostringstream cerrstr;
  if (wLength <= 0)
  {
    cerr <<
      "--Error-- MotifSamplerRun::SetMotifLength(): Trying to set motif length to value smaller than 1"
      << endl;
      cerrstr << "--Error-- MotifSamplerRun::SetMotifLength(): Trying to set motif length to value smaller than 1"
      << endl;
      _pgffio->AddComment(cerrstr.str());
      cerrstr.flush();

    // reset motif model
    if ( _localMotif != NULL )
      delete _localMotif;
    _localMotif = NULL;
  }
  
  if ( wLength != _w )
  {
    // delete present motif
    for (int i=0; i<_w; i++)
      delete[] _pLocalCounts[i];
    delete _pLocalCounts;
    
    _w = wLength;
    
    _pLocalCounts = new double*[_w];  
    for (int i=0; i<_w; i++)
    {
      _pLocalCounts[i] = new double[4];
      for (int j=0; j<4; j++)
        _pLocalCounts[i][j] = 0;
    }
  }

  // set motif length in individual sequences
  string id = "";
  MapIterator si = _pComputationMap->begin();
  while (si != _pComputationMap->end())
  {
    if ((*si).first->Length() < _w)
    { id = *(*si).first->GetID();
cerr << "--Warning::MotifSamplerRun::SetMotifLength : remove sequence " <<
id  << " as too short length (L= " 
<< (*si).first->Length() << ") for motif detection (w=" << _w << ")." << endl;
      // remove this sequence from _pSequenceList
      SeqIterator i = _pSequenceList->begin();
      for (; i != _pSequenceList->end(); i++)
      { if ( *(*i)->GetID() == id)
        { _pSequenceList->erase(i);
          break;
        }
      }
      // delete sequenceObject and SequenceComputation
      delete (*si).first; delete (*si).second;
      // remove this element from _pComputationMap
      _pComputationMap->erase(si); 
      _nbrSequences--;
      // restart (avoid iterator jump issues)
      si = _pComputationMap->begin();
    }
    else 
    { (*si).second->SetMotifLength(_w); si++;}
  }

  //cerr << "DEBUG: MotifSamplerRun::SetMotifLength" << endl;
  return;
}
/******************************************************************************
  Description:  LoadNbrInstInfo : 
                read nbr inst/seq prior input and link in each 
                SequenceComputation object
                Set the bSampling or bEstimation mode
  Author_Date:  MC_2013/12/30
******************************************************************************/
bool 
MotifSamplerRun::LoadNbrInstInfo(int maxM, string * priorInput, bool sample)
{
//cerr << "debug--MotifSamplerRun::LoadNbrInstInfo - begin" << endl;

  // iterate through _pComputationMap to set the sampling mode 
  MapIterator si = _pComputationMap->begin();
  while (si != _pComputationMap->end())
  { (*si).second->SetNbrInstSampling(sample);si++;}
  si = _pComputationMap->begin();

  ScoreVector * pProbs = NULL;
  string::size_type pos;

  // case f (the prior should be used as final  _pCopyProbDistr)
	// compute values and store in SeqComp::_p(Rev)CopyProb
  if ((*priorInput)[0] == 'f')
  {
    priorInput->erase(0,1); // erase 'f'
    pProbs = new ScoreVector;
    double p;
    while (!priorInput->empty() && (int)pProbs->size() < maxM+1)
    { 
      pos = priorInput->find_first_not_of("0.123456789");
      istringstream istr(priorInput->substr(0,pos)); // convert to double
      istr >> p;
      if (p < 0 || p > 1) // invalid prob
      { delete pProbs; return false;} // invoke to exit
      priorInput->erase(0,pos);
      while ( priorInput->find_first_of("_ \t") == 0) priorInput->erase(0,1); 
      pProbs->push_back(p);
    }
    // cut low values at the end
    for(int i = (int)pProbs->size()-1; i >= 0; i--) 
    { if ( (*pProbs)[i] <= 0.001) pProbs->pop_back();
      else break;
    }
    // make a check on the total entry
    double sum = 0;
    for (int i = 0; i < (int)pProbs->size(); i++)
    { sum += (*pProbs)[i];}
    if (sum < 0.001)
    { delete pProbs; return false;} // invoke to exit
    // copy scorevector into SeqComp
    while (si != _pComputationMap->end())
    { (*si).second->LoadNbrInstFixed(pProbs); si++;}
    // cleanup
    delete pProbs; pProbs = NULL;
  }
  else
  { // create all prior-copy-distributions for c=1,2,...,-M
    // or reset -M in case Pr(c) becomes too small

    // create next distributions
    char c = (*priorInput)[0]; // 0,u,b,e
    switch (c)
    {
      case '0' :
      { // read prior p
        double p; // first entry of X_X
        pos = priorInput->find_first_not_of("0.123456789");
        istringstream istr(priorInput->substr(0,pos)); // convert to double
        istr >> p;
        if (p <= 0 || p >= 1) // invalid prior
        { return false;} // invoke to exit
        priorInput->erase(0,pos);
        while ( priorInput->find_first_of("_ \t") == 0) priorInput->erase(0,1); 
        // read kappa 
        double kappa; 
        if ( !priorInput->empty()) 
        { 
          pos = priorInput->find_first_not_of("0.123456789");
          istringstream istr(priorInput->substr(0,pos)); // convert to double
          istr >> kappa;
          if (kappa < 0 || kappa > 1) // invalid kappa
          { return false;} // invoke to exit
          priorInput->erase(0,pos);
          while ( priorInput->find_first_of("_ \t") == 0) priorInput->erase(0,1); 
        }
        if (!priorInput->empty()) 
        { return false;} // invoke to exit 
        // construct distributions
        pProbs = new ScoreVector;
        double value, sum(1);
        pProbs->push_back(1 - p); // Pr(0)
        pProbs->push_back(p); // Pr(1)
        for(int i = 2; i <= maxM; i++) 
        { value = kappa * (*pProbs)[i-1];
          sum += value;
          if ( value/sum <= 0.001) break; // untill max or <0.001
          pProbs->push_back(value);
        }
        // create local vector to store distributions 
        _pPriorDistributions = new vector<Distribution *>;
        for (int n = 1; n <= (int)pProbs->size() -1; n++)
          _pPriorDistributions->push_back(new Distribution(pProbs, 0, n+1));
        delete pProbs; pProbs = NULL;
        break;
      }
      case 'u' : 
      {
        priorInput->erase(0,1); // erase 'u'
        if ( !priorInput->empty()) break;
        // create and store the distributions
        pProbs = new ScoreVector; pProbs->push_back(1);
        // create local vector to store distributions 
        _pPriorDistributions = new vector<Distribution *>;
        for (int n = 1; n <= maxM; n++) 
        { pProbs->push_back(1);
          _pPriorDistributions->push_back(new Distribution(pProbs, 0, n+1));
        }
        delete pProbs; pProbs = NULL;
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
        { return false;} // invoke to exit 
        priorInput->erase(0,pos);
        if ( !priorInput->empty()) break;
        // create local vector to store distributions 
        _pPriorDistributions = new vector<Distribution *>;
        for (int n = 1; n <= maxM; n++)
        { pProbs = new ScoreVector; 
          for(int i = 0; i <= n; i++)
            pProbs->push_back(pow(p,i)*pow((1-p),(n-i))*
                   INCLUSIVE::fac(n)/INCLUSIVE::fac(n-i)/INCLUSIVE::fac(i));
          for(int i = n; i >= 0; i--) // cut low values at the end
          { if ( (*pProbs)[i] <= 0.001) pProbs->pop_back();
            else break;
          }
          if ( (int)pProbs->size() == n + 1 )
          { _pPriorDistributions->push_back(new Distribution(pProbs, 0, n+1));
            delete pProbs; pProbs = NULL; 
          }
          else 
          { delete pProbs; pProbs = NULL;
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
        while (!priorInput->empty() && (int)pProbs->size() < maxM+1)
        { 
          pos = priorInput->find_first_not_of("0.123456789");
          istringstream istr(priorInput->substr(0,pos)); // convert to double
          istr >> p;
          if (p < 0 || p >= 1) // invalid prob
          { return false;} // invoke to exit 
          priorInput->erase(0,pos);
          while ( priorInput->find_first_of("_ \t") == 0) priorInput->erase(0,1); 
          pProbs->push_back(p);
        }
        // cut low values at the end
        for(int i = (int)pProbs->size()-1; i >= 0; i--) 
        { if ( (*pProbs)[i] <= 0.001) { pProbs->pop_back();}
          else break;
        }
        // make a check on the total entry
        double sum = 0;
        for (int i = 0; i < (int)pProbs->size(); i++)
        { sum += (*pProbs)[i];}
        if (sum < 0.001)
        { delete pProbs; return false;} // invoke to exit
        // create local vector to store distributions 
        _pPriorDistributions = new vector<Distribution *>;
        for (int n = 1; n <= (int)pProbs->size() - 1; n++)
        { 
          _pPriorDistributions->push_back(new Distribution(pProbs, 0, n+1)); 
        }
        delete pProbs; pProbs = NULL;
        break;
      }
    }

    // link to _priorDistrs in SequenceComputation
    while (si != _pComputationMap->end())
    { (*si).second->LinkNbrInstPrior(_pPriorDistributions); si++;}

  }

  //cerr << "debug--MotifSamplerRun::LoadNbrInstInfo - end" << endl;
  return true;
}

/****************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
****************************************************************************/
void
  MotifSamplerRun::UpdateCopyProbabilityDistribution()
{
  // cerr << "DEBUG: MotifSamplerRun::UpdateCopyProbabilityDistribution" << endl;
  MapIterator mi = _pComputationMap->begin();
  while (mi != _pComputationMap->end())
  {

    // update 
    //(*mi).second->UpdateCopyProbability(_prior, _strand);
    (*mi).second->UpdateCopyProbability(_strand);//

    // go to the next iteration
    mi++;
  }
  // cerr << "DEBUG: MotifSamplerRun::UpdateCopyProbabilityDistribution" << endl;
  return;
}



/****************************************************************************
  Method:       InitializationStep
  Class:        MotifSamplerRun
  Arguments:    int iterations
  
  Description:  initialization step of the MotifSampler loop
								in this step exactly one instance per sequence is selected
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
****************************************************************************/
void
  MotifSamplerRun::InitializationStep(int iterations)
{
  //cerr << "DEBUG: MotifSamplerRun::InitializationStep - begin" << endl;
  //cerr << "Initialize Instance Map" << endl;
  InitFixedSizeInstanceMap(1);
  //StderrPrintInstanceMap();
  
  //cerr << "Starting initial 1 copy step" << endl;
  for (int i = 0; i < iterations; i++)
  {
    UpdateMotifScoresFromReducedInstanceMap();

    // sample exactly one instance in each sequence 
    SampleFixedSizeInstanceMap(1);

  }
  //cerr << "DEBUG: MotifSamplerRun::InitializationStep -end" << endl;
  return;
}



/****************************************************************************
  Method:       CoreSamplingStep
  Class:        MotifSamplerRun
  Arguments:    int iterations, int shiftTime, int maxShift
  
  Description:  core sampling loop of the basic MotifSampler.
 This includes a step in which instances are shifted to 
                find better instances
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
****************************************************************************/
void
  MotifSamplerRun::CoreSamplingStep(int iterations, int shiftTime,
                                    int maxShift)
{
  MapIterator mi = _pComputationMap->begin();
  SequenceObject *pSeq = NULL;
  SequenceComputation *pComp = NULL;
  ostringstream cerrstr;

  // cerr << "DEBUG: MotifSamplerRun::CoreSamplingStep" << endl;
  for (int i = 0; i < iterations; i++)
  {
    if ( maxShift )
    {
      if ((i % shiftTime) == 0)
      {
        //cerr << "Shift Instances: " << endl;
        ShiftInstanceMap(maxShift);
      }
    }
    // loop over all sequences and update scores
    mi = _pComputationMap->begin();
    while (mi != _pComputationMap->end())
    {
      pSeq = (*mi).first;
      pComp = (*mi).second;

      // build motif model from instances in instance map excluding current sequence 
      BuildMotifFromReducedInstanceMap(pSeq);
      // update motif scores
      pComp->UpdateInstanceMotifScore(_localMotif, _strand);
      pComp->UpdateInstanceExpWx(_strand);
      // compute distribution to estimate number of instances
      //pComp->UpdateCopyProbability(_prior, _strand);
      pComp->UpdateCopyProbability(_strand);//
      // next sequence
      mi++;
    }

    // sample a new instance map
    SampleInstanceMap();
    // check number of instances
    if (_nbrInstances < 2 || _nbrSequencesWithInstances < 2)
    {
      //cerr << "--Warning-- not enough instances to proceed procedure." <<
       // endl;
      pSeq = NULL;
      pComp = NULL;
      break;
    }

    pSeq = NULL;
    pComp = NULL;
  }
  //cerr << "DEBUG: MotifSamplerRun::CoreSamplingStep - end" << endl;
  return;
}



/****************************************************************************
  Method:       ConvergenceStep
  Class:        MotifSamplerRun
  Arguments:    int iterations
  
  Description:  convergence loop of the basic MotifSampler algorithm
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
****************************************************************************/
void
  MotifSamplerRun::ConvergenceStep(int iterations)
{
  //cerr << "DEBUG: MotifSamplerRun::ConvergenceStep-begin" << endl;
  MapIterator
    mi = _pComputationMap->begin();
  SequenceObject *
    pSeq = NULL;
  SequenceComputation *
    pComp = NULL;
  ostringstream cerrstr;
  SequenceObject * pdummy = NULL;

  // cerr << "DEBUG: MotifSamplerRun::ConvergenceStep" << endl;
  for (int i = 0; i < iterations; i++)
  {
    // loop over all sequences
    mi = _pComputationMap->begin();
    // build motif model from all instances
    BuildMotifFromReducedInstanceMap(pdummy); // gibbs iteration step
    while (mi != _pComputationMap->end())
    {
      pSeq = (*mi).first;
      pComp = (*mi).second;

      // update motif scores
      pComp->UpdateInstanceMotifScore(_localMotif, _strand);
      pComp->UpdateInstanceExpWx(_strand);
      // compute distribution to estimate number of instances
      //pComp->UpdateCopyProbability(_prior, _strand);
      pComp->UpdateCopyProbability(_strand); //
      // next sequence
      mi++;
    }
    // create new instance map from curent scores
    SelectBestInstanceMap();

    // check number of instances
    if (_nbrInstances < 2 || _nbrSequencesWithInstances < 2)
    {
      //cerr << "--Warning-- not enough instances to proceed procedure." <<
       // endl;
      pSeq = NULL;
      pComp = NULL;
      break;
    }
    pSeq = NULL;
    pComp = NULL;
  }
  BuildMotifFromInstanceMap(); // no PSP applicable
  //StderrPrintInstanceMap();
  
  //cerr << "DEBUG: MotifSamplerRun::ConvergenceStep-end" << endl;
  return;
}

/****************************************************************************
  Method:       PrintInstanceMap
  Class:        MotifSamplerRun
  Arguments:    GFFWriter * pGFFout, string * pSource
  
  Description:  print the current instances in GFF format
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
****************************************************************************/
void
  MotifSamplerRun::PrintInstanceMap(GFFWriter * pGFFout, string * pSource)
{
  // cerr << "DEBUG: MotifSamplerRun::PrintInstanceMap" << endl;
  list < Instance * >::iterator ii = _pInstanceMap->begin();
  while (ii != _pInstanceMap->end())
  {
    if ((*ii) != NULL)
    {
      (*ii)->SetID(_pId);
      pGFFout->WriteInstance(*ii, pSource);
    }
    ii++;
  }
  // cerr << "DEBUG: MotifSamplerRun::PrintInstanceMap" << endl;
  return;
}


/******************************************************************************
  Method:       StderrPrintInstanceMap
  Class:        MotifSamplerRun
  Arguments:    none
  
  Description:  pretty print the current list of instances on STDERR 
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
void 
  MotifSamplerRun::StderrPrintInstanceMap()
{
  cerr << "\n-----------------------------------------------------------\n";
  list < Instance * >::iterator ii = _pInstanceMap->begin();
  while (ii != _pInstanceMap->end())
  {
    if ((*ii) != NULL)
    {
      cerr << "\t" << *((*ii)->PrintSeqName())
      << "\t" << *((*ii)->PrintSite())
      << "\t" << (*ii)->Start()
      << "\t" << (*ii)->Strand() << endl;
    }
    ii++;
  }
  cerr << "-----------------------------------------------------------\n\n";
  return;
}

/******************************************************************************
  Method:       ExtendLeftInstanceMap
  Class:        MotifSamplerRun
  Arguments:    int pos
  
  Description:  extend all instances pos positions to the left
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
void
MotifSamplerRun::ExtendLeftInstanceMap(int pos)
{
  if ( pos > 0 )
  {
    _pInstanceMap->ExtendLeft(pos);
  }
  return;
}

  
/******************************************************************************
  Method:       ExtendRightInstanceMap
  Class:        MotifSamplerRun
  Arguments:    int pos
  
  Description:  extend all instances pos positions to the right
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
void
MotifSamplerRun::ExtendRightInstanceMap(int pos)
{
  if ( pos > 0 )
  {
    _pInstanceMap->ExtendRight(pos);
  }
  return;
}


/******************************************************************************
  Method:       CheckLeftExtension
  Class:        MotifSamplerRun
  Arguments:    double threshold
  
  Description:  check if the motif can be extended one position to the left
								the consensus score of the position should be greater 
                than the threshold
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
int
MotifSamplerRun::CheckLeftExtension(double threshold)
{
  double score = 2.0;
  int pos = 1,
    b = -1;
  list < Instance * >::iterator lIter;
  bool getOut = false;
  vector<double> vACGT(4);

  while ( score > threshold ){
    // initialize counts
    for (int i=0; i<4; i++)
      vACGT[i] = 0.05;

    // loop through all sites and select next positions
    lIter = _pInstanceMap->begin();
    while ( lIter != _pInstanceMap->end() )
    {
      // get start, length and strand of this instance
      int sstart = (*lIter)->Start();
      strand_modes sstrand =  (*lIter)->Strand();

      // int seqLength = ((*lIter)->ParentSequence())->Length();
      if ( sstart - pos >= 0 )
      {
        // add letter to counts
        b = ((*lIter)->ParentSequence())->GetNucleotideAt(sstrand, sstart - pos);
        // cerr << b;
        if ( b != -1 )
        {
          vACGT[b]++;
        }
        else
        {
          // this site false outside the range
          getOut = true;
          break;
        }

      }else{
        // this site false outside the range
        getOut = true;
        break;
      }
      
      lIter++;
    }
    // check if the iterations were stopped because the instances fall outside the range
    if ( getOut )
    {
      pos--;
      return pos;
    }

    // compute score from counts
    double sum = 0;
    for (int i=0; i<4; i++)
      sum += vACGT[i];
    for (int i=0; i<4; i++)
      vACGT[i] = vACGT[i]/sum;
    score = 2.0;
    for (int i=0; i<4; i++)
      score += vACGT[i] * (log(vACGT[i]) / log(2.0));
    
    // cerr << "-" << pos << " --> CS: " << score << endl;
      
    if ( score > threshold ){
      pos++;
    }else{
      pos--;
      return pos;
    }      
  }
  
  return pos;
}


/******************************************************************************
  Method:       CheckRightExtension
  Class:        MotifSamplerRun
  Arguments:    double threshold
  
  Description:  check if the motif can be extended one position to the right
								the consensus score of the position should be greater 
                than the threshold
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
int
MotifSamplerRun::CheckRightExtension(double threshold)
{
  double score = 2.0;
  int pos = 1,
    b = -1;
  list < Instance * >::iterator lIter;
  bool getOut = false;
  vector<double> vACGT(4);

  while ( score > threshold ){
    // initialize counts
    for (int i=0; i<4; i++)
      vACGT[i] = 0.05;

    // loop through all sites and select next positions
    lIter = _pInstanceMap->begin();
    while ( lIter != _pInstanceMap->end() )
    {
      // get start, length and strand of this instance
      int sstart = (*lIter)->Start();
      int slength = (*lIter)->Length();
      strand_modes sstrand =  (*lIter)->Strand();

      int seqLength = ((*lIter)->ParentSequence())->Length();
      if ( sstart + slength + pos < seqLength )
      {
        // add letter to counts
        b = ((*lIter)->ParentSequence())->GetNucleotideAt(sstrand, sstart + slength + pos);
        // cerr << b;
        if ( b != -1 )
        {
          vACGT[b]++;
        }
        else
        {
          // this site false outside the range
          getOut = true;
          break;
        }
      }else{
        // this site false outside the range
        getOut = true;
        break;
      }
      
      lIter++;
    }
    // check if the iterations were stopped because the instances fall outside the range
    if ( getOut )
    {
      pos--;
      return pos;
    }

    // compute score from counts
    // normalize
    double sum = 0;
    for (int i=0; i<4; i++)
      sum += vACGT[i];   
    for (int i=0; i<4; i++)
      vACGT[i] = vACGT[i]/sum;
    score = 2.0;
    for (int i=0; i<4; i++)
      score += vACGT[i] * (log(vACGT[i]) / log(2.0));
    
    // cerr << "+" << pos << " --> CS: " << score << endl;
      
    if ( score > threshold ){
      pos++;
    }else{
      pos--;
      return pos;
    }      
  }
  
  return pos;
}



/******************************************************************************
  Method:       CheckLeftExtension
  Class:        MotifSamplerRun
  Arguments:    double threshold, int window
  
  Description:  check if the motif can be extended over a given window
								the average consensus score of the window should be 
                greater than the threshold
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
int
MotifSamplerRun::CheckLeftExtension(double threshold, int window)
{
  double score = 2.0;
  int pos = 1,
    b = -1,
    w = 0,
    i = 0;
  list < Instance * >::iterator lIter;
  bool getOut = false;
  double **vACGT;

  // initialize counts
  vACGT = new double*[window];
  for (w=0; w<window; w++)
  {
    vACGT[w] = new double[4];
    for (i=0; i<4; i++)
      vACGT[w][i] = 0.05;
  }

  while ( !getOut && score > threshold ){
    // reinitialize counts
    for (w=0; w<window; w++)
      for (i=0; i<4; i++)
        vACGT[w][i] = 0.05;
      
    // loop through all sites and select next positions
    lIter = _pInstanceMap->begin();
    while ( !getOut && lIter != _pInstanceMap->end() )
    {
      // get start, length and strand of this instance
      int sstart = (*lIter)->Start();
      strand_modes sstrand =  (*lIter)->Strand();

      for (w=0; w<window; w++)
      {
        // int seqLength = ((*lIter)->ParentSequence())->Length();
        if ( sstart - pos - w >= 0 )
        {
          // add letter to counts
          b = ((*lIter)->ParentSequence())->GetNucleotideAt(sstrand, sstart - pos - w);
          // cerr << b;
          if ( b != -1 )
          {
            vACGT[w][b]++;
          }
          else
          {
            // this site falls outside the range
            getOut = true;
            break;
          }
        }else{
          // this site falls outside the range
          getOut = true;
          break;
        }
      }
      
      lIter++;
    }
    
    // check if the iterations were stopped because the instances fall outside the range
    if ( getOut )
    {
      pos--;
    }
    else
    {
      // compute score from counts
      double sum = 0;
      score = 0;
      for ( int w=0; w<window; w++)
      {
        // normalize scores
        sum = 0;
        for (int i=0; i<4; i++)
          sum += vACGT[w][i];
        score += 2.0;
        for (int i=0; i<4; i++)
          score += (vACGT[w][i]/sum) * ((log(vACGT[w][i]) - log(sum)) / log(2.0));
      }
      // average over window
      score = score / window;
      
      // cerr << "-" << pos << " --> CS: " << score << endl;
        
      if ( score > threshold ){
        pos++;
      }else{
        pos--;
        getOut = true;
      }
    }    
  }

  // clean up local variables
  for (w=0; w<window; w++)
    delete[] vACGT[w];
  delete[] vACGT;
  
  return pos;
}


/******************************************************************************
  Method:       CheckRightExtension
  Class:        MotifSamplerRun 
  Arguments:    double threshold, int window
  
  Description:  check if the motif can be extended over a given window
                the average consensus score of the window should be 
                greater than the threshold
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
int
MotifSamplerRun::CheckRightExtension(double threshold, int window)
{
  double score = 2.0;
  int pos = 1,
    b = -1,
    i = 0,
    w = 0;
  list < Instance * >::iterator lIter;
  bool getOut = false;
  double **vACGT;

  // initialize counts
  vACGT = new double*[window];
  for (w=0; w<window; w++)
  {
    vACGT[w] = new double[4];
    for (i=0; i<4; i++)
      vACGT[w][i] = 0.05;
  }

  while ( !getOut && score > threshold ){
  
    // reinitialize counts
  for (w=0; w<window; w++)
    for (int i=0; i<4; i++)
      vACGT[w][i] = 0.05;

    // loop through all sites and select next positions
    lIter = _pInstanceMap->begin();
    while ( !getOut && lIter != _pInstanceMap->end() )
    {
      // get start, length and strand of this instance
      int sstart = (*lIter)->Start();
      int slength = (*lIter)->Length();
      strand_modes sstrand =  (*lIter)->Strand();
      int seqLength = ((*lIter)->ParentSequence())->Length();
      
      for (w=0; w<window; w++)
      {
        
        if ( sstart + slength + pos + w < seqLength )
        {
          // add letter to counts
          b = ((*lIter)->ParentSequence())->GetNucleotideAt(sstrand, sstart + slength + pos + w);
          // cerr << b;
          if ( b != -1 )
          {
            vACGT[w][b]++;
          }
          else
          {
            // this site false outside the range
            getOut = true;
            break;
          }
        }else{
          // this site false outside the range
          getOut = true;
          break;
        }
      }
      lIter++;
    }
    
    // check if the iterations were stopped because the instances fall outside the range
    if ( getOut )
    {
      pos--;
    }
    else
    {
      // compute score from counts
      double sum = 0;
      score = 0;
      for ( int w=0; w<window; w++)
      {
        // normalize scores
        sum = 0;
        for (int i=0; i<4; i++)
          sum += vACGT[w][i];
        score += 2.0;
        for (int i=0; i<4; i++)
          score += (vACGT[w][i]/sum) * ((log(vACGT[w][i]) - log(sum)) / log(2.0));
      }
      // average over window
      score = score / window;
      
      // cerr << "-" << pos << " --> CS: " << score << endl;
        
      if ( score > threshold ){
        pos++;
      }else{
        pos--;
        getOut = true;
      }
    }
  }
  
  // clean up local variables
  for (w=0; w<window; w++)
    delete[] vACGT[w];
  delete[] vACGT;
  
  return pos;
}
/******************************************************************************
  Description:  check if 'id' is present in SeqObj of _pComputationMap 
                if so return SequenceComputation
  Date:         2012/09/19
  Author:       Marleen Claeys <mclaeys@qatar.net.qa>
  
******************************************************************************/
SequenceComputation *
  MotifSamplerRun::FindSequence(string * id)
{
  // cerr << "DEBUG: MotifSamplerRun::FindSequence()" << endl;
  SequenceComputation * comp = NULL;

  // iterate through sequences 
  MapIterator seqIter = _pComputationMap->begin();
  while (seqIter != _pComputationMap->end())
  {
    if ((*id) == (*(*seqIter).first->GetID()))
    {
      // found, and assign comp-object
      comp = (*seqIter).second;
      break;
    }
    seqIter++;
  }
  // cerr << "DEBUG: MotifSamplerRun::FindSequence() -- end" << endl;
  return comp;
}
