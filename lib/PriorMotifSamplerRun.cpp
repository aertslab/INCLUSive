#include "PriorMotifSamplerRun.h" 
#include "inclusive.h"
#include "utilities.h"

using namespace INCLUSIVE;

PriorMotifSamplerRun::PriorMotifSamplerRun(string* pFasta, strand_modes strand, PWM* pPriorMotif)
  : MotifSamplerRun(pFasta, strand) 
{

  // set length of the search equal to the length of the prior motif
  // and let SetMotifLength initiate the local variables
  SetMotifLength(pPriorMotif->Length());  

  // define mixture coefficients with initial value set to 0.5
  _pMixtureCoefficients = new vector<double>(_w, 0.5);

  // set up default scaling factor of the pseudo counts
  _scale = sqrt((double)_nbrSequences);

  // initialize prior motif model
  _pPriorMotifModel = new double[_w][4];
  for (int i=0; i<_w; i++)
    for (int j=0; j<4; j++)
      _pPriorMotifModel[i][j] = pPriorMotif->GetValueAt(i,j);
    
}

// destructor
PriorMotifSamplerRun::~PriorMotifSamplerRun() 
{
  delete[] _pPriorMotifModel;
    _pPriorMotifModel= NULL;
  
  if ( _pMixtureCoefficients != NULL )
    delete _pMixtureCoefficients;
  _pMixtureCoefficients = NULL;
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
  PriorMotifSamplerRun::ResetMixtureCoefficients()
{
  for ( uint i=0; i<_pMixtureCoefficients->size(); i++)
    (*_pMixtureCoefficients)[i] = 0.5;
  
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
  PriorMotifSamplerRun::UpdateMixtureCoefficients()
{
  double 
    z1 = 0,
    z2 = 0,
    q = 0;
  // local variables to store the arrays
  vector<double> nainfo(4,0); 
  vector<double> ainfo(4,0); 
  vector<double> naun(4,0); 
  vector<double> aun(4,0); 

  for (int i=0; i<_w; i++)
  {
    // get current value of the mixture coefficients
    q = (*_pMixtureCoefficients)[i];
    // define the correct arrays
    for (int j=0; j<4; j++)
    {
      nainfo[j] = _pLocalCounts[i][j] + (_scale * _pPriorMotifModel[i][j]);
      ainfo[j] = _scale * _pPriorMotifModel[i][j];
      naun[j] = _pLocalCounts[i][j] + (_scale * _pPseudoCounts[j]);
      aun[j] = _scale * _pPseudoCounts[j];
    }
    // get the scores
    z1 = LogNormDirichlet(&nainfo) - LogNormDirichlet(&ainfo);
    z2 = LogNormDirichlet(&naun) - LogNormDirichlet(&aun);
  
    // update _pMixtureCoefficients
    q = (q * exp(z1))/((q * exp(z1)) + ((1 - q) * exp(z2)));
    if ( q > 0 && q < 1 )
      (*_pMixtureCoefficients)[i] = q;
    else if ( q <= 0 )
      (*_pMixtureCoefficients)[i] = 0;
    else 
      (*_pMixtureCoefficients)[i] = 1;
    
  }
  
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
  PriorMotifSamplerRun::SetNewPriorMotifModel(PWM* pPriorMotif)
{
  // check length and change _w if necessary
  if ( pPriorMotif->Length() != _w ){
    // use SetMotifLength to re-initiate local variables
    SetMotifLength(pPriorMotif->Length());
    
    // recreate _pPriorMotifModel
    if ( _pPriorMotifModel != NULL )
      delete[] _pPriorMotifModel;
    _pPriorMotifModel = new double[_w][4];

    if ( _pMixtureCoefficients != NULL )
      _pMixtureCoefficients->resize(_w);
    else
      _pMixtureCoefficients = new vector<double>(_w, 0.5);

    // update background scores
    UpdateBackgroundScores();
    
  }
  
  // fill in new values
  for (int i=0; i<_w; i++)
  {
    for (int j=0; j<4; j++)
      _pPriorMotifModel[i][j] = pPriorMotif->GetValueAt(i,j);

    (*_pMixtureCoefficients)[i] = 0.5;    
  }
  
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
  PriorMotifSamplerRun::BuildMotifFromInstanceMap()
{
  // cerr << "DEBUG: MotifSamplerRun::BuildMotifFromInstanceMap" << endl;

  // check instance map
  if (_pInstanceMap == NULL)
  {
    cerr <<
      "--Error-- MotifSamplerRun::MotifFromInstanceMap  Trying to make motif model from empty instance map, returning NULL Pointer."
      << endl;

    // clear local motif
    if (_localMotif != NULL)
      delete _localMotif;
    _localMotif = NULL;
    return;
  }

  // create new motif count 
  for (int i = 0; i < _w; i++)
    for (int j = 0; j < 4; j++)
      _pLocalCounts[i][j] = 0;

    // loop over all instances in the InstanceMap and add them to the the count matrix
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

  // update mixture coefficients
  UpdateMixtureCoefficients();
  
  // add _pseudoMotif
  Matrix pC = new double[_w][4];
  for (int i = 0; i < _w; i++)
  {
    // cerr << "Pos " << i << ": ";
    for (int j = 0; j < 4; j++)
    {
      pC[i][j] = ((*_pMixtureCoefficients)[i] * ( _pLocalCounts[i][j] + _scale * _pPriorMotifModel[i][j])) 
        + ((1 - (*_pMixtureCoefficients)[i]) * ( _pLocalCounts[i][j] + _scale * _pPseudoCounts[j]));
  
      // cerr << _pLocalCounts[i][j] << "+" << _scale * _pPriorMotifModel[i][j] << "+" << _scale * _pPseudoCounts[j] << "=" << pC[i][j] << "|  ";
    }
    // cerr << endl;
  }
  
  // rebuild motif model from new counts
  if (_localMotif != NULL)
  {
    _localMotif->RebuildMatrix(_w, pC);
  }
  else
  {
    _localMotif = new PWM(_w, pC);
  }

  // _localMotif->StderrPrintMatrix();
  
  // clear local variable
  delete[]pC;

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
  PriorMotifSamplerRun::BuildMotifFromReducedInstanceMap(SequenceObject * pSeq)
{
  // cerr << "DEBUG: MotifSamplerRun::BuildMotifFromReducedInstanceMap" << endl;

  // get id of sequence
  string *
    pSeqId = pSeq->GetID();
  if (_pInstanceMap == NULL)
  {
    cerr <<
      "--Error-- MotifSamplerRun::MotifFromInstanceMap():  Trying to make motif model from empty instance map."
      << endl;
    if (_localMotif != NULL)
      delete _localMotif;
    _localMotif = NULL;
    return;
  }

  // create empty motif count
  for (int i = 0; i < _w; i++)
    for (int j = 0; j < 4; j++)
      _pLocalCounts[i][j] = 0;

  // loop over all instances in the InstanceMap and add them to the the count matrix
  // cerr << "Loop over all instances." << endl;
  list < Instance * >::iterator iter = _pInstanceMap->begin();
  int
    nt;
  while (iter != _pInstanceMap->end())
  {
    if ((*iter) == NULL)
    {
      iter++;
      continue;
    }

    // get ID of ParentSequence
    string *
      pInstanceID = (*iter)->ParentSequence()->GetID();
    if (pInstanceID->compare(0, pInstanceID->size(), *pSeqId))
    {                           // is 0 if both strings are equal
      // cerr << "DEBUG next instance " << *((*iter)->PrintSite())  << " from " << *pInstanceID << " - " << (*iter)->Start() << endl;
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
    }
    iter++;                     // next
  }

  // update mixture coefficients
  UpdateMixtureCoefficients();
  
  // add _pseudoMotif 
  Matrix pC = new double[_w][4];
  for (int i = 0; i < _w; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      pC[i][j] = ((*_pMixtureCoefficients)[i] * (_pLocalCounts[i][j] + _scale * _pPriorMotifModel[i][j])) 
        + ((1 - (*_pMixtureCoefficients)[i]) * (_pLocalCounts[i][j] + _scale * _pPseudoCounts[j]));
      // cerr << pC[i][j];
    }
    // cerr << endl;
  }
  // cerr << endl;
  
  // rebuild motif model from new counts
  // cerr << "Rebuild _localMotif from pC." << endl;
  if (_localMotif != NULL)
  {
    _localMotif->RebuildMatrix(_w, pC);
  }
  else
  {
    _localMotif = new PWM(_w, pC);
  }

  // cerr << "DEBUG: motif model (reduced) = " << endl;
  // _localMotif->StderrPrintMatrix();
  
  // clear local variable
  delete[]pC;
  
  // cerr << "DEBUG: MotifSamplerRun::BuildMotifFromReducedInstanceMap" << endl;
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
void
  PriorMotifSamplerRun::InitializationStep(int iterations)
{
  MapIterator mi = _pComputationMap->begin();
  SequenceComputation *pComp = NULL;

  cerr << "Initialize motif model " << endl;
  if ( _localMotif != NULL )
    delete _localMotif;
  _localMotif = new PWM(_w, _pPriorMotifModel);

  cerr << "Score all sequences with initial motif model." << endl;
  UpdateMotifScores();

  cerr << "Create initial instance map" << endl;
  SampleFixedSizeInstanceMap(1);
  
  // StderrPrintInstanceMap();
  
  cerr << "Starting initial 1 copy step" << endl;
  for (int i = 0; i < iterations; i++)
  {

    // loop over all sequences and update scores
    mi = _pComputationMap->begin();
    while (mi != _pComputationMap->end())
    {
      pComp = (*mi).second;

      // build motif model from instances in instance map excluding current sequence 
      BuildMotifFromReducedInstanceMap((*mi).first);

      // update motif scores
      pComp->UpdateInstanceMotifScore(_localMotif, _strand);
      pComp->UpdateInstanceExpWx(_strand);

      // next sequence
      mi++;
    }

    // sample exactly one instance in each sequence 
    SampleFixedSizeInstanceMap(1);

    // print information on screen
    cerr << i << "\t";
    StderrPrintMotifInfo();
  }

  // cerr << "DEBUG: MotifSamplerRun::InitializationStep" << endl;
  pComp = NULL;
  return;
}


void
  PriorMotifSamplerRun::CoreSamplingStep(int iterations)
{
  MapIterator mi = _pComputationMap->begin();
  SequenceComputation *pComp = NULL;

  for (int i = 0; i < iterations; i++)
  {

    // loop over all sequences and update scores
    mi = _pComputationMap->begin();
    while (mi != _pComputationMap->end())
    {
      pComp = (*mi).second;

      // build motif model from instances in instance map excluding current sequence 
      BuildMotifFromReducedInstanceMap((*mi).first);

      // update motif scores
      pComp->UpdateInstanceMotifScore(_localMotif, _strand);
      pComp->UpdateInstanceExpWx(_strand);

      // compute distribution to estimate number of instances
      pComp->UpdateCopyProbability(_prior, _strand);

      // next sequence
      mi++;
    }

    // sample a new instance map
    SampleInstanceMap();

    // check number of instances
    if (_nbrInstances < 2 || _nbrSequencesWithInstances < 2)
    {
      cerr << "--Warning-- not enough instances to proceed procedure." <<
        endl;
      pComp = NULL;
      break;
    }

    // print information on current motif
    cerr << i << "\t";
    StderrPrintMotifInfo();

    pComp = NULL;
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
void
  PriorMotifSamplerRun::ConvergenceStep(int iterations)
{
  MapIterator
    mi = _pComputationMap->begin();
  SequenceObject *
    pSeq = NULL;
  SequenceComputation *
    pComp = NULL;
  // cerr << "DEBUG: MotifSamplerRun::ConvergenceStep" << endl;
  for (int i = 0; i < iterations; i++)
  {
    // build motif model from all instances
    BuildMotifFromInstanceMap();

    // loop over all sequences
    mi = _pComputationMap->begin();
    while (mi != _pComputationMap->end())
    {
      pSeq = (*mi).first;
      pComp = (*mi).second;

      // update motif scores
      pComp->UpdateInstanceMotifScore(_localMotif, _strand);
      pComp->UpdateInstanceExpWx(_strand);

      // compute distribution to estimate number of instances
      pComp->UpdateCopyProbability(_prior, _strand);

      // next sequence
      mi++;
    }

    // create new instance map from curent scores
    SelectBestInstanceMap();

    // check number of instances
    if (_nbrInstances < 2 || _nbrSequencesWithInstances < 2)
    {
      cerr << "--Warning-- not enough instances to proceed procedure." <<
        endl;
      pSeq = NULL;
      pComp = NULL;
      break;
    }

    cerr << i << "\t";
    StderrPrintMotifInfo();

    pSeq = NULL;
    pComp = NULL;
  }

  // final motif model
  BuildMotifFromInstanceMap();
  // StderrPrintInstanceMap();
  
  // cerr << "DEBUG: MotifSamplerRun::ConvergenceStep" << endl;
  return;
}


void 
  PriorMotifSamplerRun::StderrPrintMixtureCoefficients()
{
  cerr << "Mixture Coefficients: |";
  for ( int i=0; i<_w; i++ ) 
    cerr << (*_pMixtureCoefficients)[i] << "|";
  cerr << endl;

  return;
}
