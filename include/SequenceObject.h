#ifndef sequence_object_declared
#define sequence_object_declared

#include "inclusive.h"
vector<int> *_TranslateSequenceString(const string *pSeqStr);

class SequenceObject{
 private: 
  // internal variables
  int _length;
  string _id;
  string _description;
  Sequence *_sequenceVector;
  Sequence *_revSequenceVector;

  // private helper functions
  void _CreateSequenceVector(const string *pSeqStr);
  void _SetReverseComplement();
  void _ParseDescriptionHeader(string *pDescStr);
  string *_TranslateSequenceVector(strand_modes s, int start, int length);

 public:
  // constructors
  SequenceObject(const string * pSeqStr);
  SequenceObject(Sequence * pSeqVector);

  // destructor
  ~SequenceObject();

  // adaptors
  SequenceObject * SetID(string *idname);
  SequenceObject * SetDescription(string *desc);

  // inspectors
  string * GetID(){return &_id;};
  string * GetDescription(){return &_description;};
  int Length(){return _length;};
  bool CheckSequence();
  string * GetSequenceString(strand_modes s);
  string * GetSubSequenceString(strand_modes s, int start, int length);
  int GetNucleotideAt(strand_modes s, int pos);
  
  Sequence::iterator begin(strand_modes s);
  Sequence::iterator end(strand_modes s);
 
  Sequence * GetSequenceVector(){ return _sequenceVector;};
  Sequence * GetReverseSequenceVector(){ return _revSequenceVector;};
  Sequence * GetSubSequenceVector(strand_modes s, int start, int length);

};

#endif
