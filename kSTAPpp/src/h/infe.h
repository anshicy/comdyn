#pragma once
#pragma once

#include "Element.h"

using namespace std;

//! shell element class
class CInf : public CElement
{
public:

	//!	Constructor
	CInf();

	//!	Desconstructor
	~CInf();

	//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

	//!	Write element data to stream
	virtual void Write(COutputter& output);

	//!	Calculate element stiffness matrix
	virtual void ElementStiffness(double* StfMatrix);

	//!	Calculate element stress
	virtual void ElementStress(double* stress, double* Displacement);


};
