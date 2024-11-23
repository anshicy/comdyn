
/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     plate by Xiny Tang                                     */
/*                                                                           */
/*                                                    */
/*****************************************************************************/

#pragma once

#include "Element.h"

using namespace std;

//! shell element class
class CShell : public CElement
{
public:

	//!	Constructor
	CShell();

	//!	Desconstructor
	~CShell();

	//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

	//!	Write element data to stream
	virtual void Write(COutputter& output);

	//!	Calculate element stiffness matrix
	virtual void ElementStiffness(double* StfMatrix);

	//!	Calculate element stress
	virtual void ElementStress(double* stress, double* Displacement);


};