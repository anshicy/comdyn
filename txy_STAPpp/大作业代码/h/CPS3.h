#pragma once
/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include "Element.h"

using namespace std;

//! CPS3 element class
class CPS3 : public CElement
{
public:

	//!	Constructor
	CPS3();

	//!	Desconstructor
	~CPS3();

	//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

	//!	Write element data to stream
	virtual void Write(COutputter& output);

	//!	Calculate element stiffness matrix
	virtual void ElementStiffness(double* StfMatrix);

	//!	Calculate element stress
	virtual void ElementStress(double* stress, double* Displacement);
};
