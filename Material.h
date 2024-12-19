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

#include "Outputter.h"

using namespace std;

//!	Material base class which only define one data member
/*!	All type of material classes should be derived from this base class */
class CMaterial
{
public:

	unsigned int nset;	//!< Number of set
	
	double E;  //!< Young's modulus

public:

//! Virtual deconstructor
    virtual ~CMaterial() {};

//!	Read material data from stream Input
	virtual bool Read(ifstream& Input) = 0;

//!	Write material data to Stream
    virtual void Write(COutputter& output) = 0;

};

//!	Material class for bar element
class CBarMaterial : public CMaterial
{
public:

	double Area;	//!< Sectional area of a bar element

public:
	
//!	Read material data from stream Input
	virtual bool Read(ifstream& Input);

//!	Write material data to Stream
	virtual void Write(COutputter& output);
};
//!	Material class for Plate element
class CShellMaterial : public CMaterial
{
public:

	double h;	//!the width
	double niu;	//!the possion
	double rho;	//!the rho
	//double k;	//!the fix of tau
	//double p;	//!the Ewithin
	//double a;	//!the length
	//double b;	//!the kuan

public:

	//!	Read material data from stream Input
	virtual bool Read(ifstream& Input);

	//!	Write material data to Stream
	virtual void Write(COutputter& output);
};
//! Material class for beam element
class CBeamMaterial : public CMaterial
{
public:

	double Area; //!< Sectional area of a beam element
	double l; //!< Sectional length of a beam element
	double Iz; //!< Sectional inertia moment of a beam element (x dir)
	double Iy; //!< Sectional inertia moment of a beam element (y dir)
	double v;   //!< Sectional possion ratio of a beam element
	double rho;   //!< Sectional possion ratio of a beam element

public:

	//! Read material data from stream Input
	virtual bool Read(ifstream& Input);

	//! Write material data to Stream
	virtual void Write(COutputter& output);
};
class C3D8Material : public CMaterial
{
public:

	double niu; //!< Poisson Ratio of a bar element

public:

	//! Read material data from stream Input
	virtual bool Read(ifstream& Input);

	//! Write material data to Stream
	virtual void Write(COutputter& output);
};
class CPS4Material : public CMaterial
{
public:

	double niu;	//!< Poisson Ratio of a Q4 element

public:

	//!	Read material data from stream Input
	virtual bool Read(ifstream& Input);

	//!	Write material data to Stream
	virtual void Write(COutputter& output);
};
class C3D20Material : public CMaterial
{
public:

	double niu;	//!< Poisson Ratio of a C3D20 element
	double rho;
public:

	//!	Read material data from stream Input
	virtual bool Read(ifstream& Input);

	//!	Write material data to Stream
	virtual void Write(COutputter& output);
};
class CPS3Material : public CMaterial
{
public:

	double niu;	//!< Poisson Ratio of a T3 element

public:

	//!	Read material data from stream Input
	virtual bool Read(ifstream& Input);

	//!	Write material data to Stream
	virtual void Write(COutputter& output);
};

class CInfMaterial : public CMaterial
{
public:

	double niu;	//!< Poisson Ratio of a Q4 element

public:

	//!	Read material data from stream Input
	virtual bool Read(ifstream& Input);

	//!	Write material data to Stream
	virtual void Write(COutputter& output);
};
