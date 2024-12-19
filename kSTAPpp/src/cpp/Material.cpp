/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Material.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

//	Read material data from stream Input
bool CBarMaterial::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E >> Area;	// Young's modulus and section area

	return true;
}

//	Write material data to Stream
void CBarMaterial::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << Area << endl;
}




//	Read material data from stream Input
bool CShellMaterial::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E  >> h >> niu >> rho /*>> k >> p >> a >> b*/;	// Young's modulus and length and width

	return true;
}

//	Write material data to Stream
void CShellMaterial::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << h << setw(16) << niu << setw(16) << rho << setw(16) /*<< k << setw(16) << p << setw(16) << a << setw(16) << b */<< endl;
}



bool CBeamMaterial::Read(ifstream& Input)
{
	Input >> nset; // Number of property set

	Input >> E >> Area >> Iz >> Iy >> v >> rho; // Young's modulus, element length, inertial moment, Possion's ratio and density.

	return true;
}

// Write material data to Stream
void CBeamMaterial::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << Area << setw(16) << Iz << setw(16) << Iy << setw(16) << v << setw(16) << rho << endl;
}
// Read material data from stream Input

bool C3D8Material::Read(ifstream& Input)
{
	Input >> nset; // Number of property set

	Input >> E >> niu; // Young's modulus and Poisson ratio

	return true;
}

// Write material data to Stream
void C3D8Material::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << niu << endl;
}
bool C3D20Material::Read(ifstream& Input)
{
	Input >> nset; // Number of property set

	Input >> E >> niu >> rho; // Young's modulus and Poisson ratio

	return true;
}

// Write material data to Stream
void C3D20Material::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << niu << endl;
}
bool CPS4Material::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E >> niu;	// Young's modulus and Poisson ratio

	return true;
}

//	Write material data to Stream
void CPS4Material::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << niu << endl;
}
bool CInfMaterial::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E >> niu;	// Young's modulus and Poisson ratio

	return true;
}

//	Write material data to Stream
void CInfMaterial::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << niu << endl;
}

bool CPS3Material::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E >> niu;	// Young's modulus and Poisson ratio

	return true;
}

void CPS3Material::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << niu << endl;
}