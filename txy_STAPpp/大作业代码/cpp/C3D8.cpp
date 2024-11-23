/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "C3D8.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

Vector2d s;

Matrix<double, 6, 6> C;

//	Constructor
C3D8::C3D8()
{
	NEN_ = 8;	// Each element has 8 nodes
	nodes_ = new CNode*[NEN_];
    
    ND_ = 24; 
    LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;

	ElementNDF = 3;
}

//	Desconstructor
C3D8::~C3D8()
{
}

//	Read element data from stream Input
bool C3D8::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet;	// Material property set number
	unsigned int N1, N2, N3, N4, N5, N6, N7, N8;	// The four node numbers

	Input >> N1 >> N2 >> N3 >> N4 >> N5 >> N6 >> N7 >> N8 >> MSet;
    ElementMaterial_ = dynamic_cast<C3D8Material*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];
	nodes_[3] = &NodeList[N4 - 1];
	nodes_[4] = &NodeList[N5 - 1];
	nodes_[5] = &NodeList[N6 - 1];
	nodes_[6] = &NodeList[N7 - 1];
	nodes_[7] = &NodeList[N8 - 1];

	return true;
}

//	Write element data to stream
void C3D8::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber 
		   << setw(9) << nodes_[2]->NodeNumber 
		   << setw(9) << nodes_[3]->NodeNumber
		   << setw(9) << nodes_[4]->NodeNumber
		   << setw(9) << nodes_[5]->NodeNumber
		   << setw(9) << nodes_[6]->NodeNumber
		   << setw(9) << nodes_[7]->NodeNumber     
		   << setw(12) << ElementMaterial_->nset << endl;
}

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void C3D8::ElementStiffness(double* stfMatrix)
{
	clear(stfMatrix, SizeOfStiffnessMatrix());

//	Calculate bar length
	double X[8];
	double Y[8];
	double Z[8];		//	dx = X[1]-X[0], dy = Y[1]-Y[0], dz = z2-z1
	for (unsigned int i = 0; i <= 7; i++)
	  {
		X[i] = nodes_[i]->XYZ[0];
		Y[i] = nodes_[i]->XYZ[1];
		Z[i] = nodes_[i]->XYZ[2];
	  }

//	Calculate element stiffness matrix

	C3D8Material* material_ = dynamic_cast<C3D8Material*>(ElementMaterial_);	// Pointer to material of the element

	double E = material_->E; 
	double niu = material_->niu;


	//积分点
	Vector2d s;
	s << sqrt(3.)/3., -sqrt(3.)/3.;

	Vector2d  r;
	r << sqrt(3.)/3., -sqrt(3.)/3.;

	Vector2d  t;
	t << sqrt(3.)/3., -sqrt(3.)/3.;

	//单元刚度矩阵
	MatrixXd Ke;
	Ke = MatrixXd::Zero(24, 24); //这里没问题

	for (int ids = 0; ids < 2; ids++)
	{
		for (int idr = 0; idr < 2; idr++)
		{
			for (int idt = 0; idt < 2; idt++){

				Matrix<double, 3, 8> G_H8;
				G_H8 << -(1-r(idr))*(1-t(idt)),  (1-r(idr))*(1-t(idt)),  (1+r(idr))*(1-t(idt)), -(1+r(idr))*(1-t(idt)), -(1-r(idr))*(1+t(idt)),  (1-r(idr))*(1+t(idt)), (1+r(idr))*(1+t(idt)), -(1+r(idr))*(1+t(idt)), 
					    -(1-s(ids))*(1-t(idt)), -(1+s(ids))*(1-t(idt)),  (1+s(ids))*(1-t(idt)),  (1-s(ids))*(1-t(idt)), -(1-s(ids))*(1+t(idt)), -(1+s(ids))*(1+t(idt)), (1+s(ids))*(1+t(idt)),  (1-s(ids))*(1+t(idt)),
					    -(1-s(ids))*(1-r(idr)), -(1+s(ids))*(1-r(idr)), -(1+s(ids))*(1+r(idr)), -(1-s(ids))*(1+r(idr)),  (1-s(ids))*(1-r(idr)),  (1+s(ids))*(1-r(idr)), (1+s(ids))*(1+r(idr)),  (1-s(ids))*(1+r(idr)); 

				G_H8 = 1./ 8. * G_H8;

				Matrix<double, 8, 3> XYZe;
				XYZe << X[0], Y[0], Z[0],
					X[1], Y[1], Z[1],
					X[2], Y[2], Z[2],
					X[3], Y[3], Z[3],
					X[4], Y[4], Z[4],
					X[5], Y[5], Z[5],
					X[6], Y[6], Z[6],
					X[7], Y[7], Z[7];

				Matrix<double, 3, 3> Jacobi;
				Jacobi = G_H8 * XYZe;

				Matrix<double, 3, 3> Jacobi_inv;
				Jacobi_inv = Jacobi.inverse();

				Matrix<double, 3, 8> N_H8;
				N_H8 = Jacobi_inv * G_H8;

				Matrix<double, 6, 24> Bij;
				for (int i = 0; i < 8; i++)
					{
						Bij(0, i * 3)     = N_H8(0, i);
						Bij(0, i * 3 + 1) = 0;
						Bij(0, i * 3 + 2) = 0;

						Bij(1, i * 3)     = 0;
						Bij(1, i * 3 + 1) = N_H8(1, i);
						Bij(1, i * 3 + 2) = 0;

						Bij(2, i * 3)     = 0;
						Bij(2, i * 3 + 1) = 0;
						Bij(2, i * 3 + 2) = N_H8(2, i);

						Bij(3, i * 3)     = N_H8(1, i);
						Bij(3, i * 3 + 1) = N_H8(0, i);
						Bij(3, i * 3 + 2) = 0;

						Bij(4, i * 3)     = N_H8(2, i);
						Bij(4, i * 3 + 1) = 0;
						Bij(4, i * 3 + 2) = N_H8(0, i);

						Bij(5, i * 3)     = 0;
						Bij(5, i * 3 + 1) = N_H8(2, i);
						Bij(5, i * 3 + 2) = N_H8(1, i);
					}

				//材料矩阵 弹性矩阵
				Matrix<double, 6, 6>C;
				C << 1 - niu, niu, niu, 0, 0, 0,
					niu, 1 - niu, niu, 0, 0, 0,
					niu, niu, 1 - niu, 0, 0, 0,
					0, 0, 0, (1 - 2 * niu) / 2, 0, 0,
					0, 0, 0, 0, (1 - 2 * niu) / 2, 0,
					0, 0, 0, 0, 0, (1 - 2 * niu) / 2;

				C = C * E / (1 + niu) / (1 - 2 * niu);

				Ke = Ke + Bij.transpose() * C * Bij * Jacobi.determinant();
			}
		}
	}
	int p = 0;
	for (int i = 0; i < 24; i++){
		for (int j = i ; j >= 0; j--){
			stfMatrix[p] = Ke(i, j);
			p++;
		}
	}

}

//	Calculate element stress 
void C3D8::ElementStress(double* stress, double* Displacement)
{
	C3D8Material* material_ = dynamic_cast<C3D8Material*>(ElementMaterial_);	// Pointer to material of the element
	
	double DX[3];	//	dx = X[1]-X[0], dy = Y[1]-Y[0], dz = z2-z1
	double L2 = 0;	//	Square of bar length (L^2)

	for (unsigned int i = 0; i < 3; i++)
	{
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];
		L2 = L2 + DX[i]*DX[i];
	}

	double S[6];
	for (unsigned int i = 0; i < 3; i++)
	{
		S[i] = -DX[i] * material_->E / L2;
		S[i+3] = -S[i];
	}
	
	*stress = 0.0;
	for (unsigned int i = 0; i < 6; i++)
	{
		if (LocationMatrix_[i])
			*stress += S[i] * Displacement[LocationMatrix_[i]-1];
	}
}
