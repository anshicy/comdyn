/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "CPS4.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Core>

using namespace Eigen;
using namespace std;



//	Constructor
CPS4::CPS4()
{
	NEN_ = 4;	// Each element has 4 nodes
	nodes_ = new CNode*[NEN_];
    
    ND_ = 12;
    LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;

	ElementNDF = 3;
}

//	Desconstructor
CPS4::~CPS4()
{
}

//	Read element data from stream Input.Starting from up-left
bool CPS4::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet;	// Material property set niumber
	unsigned int N1, N2, N3, N4;	// Left node niumber and right node niumber

	Input >> N1 >> N2 >> N3 >> N4 >> MSet;
    ElementMaterial_ = dynamic_cast<CPS4Material*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];
	nodes_[3] = &NodeList[N4 - 1];

	return true;
}

//	Write element data to stream
void CPS4::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber 
		   << setw(9) << nodes_[2]->NodeNumber 
		   << setw(9) << nodes_[3]->NodeNumber 
		   << setw(12) << ElementMaterial_->nset << endl;
}

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CPS4::ElementStiffness(double* stfMatrix)
{
	clear(stfMatrix, SizeOfStiffnessMatrix());

//	Calculate bar length
	double X[4];
	double Y[4];
	double Z[4];//	dx = X[1]-X[0], dy = Y[1]-Y[0], dz = z2-z1
	for (unsigned int i = 0; i <= 3; i++)
	  {
		X[i] = nodes_[i]->XYZ[0];
		Y[i] = nodes_[i]->XYZ[1];
		Z[i] = nodes_[i]->XYZ[2];
	  }

//	Calculate element stiffness matrix

	CPS4Material* material_ = dynamic_cast<CPS4Material*>(ElementMaterial_);	// Pointer to material of the element

	double E = material_->E; 
	double niu = material_->niu;


	//积分点
	
			Vector2d s;
			s << sqrt(1.0/3.0), -sqrt(1.0 / 3.0); // eta gaussian

			Vector2d r;
			r << sqrt(1.0 / 3.0), -sqrt(1.0 / 3.0); // ksi gaussian

			Matrix2d af;
			af << 1., 1.,
				1., 1.;

			//单元刚度矩阵
			MatrixXd Ke;
			Ke = MatrixXd::Zero(8, 8);

			for (int ids = 0; ids < 2; ids++)
			{
				for (int idr = 0; idr < 2; idr++)
				{
					Matrix<double, 2, 2> Jacobi;


					/*Jacobi(0, 0) = -(1. / 4. * (1 - s(ids)) * X[0] +
						1. / 4. * (1 - s(ids)) * X[1] +
						1. / 4. * (1 + s(ids)) * X[2] -
						1. / 4. * (1 + s(ids)) * X[3]);

					Jacobi(1, 0) = -(1. / 4. * (1 - r(idr)) * X[0] -
						1. / 4. * (1 + r(idr)) * X[1] +
						1. / 4. * (1 + r(idr)) * X[2] +
						1. / 4. * (1 - r(idr)) * X[3]);

					Jacobi(0, 1) = -(1. / 4. * (1 - s(ids)) * Y[0] +
						1. / 4. * (1 - s(ids)) * Y[1] +
						1. / 4. * (1 + s(ids)) * Y[2] -
						1. / 4. * (1 + s(ids)) * Y[3]);

					Jacobi(1, 1) = -(1. / 4. * (1 - r(idr)) * Y[0] -
						1. / 4. * (1 + r(idr)) * Y[1] +
						1. / 4. * (1 + r(idr)) * Y[2] +
						1. / 4. * (1 - r(idr)) * Y[3]);*/

					Jacobi(0, 0) = 0.25 * ((X[1] - X[0]) * (1 - s(ids)) - (X[3] - X[2]) * (1 + s(ids)));
					Jacobi(1, 0) = 0.25 * ((X[3] - X[0]) * (1 - r(idr)) + (X[2] - X[1]) * (1 + r(idr)));
					Jacobi(0, 1) = 0.25 * ((Y[1] - Y[0]) * (1 - s(ids)) - (Y[3] - Y[2]) * (1 + s(ids)));
					Jacobi(1, 1) = 0.25 * ((Y[3] - Y[0]) * (1 - r(idr)) + (Y[2] - Y[1]) * (1 + r(idr)));

					Matrix<double, 2, 2> Jacobi_inv = Jacobi.inverse();

					/*Matrix<double, 2, 8> B_u, B_v;
					B_u << 1 + s(ids), 0, -(1 + s(ids)), 0, -(1 - s(ids)), 0, 1 - s(ids), 0,
						1 + r(idr), 0, 1 - r(idr), 0, -(1 - r(idr)), 0, -(1 + r(idr)), 0;

					B_v << 0, 1 + s(ids), 0, -(1 + s(ids)), 0, -(1 - s(ids)), 0, 1 - s(ids),
						0, 1 + r(idr), 0, 1 - r(idr), 0, -(1 - r(idr)), 0, -(1 + r(idr));

				
					Matrix<double, 2, 8> u_x_y, v_x_y;
					u_x_y = 1. / 4. * Jacobi_inv * B_u;

					v_x_y = 1. / 4. * Jacobi_inv * B_v;*/

					Matrix<double, 2, 4> GQ; //广义坐标梯度
					GQ << s(ids) - 1.0, 1.0 - s(ids), 1.0 + s(ids), -s(ids) - 1.0,
						r(idr) - 1.0, -r(idr) - 1.0, 1 + r(idr), 1.0 - r(idr);
					GQ = GQ * 0.25;

					Matrix<double, 2, 4> DN; //物理坐标梯度
					DN = Jacobi_inv * GQ;

					Matrix<double, 2, 8> u_x_y, v_x_y;
					u_x_y << DN(0, 0), 0, DN(0, 1), 0, DN(0, 2), 0, DN(0, 3), 0,
						DN(1, 0), 0, DN(1, 1), 0, DN(1, 2), 0, DN(1, 3), 0;

					v_x_y << 0, DN(0, 0), 0, DN(0, 1), 0, DN(0, 2), 0, DN(0, 3),
						0, DN(1, 0), 0, DN(1, 1), 0, DN(1, 2), 0, DN(1, 3);

					Matrix<double, 3, 8>Bij;
					for (int i = 0; i < 8; i++)
					{
						Bij(0, i) = u_x_y(0, i);
						Bij(1, i) = v_x_y(1, i);
						Bij(2, i) = u_x_y(1, i) + v_x_y(0, i);
					}

					//材料矩阵 弹性矩阵

					Matrix<double, 3, 3>C1;
					C1<< 1.0, niu, 0.0,
						niu, 1.0, 0.0,
						0.0, 0.0, (1.0 - niu) / 2.0;

					C1 = C1 * E / (1.0 - niu * niu);

					Ke = Ke + Bij.transpose() * C1 * Bij * Jacobi.determinant() * af(ids, idr);

				}
			}



			double det;
			det = Ke.determinant();
	stfMatrix[0] = Ke(0,0);

	stfMatrix[1] = Ke(1,1);
	stfMatrix[2] = Ke(0,1);

	stfMatrix[3] = 0.0001;
	stfMatrix[4] = 0.0001;
	stfMatrix[5] = 0.0001;


	stfMatrix[6] = Ke(2,2);
	stfMatrix[7] = 0.0001;
	stfMatrix[8] = Ke(1,2);
	stfMatrix[9] = Ke(0,2);

	stfMatrix[10] = Ke(3,3);
	stfMatrix[11] = Ke(2,3);
	stfMatrix[12] = 0.0001;
	stfMatrix[13] = Ke(1,3);
	stfMatrix[14] = Ke(0,3);

	stfMatrix[15] = 0.0001;
	stfMatrix[16] = 0.0001;
	stfMatrix[17] = 0.0001;
	stfMatrix[18] = 0.0001;
	stfMatrix[19] = 0.0001;
	stfMatrix[20] = 0.0001;

	stfMatrix[21] = Ke(4,4);
	stfMatrix[22] = 0.0001;
	stfMatrix[23] = Ke(3,4);
	stfMatrix[24] = Ke(2,4);
	stfMatrix[25] = 0.0001;
	stfMatrix[26] = Ke(1,4);
	stfMatrix[27] = Ke(0,4);

	stfMatrix[28] = Ke(5,5);
	stfMatrix[29] = Ke(4,5);
	stfMatrix[30] = 0.0001;
	stfMatrix[31] = Ke(3,5);
	stfMatrix[32] = Ke(2,5);
	stfMatrix[33] = 0.0001;
	stfMatrix[34] = Ke(1,5);
	stfMatrix[35] = Ke(0,5);

	stfMatrix[36] = 0.0001;
	stfMatrix[37] = 0.0001;
	stfMatrix[38] = 0.0001;
	stfMatrix[39] = 0.0001;
	stfMatrix[40] = 0.0001;
	stfMatrix[41] = 0.0001;
	stfMatrix[42] = 0.0001;
	stfMatrix[43] = 0.0001;
	stfMatrix[44] = 0.0001;

	stfMatrix[45] = Ke(6,6);
	stfMatrix[46] = 0.0001;
	stfMatrix[47] = Ke(5,6);
	stfMatrix[48] = Ke(4,6);
	stfMatrix[49] = 0.0001;
	stfMatrix[50] = Ke(3,6);
	stfMatrix[51] = Ke(2,6);
	stfMatrix[52] = 0.0001;
	stfMatrix[53] = Ke(1,6);
	stfMatrix[54] = Ke(0,6);

	stfMatrix[55] = Ke(7,7);
	stfMatrix[56] = Ke(6,7);
	stfMatrix[57] = 0.0001;
	stfMatrix[58] = Ke(5,7);
	stfMatrix[59] = Ke(4,7);
	stfMatrix[60] = 0.0001;
	stfMatrix[61] = Ke(3,7);
	stfMatrix[62] = Ke(2,7);
	stfMatrix[63] = 0.0001;
	stfMatrix[64] = Ke(1,7);
	stfMatrix[65] = Ke(0,7);

	stfMatrix[66] = 0.0001;
	stfMatrix[67] = 0.0001;
	stfMatrix[68] = 0.0001;
	stfMatrix[69] = 0.0001;
	stfMatrix[70] = 0.0001;
	stfMatrix[71] = 0.0001;
	stfMatrix[72] = 0.0001;
	stfMatrix[73] = 0.0001;
	stfMatrix[74] = 0.0001;
	stfMatrix[75] = 0.0001;
	stfMatrix[76] = 0.0001;
	stfMatrix[77] = 0.0001;

	

}

//	Calculate element stress 
void CPS4::ElementStress(double* stress, double* Displacement)
{
	CPS4Material* material_ = dynamic_cast<CPS4Material*>(ElementMaterial_);	// Pointer to material of the element

	double DX[3];	//	dx = X[1]-X[0], dy = Y[1]-Y[0], dz = z2-z1
	double L2 = 0;	//	Square of bar length (L^2)
	double X[4];
	double Y[4];
	double Z[4];//	dx = X[1]-X[0], dy = Y[1]-Y[0], dz = z2-z1
	for (unsigned int i = 0; i <= 3; i++)
	{
		X[i] = nodes_[i]->XYZ[0];
		Y[i] = nodes_[i]->XYZ[1];
		Z[i] = nodes_[i]->XYZ[2];
	}

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

	double E = material_->E;
	double niu = material_->niu;
	//材料矩阵 弹性矩阵

	Matrix<double, 3, 3>C1;
	C1 << 1, niu, 0,
		niu, 1, 0,
		0, 0, (1 - niu) / 2;

	C1 = C1 * E / (1 - niu * niu);

	Vector2d s;
	s << -0.577, 0.577; // eta gaussian

	Vector2d r;
	r << -0.577, 0.577; // ksi gaussian

	Matrix2d af;
	af << 1., 1.,
		1., 1.;

	for (int ids = 0; ids < 2; ids++)
	{
		for (int idr = 0; idr < 2; idr++)
		{
			Matrix<double, 2, 2> Jacobi;


			Jacobi(0, 0) = -(1. / 4. * (1 - s(ids)) * X[0] -
				1. / 4. * (1 - s(ids)) * X[1] -
				1. / 4. * (1 + s(ids)) * X[2] +
				1. / 4. * (1 + s(ids)) * X[3]);

			Jacobi(1, 0) = -(1. / 4. * (1 - r(idr)) * X[0] +
				1. / 4. * (1 + r(idr)) * X[1] -
				1. / 4. * (1 + r(idr)) * X[2] -
				1. / 4. * (1 - r(idr)) * X[3]);

			Jacobi(0, 1) = -(1. / 4. * (1 - s(ids)) * Y[0] -
				1. / 4. * (1 - s(ids)) * Y[1] -
				1. / 4. * (1 + s(ids)) * Y[2] +
				1. / 4. * (1 + s(ids)) * Y[3]);

			Jacobi(1, 1) = -(1. / 4. * (1 - r(idr)) * Y[0] +
				1. / 4. * (1 + r(idr)) * Y[1] -
				1. / 4. * (1 + r(idr)) * Y[2] -
				1. / 4. * (1 - r(idr)) * Y[3]);

			Matrix<double, 2, 2> Jacobi_inv = Jacobi.inverse();

			/*Matrix<double, 2, 8> B_u, B_v;
			B_u << 1 + s(ids), 0, -(1 + s(ids)), 0, -(1 - s(ids)), 0, 1 - s(ids), 0,
				1 + r(idr), 0, 1 - r(idr), 0, -(1 - r(idr)), 0, -(1 + r(idr)), 0;

			B_v << 0, 1 + s(ids), 0, -(1 + s(ids)), 0, -(1 - s(ids)), 0, 1 - s(ids),
				0, 1 + r(idr), 0, 1 - r(idr), 0, -(1 - r(idr)), 0, -(1 + r(idr));


			Matrix<double, 2, 8> u_x_y, v_x_y;
			u_x_y = 1. / 4. * Jacobi_inv * B_u;

			v_x_y = 1. / 4. * Jacobi_inv * B_v;*/

			Matrix<double, 2, 4> GQ; //广义坐标梯度
			GQ << s(ids) - 1, 1 - s(ids), 1 + s(ids), -s(ids) - 1,
				r(idr) - 1, -r(idr) - 1, 1 + r(idr), 1 - r(idr);
			GQ = GQ * 0.25;

			Matrix<double, 2, 4> DN; //物理坐标梯度
			DN = Jacobi_inv * GQ;

			Matrix<double, 2, 8> u_x_y, v_x_y;
			u_x_y << DN(0, 0), 0, DN(0, 1), 0, DN(0, 2), 0, DN(0, 3), 0,
				DN(1, 0), 0, DN(1, 1), 0, DN(1, 2), 0, DN(1, 3), 0;

			v_x_y << 0, DN(0, 0), 0, DN(0, 1), 0, DN(0, 2), 0, DN(0, 3),
				0, DN(1, 0), 0, DN(1, 1), 0, DN(1, 2), 0, DN(1, 3);

			Matrix<double, 3, 8>Bij;
			for (int i = 0; i < 8; i++)
			{
				Bij(0, i) = u_x_y(0, i);
				Bij(1, i) = v_x_y(1, i);
				Bij(2, i) = u_x_y(1, i) + v_x_y(0, i);
			}

			


		}
	}
	
	*stress = 0.0;
	for (unsigned int i = 0; i < 6; i++)
	{
		if (LocationMatrix_[i])
			*stress += E * Displacement[LocationMatrix_[i]-1];
	}
}
