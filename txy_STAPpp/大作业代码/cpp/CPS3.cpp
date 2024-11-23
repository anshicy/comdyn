/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "CPS3.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Core>

using namespace Eigen;
using namespace std;

Matrix<double, 3, 3> C2;

//	Constructor
CPS3::CPS3()
{
	NEN_ = 3;	// Each element has 3 nodes
	nodes_ = new CNode*[NEN_];

	ND_ = 9;
	LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;

	ElementNDF = 3;
}

//	Desconstructor
CPS3::~CPS3()
{
}

//	Read element data from stream Input
bool CPS3::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet;	// Material property set niumber
	unsigned int N1, N2, N3;	// Left node number and right node number

	Input >> N1 >> N2 >> N3 >> MSet;
	ElementMaterial_ = dynamic_cast<CPS3Material*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];

	return true;
}

//	Write element data to stream
void CPS3::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber
		   << setw(9) << nodes_[2]->NodeNumber
		   << setw(12) << ElementMaterial_->nset << endl;
}//setw:空格

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CPS3::ElementStiffness(double* stfMatrix)
{
	clear(stfMatrix, SizeOfStiffnessMatrix());

	//	Calculate bar length
	double X[3];
	double Y[3];	
	double Z[3];//	dx = X[1]-X[0], dy = Y[1]-Y[0], dz = z2-z1
	for (unsigned int i = 0; i <= 2; i++)
	{
		X[i] = nodes_[i]->XYZ[0];
		Y[i] = nodes_[i]->XYZ[1];
		Z[i] = nodes_[i]->XYZ[2];
	}

	//	Calculate element stiffness matrix

	CPS3Material* material_ = dynamic_cast<CPS3Material*>(ElementMaterial_);	// Pointer to material of the element

	double E = material_->E;
	double niu = material_->niu;


	//积分点有3个，对应ksi123的轮换，值为2/3，1/6，1/6 变换后空间：点1（点0）（1，0），点2（点1）（0，1），点3（点2）（0，0）
	//其实吧，CPS3的B是常数，K无需求和，直接乘就好，以下只是为了形式上的一致
	Vector3d ksi1;
	ksi1<< 2.0/3.0, 1.0/6.0, 1.0/6.0;

	Vector3d ksi2;
	ksi2 << 1.0/6.0, 2.0/3.0, 1.0/6.0;  

	Vector3d ksi3;
	ksi3 << 1.0/6.0, 1.0/6.0, 2.0/3.0;

	Vector3d w;//权重
	w << 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0;

	//单元刚度矩阵
	MatrixXd Ke;
	Ke = MatrixXd::Zero(6, 6);

			Matrix<double, 2, 2> Jacobi;

			Jacobi(0, 0) = X[0] - X[2];

			Jacobi(1, 0) = X[1] - X[2];

			Jacobi(0, 1) = Y[0] - Y[2];

			Jacobi(1, 1) = Y[1] - Y[2];

			Matrix<double, 2, 2> Jacobi_inv = Jacobi.inverse();

			Matrix<double, 3, 6>Bij;
			Bij << Y[1] - Y[2], 0, Y[2] - Y[0], 0, Y[0] - Y[1], 0,
				0, X[2] - X[1], 0, X[0] - X[2], 0, X[1] - X[0],
				X[2] - X[1], Y[1] - Y[2], X[0] - X[2], Y[2] - Y[0], X[1] - X[0], Y[0] - Y[1];

			Bij = Bij / Jacobi.determinant();

			//材料矩阵 弹性矩阵
			Matrix<double, 3, 3> C2;
			C2 << 1.0, niu, 0.0,
				niu, 1.0, 0.0,
				0.0, 0.0, (1.0 - niu) / 2.0;

			C2 = C2 * E / (1.0 - niu * niu);

			Ke = Bij.transpose() * C2 * Bij * Jacobi.determinant() * 0.5;


	stfMatrix[0] = Ke(0, 0);

	stfMatrix[1] = Ke(1, 1);
	stfMatrix[2] = Ke(0, 1);

	stfMatrix[3] = 0.0001;
	stfMatrix[4] = 0.0001;
	stfMatrix[5] = 0.0001;
	
	stfMatrix[6] = Ke(2, 2);
	stfMatrix[7] = 0.0001;
	stfMatrix[8] = Ke(1, 2);
	stfMatrix[9] = Ke(0, 2);

	stfMatrix[10] = Ke(3, 3);
	stfMatrix[11] = Ke(2, 3);
	stfMatrix[12] = 0.0001;
	stfMatrix[13] = Ke(1, 3);
	stfMatrix[14] = Ke(0, 3);

	stfMatrix[15] = 0.0001;
	stfMatrix[16] = 0.0001;
	stfMatrix[17] = 0.0001;
	stfMatrix[18] = 0.0001;
	stfMatrix[19] = 0.0001;
	stfMatrix[20] = 0.0001;

	stfMatrix[21] = Ke(4, 4);
	stfMatrix[22] = 0.0001;
	stfMatrix[23] = Ke(3, 4);
	stfMatrix[24] = Ke(2, 4);
	stfMatrix[25] = 0.0001;
	stfMatrix[26] = Ke(1, 4);
	stfMatrix[27] = Ke(0, 4);

	stfMatrix[28] = Ke(5, 5);
	stfMatrix[29] = Ke(4, 5);
	stfMatrix[30] = 0.0001;
	stfMatrix[31] = Ke(3, 5);
	stfMatrix[32] = Ke(2, 5);
	stfMatrix[33] = 0.0001;
	stfMatrix[34] = Ke(1, 5);
	stfMatrix[35] = Ke(0, 5);

	stfMatrix[36] = 0.0001;
	stfMatrix[37] = 0.0001;
	stfMatrix[38] = 0.0001;
	stfMatrix[39] = 0.0001;
	stfMatrix[40] = 0.0001;
	stfMatrix[41] = 0.0001;
	stfMatrix[42] = 0.0001;
	stfMatrix[43] = 0.0001;
	stfMatrix[44] = 0.0001;


}

//	Calculate element stress 
void CPS3::ElementStress(double* stress, double* Displacement)
{
	CPS3Material* material_ = dynamic_cast<CPS3Material*>(ElementMaterial_);	// Pointer to material of the element

	double DX[3];	//	dx = X[1]-X[0], dy = Y[1]-Y[0], dz = z2-z1
	double L2 = 0;	//	Square of bar length (L^2)
	double X[4];
	double Y[4];
	double Z[4];//	dx = X[1]-X[0], dy = Y[1]-Y[0], dz = z2-z1
	for (unsigned int i = 0; i < 3; i++)
	{
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];
		L2 = L2 + DX[i] * DX[i];
	}

	double S[6];
	for (unsigned int i = 0; i < 3; i++)
	{
		S[i] = -DX[i] * material_->E / L2;
		S[i + 3] = -S[i];
	}

	*stress = 0.0;
	for (unsigned int i = 0; i < 6; i++)
	{
		if (LocationMatrix_[i])
			*stress += S[i] * Displacement[LocationMatrix_[i] - 1];
	}
}
