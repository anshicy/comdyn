/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     plate by Xiny Tang                                      */
/*                                                                           */
/*                                                    */
/*****************************************************************************/

#include "Shell.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Core>

using namespace Eigen;
using namespace std;

//	Constructor
CShell::CShell()
{
	NEN_ = 4;	// Each element has 4 nodes
	nodes_ = new CNode * [NEN_];

	ND_ = 4 * 6; 
	LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;

	ElementNDF = 6;
}

//	Desconstructor
CShell::~CShell()
{
}
Matrix<double, 3, 12> GenerateBm(double a, double b, double h, double ksai, double ita)
{
	double ap = 1.0 / (4.0 * a * b);
	Matrix<double, 3, 12> Bm;
	Matrix<double, 3, 3> Dm;
	MatrixXd KeM;
	double x = ksai;
	double y = ita;
	KeM = MatrixXd::Zero(12, 12);

	Bm(0, 0) = -3 * b / a * (-1) * x * (1 + (-1) * y);
	Bm(0, 1) = 0;
	Bm(0, 2) = -1 * b * (-1) * (1 + 3 * (-1) * x) * (1 + (-1) * y);
	Bm(1, 0) = -3 * b / a * (-1) * y * (1 + (-1) * x);
	Bm(1, 1) = a * (-1) * (1 + 3 * (-1) * y) * (1 + (-1) * x);
	Bm(1, 2) = 0;
	Bm(2, 0) = (-1) * (-1) * (4 - 3 * x * x - 3 * y * y);
	Bm(2, 1) = b * (-1) * (3 * y * y + 2 * (-1) * y - 1);
	Bm(2, 2) = a * (-1) * (1 - 2 * (-1) * x - 3 * x * x);

	Bm(0, 3) = -3 * b / a * (1) * x * (1 + (-1) * y);
	Bm(0, 4) = 0;
	Bm(0, 5) = -1 * b * (1) * (1 + 3 * (1) * x) * (1 + (-1) * y);
	Bm(1, 3) = -3 * b / a * (-1) * y * (1 + (1) * x);
	Bm(1, 4) = a * (-1) * (1 + 3 * (-1) * y) * (1 + (1) * x);
	Bm(1, 5) = 0;
	Bm(2, 3) = (1) * (-1) * (4 - 3 * x * x - 3 * y * y);
	Bm(2, 4) = b * (1) * (3 * y * y + 2 * (-1) * y - 1);
	Bm(2, 5) = a * (-1) * (1 - 2 * (1) * x - 3 * x * x);


	Bm(0, 6) = -3 * b / a * (1) * x * (1 + (1) * y);
	Bm(0, 7) = 0;
	Bm(0, 8) = -1 * b * (1) * (1 + 3 * (1) * x) * (1 + (1) * y);
	Bm(1, 6) = -3 * b / a * (1) * y * (1 + (1) * x);
	Bm(1, 7) = a * (1) * (1 + 3 * (1) * y) * (1 + (1) * x);
	Bm(1, 8) = 0;
	Bm(2, 6) = (1) * (1) * (4 - 3 * x * x - 3 * y * y);
	Bm(2, 7) = b * (1) * (3 * y * y + 2 * (1) * y - 1);
	Bm(2, 8) = a * (1) * (1 - 2 * (1) * x - 3 * x * x);

	Bm(0, 9) = -3 * b / a * (-1) * x * (1 + (1) * y);
	Bm(0, 10) = 0;
	Bm(0, 11) = -1 * b * (-1) * (1 + 3 * (-1) * x) * (1 + (1) * y);
	Bm(1, 9) = -3 * b / a * (1) * y * (1 + (-1) * x);
	Bm(1, 10) = a * (1) * (1 + 3 * (1) * y) * (1 + (-1) * x);
	Bm(1, 11) = 0;
	Bm(2, 9) = (-1) * (1) * (4 - 3 * x * x - 3 * y * y);
	Bm(2, 10) = b * (-1) * (3 * y * y + 2 * (1) * y - 1);
	Bm(2, 11) = a * (1) * (1 - 2 * (-1) * x - 3 * x * x);



	Bm = ap * Bm;
	// !!!!!!!!!!!!!!!!!!!!!

	return Bm;
}
// 生成M部分的显示B（单元矩阵），用来高斯积分
Matrix<double, 12, 12> Generatessj(double a, double b ,double h, double ksai, double ita, Matrix3d D)
{
	double ap = 1.0 / (4.0 * a * b);
	Matrix<double, 3, 12> Bm;
	Matrix<double, 3, 3> Dm;
	MatrixXd KeM;
	double x = ksai;
	double y = ita;
	KeM = MatrixXd::Zero(12, 12);
	Dm = D;

	Bm(0, 0) = -3 * b / a * (-1) * x * (1 + (-1) * y);
	Bm(0, 1) = 0;
	Bm(0, 2) = -1 * b * (-1) * (1 + 3 * (-1) * x) * (1 + (-1) * y);
	Bm(1, 0) = -3 * b / a * (-1) * y * (1 + (-1) * x);
	Bm(1, 1) = a * (-1) * (1 + 3 * (-1) * y) * (1 + (-1) * x);
	Bm(1, 2) = 0;
	Bm(2, 0) = (-1) * (-1) * (4 - 3 * x * x - 3 * y * y);
	Bm(2, 1) = b * (-1) * (3 * y * y + 2 * (-1) * y - 1);
	Bm(2, 2) = a * (-1) * (1 - 2 * (-1) * x - 3 * x * x);

	Bm(0, 3) = -3 * b / a * (1) * x * (1 + (-1) * y);
	Bm(0, 4) = 0;
	Bm(0, 5) = -1 * b * (1) * (1 + 3 * (1) * x) * (1 + (-1) * y);
	Bm(1, 3) = -3 * b / a * (-1) * y * (1 + (1) * x);
	Bm(1, 4) = a * (-1) * (1 + 3 * (-1) * y) * (1 + (1) * x);
	Bm(1, 5) = 0;
	Bm(2, 3) = (1) * (-1) * (4 - 3 * x * x - 3 * y * y);
	Bm(2, 4) = b * (1) * (3 * y * y + 2 * (-1) * y - 1);
	Bm(2, 5) = a * (-1) * (1 - 2 * (1) * x - 3 * x * x);


	Bm(0, 6) = -3 * b / a * (1) * x * (1 + (1) * y);
	Bm(0, 7) = 0;
	Bm(0, 8) = -1 * b * (1) * (1 + 3 * (1) * x) * (1 + (1) * y);
	Bm(1, 6) = -3 * b / a * (1) * y * (1 + (1) * x);
	Bm(1, 7) = a * (1) * (1 + 3 * (1) * y) * (1 + (1) * x);
	Bm(1, 8) = 0;
	Bm(2, 6) = (1) * (1) * (4 - 3 * x * x - 3 * y * y);
	Bm(2, 7) = b * (1) * (3 * y * y + 2 * (1) * y - 1);
	Bm(2, 8) = a * (1) * (1 - 2 * (1) * x - 3 * x * x);

	Bm(0, 9) = -3 * b / a * (-1) * x * (1 + (1) * y);
	Bm(0, 10) = 0;
	Bm(0, 11) = -1 * b * (-1) * (1 + 3 * (-1) * x) * (1 + (1) * y);
	Bm(1, 9) = -3 * b / a * (1) * y * (1 + (-1) * x);
	Bm(1, 10) = a * (1) * (1 + 3 * (1) * y) * (1 + (-1) * x);
	Bm(1, 11) = 0;
	Bm(2, 9) = (-1) * (1) * (4 - 3 * x * x - 3 * y * y);
	Bm(2, 10) = b * (-1) * (3 * y * y + 2 * (1) * y - 1);
	Bm(2, 11) = a * (1) * (1 - 2 * (-1) * x - 3 * x * x);



	Bm = ap * Bm;
	// !!!!!!!!!!!!!!!!!!!!!
	KeM = Bm.transpose() * Dm * Bm * a * b ;
	
	return KeM;
}

//	Read element data from stream Input.Starting from up-left
bool CShell::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet;	// Material property set niumber
	unsigned int N1, N2, N3, N4;	// Left node niumber and right node niumber

	Input >> N1 >> N2 >> N3 >> N4 >> MSet;
	ElementMaterial_ = dynamic_cast<CShellMaterial*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];
	nodes_[3] = &NodeList[N4 - 1];

	return true;
}

//	Write element data to stream
void CShell::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
		<< setw(9) << nodes_[1]->NodeNumber
		<< setw(9) << nodes_[2]->NodeNumber
		<< setw(9) << nodes_[3]->NodeNumber
		<< setw(12) << ElementMaterial_->nset << endl;
}

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CShell::ElementStiffness(double* stfMatrix)
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
	double a;
	double b;
	a = 0.5 * (X[1] - X[0]);
	b = 0.5 * (Y[2] - Y[1]);
	//	Calculate element stiffness matrix

	CShellMaterial* material_ = dynamic_cast<CShellMaterial*>(ElementMaterial_);	// Pointer to material of the element

	double E = material_->E;
	double h = material_->h;
	double niu = material_->niu;
	double rho = material_->rho;
	/*double k = material_->k;
	double p = material_->p;
	double a = material_->a;
	double b = material_->b;*/

	Vector2d s4q;
	s4q << -sqrt(1.0/3.0), sqrt(1.0 / 3.0); // eta gaussian

	Vector2d r4q;
	r4q << -sqrt(1.0 / 3.0), sqrt(1.0 / 3.0); // ksi gaussian

	Matrix2d af4q;
	af4q << 1., 1.,
		1., 1.;

	//单元刚度矩阵
	MatrixXd Ke4q;
	Ke4q = MatrixXd::Zero(8, 8);

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
			
			Jacobi(0, 0) = 0.25 * ((X[1] - X[0]) * (1 - s4q(ids)) - (X[3] - X[2]) * (1 + s4q(ids)));
			Jacobi(1, 0) = 0.25 * ((X[3] - X[0]) * (1 - r4q(idr)) + (X[2] - X[1]) * (1 + r4q(idr)));
			Jacobi(0, 1) = 0.25 * ((Y[1] - Y[0]) * (1 - s4q(ids)) - (Y[3] - Y[2]) * (1 + s4q(ids)));
			Jacobi(1, 1) = 0.25 * ((Y[3] - Y[0]) * (1 - r4q(idr)) + (Y[2] - Y[1]) * (1 + r4q(idr)));

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
			GQ << s4q(ids) - 1, 1 - s4q(ids), 1 + s4q(ids), -s4q(ids) - 1,
				r4q(idr) - 1, -r4q(idr) - 1, 1 + r4q(idr), 1 - r4q(idr);
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
			C1 << 1, niu, 0,
				niu, 1, 0,
				0, 0, (1 - niu) / 2;

			C1 = C1 * E / (1 - niu * niu);

			Ke4q = Ke4q + Bij.transpose() * C1 * Bij * Jacobi.determinant() * af4q(ids, idr);

		}
	}

	double det;
	det = Ke4q.determinant();

	////////////板部分

	//积分点

	Vector3d sm;
	sm << 0.0, -sqrt(3.0 / 5.0), sqrt(3.0 / 5.0); // eta gaussian

	Vector3d rm;
	rm << 0.0, -sqrt(3.0 / 5.0), sqrt(3.0 / 5.0); // ksi gaussian

	Matrix3d afm;
	double t1 = 8.0 / 9.0;
	double t2 = 5.0 / 9.0;
	double t3 = 5.0 / 9.0;
	afm << t1 * t1, t1* t2, t1* t2,
		t2* t1, t2* t2, t2* t3,
		t3* t1, t3* t2, t3* t3;

	MatrixXd KeM;
	Matrix3d D;
	D << 1.0, niu, 0.,
		niu, 1.0, 0.,
		0., 0., (1.0 - niu) / 2.0;
	double D0 = E * h * h * h / (12.0 * (1 - niu * niu));
	D = D0 * D;
	KeM = MatrixXd::Zero(12, 12);
	for (int ids = 0; ids < 3; ids++) {
		for (int idr = 0; idr < 3; idr++) {

			KeM = KeM + afm(ids, idr) * Generatessj(a,b,h, sm(ids), rm(idr), D);
		}
	}




	// 平板壳总体单元刚度矩阵
	Matrix<double, 24, 24> KE;
	KE = MatrixXd::Zero(24, 24);
	Matrix<double, 1, 8> LM4q;
	Matrix<double, 1, 12> LMM;
	LM4q << 1, 2, 7, 8, 13, 14, 19, 20;
		
	LMM << 3, 4, 5, 9, 10, 11, 15, 16, 17, 21, 22, 23;
	
	//分别组装4Q和M
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 8; j++) {
			int x;
			int y;
			x = int(LM4q(0, i)) - 1;
			y = int(LM4q(0, j)) - 1;
			KE(x, y) += Ke4q(i, j);
			
		}
	}
	for (int i = 0; i < 12; i++)
	{
		for (int j = 0; j < 12; j++) {
			int x;
			int y;
			x = int(LMM(0, i)) -1 ;
			y = int(LMM(0, j)) -1 ;
			KE(x, y) += KeM(i, j);

		}
	}
	
	//存储单元刚度矩阵
	int flag = 0;
	for (int j = 0; j < 24; j++)
	{
		for (int i = j; i >= 0; i--)
		{
			stfMatrix[flag] = KE(i, j);
			flag++;
		}
	}

	/*stfMatrix[0] = Ke4q(0, 0);
	stfMatrix[1] = Ke4q(1, 1);
	stfMatrix[2] = Ke4q(0, 1);
	stfMatrix[3] = Ke4q(2, 2);
	stfMatrix[4] = Ke4q(1, 2);
	stfMatrix[5] = Ke4q(0, 2);
	stfMatrix[6] = Ke4q(3, 3);
	stfMatrix[7] = Ke4q(2, 3);
	stfMatrix[8] = Ke4q(1, 3);
	stfMatrix[9] = Ke4q(0, 3);
	stfMatrix[10] = Ke4q(4, 4);
	stfMatrix[11] = Ke4q(3, 4);
	stfMatrix[12] = Ke4q(2, 4);
	stfMatrix[13] = Ke4q(1, 4);
	stfMatrix[14] = Ke4q(0, 4);
	stfMatrix[15] = Ke4q(5, 5);
	stfMatrix[16] = Ke4q(4, 5);
	stfMatrix[17] = Ke4q(3, 5);
	stfMatrix[18] = Ke4q(2, 5);
	stfMatrix[19] = Ke4q(1, 5);
	stfMatrix[20] = Ke4q(0, 5);
	stfMatrix[21] = Ke4q(6, 6);
	stfMatrix[22] = Ke4q(5, 6);
	stfMatrix[23] = Ke4q(4, 6);
	stfMatrix[24] = Ke4q(3, 6);
	stfMatrix[25] = Ke4q(2, 6);
	stfMatrix[26] = Ke4q(1, 6);
	stfMatrix[27] = Ke4q(0, 6);
	stfMatrix[28] = Ke4q(7, 7);
	stfMatrix[29] = Ke4q(6, 7);
	stfMatrix[30] = Ke4q(5, 7);
	stfMatrix[31] = Ke4q(4, 7);
	stfMatrix[32] = Ke4q(3, 7);
	stfMatrix[33] = Ke4q(2, 7);
	stfMatrix[34] = Ke4q(1, 7);
	stfMatrix[35] = Ke4q(0, 7);*/
	/*for (int i = 0; i < 250; i++) {
		cout << stfMatrix[i] << endl;
	}*/

}

//	Calculate element stress 
void CShell::ElementStress(double* stress, double* Displacement)
{
	CShellMaterial* material_ = dynamic_cast<CShellMaterial*>(ElementMaterial_);	// Pointer to material of the element
	double sigamax1;
	double X[4];
	double Y[4];
	double Z[4];//	dx = X[1]-X[0], dy = Y[1]-Y[0], dz = z2-z1
	for (unsigned int i = 0; i <= 3; i++)
	{
		X[i] = nodes_[i]->XYZ[0];
		Y[i] = nodes_[i]->XYZ[1];
		Z[i] = nodes_[i]->XYZ[2];
	}
	double a;
	double b;
	a = 0.5 * (X[1] - X[0]);
	b = 0.5 * (Y[2] - Y[1]);
	//	Calculate element stiffness matrix

		// Pointer to material of the element

	double E = material_->E;
	double h = material_->h;
	double niu = material_->niu;
	double rho = material_->rho;
	Matrix3d D;
	D << 1.0, niu, 0.,
		niu, 1.0, 0.,
		0., 0., (1.0 - niu) / 2.0;
	double D0 = E * h * h * h / (12.0 * (1 - niu * niu));
	D = D0 * D;
	Matrix<double, 3, 12> D1;
	Matrix<double, 3, 12> D2;
	Matrix<double, 3, 12> D3;
	Matrix<double, 3, 12> D4;

	Matrix<double, 1, 3> M1;
	Matrix<double, 1, 3> M2;
	Matrix<double, 1, 3> M3;
	Matrix<double, 1, 3> M4;

	Vector2d s1;
	s1 << -1,1; // 

	Vector2d r1;
	r1 << -1,1; // 

	Matrix<double, 1, 12> wtheaxy;
	wtheaxy(0, 0) = Displacement[2];
	wtheaxy(0, 1) = Displacement[3];
	wtheaxy(0, 2) = Displacement[4];

	wtheaxy(0, 3) = Displacement[8];
	wtheaxy(0, 4) = Displacement[9];
	wtheaxy(0, 5) = Displacement[10];

	wtheaxy(0, 6) = Displacement[14];
	wtheaxy(0, 7) = Displacement[15];
	wtheaxy(0, 8) = Displacement[16];

	wtheaxy(0, 9) = Displacement[20];
	wtheaxy(0, 10) = Displacement[21];
	wtheaxy(0, 11) = Displacement[22];
	



	D1 = D * GenerateBm(a, b, h, -1,-1);
	D2 = D * GenerateBm(a, b, h, 1, -1);
	D3 = D * GenerateBm(a, b, h, 1, 1);
	D4 = D * GenerateBm(a, b, h, 1, -1);
	

	M1 = D1 * wtheaxy.transpose();
	M2 = D2 * wtheaxy.transpose();
	M3 = D3 * wtheaxy.transpose();
	M4 = D4 * wtheaxy.transpose();

	double sigamax2;
	double sigamax3;
	double sigamax4;

	double sigamay1;
	double sigamay2;
	double sigamay3;
	double sigamay4;


	sigamax1 = 12.0 * M1(0, 1) * (h / 2.0) / pow(h, 3);
	sigamax2 = 12.0 * M2(0, 1) * (h / 2.0) / pow(h, 3);
	sigamax3 = 12.0 * M3(0, 1) * (h / 2.0) / pow(h, 3);
	sigamax4 = 12.0 * M4(0, 1) * (h / 2.0) / pow(h, 3);

	sigamay1 = 12.0 * M1(0, 0) * (h / 2.0) / pow(h, 3);
	sigamay2 = 12.0 * M1(0, 0) * (h / 2.0) / pow(h, 3);
	sigamay3 = 12.0 * M1(0, 0) * (h / 2.0) / pow(h, 3);
	sigamay4 = 12.0 * M4(0, 0) * (h / 2.0) / pow(h, 3);


	*stress = 0.0;
	stress[0] = sigamax1 ;
	stress[1] = sigamay1;
	stress[6] = sigamax2;
	stress[7] = sigamay2;
	stress[12] = sigamax3;
	stress[13] = sigamay3;
	stress[18] = sigamax4;
	stress[19] = sigamay4;
	

	/*for (unsigned int i = 0; i < 6; i++)
	{
		if (LocationMatrix_[i])
			*stress += E * Displacement[LocationMatrix_[i] - 1];
	}*/
}




	
	


