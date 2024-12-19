/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "C3D20.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <Eigen/Dense>



using namespace Eigen;
using namespace std;



//	Constructor
C3D20::C3D20()
{
	NEN_ = 20;	// Each element has 20 nodes
	nodes_ = new CNode*[NEN_];
    
    ND_ = 60; 
    LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;

	ElementNDF = 3;
}

//	Desconstructor
C3D20::~C3D20()
{
}

//	Read element data from stream Input
bool C3D20::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet;	// Material property set number
	unsigned int N1, N2, N3, N4, N5, N6, N7, N8, N9, N10, N11, N12, N13, N14, N15, N16, N17, N18, N19, N20;	// The four node numbers

	Input >> N1 >> N2 >> N3 >> N4 >> N5 >> N6 >> N7 >> N8 >> N9 >> N10 >> N11 >> N12 >> N13 >> N14 >> N15 >> N16 >> N17 >> N18 >> N19 >> N20 >> MSet;
    ElementMaterial_ = dynamic_cast<C3D20Material*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];
	nodes_[3] = &NodeList[N4 - 1];
	nodes_[4] = &NodeList[N5 - 1];
	nodes_[5] = &NodeList[N6 - 1];
	nodes_[6] = &NodeList[N7 - 1];
	nodes_[7] = &NodeList[N8 - 1];
    nodes_[8] = &NodeList[N9 - 1];
    nodes_[9] = &NodeList[N10 - 1];
    nodes_[10] = &NodeList[N11 - 1];
    nodes_[11] = &NodeList[N12 - 1];
    nodes_[12] = &NodeList[N13 - 1];
    nodes_[13] = &NodeList[N14 - 1];
    nodes_[14] = &NodeList[N15 - 1];
    nodes_[15] = &NodeList[N16 - 1];
    nodes_[16] = &NodeList[N17 - 1];
    nodes_[17] = &NodeList[N18 - 1];
    nodes_[18] = &NodeList[N19 - 1];
    nodes_[19] = &NodeList[N20 - 1];
	return true;
}

//	Write element data to stream
void C3D20::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber 
		   << setw(9) << nodes_[2]->NodeNumber 
		   << setw(9) << nodes_[3]->NodeNumber
		   << setw(9) << nodes_[4]->NodeNumber
		   << setw(9) << nodes_[5]->NodeNumber
		   << setw(9) << nodes_[6]->NodeNumber
		   << setw(9) << nodes_[7]->NodeNumber   
           << setw(9) << nodes_[8]->NodeNumber
           << setw(9) << nodes_[9]->NodeNumber  
           << setw(9) << nodes_[10]->NodeNumber
           << setw(9) << nodes_[11]->NodeNumber
           << setw(9) << nodes_[12]->NodeNumber
           << setw(9) << nodes_[13]->NodeNumber
           << setw(9) << nodes_[14]->NodeNumber
           << setw(9) << nodes_[15]->NodeNumber
           << setw(9) << nodes_[16]->NodeNumber
           << setw(9) << nodes_[17]->NodeNumber
           << setw(9) << nodes_[18]->NodeNumber
           << setw(9) << nodes_[19]->NodeNumber
		   << setw(12) << ElementMaterial_->nset << endl;
}

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void C3D20::ElementMass(double* concentMassMatrix)
{
	clear(concentMassMatrix, SizeOfMassMatrix());

	//	Calculate bar length
	double X[20];
	double Y[20];
	double Z[20];		//	dx = X[1]-X[0], dy = Y[1]-Y[0], dz = z2-z1
	for (unsigned int i = 0; i <= 19; i++)
	{
		X[i] = nodes_[i]->XYZ[0];
		Y[i] = nodes_[i]->XYZ[1];
		Z[i] = nodes_[i]->XYZ[2];
	}

	//	Calculate element stiffness matrix

	C3D20Material* material_ = dynamic_cast<C3D20Material*>(ElementMaterial_);	// Pointer to material of the element(gaidong)

	double E = material_->E;
	//std::cout << E << std::endl;
	double niu = material_->niu;
	double rho = material_->rho;

	Vector3d ss;
	ss << sqrt(3.) / sqrt(5.), 0, -sqrt(3.) / sqrt(5.);

	Vector3d  rr;
	rr << sqrt(3.) / sqrt(5.), 0, -sqrt(3.) / sqrt(5.);

	Vector3d  tt;
	tt << sqrt(3.) / sqrt(5.), 0, -sqrt(3.) / sqrt(5.);

	Vector3d weight;
	weight << 5. / 9., 8. / 9., 5. / 9.;
	double r;
	double s;
	double t;
	double wi;
	double xi;
	double eta;
	double zeta;
	r = 0.;
	s = 0.;
	t = 0.;
	wi = 0.;
	MatrixXd M;
	M = MatrixXd::Zero(60, 60);

	for (int ids = 0; ids < 3; ids++)
	{
		for (int idr = 0; idr < 3; idr++)
		{
			for (int idt = 0; idt < 3; idt++) {
				r = rr(idr);
				s = ss(ids);
				t = tt(idt);
				xi = ss(ids);
				eta = rr(idr);
				zeta = tt(idt);
				wi = weight(ids) * weight(idr) * weight(idt);
				Matrix<double, 3, 20> G_H8;
				G_H8 << (r + 1) * (t + 1) * (r + s + t - 2) + (r + 1) * (s + 1) * (t + 1), (r - 1)* (t + 1)* (r - s - t + 2) - (r - 1) * (s + 1) * (t + 1), -(r - 1) * (t + 1) * (r + s - t + 2) - (r - 1) * (s - 1) * (t + 1), (r + 1)* (s - 1)* (t + 1) - (r + 1) * (t + 1) * (r - s + t - 2), -(r + 1) * (t - 1) * (r + s - t - 2) - (r + 1) * (s + 1) * (t - 1), (r - 1)* (s + 1)* (t - 1) - (r - 1) * (t - 1) * (r - s + t + 2),
					(r - 1)* (t - 1)* (r + s + t + 2) + (r - 1) * (s - 1) * (t - 1), -(r + 1) * (t - 1) * (s - r + t + 2) - (r + 1) * (s - 1) * (t - 1), -(2 * r * r - 2) * (t + 1), 4 * s * (r - 1) * (t + 1), (2 * r * r - 2)* (t + 1), -4 * s * (r + 1) * (t + 1),
					(2 * r * r - 2)* (t - 1), -4 * s * (r - 1) * (t - 1), -(2 * r * r - 2) * (t - 1), 4 * s * (r + 1) * (t - 1), -(2 * t * t - 2) * (r + 1), (2 * t * t - 2)* (r - 1),
					-(2 * t * t - 2) * (r - 1), (2 * t * t - 2)* (r + 1), (s + 1)* (t + 1)* (r + s + t - 2) + (r + 1) * (s + 1) * (t + 1), (s + 1)* (t + 1)* (r - s - t + 2) + (r - 1) * (s + 1) * (t + 1), -(s - 1) * (t + 1) * (r + s - t + 2) - (r - 1) * (s - 1) * (t + 1), -(s - 1) * (t + 1) * (r - s + t - 2) - (r + 1) * (s - 1) * (t + 1),
					-(s + 1) * (t - 1) * (r + s - t - 2) - (r + 1) * (s + 1) * (t - 1), -(s + 1) * (t - 1) * (r - s + t + 2) - (r - 1) * (s + 1) * (t - 1), (s - 1)* (t - 1)* (r + s + t + 2) + (r - 1) * (s - 1) * (t - 1), (r + 1)* (s - 1)* (t - 1) - (s - 1) * (t - 1) * (s - r + t + 2), -4 * r * (s + 1) * (t + 1), (2 * s * s - 2)* (t + 1),
					4 * r * (s - 1) * (t + 1), -(2 * s * s - 2) * (t + 1), 4 * r * (s + 1) * (t - 1), -(2 * s * s - 2) * (t - 1), -4 * r * (s - 1) * (t - 1), (2 * s * s - 2)* (t - 1),
					-(2 * t * t - 2) * (s + 1), (2 * t * t - 2)* (s + 1), -(2 * t * t - 2) * (s - 1), (2 * t * t - 2)* (s - 1), (r + 1)* (s + 1)* (r + s + t - 2) + (r + 1) * (s + 1) * (t + 1), (r - 1)* (s + 1)* (r - s - t + 2) - (r - 1) * (s + 1) * (t + 1),
					(r - 1)* (s - 1)* (t + 1) - (r - 1) * (s - 1) * (r + s - t + 2), -(r + 1) * (s - 1) * (r - s + t - 2) - (r + 1) * (s - 1) * (t + 1), (r + 1)* (s + 1)* (t - 1) - (r + 1) * (s + 1) * (r + s - t - 2), -(r - 1) * (s + 1) * (r - s + t + 2) - (r - 1) * (s + 1) * (t - 1), (r - 1)* (s - 1)* (r + s + t + 2) + (r - 1) * (s - 1) * (t - 1), -(r + 1) * (s - 1) * (s - r + t + 2) - (r + 1) * (s - 1) * (t - 1),
					-(2 * r * r - 2) * (s + 1), (2 * s * s - 2)* (r - 1), (2 * r * r - 2)* (s - 1), -(2 * s * s - 2) * (r + 1), (2 * r * r - 2)* (s + 1), -(2 * s * s - 2) * (r - 1),
					-(2 * r * r - 2) * (s - 1), (2 * s * s - 2)* (r + 1), -4 * t * (r + 1) * (s + 1), 4 * t * (r - 1) * (s + 1), -4 * t * (r - 1) * (s - 1), 4 * t * (r + 1) * (s - 1);
				G_H8 = 1. / 8. * G_H8;
				Matrix<double, 20, 1> Nt;
				Nt << (1 + xi) * (1 + eta) * (1 + zeta) * (xi + eta + zeta - 2),
					(1 + xi)* (1 - eta)* (1 + zeta)* (xi - eta + zeta - 2),
					(1 - xi)* (1 - eta)* (1 + zeta)* (-xi - eta + zeta - 2),
					(1 - xi)* (1 + eta)* (1 + zeta)* (-xi + eta + zeta - 2),
					(1 + xi)* (1 + eta)* (1 - zeta)* (xi + eta - zeta - 2),
					(1 + xi)* (1 - eta)* (1 - zeta)* (xi - eta - zeta - 2),
					(1 - xi)* (1 - eta)* (1 - zeta)* (-xi - eta - zeta - 2),
					(1 - xi)* (1 + eta)* (1 - zeta)* (-xi + eta - zeta - 2),
					2 * (1 - eta * eta) * (1 + xi) * (1 + zeta),
					2 * (1 - xi * xi) * (1 - eta) * (1 + zeta),
					2 * (1 - eta * eta) * (1 - xi) * (1 + zeta),
					2 * (1 - xi * xi) * (1 + eta) * (1 + zeta),
					2 * (1 - eta * eta) * (1 + xi) * (1 - zeta),
					2 * (1 - xi * xi) * (1 - eta) * (1 - zeta),
					2 * (1 - eta * eta) * (1 - xi) * (1 - zeta),
					2 * (1 - xi * xi) * (1 + eta) * (1 - zeta),
					2 * (1 - zeta * zeta) * (1 + xi) * (1 + eta),
					2 * (1 - zeta * zeta) * (1 + xi) * (1 - eta),
					2 * (1 - zeta * zeta) * (1 - xi) * (1 - eta),
					2 * (1 - zeta * zeta) * (1 - xi) * (1 + eta);
				Nt = 1. / 8. * Nt;
				Matrix<double, 3, 60> N;
				for (int i = 0; i < 20; i++) {
					N(0, i * 3) = Nt(i);
					N(0, i * 3 + 1) = 0;
					N(0, i * 3 + 2) = 0;
					N(1, i * 3) = 0;
					N(1, i * 3 + 1) = Nt(i);
					N(1, i * 3 + 2) = 0;
					N(2, i * 3) = 0;
					N(2, i * 3 + 1) = 0;
					N(2, i * 3 + 2) = Nt(i);
				}
				Matrix<double, 20, 3> XYZe;
				XYZe << X[0], Y[0], Z[0],
					X[1], Y[1], Z[1],
					X[2], Y[2], Z[2],
					X[3], Y[3], Z[3],
					X[4], Y[4], Z[4],
					X[5], Y[5], Z[5],
					X[6], Y[6], Z[6],
					X[7], Y[7], Z[7],
					X[8], Y[8], Z[8],
					X[9], Y[9], Z[9],
					X[10], Y[10], Z[10],
					X[11], Y[11], Z[11],
					X[12], Y[12], Z[12],
					X[13], Y[13], Z[13],
					X[14], Y[14], Z[14],
					X[15], Y[15], Z[15],
					X[16], Y[16], Z[16],
					X[17], Y[17], Z[17],
					X[18], Y[18], Z[18],
					X[19], Y[19], Z[19];

				Matrix<double, 3, 3> Jacobi;
				Jacobi = G_H8 * XYZe;
				M = M + wi * rho * N.transpose() * N * Jacobi.determinant();
			}
		}
	}
	for (int i = 0; i < 60; i++) {
		for (int j = 0; j < 60; j++) {
			if (i != j) {
				M(i, i) += M(i, j);
				M(i, j) = 0;
			}
		}
	}
	int p = 0;
	for (int i = 0; i < 60; i++) {
		for (int j = i; j >= 0; j--) {
			concentMassMatrix[p] = M(i, j);
			p++;
		}
	}

}
//calculate element mass
void C3D20::ElementStiffness(double* stfMatrix)
{
	clear(stfMatrix, SizeOfStiffnessMatrix());

//	Calculate bar length
	double X[20];
	double Y[20];
	double Z[20];		//	dx = X[1]-X[0], dy = Y[1]-Y[0], dz = z2-z1
	for (unsigned int i = 0; i <= 19; i++)
	  {
		X[i] = nodes_[i]->XYZ[0];
		Y[i] = nodes_[i]->XYZ[1];
		Z[i] = nodes_[i]->XYZ[2];		
	  }
    
//	Calculate element stiffness matrix

	C3D20Material* material_ = dynamic_cast<C3D20Material*>(ElementMaterial_);	// Pointer to material of the element(gaidong)
	
	double E = material_->E; 
	//std::cout << E << std::endl;
	double niu = material_->niu;
	double rho = material_->rho;

	Vector3d  ss;
	ss << sqrt(3.)/sqrt(5.), 0, -sqrt(3.)/sqrt(5.);

	Vector3d  rr;
	rr << sqrt(3.)/sqrt(5.), 0, -sqrt(3.)/sqrt(5.);

	Vector3d  tt;
	tt << sqrt(3.)/sqrt(5.), 0, -sqrt(3.)/sqrt(5.);
    
	Vector3d weight;
	weight << 5./9., 8./9., 5./9.;
    double r ;
	double s ;
	double t ;
	double wi ;
	r = 0.;
	s = 0.;
	t = 0.;
	wi = 0.;
	//单元刚度矩阵
	MatrixXd Ke;
	Ke = MatrixXd::Zero(60, 60); //这里没问题
    // std::cout << Ke << std::endl;
	for (int ids = 0; ids < 3; ids++)
	{
		for (int idr = 0; idr < 3; idr++)		{
			for (int idt = 0; idt < 3; idt++){
				r = rr(idr);
				s = ss(ids);
				t = tt(idt);
                wi = weight(ids)*weight(idr)*weight(idt);
				Matrix<double, 3, 20> G_H8;
				G_H8 << (r + 1)*(t + 1)*(r + s + t - 2) + (r + 1)*(s + 1)*(t + 1), (r - 1)*(t + 1)*(r - s - t + 2) - (r - 1)*(s + 1)*(t + 1), - (r - 1)*(t + 1)*(r + s - t + 2) - (r - 1)*(s - 1)*(t + 1), (r + 1)*(s - 1)*(t + 1) - (r + 1)*(t + 1)*(r - s + t - 2), - (r + 1)*(t - 1)*(r + s - t - 2) - (r + 1)*(s + 1)*(t - 1), (r - 1)*(s + 1)*(t - 1) - (r - 1)*(t - 1)*(r - s + t + 2),
				        (r - 1)*(t - 1)*(r + s + t + 2) + (r - 1)*(s - 1)*(t - 1), - (r + 1)*(t - 1)*(s - r + t + 2) - (r + 1)*(s - 1)*(t - 1), -(2*r*r - 2)*(t + 1), 4*s*(r - 1)*(t + 1), (2*r*r - 2)*(t + 1), -4*s*(r + 1)*(t + 1),
						(2*r*r - 2)*(t - 1), -4*s*(r - 1)*(t - 1), -(2*r*r - 2)*(t - 1), 4*s*(r + 1)*(t - 1), -(2*t*t - 2)*(r + 1), (2*t*t - 2)*(r - 1),
						-(2*t*t - 2)*(r - 1), (2*t*t - 2)*(r + 1), (s + 1)*(t + 1)*(r + s + t - 2) + (r + 1)*(s + 1)*(t + 1), (s + 1)*(t + 1)*(r - s - t + 2) + (r - 1)*(s + 1)*(t + 1), - (s - 1)*(t + 1)*(r + s - t + 2) - (r - 1)*(s - 1)*(t + 1), - (s - 1)*(t + 1)*(r - s + t - 2) - (r + 1)*(s - 1)*(t + 1), 
						- (s + 1)*(t - 1)*(r + s - t - 2) - (r + 1)*(s + 1)*(t - 1), - (s + 1)*(t - 1)*(r - s + t + 2) - (r - 1)*(s + 1)*(t - 1), (s - 1)*(t - 1)*(r + s + t + 2) + (r - 1)*(s - 1)*(t - 1), (r + 1)*(s - 1)*(t - 1) - (s - 1)*(t - 1)*(s - r + t + 2), -4*r*(s + 1)*(t + 1), (2*s*s - 2)*(t + 1),
						4*r*(s - 1)*(t + 1), -(2*s*s - 2)*(t + 1), 4*r*(s + 1)*(t - 1), -(2*s*s - 2)*(t - 1), -4*r*(s - 1)*(t - 1), (2*s*s - 2)*(t - 1),
						-(2*t*t - 2)*(s + 1), (2*t*t - 2)*(s + 1), -(2*t*t - 2)*(s - 1), (2*t*t - 2)*(s - 1), (r + 1)*(s + 1)*(r + s + t - 2) + (r + 1)*(s + 1)*(t + 1), (r - 1)*(s + 1)*(r - s - t + 2) - (r - 1)*(s + 1)*(t + 1),
						(r - 1)*(s - 1)*(t + 1) - (r - 1)*(s - 1)*(r + s - t + 2), - (r + 1)*(s - 1)*(r - s + t - 2) - (r + 1)*(s - 1)*(t + 1), (r + 1)*(s + 1)*(t - 1) - (r + 1)*(s + 1)*(r + s - t - 2), - (r - 1)*(s + 1)*(r - s + t + 2) - (r - 1)*(s + 1)*(t - 1), (r - 1)*(s - 1)*(r + s + t + 2) + (r - 1)*(s - 1)*(t - 1), - (r + 1)*(s - 1)*(s - r + t + 2) - (r + 1)*(s - 1)*(t - 1),
						-(2*r*r - 2)*(s + 1), (2*s*s - 2)*(r - 1), (2*r*r - 2)*(s - 1), -(2*s*s - 2)*(r + 1), (2*r*r - 2)*(s + 1), -(2*s*s - 2)*(r - 1),
						-(2*r*r - 2)*(s - 1), (2*s*s - 2)*(r + 1), -4*t*(r + 1)*(s + 1), 4*t*(r - 1)*(s + 1), -4*t*(r - 1)*(s - 1), 4*t*(r + 1)*(s - 1);
				G_H8 = 1./ 8. * G_H8;

				Matrix<double, 20, 3> XYZe;
				XYZe << X[0], Y[0], Z[0],
					X[1], Y[1], Z[1],
					X[2], Y[2], Z[2],
					X[3], Y[3], Z[3],
					X[4], Y[4], Z[4],
					X[5], Y[5], Z[5],
					X[6], Y[6], Z[6],
					X[7], Y[7], Z[7],
					X[8], Y[8], Z[8],
					X[9], Y[9], Z[9],
					X[10], Y[10], Z[10],
					X[11], Y[11], Z[11],
					X[12], Y[12], Z[12],
					X[13], Y[13], Z[13],
					X[14], Y[14], Z[14],
					X[15], Y[15], Z[15],
					X[16], Y[16], Z[16],
					X[17], Y[17], Z[17],
					X[18], Y[18], Z[18],
					X[19], Y[19], Z[19];

				Matrix<double, 3, 3> Jacobi;
				Jacobi = G_H8 * XYZe;

				Matrix<double, 3, 3> Jacobi_inv;
				Jacobi_inv = Jacobi.inverse();

				Matrix<double, 3, 20> N_H8;
				N_H8 = Jacobi_inv * G_H8;

				Matrix<double, 6, 60> Bij;
				for (int i = 0; i < 20; i++)
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

						Bij(4, i * 3)     = 0;
						Bij(4, i * 3 + 1) = N_H8(2, i);
						Bij(4, i * 3 + 2) = N_H8(1, i);

						Bij(5, i * 3)     = N_H8(2, i);
						Bij(5, i * 3 + 1) = 0;
						Bij(5, i * 3 + 2) = N_H8(0, i);
					}

				
				Matrix<double, 6,6> C;
				C << 1 - niu, niu, niu, 0, 0, 0,
					niu, 1 - niu, niu, 0, 0, 0,
					niu, niu, 1 - niu, 0, 0, 0,
					0, 0, 0, (1 - 2 * niu) / 2, 0, 0,
					0, 0, 0, 0, (1 - 2 * niu) / 2, 0,
					0, 0, 0, 0, 0, (1 - 2 * niu) / 2;

				C = C * E / (1 + niu) / (1 - 2 * niu);

				Ke = Ke + wi * Bij.transpose() * C * Bij * Jacobi.determinant();
			}
		}
	}
	
	int p = 0;
	for (int i = 0; i < 60; i++){
		for (int j = i ; j >= 0; j--){
			stfMatrix[p] = Ke(i, j);
			p++;
		}
	}

}







//	Calculate element stress 
void C3D20::ElementStress(double* stress, double* Displacement)
{
	C3D20Material* material_ = dynamic_cast<C3D20Material*>(ElementMaterial_);	// Pointer to material of the element
	
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
