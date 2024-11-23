#include "Beam.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

//	Constructor
CBeam::CBeam()
{
	NEN_ = 2;	// Each element has 2 nodes
	nodes_ = new CNode*[NEN_];
    
    ND_ = 12;
    LocationMatrix_ = new unsigned int[ND_];
	ElementMaterial_ = nullptr;

	ElementNDF = 6;
}

//	Desconstructor
CBeam::~CBeam()
{
}

//	Read element data from stream Input
bool CBeam::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet;	// Material property set number
	unsigned int N1, N2;	// Left node number and right node number

	Input >> N1 >> N2 >> MSet;
    ElementMaterial_ = dynamic_cast<CBeamMaterial*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];

	return true;
}

//	Write element data to stream
void CBeam::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber << setw(12) << ElementMaterial_->nset << endl;
}

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CBeam::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

	CBeamMaterial* material_ = dynamic_cast<CBeamMaterial*>(ElementMaterial_);	// Pointer to material of the element
	double E = material_->E;
	double A = material_->Area;
	double l;
	double Iz = material_->Iz;
	double Iy = material_->Iy;
	double v = material_->v;
	double G = E / (2 * (1 + v));
	double Ij = Iz + Iy;

//  calculate element length.
	CNode** nodes = this->GetNodes();
	double * XYZ1 = nodes[0] -> XYZ;
	double * XYZ2 = nodes[1] -> XYZ;
	l = sqrt(pow(XYZ1[0] - XYZ2[0], 2) + pow(XYZ1[1] - XYZ2[1], 2) + pow(XYZ1[2] - XYZ2[2], 2));

//  calculate stiffness matrix.
	MatrixXd K;
	K = MatrixXd::Zero(12,12);
	MatrixXd R;
	R = MatrixXd::Zero(12,12);

	K(0,0) = E * A / l;
	K(1,1) = 12 * E * Iz / pow(l, 3);
	K(2,2) = 12 * E * Iy / pow(l, 3);
	K(3,3) = G * Ij / l;
	K(4,4) = 4 * E * Iy / l;
	K(2,4) = - 6 * E * Iy / pow(l, 2); K(4,2) = K(2,4);
	K(5,5) = 4 * E * Iz / l;
	K(1,5) = 6 * E * Iz / pow(l, 2); K(5,1) = K(1,5);
	K(6,6) = E * A / l;
	K(0,6) = - E * A / l; K(6,0) = K(0,6);
	K(7,7) = 12 * E * Iz / pow(l, 3);
	K(5,7) = - 6 * E * Iz / pow(l, 2); K(7,5) = K(5,7);
	K(1,7) = - 12 * E * Iz / pow(l, 3); K(7,1) = K(1,7);
	K(8,8) = 12 * E * Iy / pow(l, 3);
	K(4,8) = 6 * E * Iy / pow(l, 2); K(8,4) = K(4,8);
	K(2,8) = - 12 * E * Iy / pow(l, 3); K(8,2) = K(2,8);
	K(9,9) = G * Ij / l;
	K(3,9) = - G * Ij / l; K(9,3) = K(3,9);
	K(10,10) = 4 * E * Iz / l;
	K(8,10) = 6 * E * Iz / pow(l, 2); K(10,8) = K(8,10);
	K(4,10) = 2 * E * Iz / l; K(10,4) = K(4,10);
	K(2,10) = - 6 * E * Iz / pow(l, 2); K(10,2) = K(2,10);
	K(11,11) = 4 * E * Iz / l;
	K(7,11) = - 6 * E * Iz / pow(l,2); K(11,7) = K(7,11);
	K(5,11) = 2 * E * Iz / l; K(11,5) = K(5,11);
	K(1,11) = 6 * E * Iz / pow(l, 2); K(11,1) = K(1,11);

//	calculate transpose matrix.

	MatrixXd r; // r(i,:) means the ith direction.
	r = MatrixXd::Zero(3,3);
	
	r(0,0) = XYZ2[0] - XYZ1[0];
	r(0,1) = XYZ2[1] - XYZ1[1];
	r(0,2) = XYZ2[2] - XYZ1[2];
	r(1,0) = -r(0,1);
	r(1,1) = r(0,0);
	r(1,2) = 0;
	Vector3d u(r(0,0), r(0,1), r(0,2));
  	Vector3d w(r(1,0), r(1,1), r(0,2));
	Vector3d wv = u.cross(w);
	
	r(2,0) = wv(0);
	r(2,1) = wv(1);
	r(2,2) = wv(2);

	for(int i=0; i<12; i++) {
		R(i,i) = 1;
	}
	for(int i=0; i<3; i++) {
		for(int j=0; j<3; j++) {
			double a = r(i,j) / sqrt(pow(r(i,0),2) + pow(r(i,1),2) + pow(r(i,2),2)); // a(i,j)
			R(i,j) = a;
			R(i+6, j+6) = a;
		}
	}

//  calculate element stiffness matrix.
	MatrixXd M;
	M = R.transpose() * K * R;
	int flag = 0;
	for(int j=0; j<12; j++) {
		for(int i=j; i>=0; i--) {
			Matrix[flag] = M(i, j);
			flag++;
		}
	}
}

//	Calculate element stress 
void CBeam::ElementStress(double* stress, double* Displacement)
{
	CBeamMaterial* material_ = dynamic_cast<CBeamMaterial*>(ElementMaterial_);	// Pointer to material of the element
	
	double DX[3];	//	dx = x2-x1, dy = y2-y1, dz = z2-z1
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
