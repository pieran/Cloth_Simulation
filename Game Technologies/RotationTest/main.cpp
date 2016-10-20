
#include "console.h"
#include <nclgl\Quaternion.h>
#include <nclgl\Vector3.h>
#include <nclgl\Matrix3.h>
#include <functional>
#include <iostream>
#include <iomanip>

#include <Eigen\Eigen.h>
#include <Eigen\Eigenvalues.h>


#define SYSTEM_COLOR_BACKGROUND		black

#define SYSTEM_COLOR_DEFAULT		white
#define SYSTEM_COLOR_PASSED			green
#define SYSTEM_COLOR_FAILED			red

#define CONSOLE_CHAR_WIDTH			80
#define INDENT						"    "


const int NUM_POINTS = 4;
const double NUM_POINTS_INV  = (1.0 / 4.0);


using namespace std;

void output_error(const char* error)
{
	setcolor(red, SYSTEM_COLOR_BACKGROUND);
	cout << INDENT << error << endl;
	setcolor(white, SYSTEM_COLOR_BACKGROUND);
}

void output_result(const char* description, bool result)
{
	cout << description;

	size_t len = strlen(description);
	for (size_t i = len; i < CONSOLE_CHAR_WIDTH - 6; ++i)
		cout << ".";

	if (result)
	{
		setcolor(green, SYSTEM_COLOR_BACKGROUND);
		cout << "PASSED" << endl;
	}
	else
	{
		setcolor(red, SYSTEM_COLOR_BACKGROUND);
		cout << "FAILED" << endl;
	}
	setcolor(white, SYSTEM_COLOR_BACKGROUND);

}

bool RunUnitTest(const char* header, const std::function<bool()>& func)
{
	cout << header << endl << "-------------------" << endl;
	bool result = func();

	output_result(INDENT"Final Result", result);
	

	cout << endl;
	return result;
}




Quaternion orig_rot, new_rot;

Vector3 p0[4]{ Vector3(-1.0f, 0.0f, 0.0f), Vector3(1.0f, 0.0f, -1.0f), Vector3(0.0f, 1.0f, 0.0f), Vector3(1.0f, 1.0f, 1.0f) };
Vector3 p[4];
Vector3 d[4]; //defined as Q(p - c) - (p0 - c0)



//Partial Derivatives of quaternion rotation
void DiffQuat1(const Quaternion& q, const Vector3& v, Vector3* out_vec)
{
	/*// dV/dQi
	out_vec[0].x = 2.0f * (v.x*q.x + v.y*q.y + v.z*q.z);
	out_vec[0].y = 2.0f * (v.x*q.y - v.y*q.x - v.z*q.w);
	out_vec[0].z = 2.0f * (v.x*q.z + v.y*q.w - v.z*q.x);

	// dV/dQj
	out_vec[1].x = 2.0f * (-v.x*q.y + v.y*q.x + v.z*q.w);
	out_vec[1].y = 2.0f * ( v.x*q.x + v.y*q.y + v.z*q.z);
	out_vec[1].z = 2.0f * (-v.x*q.w + v.y*q.z - v.z*q.y);

	// dV/dQk
	out_vec[2].x = 2.0f * (-v.x*q.z - v.y*q.w + v.z*q.x);
	out_vec[2].y = 2.0f * ( v.x*q.w - v.y*q.z + v.z*q.y);
	out_vec[2].z = 2.0f * ( v.x*q.x + v.y*q.y + v.z*q.z);

	// dV/dQk
	out_vec[3].x = 2.0f * ( v.x*q.w - v.y*q.z + v.z*q.y);
	out_vec[3].y = 2.0f * ( v.x*q.z + v.y*q.w - v.z*q.x);
	out_vec[3].z = 2.0f * (-v.x*q.y + v.y*q.x + v.z*q.w);*/

	out_vec[0].x = 2.0f * (v.y * q.y + v.z * q.z);
	out_vec[0].y = 2.0f * (v.x * q.y - 2.0f * v.y * q.x - v.z * q.w);
	out_vec[0].z = 2.0f * (v.x * q.z + v.y * q.w - 2.0f * v.z * q.x);

	// dV/dQj
	out_vec[1].x = 2.0f * (-2.0f * v.x * q.y + v.y * q.x + v.z * q.w);
	out_vec[1].y = 2.0f * (v.x * q.x + v.z * q.z);
	out_vec[1].z = 2.0f * (-v.x * q.w + v.y * q.z - 2.0f * v.z * q.y);

	// dV/dQk
	out_vec[2].x = 2.0f * (-2.0f * v.x * q.z - v.y * q.w + v.z * q.x);
	out_vec[2].y = 2.0f * (v.x * q.w - 2.0f * v.y * q.z + v.z * q.y);
	out_vec[2].z = 2.0f * (v.x * q.x + v.y * q.y);

	// dV/dQw
	out_vec[3].x = 2.0f * (-v.y * q.z + v.z * q.y);
	out_vec[3].y = 2.0f * (v.x * q.z - v.z * q.x);
	out_vec[3].z = 2.0f * (-v.x * q.y + v.y * q.x);
}

//Partial Second Derivatives of quaternion rotation (Note. The quaternion is not needed at this point)
void DiffQuat2(const Vector3& v, Vector3* out_vec)	//Out vec is a 4*4 matrix of vector3's
{
	
	/*out_vec[0]  = Vector3(v.x, -v.y, -v.z)  * 2.0f;		// d2V/dQidQi
	out_vec[1]  = Vector3(v.y, v.x, 0.0f)   * 2.0f;		// d2V/dQidQj
	out_vec[2]  = Vector3(v.z, 0.0f, v.x)   * 2.0f;		// d2V/dQidQk
	out_vec[3]  = Vector3(0.0f, -v.z, v.y)  * 2.0f;		// d2V/dQidQw

	out_vec[4]  = Vector3(v.y, v.x, 0.0f)   * 2.0f;		// d2V/dQjdQi
	out_vec[5]  = Vector3(-v.x, v.y, -v.z)  * 2.0f;		// d2V/dQjdQj
	out_vec[6]  = Vector3(0.0f, v.z, v.y)   * 2.0f;		// d2V/dQjdQk
	out_vec[7]  = Vector3(v.z, 0.0f, -v.x)  * 2.0f;		// d2V/dQjdQw

	out_vec[8]  = Vector3(v.z, 0.0f, v.x)   * 2.0f;		// d2V/dQkdQi
	out_vec[9]  = Vector3(0.0f, v.z, v.y)   * 2.0f;		// d2V/dQkdQj
	out_vec[10] = Vector3(-v.x, -v.y, v.z)  * 2.0f;		// d2V/dQkdQk
	out_vec[11] = Vector3(-v.y, v.x, 0.0f)  * 2.0f;		// d2V/dQkdQw

	out_vec[12] = Vector3(0.0f, -v.z, v.y)  * 2.0f;		// d2V/dQwdQi
	out_vec[13] = Vector3(v.z, 0.0f, -v.x)  * 2.0f;		// d2V/dQwdQj
	out_vec[14] = Vector3(-v.y, v.x, 0.0f)  * 2.0f;		// d2V/dQwdQk
	out_vec[15] = Vector3(v.x, v.y, v.z)    * 2.0f;		// d2V/dQwdQw*/


	Vector3 ii = Vector3(0.0f, v.y, v.z)  * -4.0f;		// d2V/dQidQi
	Vector3 jj = Vector3(v.x, 0.0f, v.z)  * -4.0f;		// d2V/dQjdQj
	Vector3 kk = Vector3(v.x, v.y, 0.0f)  * -4.0f;		// d2V/dQkdQk
	Vector3 ww = Vector3(0.0f, 0.0f, 0.0f);				// d2V/dQwdQw

	Vector3 ij = Vector3(v.y, v.x, 0.0f)   * 2.0f;		// d2V/dQidQj
	Vector3 ik = Vector3(v.z, 0.0f, v.x)   * 2.0f;		// d2V/dQidQk
	Vector3 iw = Vector3(0.0f, -v.z, v.y)  * 2.0f;		// d2V/dQidQw

	Vector3 jk = Vector3(0.0f, v.z, v.y)   * 2.0f;		// d2V/dQjdQk
	Vector3 jw = Vector3(v.z, 0.0f, -v.x)  * 2.0f;		// d2V/dQjdQw

	Vector3 kw = Vector3(-v.y, v.x, 0.0f)  * 2.0f;		// d2V/dQkdQw

	out_vec[0] = ii;
	out_vec[1] = ij;
	out_vec[2] = ik;
	out_vec[3] = iw;

	out_vec[4] = ij;
	out_vec[5] = jj;
	out_vec[6] = jk;
	out_vec[7] = jw;

	out_vec[8] = ik;
	out_vec[9] = jk;
	out_vec[10] = kk;
	out_vec[11] = kw;

	out_vec[12] = iw;
	out_vec[13] = jw;
	out_vec[14] = kw;
	out_vec[15] = ww;


}


void pinv(const Eigen::MatrixXf& mat, Eigen::MatrixXf& pinvmat)
{
	auto svd = mat.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);

	double  pinvtoler = 1.e-6; // choose your tolerance wisely!
	auto singularValues_inv = svd.singularValues();
	for (long i = 0; i< mat.cols(); ++i) {
		if (svd.singularValues()(i) > pinvtoler)
			singularValues_inv(i) = 1.0 / svd.singularValues()(i);
		else singularValues_inv(i) = 0;
	}
	pinvmat = (svd.matrixV()*singularValues_inv.asDiagonal()*svd.matrixU().transpose());
}

void FindLeastRotation()
{	
	
	//Assume equal weighting for triangle points, thus the centre is just an avg of the positions
	Vector3 c0 = Vector3(0.0f, 0.0f, 0.0f);
	Vector3 c = Vector3(0.0f, 0.0f, 0.0f);

	for (int i = 0; i < NUM_POINTS; ++i)
	{
		p[i] = Quaternion::RotatePointByQuaternion(orig_rot, p0[i]);
		c0 += p0[i];
		c += p[i];
	}
	c0 *= NUM_POINTS_INV;
	c *= NUM_POINTS_INV;


	//Offset from center points
	Vector3 offset0[NUM_POINTS];
	Vector3 offset[NUM_POINTS];
	for (int i = 0; i < NUM_POINTS; ++i)
	{
		offset0[i] = p0[i] - c0;
		offset[i] = p[i] - c;
	}
	

	//Initial starting values
	Eigen::VectorXd state(5); state.setZero(); 
	state[0] = -0.5f;
	state[1] = -0.2f;
	state[2] = -0.5f;
	state[3] = 0.5f; 
	state[4] = 1.0f;

	/*Quaternion tmp = orig_rot.Conjugate();
	state[0] = tmp.x;
	state[1] = tmp.y;
	state[2] = tmp.z;
	state[3] = tmp.w;
	state[4] = 0.0f;
	*/

	Eigen::VectorXd new_state(5);

	Eigen::VectorXd eFD(5);
	Eigen::MatrixXd eHessian(5, 5);

	Vector3 dDdQi[4];				//Partial Derivatives of dD/dQ
	Vector3 d2DdQiQi[16];			//Partial Second Derivatives of dD/dQ


	//Newtonian solver
	float tolerance = 1e-7f; // 7 digit accuracy is desired
	float epsilon = 1e-14f; // Don't want to divide by a number smaller than this

	int maxIterations = 10000; // Don't allow the iterations to continue indefinitely
	bool haveWeFoundSolution = false; // Have not converged to a solution yet

	for (;;)// int i = 0; i < maxIterations; ++i)
	{
		//Compute D at current rotation
		Quaternion q = Quaternion(state[0], state[1], state[2], state[3]);
		//q.Normalise();
		/*state[0] = q.x;
		state[1] = q.y;
		state[2] = q.z;
		state[3] = q.w;*/
		for (int j = 0; j < NUM_POINTS; ++j)
		{

			d[j] = Quaternion::RotatePointByQuaternion(q, offset[j]) - offset0[j];
		}

		//Calculate the first derivatives 
		//----------------------------------------------------
			eFD.setZero();

			//Integrate over area - In our case, just 1/3 of each
			for (int j = 0; j < NUM_POINTS; ++j)
			{
				DiffQuat1(q, d[j], dDdQi);
				eFD[0] += Vector3::Dot(dDdQi[0], d[j]);
				eFD[1] += Vector3::Dot(dDdQi[1], d[j]);
				eFD[2] += Vector3::Dot(dDdQi[2], d[j]);
				eFD[3] += Vector3::Dot(dDdQi[3], d[j]);
			}
			
			eFD[0] = eFD[0] * NUM_POINTS_INV +state[4] * state[0];
			eFD[1] = eFD[1] * NUM_POINTS_INV +state[4] * state[1];
			eFD[2] = eFD[2] * NUM_POINTS_INV +state[4] * state[2];
			eFD[3] = eFD[3] * NUM_POINTS_INV +state[4] * state[3];
			eFD *= 2.0;

			eFD[4] = Quaternion::Dot(q, q) - 1.0;
		//


		//Calculate the second derivatives
		//----------------------------------------------------
			eHessian.setZero();

			for (int j = 0; j < NUM_POINTS; ++j)
			{
				DiffQuat1(q, d[j], dDdQi);
				DiffQuat2(d[j], d2DdQiQi);

				for (int k = 0; k < 4; ++k)
				{
					for (int l = 0; l < 4; ++l)
					{
						double val = Vector3::Dot(dDdQi[k], dDdQi[l]);
						val += Vector3::Dot(d2DdQiQi[k * 4 + l], d[j]);

						eHessian(k, l) = val;
					}
				}
			}

			eHessian *= NUM_POINTS_INV;
			
			for (int k = 0; k < 4; ++k)
			{
				eHessian(k, k) += state[4];				

				eHessian(4, k) = state[k];
				eHessian(k, 4) = state[k];
			}
			eHessian *= 2.0;
			//eHessian(4, 4) = 2.0f / state[4];			//1.0f; //?????????????????
		//


		//Update Quaternion and Lagrange constraint
		//----------------------------------------------------
		/*Eigen::MatrixXf invHessian = eHessian.inverse();
		new_state = state - invHessian * eFD;*/
		Eigen::MatrixXd invHessian = eHessian.inverse();
		//pinv(eHessian, invHessian);
		Eigen::VectorXd state_change = invHessian * eFD;
		new_state[0] = state[0] - state_change[0];// *state_change[4];
		new_state[1] = state[1] - state_change[1];// *state_change[4];
		new_state[2] = state[2] - state_change[2];// *state_change[4];
		new_state[3] = state[3] - state_change[3];// *state_change[4];
		new_state[4] = state[4] - state_change[4];
		//new_state = state - eHessian.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(eFD);
		cout << "The Old State is: " << endl;
		cout << "\t" << state[0] << ", " << state[1] << ", " << state[2] << ", " << state[3] << ", " << state[4] << endl;

		cout << "The New State is: " << endl;
		cout << "\t" << new_state[0] << ", " << new_state[1] << ", " << new_state[2] << ", " << new_state[3] << ", " << new_state[4] << endl;

		cout << "The First Derivatives are: " << endl;
		cout << "\t" << eFD[0] << ", " << eFD[1] << ", " << eFD[2] << ", " << eFD[3] << ", " << eFD[4] << endl;
		
		
		cout << endl << "The Hessian Matrix is: " << endl;
		for (int j = 0; j < 5; ++j)
			cout << "\t" << eHessian(j, 0) << ", " << eHessian(j, 1) << ", " << eHessian(j, 2) << ", " << eHessian(j, 3) << ", " << eHessian(j, 4) << endl;
		//

		cout << endl << "The Inverse Hessian Matrix is: " << endl;
		for (int j = 0; j < 5; ++j)
			cout << "\t" << invHessian(j, 0) << ", " << invHessian(j, 1) << ", " << invHessian(j, 2) << ", " << invHessian(j, 3) << ", " << invHessian(j, 4) << endl;
		//
		
		//Compute Energy
		//----------------------------------------------------
		float energy = 0.0f;
		for (int j = 0; j < NUM_POINTS; ++j)
		{
			energy += Vector3::Dot(d[j], d[j]);
		}
		energy *= NUM_POINTS_INV;
		energy += state[4] * eFD[4];
	//	energy = eFD[0] + eFD[1] + eFD[2] + eFD[3];

		cout << "Final Energy in system: " << energy << endl << endl << endl;
		if (abs(energy) < tolerance)
		{
			haveWeFoundSolution = true;
			//break;
		}

		state = new_state;
	}


	if (haveWeFoundSolution)
	{
		Quaternion q = Quaternion(state[0], state[1], state[2], state[3]);

		cout << "Final Rotation Matrix found is: " << endl;
		Matrix4 rot = q.ToMatrix4();
		for (int i = 0; i < 4; ++i)
			cout << "\t" << rot.values[i * 4] << ", " << rot.values[i * 4 + 1] << ", " << rot.values[i * 4 + 2] << ", " << rot.values[i * 4 + 3] << endl;

		cout << endl << "The Final Quaternion was: " << endl;
		cout << "\t" << q.x << ", " << q.y << ", " << q.z << ", " << q.w << endl;


		cout << endl << "The original Rotation Matrix was: " << endl;
		rot = orig_rot.ToMatrix4();
		for (int i = 0; i < 4; ++i)
			cout << "\t" << rot.values[i * 4] << ", " << rot.values[i * 4 + 1] << ", " << rot.values[i * 4 + 2] << ", " << rot.values[i * 4 + 3] << endl;

		cout << endl << "The original Quaternion was: " << endl;
		cout << "\t" << orig_rot.x << ", " << orig_rot.y << ", " << orig_rot.z << ", " << orig_rot.w << endl;

		cout << endl << endl;

	}
	else
	{
		output_error("Did not converge or not enough iterations!");
	}
}
















bool ut_QuaternionRotation()
{
	Matrix3 mat_rot = orig_rot.ToMatrix3();

	Vector3 pm[3];
	Vector3 pq[3];
	bool identical = true;
	for (int i = 0; i < 3; ++i)
	{
		pm[i] = mat_rot * p0[i];
		pq[i] = Quaternion::RotatePointByQuaternion(orig_rot, p0[i]);	//Equiv to Q . v. Q'

		if (identical)
			identical = pm[i] == pq[i];
	}

	return identical;
}

bool ut_QuaternionRotationComparison()
{
	Vector3 pm[NUM_POINTS];
	Vector3 pq[NUM_POINTS];
	bool identical = true;
	for (int i = 0; i < NUM_POINTS; ++i)
	{
		pm[i] = Quaternion::RotatePointByQuaternion(new_rot, p0[i]);
		pq[i] = Quaternion::RotatePointByQuaternion(orig_rot, p0[i]);	//Equiv to Q . v. Q'

		if (identical)
			identical = pm[i] == pq[i];
	}

	return identical;
}

bool ut_NewtonianSolver()
{
	cout << "Attempting to solve x^2 - 2" << endl;
	//These choices depend on the problem being solved
	float x0 = 1.0f, x1 = x0; // The initial value
	std::function<float(float)> f = [](float x) { return x * x - 2; };			//The function whose root we are trying to find
	std::function<float(float)> fprime = [](float x) { return 2 * x; };          //The derivative of f(x)

	std::function<float(float)> f2 = [](float x) { return x * x * x + 4; };			//The function whose root we are trying to find
	std::function<float(float)> f2prime = [](float x) { return 3 * x * x; };          //The derivative of f(x)


	float tolerance = 1e-7f; // 7 digit accuracy is desired
	float epsilon = 1e-14f; // Don't want to divide by a number smaller than this

	int maxIterations = 100; // Don't allow the iterations to continue indefinitely
	bool haveWeFoundSolution = false; // Have not converged to a solution yet

	for (int i = 0; i < maxIterations; ++i)
	{
		float y = f(x0) - f2(x0);
		float yprime = fprime(x0) - f2prime(x0);

		if (abs(yprime) < epsilon)	//Dont want to divide by too small of a number
		{
			output_error("Denominator is too small!");
			return false;
		}

		x1 = x0 - y / yprime;                                //Do Newton's computation

		if ((abs(x1 - x0) <= tolerance * abs(x1))) // If the result is within the desired tolerance
		{
			haveWeFoundSolution = true;
			break;                                        //Done, so leave the loop
		}

		x0 = x1;                                           //Update x0 to start the process again

	}

	if (haveWeFoundSolution)
	{
		cout << "Answer is x = " << x1 << endl;
		return true;
	}
	else
	{
		output_error("Did not converge or not enough iterations!");
		return false;
	}
	
}


// Slightly modified version of  Stan Melax's code for 3x3 matrix diagonalization (Thanks Stan!)
// source: http://www.melax.com/diag.html?attredirects=0
typedef double Real;
void Diagonalize(const Matrix3&  A, Matrix3& Q, Matrix3& D)
{
	// A must be a symmetric matrix.
	// returns Q and D such that 
	// Diagonal matrix D = QT * A * Q;  and  A = Q*D*QT
	const int maxsteps = 24;  // certainly wont need that many.
	int k0, k1, k2;
	Real o[3], m[3];
	Real q[4] = { 0.0,0.0,0.0,1.0 };
	Real jr[4];
	Real sqw, sqx, sqy, sqz;
	Real tmp1, tmp2, mq;
	Real AQ[3][3];
	Real thet, sgn, t, c;
	for (int i = 0; i < maxsteps; ++i)
	{
		// quat to matrix
		sqx = q[0] * q[0];
		sqy = q[1] * q[1];
		sqz = q[2] * q[2];
		sqw = q[3] * q[3];
		Q[0] = (sqx - sqy - sqz + sqw);
		Q[3+1] = (-sqx + sqy - sqz + sqw);
		Q[6 + 2] = (-sqx - sqy + sqz + sqw);
		tmp1 = q[0] * q[1];
		tmp2 = q[2] * q[3];
		Q[3 + 0] = 2.0 * (tmp1 + tmp2);
		Q[0 + 1] = 2.0 * (tmp1 - tmp2);
		tmp1 = q[0] * q[2];
		tmp2 = q[1] * q[3];
		Q[6 + 0] = 2.0 * (tmp1 - tmp2);
		Q[0 + 2] = 2.0 * (tmp1 + tmp2);
		tmp1 = q[1] * q[2];
		tmp2 = q[0] * q[3];
		Q[6 + 1] = 2.0 * (tmp1 + tmp2);
		Q[3 + 2] = 2.0 * (tmp1 - tmp2);

		// AQ = A * Q
		AQ[0][0] = Q[0+0] * A[0+0] + Q[3+0] * A[0+1] + Q[6+0] * A[0+2];
		AQ[0][1] = Q[0+1] * A[0+0] + Q[3+1] * A[0+1] + Q[6+1] * A[0+2];
		AQ[0][2] = Q[0+2] * A[0+0] + Q[3+2] * A[0+1] + Q[6+2] * A[0+2];
		AQ[1][0] = Q[0+0] * A[0+1] + Q[3+0] * A[3+1] + Q[6+0] * A[3+2];
		AQ[1][1] = Q[0+1] * A[0+1] + Q[3+1] * A[3+1] + Q[6+1] * A[3+2];
		AQ[1][2] = Q[0+2] * A[0+1] + Q[3+2] * A[3+1] + Q[6+2] * A[3+2];
		AQ[2][0] = Q[0+0] * A[0+2] + Q[3+0] * A[3+2] + Q[6+0] * A[6+2];
		AQ[2][1] = Q[0+1] * A[0+2] + Q[3+1] * A[3+2] + Q[6+1] * A[6+2];
		AQ[2][2] = Q[0+2] * A[0+2] + Q[3+2] * A[3+2] + Q[6+2] * A[6+2];
		// D = Qt * AQ
		D[0+0] = AQ[0][0] * Q[0+0] + AQ[1][0] * Q[3+0] + AQ[2][0] * Q[6+0];
		D[0+1] = AQ[0][0] * Q[0+1] + AQ[1][0] * Q[3+1] + AQ[2][0] * Q[6+1];
		D[0+2] = AQ[0][0] * Q[0+2] + AQ[1][0] * Q[3+2] + AQ[2][0] * Q[6+2];
		D[3+0] = AQ[0][1] * Q[0+0] + AQ[1][1] * Q[3+0] + AQ[2][1] * Q[6+0];
		D[3+1] = AQ[0][1] * Q[0+1] + AQ[1][1] * Q[3+1] + AQ[2][1] * Q[6+1];
		D[3+2] = AQ[0][1] * Q[0+2] + AQ[1][1] * Q[3+2] + AQ[2][1] * Q[6+2];
		D[6+0] = AQ[0][2] * Q[0+0] + AQ[1][2] * Q[3+0] + AQ[2][2] * Q[6+0];
		D[6+1] = AQ[0][2] * Q[0+1] + AQ[1][2] * Q[3+1] + AQ[2][2] * Q[6+1];
		D[6+2] = AQ[0][2] * Q[0+2] + AQ[1][2] * Q[3+2] + AQ[2][2] * Q[6+2];
		o[0] = D[3+2];
		o[1] = D[0+2];
		o[2] = D[0+1];
		m[0] = fabs(o[0]);
		m[1] = fabs(o[1]);
		m[2] = fabs(o[2]);

		k0 = (m[0] > m[1] && m[0] > m[2]) ? 0 : (m[1] > m[2]) ? 1 : 2; // index of largest element of offdiag
		k1 = (k0 + 1) % 3;
		k2 = (k0 + 2) % 3;
		if (o[k0] == 0.0)
		{
			break;  // diagonal already
		}
		thet = (D[k2 * 3 + k2] - D[k1 * 3 + k1]) / (2.0*o[k0]);
		sgn = (thet > 0.0) ? 1.0 : -1.0;
		thet *= sgn; // make it positive
		t = sgn / (thet + ((thet < 1.E6) ? sqrt(thet*thet + 1.0) : thet)); // sign(T)/(|T|+sqrt(T^2+1))
		c = 1.0 / sqrt(t*t + 1.0); //  c= 1/(t^2+1) , t=s/c 
		if (c == 1.0)
		{
			break;  // no room for improvement - reached machine precision.
		}
		jr[0] = jr[1] = jr[2] = jr[3] = 0.0;
		jr[k0] = sgn*sqrt((1.0 - c) / 2.0);  // using 1/2 angle identity sin(a/2) = sqrt((1-cos(a))/2)  
		jr[k0] *= -1.0; // since our quat-to-matrix convention was for v*M instead of M*v
		jr[3] = sqrt(1.0f - jr[k0] * jr[k0]);
		if (jr[3] == 1.0)
		{
			break; // reached limits of floating point precision
		}
		q[0] = (q[3] * jr[0] + q[0] * jr[3] + q[1] * jr[2] - q[2] * jr[1]);
		q[1] = (q[3] * jr[1] - q[0] * jr[2] + q[1] * jr[3] + q[2] * jr[0]);
		q[2] = (q[3] * jr[2] + q[0] * jr[1] - q[1] * jr[0] + q[2] * jr[3]);
		q[3] = (q[3] * jr[3] - q[0] * jr[0] - q[1] * jr[1] - q[2] * jr[2]);
		mq = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
		q[0] /= mq;
		q[1] /= mq;
		q[2] /= mq;
		q[3] /= mq;
	}
	printf("done");
}

int main()
{
	cout << std::setprecision(2);
	cout << "Starting Rotation Tests" << endl << "---------------------" << endl;

	Vector3 axis = Vector3(0.0, 0.0, 1);
	axis.Normalise();
	orig_rot = Quaternion::AxisAngleToQuaterion(axis, 90.0f);
	orig_rot.Normalise();

	RunUnitTest("Quaternion Rotation", &ut_QuaternionRotation);

	//Rotate Triangle Directly (we will reverse engineer the same rotation for validitity checking)

	//Assume equal weighting for triangle points, thus the centre is just an avg of the positions/mass
	Vector3 c0 = Vector3(0.0f, 0.0f, 0.0f);
	Vector3 c = Vector3(0.0f, 0.0f, 0.0f);

	for (int i = 0; i < NUM_POINTS; ++i)
	{
		p[i] = Quaternion::RotatePointByQuaternion(orig_rot, p0[i]);
		c0 += p0[i] * NUM_POINTS_INV;
		c += p[i] * NUM_POINTS_INV;
	}

	for (int i = 0; i < NUM_POINTS; ++i)
	{
		d[i] = Quaternion::RotatePointByQuaternion(orig_rot, p[i] - c) - p0[i] - c0;
	}



	cout << endl << "The original Quaternion was: " << endl;
	cout << "\t" << orig_rot.x << ", " << orig_rot.y << ", " << orig_rot.z << ", " << orig_rot.w << endl;

	new_rot = Quaternion(0, 0, 1.0f, 0.001);
	new_rot.Normalise();

	cout << endl << "The New Quaternion was: " << endl;
	cout << "\t" << new_rot.x << ", " << new_rot.y << ", " << new_rot.z << ", " << new_rot.w << endl;

	RunUnitTest("Quaternion Comparison", &ut_QuaternionRotationComparison);

	//RunUnitTest("Newtonian Solver", &ut_NewtonianSolver);

	Matrix3 Apq = Matrix3::ZeroMatrix;
	for (int i = 0; i < NUM_POINTS; ++i)
	{
		Matrix3 current = Matrix3::OuterProduct(p[i] - c, p0[i] - c0) * NUM_POINTS_INV;
		Apq += current;
	}

	Matrix3 S2 = Matrix3::Transpose(Apq) * Apq;

	Matrix3 Q, D;
	Diagonalize(S2, Q, D);
	
	Vector3 sqrt_rows = Vector3(sqrt(D._11), sqrt(D._22), sqrt(D._33));
	
	Matrix3 sD(sqrtf(D._11), 0.0f, 0.0f,
		0.0f, sqrtf(D._22), 0.0f,
		0.0f, 0.0f, sqrtf(D._33));

	Matrix3 S = Matrix3::Transpose(Q) * sD * Q;

	Matrix3 rot = Apq * Matrix3::Inverse(S);
	Quaternion q_new = Quaternion::FromMatrix(rot);


	

	FindLeastRotation();

	system("Pause");
	return 0;
}