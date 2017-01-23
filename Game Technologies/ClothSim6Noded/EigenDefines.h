
#pragma once

//Allows simple switching between float and double precision in one place

#include <Eigen\Dense.h>
#include <Eigen\Sparse.h>

typedef Eigen::Matrix2f Mat22;
typedef Eigen::Matrix3f Mat33;
typedef Eigen::Matrix4f Mat44;

typedef Eigen::VectorXf VecX;
typedef Eigen::Vector2f Vec2;
typedef Eigen::Vector3f Vec3;
typedef Eigen::Vector4f Vec4;

typedef Eigen::SparseMatrix<float> SparseMatrix;
typedef Eigen::SparseMatrix<Mat33> SparseMatrixMat33;

#ifndef uint
typedef unsigned int uint;
#endif