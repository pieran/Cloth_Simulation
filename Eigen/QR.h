#ifndef EIGEN_QR_MODULE_H
#define EIGEN_QR_MODULE_H

#include "Core.h"

#include "src/Core/util/DisableStupidWarnings.h"

#include "Cholesky.h"
#include "Jacobi.h"
#include "Householder.h"

/** \defgroup QR_Module QR module
  *
  *
  *
  * This module provides various QR decompositions
  * This module also provides some MatrixBase methods, including:
  *  - MatrixBase::householderQr()
  *  - MatrixBase::colPivHouseholderQr()
  *  - MatrixBase::fullPivHouseholderQr()
  *
  * \code
  * #include <Eigen/QR>
  * \endcode
  */

#include "src/QR/HouseholderQR.h"
#include "src/QR/FullPivHouseholderQR.h"
#include "src/QR/ColPivHouseholderQR.h"
#ifdef EIGEN_USE_LAPACKE
#include "src/QR/HouseholderQR_MKL.h"
#include "src/QR/ColPivHouseholderQR_MKL.h"
#endif

#include "src/Core/util/ReenableStupidWarnings.h"

#endif // EIGEN_QR_MODULE_H
/* vim: set filetype=cpp et sw=2 ts=2 ai: */
