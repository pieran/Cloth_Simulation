
#ifndef MCG_H
#define MCG_H


#include <Eigen\Dense.h>
#include <vector>
#include <map>

using Eigen::Vector3f;
using Eigen::Matrix3f;

typedef std::map<int, Matrix3f> matrix_map;
typedef matrix_map::iterator matrix_iterator;

//Modified Conjugate Gradient
class mcg
{
public:


protected:
	std::vector<matrix_map> m_AMtx_Rows;
	std::vector<Vector3f>	m_PhyxelsResidual;
	std::vector<Vector3f>	m_PhyxelsPrevious;
	std::vector<Vector3f>	m_PhyxelsUpdate;
	std::vector<Vector3f>	m_PhyxelsB;
};

#endif //MCG_H