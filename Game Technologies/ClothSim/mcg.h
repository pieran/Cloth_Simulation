
#ifndef MCG_H
#define MCG_H


#include <nclgl\Vector3.h>
#include <nclgl\Matrix3.h>
#include <vector>
#include <map>

typedef std::map<int, Matrix3> matrix_map;
typedef matrix_map::iterator matrix_iterator;

//Modified Conjugate Gradient
class mcg
{
public:


protected:
	std::vector<matrix_map> m_AMtx_Rows;
	std::vector<Vector3>	m_PhyxelsResidual;
	std::vector<Vector3>	m_PhyxelsPrevious;
	std::vector<Vector3>	m_PhyxelsUpdate;
	std::vector<Vector3>	m_PhyxelsB;
};

#endif //MCG_H