#pragma once


#include "EigenDefines.h"
#include "FECommon.h"
#include "common.h"
#include <vector>
#include <stdarg.h>


enum VertexFlags : unsigned int
{
	VERTEXFLAGS_NONE = 0x0,
	VERTEXFLAGS_IS_STATIC = 0x1,
	VERTEXFLAGS_IS_CORE = 0x2
};

struct Vertex
{
	Vec3 initial_position;
	Vec3 position;
	Vec2 texCoord;
	unsigned int flags;
};

template<uint nverts>
struct TriangleG
{
	TriangleG()
	{
		memset(verts, 0, nverts * sizeof(uint));
	}

	TriangleG(int nset, ...)
	{
		va_list vl;
		va_start(vl, nset);
		for (int i = 0; i < nset; ++i)
		{
			verts[i] = va_arg(vl, uint);
		}
		va_end(vl);

		int leftover = nverts - nset;
		if (leftover > 0) memset(&verts[nset], 0, leftover * sizeof(uint));
	}

	uint verts[nverts];
};

typedef TriangleG<3> Triangle;



struct Quad	//For Generating Smooth Normals
{
	Quad() : vert_a(0), vert_b(0), vert_c(0), vert_d(0) {}
	Quad(unsigned int a, unsigned int b, unsigned int c, unsigned int d) : vert_a(a), vert_b(b), vert_c(c), vert_d(d) {}

	unsigned int vert_a, vert_b, vert_c, vert_d;
};

template <uint nVerts> //Number of verts per triangle (currently supported values are 3/6)
class ClothDesignEntity
{
public:
	ClothDesignEntity();
	~ClothDesignEntity();

	inline const std::vector<Vertex>& Vertices() const { return m_Vertices; }
	inline const std::vector<TriangleG<nVerts>>& Triangles() const { return m_Triangles; }
	inline const std::vector<Quad>& Quads() const { return m_Quads; }

	inline std::vector<Vertex>& Vertices() { return m_Vertices; }
	inline std::vector<TriangleG<nVerts>>& Triangles() { return m_Triangles; }
	inline std::vector<Quad>& Quads() { return m_Quads; }

	void		GenerateDefaultClothVertices(unsigned int divisor, float width, float height, Mat44* transform);
	void		GenerateDefaultClothTriangles(unsigned int divisor);

	inline int GetWidth() const { return m_NumWidth; }
	inline int GetHeight() const { return m_NumHeight; }
protected:
	int m_NumWidth,	m_NumHeight;
	std::vector<Vertex>					m_Vertices;
	std::vector<TriangleG<nVerts>>		m_Triangles;
	std::vector<Quad>					m_Quads;
};

#include "ClothDesignEntity.inl"