#pragma once

#include <nclgl\Vector2.h>
#include <nclgl\Vector3.h>
#include <nclgl\Matrix4.h>
#include <vector>


enum VertexFlags : unsigned int
{
	VERTEXFLAGS_NONE = 0x0,
	VERTEXFLAGS_IS_STATIC = 0x1,
	VERTEXFLAGS_IS_CORE = 0x2
};

struct Vertex
{
	Vector3 position;
	Vector2 texCoord;
	Vector3 force;
	unsigned int flags;
};

struct Triangle
{
	Triangle() : vert_a(0), vert_b(0), vert_c(0) {}
	Triangle(unsigned int a, unsigned int b, unsigned int c) : vert_a(a), vert_b(b), vert_c(c) {}

	unsigned int vert_a, vert_b, vert_c;
};

class ClothDesignEntity
{
public:
	ClothDesignEntity();
	~ClothDesignEntity();

	inline const std::vector<Vertex>& Vertices() const { return m_Vertices; }
	inline const std::vector<Triangle>& Triangles() const { return m_Triangles; }
	
	inline std::vector<Vertex>& Vertices() { return m_Vertices; }
	inline std::vector<Triangle>& Triangles() { return m_Triangles; }

	void		GenerateDefaultClothVertices(unsigned int divisor, float width, float height, Matrix4 transform = Matrix4());
	void		GenerateDefaultClothTriangles(unsigned int divisor);

protected:
	inline void AddTriangle(const Triangle& t) { m_Triangles.push_back(t); }

protected:
	std::vector<Vertex>			m_Vertices;
	std::vector<Triangle>		m_Triangles;
};