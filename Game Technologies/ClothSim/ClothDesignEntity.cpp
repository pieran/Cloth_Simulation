#include "ClothDesignEntity.h"

ClothDesignEntity::ClothDesignEntity() : m_Vertices(0), m_Triangles(0)
{
}

ClothDesignEntity::~ClothDesignEntity()
{
	m_Vertices.clear();
	m_Triangles.clear();
}

void ClothDesignEntity::GenerateDefaultClothVertices(unsigned int divisor, float width, float height, Matrix4 transform)
{
	unsigned int num_total = divisor * divisor;
	m_Vertices.resize(num_total);

	const float divisor_scalar = 1.f / (divisor - 1);

	unsigned int count = 0;
	for (unsigned int y = 0; y < divisor; ++y)
	{
		for (unsigned int x = 0; x < divisor; ++x)
		{
			m_Vertices[count].position = transform * Vector3(
				((float(x) * divisor_scalar) - 0.5f) *  width,
				((float(y) * divisor_scalar) - 0.5f)  * -height,
				0.0f
				);

			m_Vertices[count].texCoord = Vector2(
				(float(x) * divisor_scalar),
				(float(y) * divisor_scalar)
				);

			m_Vertices[count].force = Vector3(0, 0, 0);

			if (false)//y == 0)
				m_Vertices[count].flags = VERTEXFLAGS_IS_STATIC;
			//else if (y == divisor - 1)
			//	m_Vertices[count].flags = VERTEXFLAGS_IS_STATIC | VERTEXFLAGS_IS_CORE;
			else
				m_Vertices[count].flags = VERTEXFLAGS_NONE;
			count++;
		}
	}
}

void ClothDesignEntity::GenerateDefaultClothTriangles(unsigned int divisor)
{
	unsigned int num_tris = (divisor - 1) * (divisor - 1) * 2;
	m_Triangles.resize(num_tris);

	unsigned int count = 0;

	for (unsigned int x = 1; x < divisor; ++x)
	{
		for (unsigned int y = 1; y < divisor; ++y)
		{
			unsigned int a = divisor*(y - 1) + x - 1;
			unsigned int b = divisor*(y - 1) + x;
			unsigned int c = divisor*(y)+x;
			unsigned int d = divisor*(y)+x - 1;
			//Alternate triangle orientation to allow FEM to better represent the forces
			/*if (((y ^ x) & 0x1) == 0)
			{
			m_Triangles[count++] = Triangle(a, c, b);
			m_Triangles[count++] = Triangle(a, d, c);
			}
			else
			{
			m_Triangles[count++] = Triangle(a, d, b);
			m_Triangles[count++] = Triangle(b, d, c);
			}*/

			//Alternate triangle orientation to allow FEM to better represent the forces
			if (((y ^ x) & 0x1) == 0)
			{
				m_Triangles[count++] = Triangle(c, a, b);
				m_Triangles[count++] = Triangle(d, a, c);
			}
			else
			{
				m_Triangles[count++] = Triangle(d, a, b);
				m_Triangles[count++] = Triangle(d, b, c);
			}
		}
	}
}