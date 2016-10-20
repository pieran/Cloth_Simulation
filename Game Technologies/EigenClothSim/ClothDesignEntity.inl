#pragma once

#include "ClothDesignEntity.h"

template <uint nVerts>
ClothDesignEntity<nVerts>::ClothDesignEntity() : m_Vertices(0), m_Triangles(0)
{
}

template <uint nVerts>
ClothDesignEntity<nVerts>::~ClothDesignEntity()
{
	m_Vertices.clear();
	m_Triangles.clear();
}

template <uint nVerts>
void ClothDesignEntity<nVerts>::GenerateDefaultClothVertices(unsigned int divisor, float width, float height, Mat44 *transform)
{
	m_NumWidth = divisor;
	m_NumHeight = divisor;

	unsigned int num_total = divisor * divisor;
	m_Vertices.resize(num_total);

	const float divisor_scalar = 1.f / (divisor - 1);

	unsigned int count = 0;
	for (unsigned int y = 0; y < divisor; ++y)
	{
		for (unsigned int x = 0; x < divisor; ++x)
		{
			float z = float(rand() % 2000) / 1000.0f - 1.0f;
			Vec4 pos = *transform * Vec4(
				((float(x) * divisor_scalar) - 0.5f) *  width,
				((float(y) * divisor_scalar) - 0.5f)  * -height,
				0.0f,
				1.0f
				);

			m_Vertices[count].position = Vec3(pos[0], pos[1], pos[2]) / pos[3];

			m_Vertices[count].initial_position = Vec3(
				((float(x) * divisor_scalar) - 0.5f) *  width,
				((float(y) * divisor_scalar) - 0.5f)  * -height,
				0.0f);


			m_Vertices[count].texCoord = Vec2(
				(float(x) * divisor_scalar),
				(float(y) * divisor_scalar)
				);

			m_Vertices[count].flags = VERTEXFLAGS_NONE;
			count++;
		}
	}
}

template <uint nVerts>
void ClothDesignEntity<nVerts>::GenerateDefaultClothTriangles(unsigned int divisor)
{
	if (nVerts == 3)
	{
		unsigned int num_quads = (divisor - 1) * (divisor - 1);
		unsigned int num_tris = num_quads * 2;
		m_Triangles.resize(num_tris);
		m_Quads.resize(num_quads);


		unsigned int tcount = 0;
		unsigned int qcount = 0;
		for (unsigned int x = 1; x < divisor; ++x)
		{
			for (unsigned int y = 1; y < divisor; ++y)
			{
				unsigned int a = divisor*(y - 1) + x - 1;
				unsigned int b = divisor*(y - 1) + x;
				unsigned int c = divisor*(y)+x;
				unsigned int d = divisor*(y)+x - 1;

				m_Quads[qcount++] = Quad(a, b, c, d);

				//Alternate triangle orientation to allow FEM to better represent the forces
				if (((y ^ x) & 0x1) == 0)
				{
					m_Triangles[tcount++] = TriangleG<nVerts>(3, a, b, c );
					m_Triangles[tcount++] = TriangleG<nVerts>(3, a, c, d );

				//	m_Triangles[tcount++] = TriangleG<nVerts>(3, b, a, c);
				//	m_Triangles[tcount++] = TriangleG<nVerts>(3, c, a, d);
				}
				else
				{
					m_Triangles[tcount++] = TriangleG<nVerts>(3, b, d, a );
					m_Triangles[tcount++] = TriangleG<nVerts>(3, b, c, d );
				}
			}
		}
	}
	else if (nVerts == 6)
	{
		uint num_quads = (divisor - 1) * (divisor - 1) / 4;
		uint num_tris = num_quads * 2;
		m_Triangles.resize(num_tris);

		uint tcount = 0;
		uint qcount = 0;
		for (uint x = 2; x < divisor; x+=2)
		{
			for (uint y = 2; y < divisor; y+=2)
			{
				uint a = divisor*(y - 2) + x - 2;
				uint b = divisor*(y - 2) + x;
				uint c = divisor*(y)+x;
				uint d = divisor*(y)+x - 2;

				uint ab = divisor*(y - 2) + x - 1;
				uint dc = divisor*(y)+x - 1;
				uint ad = divisor*(y - 1) + x - 2;
				uint bc = divisor*(y - 1) + x;
				
				uint mid = divisor*(y - 1) + x - 1;



				//Alternate triangle orientation to allow FEM to better represent the forces
				if (true)// ((y ^ x) & 0x2) == 0)
				{
					m_Triangles[tcount++] = TriangleG<nVerts>(6, a, d, c, ad, dc, mid);
					m_Triangles[tcount++] = TriangleG<nVerts>(6, a, c, b, mid, bc, ab);
				}
				else
				{
					m_Triangles[tcount++] = TriangleG<nVerts>(6, d, b, a, mid, ab, ad);
					m_Triangles[tcount++] = TriangleG<nVerts>(6, d, c, b, dc, bc, mid);
				}
			}
		}
	}
	else
	{
		assert(false);
	}
}