#pragma once

//Custom Array Class to abstract and allow memory optimisation at a later date
#include "SimulationDefines.h"
#include <iostream>
#include <assert.h>

template <class T>
class PArray
{
public:
	PArray()
	{
		m_len = 0;
		m_arr = NULL;
	}

	~PArray()
	{
		if (m_len > 0)
		{
			delete[] m_arr;
			m_len = 0;
		}
	}

	void resize(uint size) {
		if (size != m_len)
		{
			T* tmp = m_arr;
			m_arr = new T[size];
			if (tmp != NULL && m_len > 0)
			{
				uint cpy_len = (m_len < size) ? m_len : size;
				memcpy(m_arr, tmp, cpy_len * sizeof(T));
				delete[] tmp;
			}
			m_len = size;
		}
	}

	inline const uint size() const { return m_len; }

	inline T& operator[](uint idx) {
#if _DEBUG
		assert(idx >= 0 && idx < m_len);
#endif	
		return m_arr[idx];
	}

	inline const T& operator[](uint idx) const {
#if _DEBUG
		assert(idx >= 0 && idx < m_len);
#endif	
		return m_arr[idx];
	}

protected:
	uint m_len;
	T*   m_arr;
};