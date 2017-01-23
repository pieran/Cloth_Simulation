#include "ProfilingTimer.h"


ProfilingTimer::ProfilingTimer()
{
	QueryPerformanceFrequency((LARGE_INTEGER *)&m_CPUFrequency);
	QueryPerformanceCounter((LARGE_INTEGER *)&m_CPUStart);

	m_StartTime = 0.0f;
	m_TotalMs = 0.0f;
}

ProfilingTimer::~ProfilingTimer()
{

}

void ProfilingTimer::BeginTiming()
{
	LARGE_INTEGER t;
	QueryPerformanceCounter(&t);
	m_StartTime = (float)((t.QuadPart - m_CPUStart.QuadPart) * 1000.0 / m_CPUFrequency.QuadPart);
}

void ProfilingTimer::EndTiming()
{
	LARGE_INTEGER t;
	QueryPerformanceCounter(&t);
	float end_time = (float)((t.QuadPart - m_CPUStart.QuadPart) * 1000.0 / m_CPUFrequency.QuadPart);

	m_TotalMs = end_time - m_StartTime;
}

void ProfilingTimer::EndTimingAdditive()
{
	LARGE_INTEGER t;
	QueryPerformanceCounter(&t);
	float end_time = (float)((t.QuadPart - m_CPUStart.QuadPart) * 1000.0 / m_CPUFrequency.QuadPart);

	m_TotalMs += end_time - m_StartTime;
}