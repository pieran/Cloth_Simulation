#include "SimulationProfiler.h"

#if PROFILING_ENABLED == TRUE
SimulationProfiler::SimulationProfiler()
{
	QueryPerformanceFrequency((LARGE_INTEGER *)&m_CPUFrequency);
	QueryPerformanceCounter((LARGE_INTEGER *)&m_CPUStart);

	memset(m_StartTime, 0, PROFILERID_MAX * sizeof(double));
	memset(m_TotalMS, 0, PROFILERID_MAX * sizeof(double));
}

SimulationProfiler::~SimulationProfiler()
{

}

void SimulationProfiler::ResetAllTiming()
{
	memset(m_TotalMS, 0, PROFILERID_MAX * sizeof(double));
}

void SimulationProfiler::ResetTiming(ProfilerID timer_idx)
{
	m_TotalMS[timer_idx] = 0.0;
}

void SimulationProfiler::BeginTiming(ProfilerID timer_idx)
{
	LARGE_INTEGER t;
	QueryPerformanceCounter(&t);
	m_StartTime[timer_idx] = (t.QuadPart - m_CPUStart.QuadPart) * 1000.0 / m_CPUFrequency.QuadPart;
}

void SimulationProfiler::EndTiming(ProfilerID timer_idx)
{
	LARGE_INTEGER t;
	QueryPerformanceCounter(&t);
	double end_time = (t.QuadPart - m_CPUStart.QuadPart) * 1000.0 / m_CPUFrequency.QuadPart;
	m_TotalMS[timer_idx] = end_time - m_StartTime[timer_idx];
}

void SimulationProfiler::EndTimingAccumulative(ProfilerID timer_idx)
{
	LARGE_INTEGER t;
	QueryPerformanceCounter(&t);
	double end_time = (t.QuadPart - m_CPUStart.QuadPart) * 1000.0 / m_CPUFrequency.QuadPart;
	m_TotalMS[timer_idx] += end_time - m_StartTime[timer_idx];
}

double SimulationProfiler::GetTimingMS(ProfilerID timer_idx)
{
	return m_TotalMS[timer_idx];
}
#endif