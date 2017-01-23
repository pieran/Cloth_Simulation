#pragma once

#define PROFILING_ENABLED TRUE

enum ProfilerID
{
	PROFILERID_CLOTH_TOTALTIME = 0,
	PROFILERID_CLOTH_COMPUTEFORCES,
	PROFILERID_CLOTH_SOLVER,
	PROFILERID_CLOTH_BUILDROTATIONANDNORMALS,
	PROFILERID_CLOTH_EXTERNALCOLLISIONS,
	PROFILERID_SOLVER_INITIALIZATION,
	PROFILERID_SOLVER_UPPERAMATRIX,
	PROFILERID_SOLVER_LOWERAMATRIX,
	PROFILERID_SOLVER_AVERAGEITERATIONS,

	PROFILERID_MAX
};

#if PROFILING_ENABLED == TRUE

#include <Windows.h>

class SimulationProfiler
{
public:
	SimulationProfiler();
	~SimulationProfiler();

	void ResetAllTiming();
	void ResetTiming(ProfilerID timer_idx);
	void BeginTiming(ProfilerID timer_idx);
	void EndTiming(ProfilerID timer_idx);
	void EndTimingAccumulative(ProfilerID timer_idx);


	double GetTimingMS(ProfilerID timer_idx);

protected:
	double m_StartTime[PROFILERID_MAX];
	double m_TotalMS[PROFILERID_MAX];

	LARGE_INTEGER	m_CPUStart;		//Start of timer
	LARGE_INTEGER	m_CPUFrequency;	//Ticks Per Second
};
#else
class SimulationProfiler
{
public:
	SimulationProfiler() {};
	~SimulationProfiler() {};

	void ResetAllTiming() {};
	void ResetTiming(ProfilerID timer_idx) {};
	void BeginTiming(ProfilerID timer_idx) {};
	void EndTiming(ProfilerID timer_idx) {};
	void EndTimingAccumulative(ProfilerID timer_idx) {};


	double GetTimingMS(ProfilerID timer_idx) { return (double)NAN; }
};
#endif