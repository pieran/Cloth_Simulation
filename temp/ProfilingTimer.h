#pragma once

#include "Windows.h"

class ProfilingTimer
{
public:
	ProfilingTimer();
	virtual ~ProfilingTimer();

	void BeginTiming();
	void EndTiming();
	void EndTimingAdditive();

	inline void ResetTotalMs() { m_TotalMs = 0.0f; } //For Additive Profiling
	inline float GetTimedMilliSeconds() { return m_TotalMs; }

protected:
	float m_StartTime;
	float m_TotalMs;

	LARGE_INTEGER	m_CPUStart;		//Start of timer
	LARGE_INTEGER	m_CPUFrequency;	//Ticks Per Second
};