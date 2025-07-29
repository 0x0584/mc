#ifndef PROFILER_HPP
#define PROFILER_HPP
#ifndef NDEBUG
#include <gperftools/profiler.h>
#define profiler_start(prof) ProfilerStart(prof)
#define profiler_stop() ProfilerStop()
#else
#define profiler_start(prof)
#define profiler_stop()
#endif
#endif
