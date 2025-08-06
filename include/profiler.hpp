#ifndef PROFILER_HPP
#define PROFILER_HPP
#ifdef ENABLE_PROFILING
#include <gperftools/profiler.h>
#define profiler_start(prof) ProfilerStart(prof)
#define profiler_stop() ProfilerStop()
#else
#define profiler_start(prof)                                                   \
  do {                                                                         \
  } while (0)
#define profiler_stop()                                                        \
  do {                                                                         \
  } while (0)
#endif
#endif
