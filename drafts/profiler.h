#pragma once

#include <algorithm>
#include <chrono>
#include <fstream>
#include <mutex>
#include <string>
#include <thread>

namespace Profiler {

#define PROFILING 1
#ifdef PROFILING
#include <string>

#define PROFILE_SESSION(sessionName)                                           \
  Profiler::Instrumentor::Instance().beginSession(sessionName)

#define PROFILE_SESSION_SAVE(sessionName, file)                                \
  Profiler::Instrumentor::Instance().beginSession(sessionName, file)

#define PROFILE_END_SESSION(sessionName)                                       \
  Profiler::Instrumentor::Instance().endSession()

#if defined(__GNUC__) || defined(__GNUG__)
#define FUNCTION_SIGNATURE __PRETTY_FUNCTION__ // GCC only
#else
#define FUNCTION_SIGNATURE __FUNCSIG__
#endif

// #   if !_HAS_CXX17
//         std::string eraseSubStr(std::string&& mainStr, const std::string&
//         toErase)
// #   else
#include <string_view>
std::string eraseSubStr(std::string &&mainStr, const std::string_view &toErase)
//#   endif // _HAS_CXX17
{
  // Search for the substring in string
  size_t pos = mainStr.find(toErase);
  if (pos != std::string::npos) {
    // If found then erase it from string
    mainStr.erase(pos, toErase.length());
  }
  return mainStr;
}
#define CLEAN_FUNCTION_SIGNATURE()                                             \
  Profiler::eraseSubStr(FUNCTION_SIGNATURE, "__cdecl ")

#define PROFILE()                                                              \
  Profiler::InstrumentationTimer timer##__LINE__(CLEAN_FUNCTION_SIGNATURE())

#define PROFILE_SCOPE_BEGIN()                                                  \
  {                                                                            \
    PROFILE()

#define PROFILE_SCOPE_END() }
#else
#define PROFILE_SCOPE_BEGIN()
#define PROFILE_SCOPE_END()
#define PROFILE_END_SESSION(sessionName)
#define PROFILE_SESSION_SAVE(sessionName, file)
#define PROFILE_SESSION(name)
#define PROFILE()
#endif

struct ProfileResult {
  const std::string name;
  long long start, end;
  uint32_t threadID;
};

class Instrumentor {
  std::string m_sessionName = "None";
  std::ofstream m_outputStream;
  int m_profileCount = 0;
  std::mutex m_lock;
  bool m_activeSession = false;

  Instrumentor() {}

public:
  static Instrumentor &Instance() {
    static Instrumentor instance;
    return instance;
  }

  ~Instrumentor() { endSession(); }

  void beginSession(const std::string &name, const std::string &filepath) {
    if (m_activeSession)
      endSession();
    m_activeSession = true;
    m_outputStream.open(filepath);
    writeHeader();
    m_sessionName = name;
  }

  void endSession() {
    if (!m_activeSession)
      return;
    m_activeSession = false;
    writeFooter();
    m_outputStream.close();
    m_profileCount = 0;
  }

  void writeProfile(const ProfileResult &result) {
    std::lock_guard<std::mutex> lock(m_lock);

    if (m_profileCount++ > 0)
      m_outputStream << ",";

    std::string name = result.name;
    std::replace(name.begin(), name.end(), '"', '\'');

    m_outputStream << "{";
    m_outputStream << "\"cat\":\"function\",";
    m_outputStream << "\"dur\":" << (result.end - result.start) << ',';
    m_outputStream << "\"name\":\"" << name << "\",";
    m_outputStream << "\"ph\":\"X\",";
    m_outputStream << "\"pid\":0,";
    m_outputStream << "\"tid\":" << result.threadID << ",";
    m_outputStream << "\"ts\":" << result.start;
    m_outputStream << "}";
  }

  void writeHeader() {
    m_outputStream << "{\"otherData\": {},\"traceEvents\":[";
  }

  void writeFooter() { m_outputStream << "]}"; }
};

void beginSession(const std::string &sessionName,
                  const std::string &filepath = "results.json") {
  Instrumentor::Instance().beginSession(sessionName, filepath);
}

void endSession() { Instrumentor::Instance().endSession(); }

class InstrumentationTimer {
  ProfileResult m_result;

  std::chrono::time_point<std::chrono::high_resolution_clock> m_startTimepoint;
  bool m_stopped;

public:
  InstrumentationTimer(const std::string &name)
      : m_result({name, 0, 0, 0}), m_stopped(false) {
    m_startTimepoint = std::chrono::high_resolution_clock::now();
  }

  ~InstrumentationTimer() {
    if (!m_stopped)
      stop();
  }

  void stop() {
    auto endTimepoint = std::chrono::high_resolution_clock::now();

    m_result.start = std::chrono::time_point_cast<std::chrono::microseconds>(
                         m_startTimepoint)
                         .time_since_epoch()
                         .count();
    m_result.end =
        std::chrono::time_point_cast<std::chrono::microseconds>(endTimepoint)
            .time_since_epoch()
            .count();
    m_result.threadID =
        std::hash<std::thread::id>{}(std::this_thread::get_id());
    Instrumentor::Instance().writeProfile(m_result);

    m_stopped = true;
  }
};
} // namespace Profiler
//
// Basic instrumentation profiler by Cherno
// Source: Yan Chernikov TheCherno
// https://gist.github.com/TheCherno/31f135eea6ee729ab5f26a6908eb3a5e Modified
// by: davechurchill https://pastebin.com/qw5Neq4U Modified again by: Myself
// (Callum Poole)
//
// URL: chrome://tracing/
