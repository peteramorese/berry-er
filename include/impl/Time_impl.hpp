#include "Time.h"

#include "berry/Logging.h"

namespace _BRY {

double convert(const std::chrono::microseconds& duration, BRY::TimeUnit unit) {
    switch (unit) {
        case BRY::TimeUnit::us:
            return static_cast<double>(duration.count());
        case BRY::TimeUnit::ms:
            return static_cast<double>(duration.count()) / 1000;
        case BRY::TimeUnit::s:
            return static_cast<double>(duration.count()) / 1000000;
    }
    ERROR("Unrecognized time unit");
    return 0.0;
}

}

double BRY::Profiler::getMostRecentProfile(const std::string& key, TimeUnit unit) {
    auto it = s_profiles.find(key);
    if (it == s_profiles.end()) {
        ERROR("Profile key '" << key << "' not found, has a Timer object been created with that key?");
        return 0.0;
    }
    return _BRY::convert(it->second.most_recent_duration, unit);
}

double BRY::Profiler::getTotalProfile(const std::string& key, TimeUnit unit) {
    auto it = s_profiles.find(key);
    if (it == s_profiles.end()) {
        ERROR("Profile key '" << key << "' not found, has a Timer object been created with that key?");
        return 0.0;
    }
    return _BRY::convert(it->second.total_duration, unit);
}

BRY::Profiler::Data& BRY::Profiler::getData(const std::string& key) {
    auto it = s_profiles.find(key);
    if (it == s_profiles.end()) {
        BRY::Profiler::Data data{std::chrono::microseconds(0), std::chrono::microseconds(0)};
        it = s_profiles.insert(std::make_pair(key, data)).first;
    }
    return it->second;
}

BRY::Timer::Timer(const std::string& key) 
    : m_key(key)
    , m_start(std::chrono::system_clock::now()) {}

void BRY::Timer::stop() {
    BRY::Profiler::Data& data = Profiler::getData(m_key);
    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    std::chrono::microseconds duration = std::chrono::duration_cast<std::chrono::microseconds>(end - m_start);
    data.most_recent_duration = duration;
    data.total_duration += duration;
    m_key.clear();
}

double BRY::Timer::now(TimeUnit unit) {
    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    std::chrono::microseconds duration = std::chrono::duration_cast<std::chrono::microseconds>(end - m_start);
    return _BRY::convert(duration, unit);
}

BRY::Timer::~Timer() {
    if (!m_key.empty()) {
        stop();
    }
}