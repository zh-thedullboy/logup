// // Timer.h
// #pragma once
// #include <iostream>
// #include <chrono>
// #include <unordered_map>
// #include <string>
// #include <stdexcept>

// class Timer {
// public:
//     static Timer& getInstance() {
//         static Timer instance;
//         return instance;
//     }

//     void start(const std::string& label) {
//         if (timestamps.count(label)) {
//             throw std::runtime_error("Timer '" + label + "' already exists");
//         }
//         timestamps[label] = Clock::now();
//     }

//     void stop(const std::string& label, const bool& add_to_total = true) {
//         auto end = Clock::now();
//         if (timestamps.count(label)) {
//             auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - timestamps[label]).count();
//             std::cout << label << " cost: " << duration << " ms" << std::endl;
//             if(add_to_total){
//                 totalTime += duration;
//             }
//             timestamps.erase(label);
//         } else {
//             std::cerr << "label not found: " << label << std::endl;
//         }
//     }

//     void printAll() {
//         std::cout << "total cost: " << totalTime << " ms" << std::endl;
//         totalTime = 0;
//     }


// private:
//     using Clock = std::chrono::high_resolution_clock;
//     std::unordered_map<std::string, Clock::time_point> timestamps;

//     Timer() = default;
//     Timer(const Timer&) = delete;
//     Timer& operator=(const Timer&) = delete;
//     uint64_t totalTime = 0;
// };

#pragma once
#include <iostream>
#include <chrono>
#include <unordered_map>
#include <string>
#include <stdexcept>

class Timer {
public:
    static Timer& getInstance() {
        static Timer instance;
        return instance;
    }

    void start(const std::string& label) {
        if (timers.count(label)) {
            throw std::runtime_error("Timer '" + label + "' already exists");
        }
        timers[label] = TimerData{Clock::now(), 0, false};
    }

    void pause(const std::string& label) {
        auto it = timers.find(label);
        if (it == timers.end()) {
            std::cerr << "pause failed, label not found: " << label << std::endl;
            return;
        }
        if (it->second.paused) {
            std::cerr << "Timer '" << label << "' is already paused" << std::endl;
            return;
        }
        auto now = Clock::now();
        it->second.accumulated += std::chrono::duration_cast<std::chrono::milliseconds>(now - it->second.startTime).count();
        it->second.paused = true;
    }

    void resume(const std::string& label) {
        auto it = timers.find(label);
        if (it == timers.end()) {
            std::cerr << "resume failed, label not found: " << label << std::endl;
            return;
        }
        if (!it->second.paused) {
            std::cerr << "Timer '" << label << "' is not paused" << std::endl;
            return;
        }
        it->second.startTime = Clock::now();
        it->second.paused = false;
    }

    void stop(const std::string& label, const bool& add_to_total = true) {
        auto it = timers.find(label);
        if (it == timers.end()) {
            std::cerr << "stop failed, label not found: " << label << std::endl;
            return;
        }
        auto now = Clock::now();
        uint64_t duration = it->second.accumulated;
        if (!it->second.paused) {
            duration += std::chrono::duration_cast<std::chrono::milliseconds>(now - it->second.startTime).count();
        }
        std::cout << label << " cost: " << duration << " ms" << std::endl;
        if (add_to_total) {
            totalTime += duration;
        }
        timers.erase(it);
    }

    void printAll() {
        std::cout << "total cost: " << totalTime << " ms" << std::endl;
        totalTime = 0;
    }

private:
    using Clock = std::chrono::high_resolution_clock;

    struct TimerData {
        Clock::time_point startTime;
        uint64_t accumulated; // 累积的时间(ms)
        bool paused;
    };

    std::unordered_map<std::string, TimerData> timers;

    Timer() = default;
    Timer(const Timer&) = delete;
    Timer& operator=(const Timer&) = delete;
    uint64_t totalTime = 0;
};


inline void pause_timer(const std::string& label) {
    Timer::getInstance().pause(label);
}

inline void resume_timer(const std::string& label) {
    Timer::getInstance().resume(label);
}

inline void set_timer(const std::string& label) {
    Timer::getInstance().start(label);
}

inline void end_timer(const std::string& label, bool add_to_total = true) {
    Timer::getInstance().stop(label, add_to_total);
}

inline void totaltime() {
    Timer::getInstance().printAll();
}
