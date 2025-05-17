// Timer.h
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
        if (timestamps.count(label)) {
            throw std::runtime_error("Timer '" + label + "' already exists");
        }
        timestamps[label] = Clock::now();
    }

    void stop(const std::string& label, const bool& add_to_total = true) {
        auto end = Clock::now();
        if (timestamps.count(label)) {
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - timestamps[label]).count();
            std::cout << label << " cost: " << duration << " ms" << std::endl;
            if(add_to_total){
                totalTime += duration;
            }
            timestamps.erase(label);
        } else {
            std::cerr << "label not found: " << label << std::endl;
        }
    }

    void printAll() {
        std::cout << "total cost: " << totalTime << " ms" << std::endl;
        totalTime = 0;
    }


private:
    using Clock = std::chrono::high_resolution_clock;
    std::unordered_map<std::string, Clock::time_point> timestamps;

    Timer() = default;
    Timer(const Timer&) = delete;
    Timer& operator=(const Timer&) = delete;
    uint64_t totalTime = 0;
};

inline void set_timer(const std::string& label) {
    Timer::getInstance().start(label);
}

inline void end_timer(const std::string& label, bool add_to_total = true) {
    Timer::getInstance().stop(label, add_to_total);
}

inline void totaltime() {
    Timer::getInstance().printAll();
}
