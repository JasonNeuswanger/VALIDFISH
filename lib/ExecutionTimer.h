//
// Taken from https://stackoverflow.com/questions/22387586/measuring-execution-time-of-a-function-in-c
// Heavily modified

#ifndef DRIFTMODELC_EXECUTIONTIMER_H
#define DRIFTMODELC_EXECUTIONTIMER_H

#include <chrono>
#include <iostream>
#include <sstream>

template<class Resolution = std::chrono::milliseconds>
class ExecutionTimer {
public:
    using Clock = std::conditional_t<std::chrono::high_resolution_clock::is_steady,
            std::chrono::high_resolution_clock,
            std::chrono::steady_clock>;
private:
    const Clock::time_point mStart = Clock::now();
    std::string message;

public:
    explicit ExecutionTimer(std::string _message) {
        message = _message;
    }

    inline void stop() {
        const auto end = Clock::now();
        std::ostringstream strStream;
        strStream << "Stop Elapsed: "
                  << std::chrono::duration_cast<Resolution>(end - mStart).count()
                  << " (" << message << ")";
        std::cout << strStream.str() << std::endl;
    }

}; // ExecutionTimer

#endif //DRIFTMODELC_EXECUTIONTIMER_H
