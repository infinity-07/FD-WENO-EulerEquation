#ifndef TIMER_HPP
#define TIMER_HPP

//////////////////////////////////////////////////////////////////////////
/// @file      Timer.hpp
///
/// @brief     This file defines a hierarchical timer class for profiling
///            and analyzing the runtime performance of nested code blocks.
///
/// @note      This timer supports nested timing, progress update reporting,
///            and human-readable formatting of elapsed time and ETA.
///
/// @author    Yang Zhang
/// @date      May 19, 2025
/// @version   1.0
///
/// @copyright (c) 2025 USTB.
///            All rights reserved.
///
/// @license   This code is released under the MIT License.
///
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <chrono>
#include <vector>
#include <string>
#include <iomanip>
#include <sstream>
#include <ctime>
#include <algorithm>

/// Convert seconds to a human-readable string format (e.g., "2h 15m")
inline std::string formatTime(double seconds)
{
    int days = static_cast<int>(seconds) / 86400;
    int hours = (static_cast<int>(seconds) % 86400) / 3600;
    int minutes = (static_cast<int>(seconds) % 3600) / 60;
    double secs = std::fmod(seconds, 60.0); // 保留小数

    std::stringstream ss;
    ss << std::fixed << std::setprecision(2); // 设置两位小数

    if (days > 0)
    {
        ss << days << "d";
        if (hours > 0)
            ss << " " << hours << "h";
    }
    else if (hours > 0)
    {
        ss << hours << "h";
        if (minutes > 0)
            ss << " " << minutes << "m";
    }
    else if (minutes > 0)
    {
        ss << minutes << "m";
        if (secs > 0)
            ss << " " << secs << "s";
    }
    else
    {
        ss << secs << "s";
    }
    return ss.str();
}

/// Get the estimated completion time as a formatted string
inline std::string getCompletionTime(double eta_seconds)
{
    std::time_t now = std::time(nullptr);
    std::time_t complete_time = now + static_cast<std::time_t>(eta_seconds);

    std::tm *lt = std::localtime(&complete_time);
    char buffer[64];
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", lt); // e.g., "17:45:30"
    return std::string(buffer);
}

inline std::string getTimestamp()
{
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&in_time_t), "%Y%m%d_%H%M%S");
    return ss.str();
}

/// Timer class for nested timing and progress estimation
class Timer
{
public:
    /**
     * @brief Constructor for Timer
     * @param name Name of the timer
     * @param parent Optional parent timer for hierarchy
     */

    // 默认构造函数
    Timer()
        : name("main"), running(false), terminated(false), elapsed_us(0), parent(nullptr)
    {
    }

    Timer(const std::string &name, Timer *parent = nullptr)
        : name(name), running(false), terminated(false), elapsed_us(0), parent(parent)
    {
        if (parent)
        {
            parent->children.push_back(this);
        }
    }

    ~Timer() {}

    void reset(const std::string &newName, Timer *newParent = nullptr)
    {
        // Remove self from previous parent's children list if needed
        if (parent)
        {
            auto &siblings = parent->children;
            siblings.erase(std::remove(siblings.begin(), siblings.end(), this), siblings.end());
        }

        // Update name and parent
        name = newName;
        parent = newParent;

        // Add to new parent's children list if applicable
        if (parent)
        {
            parent->children.push_back(this);
        }

        // Reset internal states
        running = false;
        terminated = false;
        elapsed_us = 0;
        startTime = std::chrono::high_resolution_clock::time_point(); // Reset start time
        children.clear();                                             // Remove all nested children (optional based on use case)
    }

    /// Start the timer
    void start()
    {
        if (terminated)
        {
            std::cerr << "[Timer Warning] Attempting to start a terminated timer \"" << name << "\".\n";
            return;
        }

        if (running)
        {
            std::cerr << "[Timer Warning] Timer \"" << name << "\" is already running. start() ignored.\n";
            return;
        }

        // 新增检查：如果该计时器有父计时器，则父计时器必须处于运行状态
        if (parent)
        {
            if (!parent->running)
            {
                std::cerr << "[Timer Error] Cannot start timer \"" << name
                          << "\" because its parent timer \"" << parent->name
                          << "\" is not running.\n";
                return;
            }

            // 新增检查：同一个父计时器下不允许有其他兄弟计时器正在运行
            for (auto sibling : parent->children)
            {
                if (sibling != this && sibling->running)
                {
                    std::cerr << "[Timer Error] Cannot start timer \"" << name
                              << "\" because sibling timer \"" << sibling->name
                              << "\" is already running.\n";
                    return;
                }
            }
        }

        startTime = std::chrono::high_resolution_clock::now();
        running = true;
    }

    /// Pause the timer (accumulates elapsed time)
    void pause()
    {
        if (terminated)
        {
            std::cerr << "[Timer Warning] Timer \"" << name << "\" is terminated. Cannot pause.\n";
            return;
        }

        if (!running)
        {
            std::cerr << "[Timer Warning] Timer \"" << name << "\" is not running. pause() ignored.\n";
            return;
        }

        for (auto child : children)
        {
            if (child->running)
            {
                std::cerr << "[Timer Error] Cannot pause timer \"" << name
                          << "\" because child timer \"" << child->name
                          << "\" is still running.\n";
                return;
            }
        }

        auto now = std::chrono::high_resolution_clock::now();
        elapsed_us += std::chrono::duration_cast<std::chrono::microseconds>(now - startTime).count();
        running = false;
    }

    /// Stop the timer (finalize timing and mark as terminated)
    void stop()
    {
        if (terminated)
        {
            std::cerr << "[Timer Warning] Timer \"" << name << "\" is already terminated. stop() ignored.\n";
            return;
        }

        pause();           // Stop and accumulate time
        terminated = true; // Mark as done
    }

    /// Get the elapsed time in seconds
    double getElapsed() const
    {
        if (running)
        {
            auto now = std::chrono::high_resolution_clock::now();
            return (elapsed_us + std::chrono::duration_cast<std::chrono::microseconds>(now - startTime).count()) / 1e6;
        }
        return elapsed_us / 1e6;
    }

    /**
     * @brief Update and print progress estimation info to console
     * @param progress A value between 0.0 and 1.0
     */
    void updateProgress(double progress)
    {
        if (!running)
        {
            std::cerr << "[Timer Warning] updateProgress() called while timer \"" << name << "\" is not running.\n";
            return;
        }

        double elapsedTime = getElapsed();
        double eta = (progress > 0.0 && progress <= 1.0) ? (elapsedTime / progress) * (1.0 - progress) : 0.0;

        std::string elapsedStr = formatTime(elapsedTime);
        std::string etaStr = formatTime(eta);
        std::string completionTime = getCompletionTime(eta);

        std::cout << "\r"
                  << "Progress: \033[1;32m" << std::fixed << std::setprecision(2) << (progress * 100.0) << "%\033[0m "
                  << "| Elapsed: \033[1;34m" << elapsedStr << "\033[0m "
                  << "| ETA: \033[1;33m" << etaStr << "\033[0m "
                  << "| Complete at: \033[1;35m" << completionTime << "\033[0m"
                  << std::flush;
    }

    /// @brief 递归打印计时器及其子计时器的时间使用情况
    /// @param level 缩进层级，用于格式化输出
    /// @param showCumulative 是否显示累计耗时占比，默认为false
    /// @param rootElapsed 根计时器总耗时，用于计算百分比
    void printHierarchy(bool showCumulative = false, int level = 0, double rootElapsed = -1.0) const
    {
        std::string indent(level * 4, ' ');
        double timeSec = getElapsed();

        if (level == 0)
        {
            // 初始化根计时器耗时
            rootElapsed = timeSec;
            std::cout << name << ": " << formatTime(timeSec) << std::endl;
        }

        // 按耗时降序排序子计时器
        std::vector<Timer *> sortedChildren = children;
        std::sort(sortedChildren.begin(), sortedChildren.end(),
                  [](Timer *a, Timer *b)
                  {
                      return a->getElapsed() > b->getElapsed();
                  });

        double cumulativePercent = 0.0;
        for (auto child : sortedChildren)
        {
            double childElapsed = child->getElapsed();
            double childPercent = (rootElapsed > 0.0) ? (childElapsed / rootElapsed * 100.0) : 0.0;
            cumulativePercent += childPercent;

            std::cout << indent << "    " << child->name << ": " << formatTime(childElapsed);
            if (level >= 0) // 从第一层级开始显示百分比
            {
                std::cout << ", Percent: " << std::fixed << std::setprecision(2) << childPercent << "%";
                if (showCumulative)
                {
                    std::cout << " " << std::setprecision(2) << cumulativePercent << "%";
                }
            }
            std::cout << std::endl;

            // 递归打印子节点，保持参数传递一致性
            child->printHierarchy(showCumulative, level + 1, rootElapsed);
        }
    }

    /// Debug打印：显示当前计时器树的状态，包括运行状态和时间
    void debug(int level = 0) const
    {
        std::string indent(level * 4, ' ');
        std::cout << indent << "- " << name << " | ";

        if (terminated)
            std::cout << "Stopped";
        else if (running)
            std::cout << "Running";
        else
            std::cout << "Paused";

        std::cout << " | Elapsed: " << std::fixed << std::setprecision(2) << getElapsed() << "s\n";

        for (const auto &child : children)
        {
            child->debug(level + 1);
        }
    }

private:
    std::string name;                                         ///< Timer name
    bool running;                                             ///< Whether the timer is currently running
    bool terminated;                                          ///< Whether the timer has been stopped
    long long elapsed_us;                                     ///< Total elapsed time in microseconds
    Timer *parent;                                            ///< Pointer to parent timer (for hierarchy)
    std::chrono::high_resolution_clock::time_point startTime; ///< Start time
    std::vector<Timer *> children;                            ///< Child timers
};

#endif // TIMER_HPP
