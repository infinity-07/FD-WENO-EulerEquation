#ifndef _NEWCOMMON_H_
#define _NEWCOMMON_H_
#include <iostream>
#include <chrono>
#include <thread>
#include <iomanip>
#include <fstream>
#include <string>
#include <ctime>

inline void showProgressBar(double total, double progress, double elapsedTime)
{
    const int barWidth = 50; // 进度条宽度
    double progressFraction = progress / total;

    // 计算预计剩余时间
    double estimatedTime = elapsedTime / progressFraction - elapsedTime;

    // 时间转换函数
    auto formatTime = [](double timeInSeconds)
    {
        int days = static_cast<int>(timeInSeconds) / 86400;
        timeInSeconds -= days * 86400;
        int hours = static_cast<int>(timeInSeconds) / 3600;
        timeInSeconds -= hours * 3600;
        int minutes = static_cast<int>(timeInSeconds) / 60;
        timeInSeconds -= minutes * 60;
        int seconds = static_cast<int>(timeInSeconds);

        return std::make_tuple(days, hours, minutes, seconds);
    };

    auto [elapsedDays, elapsedHours, elapsedMinutes, elapsedSeconds] = formatTime(elapsedTime);
    auto [estimatedDays, estimatedHours, estimatedMinutes, estimatedSeconds] = formatTime(estimatedTime);

    std::cout << "[";
    int pos = barWidth * progressFraction;

    // 输出进度条
    for (int i = 0; i < barWidth; ++i)
    {
        if (i < pos)
            std::cout << "=";
        else if (i == pos)
            std::cout << ">";
        else
            std::cout << "-"; // 用 '-' 补全未完成的部分
    }

    // 输出进度百分比、已用时间和预计剩余时间
    std::cout << "] " << std::fixed << std::setprecision(2) << progressFraction * 100.0 << " % | "
              << "已用时: " << elapsedDays << "d " << elapsedHours << "h "
              << elapsedMinutes << "m " << elapsedSeconds << "s | "
              << "预计剩余: " << estimatedDays << "d " << estimatedHours << "h "
              << estimatedMinutes << "m " << estimatedSeconds << "s\r";
    std::cout.flush(); // 刷新输出缓冲区
}

inline void logMessage(const std::string &message)
{
    // 打开日志文件，使用 std::ios::app 以追加的方式打开文件
    std::ofstream logFile("logfile.log", std::ios::app);

    // 检查文件是否成功打开
    if (logFile.is_open())
    {
        // 获取当前时间
        std::time_t now = std::time(nullptr);
        std::tm *localTime = std::localtime(&now);

        // 格式化时间字符串
        char timeBuffer[80];
        std::strftime(timeBuffer, sizeof(timeBuffer), "%Y-%m-%d %H:%M:%S", localTime);

        // 写入时间戳和消息到日志文件
        logFile << "[" << timeBuffer << "] " << message << std::endl;

        // 关闭文件
        logFile.close();
    }
    else
    {
        std::cerr << "无法打开日志文件!" << std::endl;
    }
}
#endif