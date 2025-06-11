#ifndef READCONFIG_HPP
#define READCONFIG_HPP

#include <map>
#include <string>
#include <iostream>
#include <fstream>

/**
 *	这个函数的功能是去除字符串两端的空格，并将更新后的字符串保存在传入的参数 str 中。
 */
inline void trim(std::string &str)
{
    // 找到第一个非空格字符的位置
    size_t first = str.find_first_not_of(' ');

    // 找到最后一个非空格字符的位置
    size_t last = str.find_last_not_of(' ');

    // 如果最后一个字符为0x0D（回车符），那么 last=last-1
    if (str[last] == '\r')
    {
        last--;
    }

    // 更新字符串，使其只包含非空格字符部分
    str = str.substr(first, (last - first + 1));
}

class CConfig
{

public:
    std::string m_inputFileName;

public:
    CConfig(const std::string &str) : m_inputFileName(str) {}

    void read(std::map<std::string, std::string> &option)
    {
        std::ifstream input(m_inputFileName);

        // 如果无法打开文件，输出错误信息并退出程序
        if (!input.is_open())
        {
            std::cerr << "Failed to open file" << m_inputFileName << std::endl;
            exit(1);
        }

        std::string line;
        while (std::getline(input, line))
        {
            // 跳过空行、以#开头的注释行以及回车符行
            if (line.empty() || line[0] == '#' || line[0] == '\r')
            {
                continue;
            }

            // 将行按等号分割为键和值
            const size_t equal_pos = line.find('=');

            // 如果找不到等号，输出错误信息并继续下一行
            if (equal_pos == std::string::npos)
            {
                std::cerr << "Invalid line: " << line << std::endl;
                continue;
            }

            // 提取键和值
            std::string key = line.substr(0, equal_pos);
            std::string value = line.substr(equal_pos + 1);

            // 去除键和值两端的空格
            trim(key);
            trim(value);

            // 将键值对添加到配置map中
            option[key] = value;
        }
    }

    void output(std::map<std::string, std::string> option)
    {
        // // 输出配置信息
        std::map<std::string, std::string>::iterator it;

        // 输出一份配置信息到文件夹中，方便日后查看
        // 打开文件流，将输出定向到 /output/config.txt 文件
        std::ofstream file("../output/config.txt");
        if (file.is_open())
        {
            for (it = option.begin(); it != option.end(); ++it)
            {
                file << it->first << ": " << it->second << std::endl;
            }
        }
        else
        {
            file << "无法打开文件。" << std::endl;
        }
        // 关闭文件流
        file.close();
    }
};
#endif