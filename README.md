# WENO Finite Difference Solver (WENOFD)

## 🔍 简介

`WENOFD` 是一个基于 **加权本质无振荡有限差分方法（WENO-FD）** 的二维双曲守恒律数值求解器。该程序适用于多种测试算例，支持高阶时间推进、边界条件处理及 MPI 并行计算，适合科研与教学使用。

## ✨ 特点
- [x] 支持二维带源项欧拉方程组求解并且可以很便捷地推广到其他二维守恒型方程
- [x] 支持 WENO/WENO-Z 重构（五阶）
- [x] 支持 RK1/RK2/RK3 时间推进
- [x] 多种边界条件：周期、对称、滑移、特殊
- [x] 支持 MPI 并行计算，支持快速切换串行和并行
- [x] Tecplot 可视化输出
- [x] 自动误差评估与精度验证

## 🚧 开发计划
- [ ] 增加二进制格式的 tecplot 文件输出方法
- [ ] 增加 paraview 文件输出方法
- [ ] 增加保正限制器
- [ ] 给通用文件前增加文件信息
- [ ] 每次运行时能够将数据输出到不同的文件夹，并且备份一份代码
- [ ] 给 readme 文件增加算例计算结果图
- [ ] 给出 MPI 加速比，提升并行效率，一个可行的路线是在 MPI 并行的同时完成计算
- [ ] 现在只支持五阶格式重构，三点格式七点格式很不方便
- [ ] 删掉大部分 prolog 代码
- [ ] 储存一些结果便于对照

## 📁 文件结构
```
├── WENOFV.cpp             # 主程序实现
├── WENOFV.hpp             # 类声明
├── Equations.hpp     # 方程定义与参数设置
├── Utility.hpp       # 工具函数
├── input/
│   └── config.cfg         # 配置文件
├── output/
│   ├── average\_*.plt      # 平均变量输出
│   ├── error\_*.plt        # 解误差图
│   └── accuracy\_\*.csv     # 误差范数表格
└── Makefile               # 编译脚本
````




## 🚀 运行方式
以下运行方式说明均基于 linux 系统

1. 安装 支持 C++17 的编译器（如 `g++`）
- MPI 库（如 OpenMPI 或 MPICH）

2. 修改 `input/config.cfg` 文件，设置模拟参数：

| 参数名      | 含义                             |
| -------- | ------------------------------ |
| TESTCASE | 测试算例编号（                       |
| SCHEME   | 重构方案（0=WENO, 1=WENOZ）          |
| CFL      | CFL 数                          |
| RKMETH   | Runge-Kutta 方法（0,1,2 对应 RK1-3） |
| ELEMNUMX | x 方向网格数                        |
| ELEMNUMY | y 方向网格数                        |

3. 编译程序
```bash
make clean
make
```
 
4. 运行程序（以 4 个线程为例）：

```bash
mpirun -np 4 ./wenofv
```

输出结果储存在 output 文件夹中


## 📤 输出说明

| 文件名格式            | 含义                 |
| ---------------- | ------------------ |
| `average_*.plt`  | Tecplot 平均变量输出     |
| `error_*.plt`    | 与解析解的误差分布          |
| `accuracy_*.csv` | 误差范数表格（L1, L2, L∞） |


## 🔧 主要函数说明

| 函数名                  | 功能描述                      |
| -------------------- | ------------------------- |
| `initializeSolver()` | 初始化网格、变量、MPI分区等           |
| `run()`              | 控制整体模拟流程                  |
| `calculateDeltaT()`  | 根据 CFL 条件计算时间步长           |
| `RunRK1/2/3()`       | 一阶/二阶/三阶 Runge-Kutta 时间推进 |
| `getFlux()`          | 使用 WENO 特征分解计算数值通量        |
| `MPICommunication()` | 进行左右子域数据通信                |
| `outputAve()`        | 输出平均变量到 Tecplot           |
| `outputError()`      | 输出解误差图                    |
| `outputAccuracy()`   | 误差范数评估并写入 csv 文件          |

## 🙋 联系方式

如需协助或合作，请联系作者：

邮箱: cauzhangyang@outlook.com

