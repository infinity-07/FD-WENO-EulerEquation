# WENO Finite Volume Solver (WENOFV)

## 🔍 简介

`WENOFD` 是一个基于 **加权本质无振荡有限差分方法（WENO-FD）** 的二维双曲守恒律数值求解器。该程序适用于多种测试算例，支持高阶时间推进、边界条件处理及 MPI 并行计算，适合科研与教学使用。

## ✨ 特点
- 支持二维带源项欧拉方程组求解并且可以很便捷地推广到其他二维守恒型方程
- 支持 WENO/WENO-Z 重构（五阶）
- 支持 RK1/RK2/RK3 时间推进
- 多种边界条件：周期、对称、滑移、特殊
- 支持 MPI 并行计算
- Tecplot 可视化输出
- 自动误差评估与精度验证

## 📁 文件结构


```

├── WENOFV.cpp             # 主程序实现
├── WENOFV.hpp             # 类声明（需要提供）
├── Equations.hpp/.cpp     # 方程定义与参数设置（需要提供）
├── Utility.hpp/.cpp       # 工具函数与定时器（需要提供）
├── input/
│   └── config.cfg         # 配置文件
├── output/
│   ├── average\_*.plt      # 平均变量输出
│   ├── error\_*.plt        # 解误差图
│   └── accuracy\_\*.csv     # 误差范数表格
└── Makefile               # 编译脚本（推荐添加）

````


## ⚙️ 编译说明

依赖：

- 支持 C++17 的编译器（如 `g++`）
- MPI 库（如 OpenMPI 或 MPICH）

编译命令示例：

```bash
mpic++ -std=c++17 -O3 WENOFV.cpp Utility.cpp Equations.cpp -o wenofv
````

或使用 `Makefile`：

```bash
make
```

## 🚀 运行方式

1. 修改 `input/config.cfg` 文件，设置模拟参数：

| 参数名      | 含义                             |
| -------- | ------------------------------ |
| TESTCASE | 测试算例编号                         |
| SCHEME   | 重构方案（0=WENO, 1=WENOZ）          |
| CFL      | CFL 数                          |
| RKMETH   | Runge-Kutta 方法（0,1,2 对应 RK1-3） |
| ELEMNUMX | x 方向网格数                        |
| ELEMNUMY | y 方向网格数                        |

2. 运行程序（以 4 个线程为例）：

```bash
mpirun -np 4 ./wenofv
```


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


## 📌 注意事项

* 请确保 `input/config.cfg` 文件存在，且格式正确
* 所有输出文件默认保存至 `./output/` 文件夹
* 若使用并行，请确保节点间文件系统一致（如共享目录）

## 🙋 联系方式

如需协助或合作，请联系作者：

邮箱: cauzhangyang@outlook.com

