# 「极限阶乘」高性能计算挑战

## 1. 挑战内容

给定：

$$
k\in[4,8]
$$

计算：

$$
n=10^k+\text{student\_id}\times 1000
$$

本项目中：

```text
student_id = 5
```

因此需要计算：

|  k |           n |
| -: | ----------: |
|  4 |      15,000 |
|  5 |     105,000 |
|  6 |   1,005,000 |
|  7 |  10,005,000 |
|  8 | 100,005,000 |

要求程序在 **3 GB 内存、600 s 时间限制**内完成计算，且不能调用第三方大整数库完成主体计算。

## 2. 项目简介

本项目采用 Product Tree 计算大规模阶乘，并针对超大整数乘法进行了专门优化。

核心方法包括：

* PrimeSwing + Product Tree 分治计算阶乘
* `10^14` 进制大整数压位
* 双模数 NTT
* CRT 重构
* Montgomery Reduction
* OpenMP 并行 NTT
* 针对不同规模选择不同的乘法路径
* 大规模 NTT 中复用工作空间以控制内存占用

大整数采用 `10^14` 作为进制，每个 limb 使用 `uint64_t` 存储，从而减少超大整数所需的 limb 数量。

目录结构：

```
FastFactory/
├── src/
│   ├── bigint.h          # 大整数实现及乘法
│   ├── common.h          # 公共类型、常量及辅助函数
│   ├── mulmod.h          # 模乘及 Montgomery Reduction
│   ├── primeswing.h      # 实现 PrimeSwing + Product Tree
│   ├── solution.cpp      # 主程序
│   ├── speed.h           # 性能测试及计时相关代码
│   └── verify.cpp        # 独立正确性验证程序
│
├── benchmark.sh          # 自动编译并测试 k=4~8，生成 fingerprint.txt
├── extreme_test.sh       # 自动极限测试 n=3e8
├── extreme_test.txt      # 极限测试 n=3e8 的结果
├── verify.sh             # 自动执行独立正确性验证
├── Makefile              # Linux 下的编译入口
├── Dockerfile            # 统一测试环境
├── fingerprint.txt       # 各个 k 的结果指纹、用时及峰值内存
└── README.md             # 项目说明、算法介绍及测试结果
```

## 3. 编译与运行

### 3.1 直接编译运行

直接编译 `src/solution.cpp`：

```bash
g++ src/solution.cpp -o solution -Ofast -fopenmp -march=native
```

运行时传入 `n`：

```bash
./solution 15000
```

### 3.2 Make

项目提供 Makefile 时，可以直接在项目根目录执行：

```bash
make
```

然后运行：

```bash
./solution 15000
```

### 3.3 Docker 自动测试

项目提供 Docker 环境用于统一测试环境。

构建镜像：

```bash
docker build -t judge-cpp:13 .
```

然后运行：

```bash
bash benchmark.sh
```

`benchmark.sh` 会自动完成编译，并依次测试：

```text
k = 4
k = 5
k = 6
k = 7
k = 8
```

测试环境固定使用：

```text
Ubuntu 24.04
OpenMP = 8 threads
CPU affinity = CPU 0-7
Memory limit = 3 GiB
Network = disabled
Compiler optimization = -Ofast
```

测试结果自动保存到：

```text
fingerprint.txt
```

**【补充】** 如果在已经进入 Docker 容器的环境中运行 `benchmark.sh`，脚本会直接使用当前环境；在宿主机运行时，如果检测到 Docker，则会自动启动指定镜像执行测试。

## 4. 自测数据

### 4.1 测量方法

性能测试使用 Linux 的 `/usr/bin/time` 测量。

示例：

```bash
/usr/bin/time -f 'Time: %es\nMemory: %M' ./solution 15000
```

其中：

* `Time`：程序实际 wall-clock elapsed time，单位为秒
* `Memory`：程序运行期间的最大 Resident Set Size（RSS），单位为 KiB

最终将 RSS 转换为 MB：

```text
Peak Memory (MB) = Peak RSS (KiB) / 1024
```

Docker 中使用：

```text
--memory=3g
--memory-swap=3g
```

作为内存上限。

需要注意，3 GB 是运行环境的**最大内存限制**，而 `Memory` 是程序实际测得的**峰值 RSS**，两者并不相同。

该自测已经在 `benchmark.sh` 中实现，具体见上面。

### 4.2 自测结果

在上述固定测试环境下得到：

| k   | n           | 用时      | 峰值内存       |
| ---:| -----------:| -------:| ----------:|
| 4   | 15,000      | 0.02 s  | 5.25 MB    |
| 5   | 105,000     | 0.06 s  | 7.36 MB    |
| 6   | 1,005,000   | 0.48 s  | 24.96 MB   |
| 7   | 10,005,000  | 5.89 s  | 312.79 MB  |
| 8   | 100,005,000 | 62.99 s | 2404.18 MB |

其中 `k=8`：

```text
n = 100005000
```

计算得到：

```text
DIGITS = 756610557
```

即结果包含约 **7.57 亿个十进制数字**。

最大测试规模的峰值 RSS 为：

```text
2404.18 MB
```

低于 3 GiB 的 Docker 内存限制，并在 600 s 时间限制内完成。

**【补充】** 以上数据来自固定的 Docker 环境和固定的 8 个 OpenMP 线程。实际运行时间会受到 CPU 型号、系统负载等因素影响，因此这里的时间主要用于记录本次自测环境下的结果，而不是声称所有机器都能达到相同时间。

## 5. 正确性验证

由于 `k=8` 时结果包含约 **7.57 亿位十进制数字**，不可能通过人工逐位检查整个结果。

因此，本项目没有仅依赖最终输出是否“看起来正确”，而是采用多种相互独立的方式进行交叉验证。

### 5.1 GMP 交叉验证

对于 `k=4~8`，将程序输出的以下五项 fingerprint 与独立 GMP 实现的计算结果进行比较：

```text
DIGITS
ZEROS
DIGITSUM
HEAD50
TAIL50
```

目前 `k=4~8` 的全部 fingerprint 均与 GMP 结果一致，没有发现差异。

其中 `k=8` 的结果为：

```text
DIGITS   756610557
ZEROS    25001248
DIGITSUM 3292323588
HEAD50   18325739734648128139785673730312153899158313606708
TAIL50   90819894298634815510347177387550939212635881603072
```

这一步能够从整体上检查最终结果。

### 5.2 末尾 0 的数学验证

对于阶乘：

$$
n! = 1\times2\times\cdots\times n
$$

十进制末尾 `0` 的数量为：

$$
v_{10}(n!)=\min(v_2(n!),v_5(n!))
$$

由于 `n!` 中 2 的因子数量远多于 5，因此：

$$
v_{10}(n!)=v_5(n!)
$$

根据 Legendre 公式：

$$
v_5(n!)=
\left\lfloor\frac n5\right\rfloor+
\left\lfloor\frac n{25}\right\rfloor+
\left\lfloor\frac n{125}\right\rfloor+\cdots
$$

因此，`ZEROS` 可以完全不依赖大整数乘法重新计算。

程序一方面从最终 `BigInt` 的 `vector<ull> w` 中直接统计末尾 `0`，另一方面使用 Legendre 公式独立计算 `v5(n!)`，然后比较两者。

例如 `k=8`：

```text
BigInt trailing zeros = 25001248
Legendre v5(n!)       = 25001248
```

两者一致。

这种验证不依赖 NTT、CRT 或 Product Tree。

### 5.3 多模数独立验证

为了进一步验证最终大整数结果，使用多个模数进行检查：

```text
1000000007
1000000009
998244353
100005001
100005023
```

对于每一个模数 `p`，采用两条不同的计算路径。

第一条路径从最终 `BigInt` 计算：

```text
BigInt mod p
```

由于本项目使用 `10^14` 作为进制，可以直接从最高 limb 到最低 limb 计算：

```text
r = (r * BASE + limb) mod p
```

第二条路径完全不使用 BigInt、Product Tree 或 NTT，而是直接计算：

```text
1 × 2 × 3 × ... × n mod p
```

然后比较：

```text
BigInt mod p
=
independent n! mod p
```

所有测试模数的结果均一致。

因此这一步提供了一条与主算法明显不同的独立验证路径。

**【补充】** 后两种验证方法可以通过运行 `verify.sh` 自动完成。

### 5.4 为什么可以相信 7.6 亿位结果？

本项目的正确性并不是建立在“把 7.6 亿位数字全部看一遍”的基础上，而是建立在多种独立检查之上，其中最重要的是多模数验证。

该两条路径没有复用 Product Tree、NTT 或 CRT 的计算过程。

因此，如果主程序的 NTT、CRT、Product Tree 或大整数进位出现错误，最终结果的模值通常会与独立计算得到的 `n! mod p` 不一致。

同时，GMP 对 fingerprint 的交叉验证以及 Legendre 对末尾 0 的验证，又从另外两个角度提供了独立检查。

因此，对于约 7.6 亿位的最终结果，可以通过**独立实现、数学恒等式和模运算验证的组合结果**建立对正确性的信心，而不是依赖人工逐位检查。

## 6. 遇到的问题与优化

### 6.1 OpenMP 在小规模 NTT 中的额外开销

OpenMP 并行并不意味着所有规模都能获得加速。

在较小长度的 NTT 中，线程创建、同步以及任务划分产生的开销可能超过并行计算本身，从而导致运行速度下降。

因此代码中设置了 `OMP_PIVOT`，只在达到一定规模后启用 OpenMP，以降低小规模计算中的并行开销。

### 6.2 大整数内存占用

> 这一点是在加入 PrimeSwing 之前遇到的，算是在尝试过程中试过的解决方案，在最终版本中已经不需要使用了。

`k=8` 时最终结果约有 7.57 亿位，即使采用 `10^14` 压位，仍需要大量 limb 存储。

因此大规模乘法中的临时 NTT 数组会成为重要的内存开销来源。

实测发现，在 $k=8$ 时最高可能会达到 $\text{lim}=2^{26}$，此时若同时存在 $4$ 个该规模的 $64$ 位整数数组，那么内存开销会直逼 3 GB，不可接受。

开始我尝试将 $4$ 个数组压到 $3$ 个，但是由于双模数的需求，$4$ 个数组似乎无法避免。

此时其实可以选择分块 NTT，但同时会增加：

- 数据搬运
- block 管理
- 多次变换
- 额外的边界处理
- 实现复杂度

等多个问题。因此我选择了另一种方法：仅在 $\text{lim}=2^{26}$ 时，通过多一次 $iNTT$ 来节省一个数组的内存开销。具体地：

```
input x, y
z = x
NTT1(z), NTT1(y)
z *= y
iNTT1(z) # 得到一个模数的答案 z
# 此时 y 的原始数据已经被破坏，所以这里用 iNTT1 来还原（此处会多一次计算开销，但减少一个数组的内存）
iNTT1(y)
NTT2(x), NTT2(y)
x *= y
iNTT2(x)
CRT(x, z) => Answer
```

最终在当前测试环境下：

```text
k=8
Peak RSS = 2052.90 MB
```

在 3 GB 限制内完成计算。

### 6.3 DIF 与 DIT

由于在 $lim$ 较大时，蝴蝶变换会面临访问不连续、占用大数组等诸多不便，所以这里用 DIF 与 DIT 规避掉了蝴蝶变换，同时加快了运行效率。

### PrimeSwing 与自乘

从普通 Product Tree 改成 PrimeSwing 的一大优点便是可以用自乘代替部分乘法，这时只需要一次 NTT/DIF，能够做到显著优化。

## 7. 关于 AI 使用

本项目使用 AI 主要作为辅助工具，而非直接生成整个算法实现。

### AI 主要参与的部分

* `verify.sh`、`benchmark.sh`、Dockerfile、Makefile 等辅助工具
* `speed.h` 中与性能测试相关的部分
* `verify.cpp` 中与正确性测试相关的部分
* Montgomery, PrimeSwing 等算法的辅助理解
* README 的结构与措辞润色
* 算法优化方向的讨论与参考
* 部分代码 Debug
* 正确性验证方案的讨论

### 由本人完成的主要部分

* `bigint.h`
* `common.h`
* `mulmod.h`
* `primeswing.h`
* `solution.cpp`

以及上述算法相关代码的主体设计与主要编写工作。

最终采用的算法路线、数据结构和优化方案由本人分析并决定。

**【补充】** 在使用 AI 辅助 Debug 和优化时，最终代码均经过本人理解、修改和实际测试，而不是直接将 AI 输出作为未经验证的程序提交。
