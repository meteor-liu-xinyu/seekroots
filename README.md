# README.md

# 方程求根程序

## 程序简介

这是一个用 C++编写的方程求根程序，适用于一次至四次方程。对于这些方程，程序使用标准的求根公式来计算复数域上的所有根。对于五次及以上的方程，程序采用牛顿迭代法来求解，并尽可能找到所有的根。需要注意的是，当根的值较大或分布较为分散时，求得的根可能不够精确。

## 功能特点

- **一次至四次方程**：直接使用求根公式计算所有根。
- **五次及以上方程**：使用牛顿迭代法求解。
- **复数域**：支持计算复数根。
- **用户友好**：简单的命令行界面，易于使用。

## 使用方法
1. 运行bin_static文件夹中的main_static.exe程序
2. 将方程的系数以空格分隔的形式输入。
3. 运行程序，程序将输出方程的所有根。

## 示例

对于方程  $ ax^4 + bx^3 + cx^2 + dx + e = 0 $ ，输入系数 `a b c d e`，程序将输出所有根。

## 许可证

本程序采用 MIT 许可证，详情见 LICENSE 文件。

## Example
For the equation  $ ax^4 + bx^3 + cx^2 + dx + e = 0 $ , enter the coefficients  $ a \ b \ c \ d \ e $ , and the program will output all roots.

## License
This program is licensed under the MIT License. See the LICENSE file for details.
