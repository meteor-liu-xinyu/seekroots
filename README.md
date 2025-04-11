# README.md

# 方程求根程序

<!-- HTML注释实现中英文切换 -->
<div>
    <button onclick="toggleLanguage('zh')">中文</button>
    <button onclick="toggleLanguage('en')">English</button>
</div>

<script>
    function toggleLanguage(lang) {
        var elements = document.querySelectorAll('.lang');
        elements.forEach(function(element) {
            if (element.classList.contains(lang)) {
                element.style.display = 'block';
            } else {
                element.style.display = 'none';
            }
        });
    }
</script>

## <span class="lang zh">程序介绍</span>
## <span class="lang en">Program Introduction</span>

<span class="lang zh">
  
  ## 程序简介
  
  这是一个用 C++编写的方程求根程序，适用于一次至四次方程。对于这些方程，程序使用标准的求根公式来计算复数域上的所有根。对于五次及以上的方程，程序采用牛顿迭代法来求解，并尽可能找到所有的根。需要注意的是，当根的值较大或分布较为分散时，求得的根可能不够精确。
  
  ## 功能特点
  
  - **一次至四次方程**：直接使用求根公式计算所有根。
  - **五次及以上方程**：使用牛顿迭代法求解。
  - **复数域**：支持计算复数根。
  - **用户友好**：简单的命令行界面，易于使用。
  
  ## 使用方法
  
  1. 将方程的系数以空格分隔的形式输入。
  2. 运行程序，程序将输出方程的所有根。
  
  ## 示例
  
  对于方程  $ ax^4 + bx^3 + cx^2 + dx + e = 0 $ ，输入系数 `a b c d e`，程序将输出所有根。
  
  ## 许可证
  
  本程序采用 MIT 许可证，详情见 LICENSE 文件。

</span>
<span class="lang en">
  ## Program Introduction
  This is a root-finding program written in C++, suitable for equations from first to fourth degree. For these equations, the program uses standard root-finding formulas to calculate all roots in the complex domain. For equations of fifth degree and higher, the program employs the Newton-Raphson iteration method to find as many roots as possible. Note that the roots may not be precise when they are large or widely scattered.
  
  ## Features
  - **First to Fourth Degree Equations:** Directly uses root-finding formulas to calculate all roots.
  - **Fifth Degree and Higher Equations:** Uses Newton-Raphson iteration method for solving.
  - **Complex Domain:** Supports calculation of complex roots.
  - **User-Friendly:** Simple command-line interface, easy to use.
  
  ## Usage
  1. Enter the coefficients of the equation separated by spaces.
  2. Run the program, and it will output all roots of the equation.
  
  ## Example
  For the equation  $ ax^4 + bx^3 + cx^2 + dx + e = 0 $ , enter the coefficients  $ a \ b \ c \ d \ e $ , and the program will output all roots.

  ## License
  This program is licensed under the MIT License. See the LICENSE file for details.
  
</span>
