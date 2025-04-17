#include "seekroots.h"
#include <windows.h>
#include <conio.h> // 用于 _getch() 函数

int main()
{
    // 设置控制台输入和输出编码为 UTF-8
    SetConsoleCP(65001); // 设置控制台输入代码页为 UTF-8
    SetConsoleOutputCP(65001); // 设置控制台输出代码页为 UTF-8

    while (true)
    {
        vector<double> coefficients;
        cout << "请输入多项式的系数（从高到低），以空格分隔： ";
        string input;
        getline(cin, input);
        stringstream ss(input);
        double coefficient;
        while (ss >> coefficient)
        {
            coefficients.push_back(coefficient);
        }

        SeekRoots function(coefficients);
        function.FindRoots();
        function.PrintFunction();
        function.PrintRoots();

        cout << "是否继续？(y/n)  ";
        char choice;
        cin >> choice;
        while(getchar()!='\n');
        if (choice != 'y' && choice != 'Y')
        {
            break;
        }
    }

    // SeekRoots function1({1,-4,0,16,-16}); // !test
    // function1.FindRoots();
    // function1.PrintFunction();
    // function1.PrintRoots();

    // SeekRoots function2({1,0,-5,0,4,0}); // !test
    // function2.FindRoots();
    // function2.PrintFunction();
    // function2.PrintRoots();

    // SeekRoots function3({1,-10,-3400,34000,2250000,-22500000}); // !test
    // function3.FindRoots();
    // function3.PrintFunction();
    // function3.PrintRoots();

    // SeekRoots function4({1,-265,25885,-1138175,21591794,-130945680}); // !test
    // function4.FindRoots();
    // function4.PrintFunction();
    // function4.PrintRoots();

    cout << "按任意键退出...";
    _getch();
    return 0;
}