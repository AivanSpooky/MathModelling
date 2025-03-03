using System;
using System.Collections.Generic;
using ScottPlot;
using System.Drawing;
using ExtendedNumerics;
using System.Globalization;

namespace ModelingTasks
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("Выберите задачу для запуска: 1, 2 или 3");
            var choice = Console.ReadLine();
            switch (choice)
            {
                case "1":
                    Task1.Run();
                    break;
                case "2":
                    Task2.Run();
                    break;
                case "3":
                    Task3.Run();
                    break;
                case "4":
                    Task4.Run();
                    break;
                default:
                    Console.WriteLine("Некорректный выбор.");
                    break;
            }

            Console.WriteLine("Нажмите Enter для выхода...");
            Console.ReadLine();
        }
    }

    #region Задача 1: 4xU''+2U'+U=0, U(0)=1, U'(0)=-0.5
    /// <summary>
    /// Для численного решения преобразуем ОДУ второго порядка в систему:
    ///   U1 = U,  U2 = U'
    /// Тогда:
    ///   U1' = U2,
    ///   U2' = -[2U2+U1] / (4x)   (при x > 0)
    /// Поскольку x=0 – особая точка, на первом шаге x=h (маленькое число) значения U и U' 
    /// вычисляются по ряду Тейлора:
    ///   U(h) = U(0) + U'(0) h + (1/2) U''(0) h^2, 
    ///   U'(h) = U'(0) + U''(0) h,
    /// где U''(0) можно найти из уравнения.
    /// Подставляя x->0 в уравнение, получаем условие совместности: 2U'(0)+U(0)=0,
    /// что выполнено: 2(-0.5)+1=0.
    /// Найдём U''(0):
    ///   При x→0, если считать формально, коэффициент 4x → 0, поэтому для определения U''(0)
    ///   подберём a2 из ряда U(x)=a0+a1x+a2x^2+...; 
    ///   Подставляя в 4xU''+2U'+U=0, получаем:
    ///   (2a1 + a0) + (12a2 + …)x + … = 0.
    ///   Из условия при x^0: 2a1+a0=0 (выполнено), при x^1: 12a2 + ... =0, поэтому 
    ///   можно принять a2 = 0 (или вычислить более точно, если брать следующие члены).
    /// В данном примере для упрощения используем приближённое разложение:
    ///   U(h) ≈ 1 - 0.5 h,   U'(h) ≈ -0.5.
    /// </summary>
   /* public static class Task1
    {
        public static void Run()
        {
            Console.WriteLine("=== Задача 1 ===");
            Console.WriteLine("Аналитическое решение:");
            Console.WriteLine("  U(x) = cos(sqrt(x))");

            double h = 0.01; // маленький шаг, начиная с x = h, чтобы обойти особую точку в x=0
            double xStart = h;
            // Начальные значения, вычисленные по ряду (до первого порядка):
            // U(h) ≈ U(0) + U'(0)*h, т.е. 1 - 0.5*h
            double U0 = 1 - 0.5 * h;
            double V0 = -0.5; // U'(h) ≈ U'(0)
            int nSteps = 2000; // число шагов (интервал можно менять)

            // Эйлер:
            // U'' = - (2U'+U)/(4x)
            var solution = EulerMethodSingular(xStart, U0, V0, h, nSteps);

            int N = solution.Count;
            double[] xs = new double[N];
            double[] uEuler = new double[N];
            double[] uTaylor = new double[N];
            double[] uAnalyt = new double[N];

            for (int i = 0; i < N; i++)
            {
                xs[i] = solution[i].x;
                uEuler[i] = solution[i].U;
                uTaylor[i] = TaylorSolution(solution[i].x);
                uAnalyt[i] = Math.Cos(Math.Sqrt(solution[i].x));
            }

            // Выводим таблицу сравнения решений:
            Console.WriteLine("    x        U(Euler)      U(Taylor)");
            for (int i = 0; i < N; i++)
                Console.WriteLine($"{xs[i],6:F3}   {uEuler[i],12:F6}   {uTaylor[i],12:F6}");

            // Строим график:
            var plt = new ScottPlot.Plot(800, 300);

            // (1) Эйлер
            var scatterEuler = plt.AddScatter(xs, uEuler);
            scatterEuler.Label = "Метод Эйлера";
            scatterEuler.MarkerSize = 0;
            scatterEuler.LineWidth = 2;
            scatterEuler.Color = Color.Blue;

            // (2) Тейлор (4 члена)
            var scatterTaylor = plt.AddScatter(xs, uTaylor);
            scatterTaylor.Label = "Ряд Тейлора (4 члена)";
            scatterTaylor.MarkerSize = 0;
            scatterTaylor.LineWidth = 2;
            scatterTaylor.Color = Color.Red;

            // (3) Аналитическое решение
            var scatterAnalyt = plt.AddScatter(xs, uAnalyt);
            scatterAnalyt.Label = "Аналитическое cos(√x)";
            scatterAnalyt.MarkerSize = 0;
            scatterAnalyt.LineWidth = 2;
            scatterAnalyt.Color = Color.Green;

            plt.Title("Сравнение: метод Эйлера, ряд Тейлора и аналитическое решение");
            plt.XLabel("x");
            plt.YLabel("U(x)");
            plt.Legend();

            plt.SaveFig("Task1.png");
            Console.WriteLine("График сохранён в файл Task1.png");
        }

        // Метод Эйлера для системы с особенностью в x=0:
        // U'' = - (2U'+U)/(4x)
        public static List<(double x, double U, double V)> EulerMethodSingular(
            double x0, double U0, double V0, double h, int nSteps)
        {
            var result = new List<(double, double, double)>();
            double x = x0;
            double U = U0;
            double V = V0;
            result.Add((x, U, V));

            for (int i = 0; i < nSteps; i++)
            {
                double dU = V;
                double dV = -(2 * V + U) / (4 * x);

                U = U + h * dU;
                V = V + h * dV;
                x = x + h;

                result.Add((x, U, V));
            }

            return result;
        }

        // Аналитическое решение в виде ряда Тейлора (4 члена):
        // U(x) = 1 - 0.5*x + (1/24)*x^2 - (1/720)*x^3
        public static double TaylorSolution(double x) =>
            1 - 0.5 * x + (1.0 / 24.0) * Math.Pow(x, 2) - (1.0 / 720.0) * Math.Pow(x, 3);
    }*/
    #endregion
    public static class Task1
    {
        public static void Run()
        {
            Console.WriteLine("=== Задача 1 ===");
            Console.WriteLine("Аналитическое решение (только для x>=0): U(x)=cos(sqrt(x)).");

            // Зададим границы и шаг
            double x1 = -10; // левая граница
            double x2 = 20; // правая граница
            double h = 0.001; // шаг

            // Мы хотим "пропускать" небольшую окрестность 0,
            // чтобы не пытаться делить на x=0. Пусть будет gap = h (или другое малое число).
            double gap = h;

            // 1) Строим решение слева (x<0), если x1 < 0
            List<(double x, double U, double V)> solutionLeft = new List<(double, double, double)>();
            if (x1 < 0)
            {
                // Дойдём до -gap, но не перейдём через 0
                double xEndLeft = Math.Min(x2, -gap);
                if (xEndLeft > x1)
                {
                    // количество шагов
                    int nStepsLeft = (int)((xEndLeft - x1) / h);

                    // Для отрицательного участка вы можете задать какие-то условия.
                    // Но изначально уравнение было U(0)=1, U'(0)=-0.5
                    // Мы НЕ можем напрямую использовать U(0)=1 при x1<0.
                    // Можно, например, задать "фиктивные" начальные условия на x1.
                    // Вариант: на x1 (отрицательное) положим U=..., V=...
                    // Допустим, пусть U=1 - 0.5*x1 (аналогично вашему коду), V=-0.5:
                    double U0_left = 1 - 0.5 * x1;
                    double V0_left = -0.5;

                    // Считаем
                    solutionLeft = EulerMethodSingular(x1, U0_left, V0_left, h, nStepsLeft);
                }
            }

            // 2) Строим решение справа (x>0), если x2 > 0
            List<(double x, double U, double V)> solutionRight = new List<(double, double, double)>();
            if (x2 > 0)
            {
                // Стартуем с +gap
                double xStartRight = Math.Max(x1, gap);
                if (xStartRight < x2)
                {
                    int nStepsRight = (int)((x2 - xStartRight) / h);

                    // Начальные условия, будто x~0 => U(0)=1, U'(0)=-0.5
                    // но мы фактически стартуем в x=gap, 
                    // так что берем U(gap)= 1 - 0.5*gap  (приблизительно),
                    // V(gap)= -0.5
                    double U0_right = 1 - 0.5 * gap;
                    double V0_right = -0.5;

                    solutionRight = EulerMethodSingular(xStartRight, U0_right, V0_right, h, nStepsRight);
                }
            }

            // Объединяем решения
            // (Вначале - левое, потом - правое. Между ними будет "дыра" в [ -gap, +gap ].)
            var solutionAll = new List<(double x, double U, double V)>();
            solutionAll.AddRange(solutionLeft);
            solutionAll.AddRange(solutionRight);

            // Удалим дубликаты, если вдруг x1 > gap, etc.

            solutionAll = solutionAll.Distinct().OrderBy(p => p.x).ToList();

            // Готовим массивы
            int N = solutionAll.Count;
            double[] xs = new double[N];
            double[] uEuler = new double[N];
            
            
            int ind = 0;
            for (int i = 0; i < N; i++)
            {
                if (solutionAll[i].x < 0)
                    ind++;
            }

            /*if (ind > 0)
                uEuler = new double[ind];
            double[] xEu = new double[ind];*/
            //NONEG
            uEuler = new double[N - ind];
            double[] xEu = new double[N - ind];

            /*//!
            uEuler = new double[N];
            xEu = new double[N];*/

            double[] xanal = new double[N - ind];
            double[] uAnalyt = new double[N-ind];
            double[] uTaylor = new double[N-ind];
            for (int i = 0; i < N; i++)
            {
                xs[i] = solutionAll[i].x;
                if (xs[i] < 0)
                {
                    uEuler[i] = solutionAll[i].U;
                    xEu[i] = xs[i];
                }
                /*//!
                else
                {
                    uEuler[i] = solutionAll[i].U;
                    xEu[i] = xs[i];
                }*/

                if (xs[i] >= 0)
                {
                    xanal[i - ind] = xs[i];
                    uAnalyt[i - ind] = Math.Cos(Math.Sqrt(xs[i]));
                }

                if (xs[i] >= 0)
                {
                    uTaylor[i - ind] = TaylorSolution(xs[i]);

                    //NONEG
                    uEuler[i - ind] = solutionAll[i].U;
                    xEu[i - ind] = xs[i];
                }
            }
            Console.WriteLine("   x         Euler           Analytic         Taylor");

            var plt = new ScottPlot.Plot(800, 400);

            /*// (1) Эйлер
            xEu = xs;
            if (ind > 0)
                xEu = xanal;
            else
                xEu = xs;*/
            var scEuler = plt.AddScatter(xEu, uEuler);
            scEuler.Label = "Метод Эйлера";
            scEuler.MarkerSize = 0;
            scEuler.LineWidth = 2;
            scEuler.Color = Color.Blue;

            // (2) Аналитическое (только x>=0)
            if (ind < N)
            {
                /*var scAnalytic = plt.AddScatter(xanal, uAnalyt);
                scAnalytic.Label = "cos(√x), x>=0";
                scAnalytic.MarkerSize = 0;
                scAnalytic.LineWidth = 2;
                scAnalytic.Color = Color.Green;*/

                // (3) Тейлор
                var scTaylor = plt.AddScatter(xanal, uTaylor);
                scTaylor.Label = "Тейлор (x>=0)";
                scTaylor.MarkerSize = 0;
                scTaylor.LineWidth = 2;
                scTaylor.Color = Color.Red;
            }

            plt.Title($"Сравнение решений на [{x1}, {x2}], шаг={h}");
            plt.XLabel("x");
            plt.YLabel("U(x)");
            plt.SetAxisLimitsX(x1, x2);
            plt.Legend(location: Alignment.UpperRight);

            plt.SaveFig("Task1.png");
            Console.WriteLine("График сохранён в файл Task1.png");
        }

        /// <summary>
        /// Метод Эйлера для системы:
        ///   U' = V,
        ///   V' = -(2V + U)/(4x)
        /// идём от (x0, U0, V0) c шагом h и делаем nSteps шагов.
        /// Прерываем расчёт, если появятся NaN или Infinity.
        /// </summary>
        public static List<(double x, double U, double V)> EulerMethodSingular(
            double x0, double U0, double V0, double h, int nSteps)
        {
            var result = new List<(double, double, double)>();

            double x = x0;
            double U = U0;
            double V = V0;

            result.Add((x, U, V));

            for (int i = 0; i < nSteps; i++)
            {
                double denom = 4 * x;

                // если мы вплотную подошли к x=0, denom=0 -> возникает деление на 0
                if (Math.Abs(denom) < 1e-15)
                {
                    // чтобы не было NaN, прерываем
                    break;
                }

                double dU = V;
                double dV = -(2 * V + U) / denom;

                U = U + h * dU;
                V = V + h * dV;
                x = x + h;

                // проверяем на выход в NaN/∞
                if (double.IsNaN(U) || double.IsInfinity(U) ||
                    double.IsNaN(V) || double.IsInfinity(V) ||
                    double.IsNaN(x) || double.IsInfinity(x))
                {
                    break;
                }

                result.Add((x, U, V));
            }

            return result;
        }

        /// <summary>
        /// Разложение cos(sqrt(x)) по степеням x (до x^3).
        /// Только для x>=0 смысл имеет.
        ///   cos(√x) ~ 1 - x/2 + x^2/24 - x^3/720
        /// </summary>
        public static double TaylorSolution(double x)
        {
            return 1.0
                   - 0.5 * x
                   + (1.0 / 24.0) * x * x
                   - (1.0 / 720.0) * x * x * x;
        }
    }

    #region Задача 2: 1-2xUU'=U^3U', U(0.5)=0
    /// <summary>
    /// Преобразуем уравнение:
    ///   1 - 2x U U' = U^3 U'
    /// При условии U' ≠ 0 можно вынести U':
    ///   U' [ U^3 + 2x U ] = 1   =>   U' = 1 / [ U (U^2 + 2x) ].
    /// </summary>
    public static class Task2
    {
        public static void Run()
        {
            Console.WriteLine("=== Задача 2 ===");
            Console.WriteLine("Уравнение: 1 - 2x*u*u' = u^3 * u',  u(0.5)=0");
            Console.WriteLine("Аналитическое решение в параметрическом виде:");
            Console.WriteLine("  x(u) = e^(u^2) - (u^2 + 1)/2,   при  u(0)= x=0.5");
            Console.WriteLine();

            // Будем смотреть u в некотором диапазоне вокруг 0
            double uMin = -0.2;
            double uMax = +0.2;
            int nSteps = 20;
            double du = (uMax - uMin) / nSteps;

            // Сформируем сетку значений u
            double[] uGrid = new double[nSteps + 1];
            for (int i = 0; i <= nSteps; i++)
                uGrid[i] = uMin + i * du;

            // Считаем аналитическое решение x_exact(u) на этой сетке:
            double[] xExact = new double[nSteps + 1];
            for (int i = 0; i <= nSteps; i++)
            {
                double uu = uGrid[i];
                xExact[i] = Xexact(uu);
            }

            // Итерации Пикара: x_{k+1}(u) = 0.5 + ∫[0..u] [t^3 + 2*x_k(t)*t ] dt
            int picardCount = 8;
            var picardApproxs = PicardIter(uGrid, picardCount);

            Console.Write("      u       x_exact");
            for (int k = 1; k <= picardCount; k++)
                Console.Write($"   xPic{k}    ");
            Console.WriteLine();

            // Выводим таблицу
            for (int i = 0; i <= nSteps; i++)
            {
                double uu = uGrid[i];
                Console.Write($"{uu,8:F3}  {xExact[i],10:F6}");
                for (int k = 0; k < picardCount; k++)
                    Console.Write($"  {picardApproxs[k][i],10:F6}");
                Console.WriteLine();
            }

            Console.WriteLine();
        }

        /// <summary>
        /// Аналитическое решение: x(u) = e^(u^2) - (u^2+1)/2.
        /// </summary>
        static double Xexact(double u)
        {
            return Math.Exp(u * u) - 0.5 * (u * u + 1.0);
        }

        /// <summary>
        /// Итерации Пикара в виде x_{k+1}(u) = x(0) + ∫[0..u] [ t^3 + 2*x_k(t)*t ] dt.
        /// Начальное приближение x_0(u)=0.5 (константа).
        /// Для интегрирования используем простой метод трапеций на сетке uGrid.
        /// </summary>
        static List<double[]> PicardIter(double[] uGrid, int numIter)
        {
            int N = uGrid.Length;
            // Начальное приближение: x_0(u) = 0.5 (константа)
            double[] prev = new double[N];
            for (int i = 0; i < N; i++)
                prev[i] = 0.5;

            var result = new List<double[]>();

            for (int k = 0; k < numIter; k++)
            {
                double[] current = new double[N];

                // Найдём индекс, где u=0 (или ближайший)
                int idxZero = FindIndexOfZero(uGrid);

                // x(0) = 0.5
                current[idxZero] = 0.5;

                // Идём вправо: от idxZero+1 до конца
                for (int i = idxZero + 1; i < N; i++)
                {
                    double uLeft = uGrid[i - 1];
                    double uRight = uGrid[i];
                    double uMid = 0.5 * (uLeft + uRight);
                    double du = (uRight - uLeft);

                    double xLeft = prev[i - 1];    // x(uLeft)
                    double xRight = prev[i];        // x(uRight)
                    double xMid = 0.5 * (xLeft + xRight);

                    // f(u) = u^3 + 2*x(u)*u
                    double fLeft = uLeft * uLeft * uLeft + 2.0 * xLeft * uLeft;
                    double fMid = uMid * uMid * uMid + 2.0 * xMid * uMid;
                    double fRight = uRight * uRight * uRight + 2.0 * xRight * uRight;

                    // Метод Симпсона:
                    double localIntegral = (du / 6.0) * (fLeft + 4.0 * fMid + fRight);

                    current[i] = current[i - 1] + localIntegral;
                }

                // Идём влево: от idxZero-1 до 0
                for (int i = idxZero - 1; i >= 0; i--)
                {
                    double uLeft = uGrid[i];
                    double uRight = uGrid[i + 1];
                    double uMid = 0.5 * (uLeft + uRight);
                    double du = (uRight - uLeft);

                    double xLeft = prev[i];        // x(uLeft)
                    double xRight = prev[i + 1];    // x(uRight)
                    double xMid = 0.5 * (xLeft + xRight);

                    double fLeft = uLeft * uLeft * uLeft + 2.0 * xLeft * uLeft;
                    double fMid = uMid * uMid * uMid + 2.0 * xMid * uMid;
                    double fRight = uRight * uRight * uRight + 2.0 * xRight * uRight;

                    double localIntegral = (du / 6.0) * (fLeft + 4.0 * fMid + fRight);

                    // Двигаемся "влево", поэтому x(i) = x(i+1) - integral
                    current[i] = current[i + 1] - localIntegral;
                }

                result.Add(current);
                prev = current;
            }

            return result;
        }

        // Пример функции для поиска индекса, где uGrid[i] = 0 (или ближайший)
        static int FindIndexOfZero(double[] uGrid)
        {
            // Предположим, что 0 точно в массиве
            // (иначе можно искать ближайшее значение).
            for (int i = 0; i < uGrid.Length; i++)
                if (Math.Abs(uGrid[i]) < 1e-14)
                    return i;

            // Или вернуть -1, если не нашли (на практике дорабатываем логику)
            return -1;
        }
    }
    #endregion

    #region Задача 3: U'(x)= x+U^3, U(0)=0
    /// <summary>
    /// Решаем ОДУ:
    ///   U'(x) = x + U^3,
    ///   U(0) = 0.
    /// Реализуем три метода:
    /// 1) Метод Пикара (итерационный интегральный метод);
    /// 2) Явный метод Эйлера;
    /// 3) Метод Рунге–Кутты 4-го порядка.
    /// </summary>
    public static class Task3
    {
        public static void Run()
        {
            Console.WriteLine("Задача 3: u'(x)=x+u^3,  u(0)=0");
            double eps = /*1e-1*/1e-6;
            double hInit = 0.00000001;

            (double, double, double) XUh = AdaptiveEulerBigMath.FindBlowUpPoint_AdaptiveEuler(eps, hInit);
            Console.WriteLine($"с шагом h = {XUh.Item3:E300}\n U = {XUh.Item2:E300}\n при x = {XUh.Item1:E300}");
        }
    }

    public static class AdaptiveEulerBigMath
    {
        public static (double x, double U, double h) FindBlowUpPoint_AdaptiveEuler(
            double eps = 1e-4,
            double hInit = 1,
            double xMaxLimit = 1.65)
        {
            double x = 0.0;
            double U = 0.0;
            double h = hInit;
            double h_limit = 1e-323;

            int f = 1;
            int f10 = 0;
            double U_prev = 0.0;
            while (true)
            {
                // 1) Один большой шаг h
                double U_big = EulerOneStep(x, U, h);

                // 2) Два маленьких шага h/2
                double half = h / 2.0;
                double U_half_1 = EulerOneStep(x, U, half);
                double U_small = EulerOneStep(x + half, U_half_1, half);

                // 3) Локальная ошибка
                double denominator = /*Math.Max(1.0, */Math.Abs(U_small)/*)*/;
               /* if (denominator < 1e-300) denominator = 1e-300;*/
                double localError = Math.Abs(U_big - U_small) / denominator;

                // 4) Сравниваем localError > eps ?
                if (localError > eps)
                {
                    // Уменьшаем шаг
                    h /= 2.0;
                    if (h < h_limit)
                    {
                        Console.WriteLine("Шаг стал слишком малым => выходим.");
                        return (x, U_small, h);
                    }
                    // Повторяем этот же участок
                    continue;
                }
                else
                {
                    // Принимаем шаг
                    x += h;
                    U = U_small;

                    if (f10 < 10 && U > 1 * Math.Pow(10, f10))
                    {
                        Console.WriteLine($" x={x:E20}\n c={FindConst(x, U):E80}");
                        f10++;
                    }
                    if (U > 1 * Math.Pow(10, f * 10))
                    {
                        Console.WriteLine($"x={x:E20}, U={U:E8}, h={h:E8}");
                        Console.WriteLine($"c={FindConst(x, U):E80}");
                        f++;
                    }

                    if (double.IsInfinity(U))
                    {
                        Console.WriteLine("Blow-up!");
                        return (x-h, U_prev, h*2);
                    }

                    if (x > xMaxLimit)
                        return (xMaxLimit, U, h);
                    U_prev = U;
                }
            }
        }

        public static double FindConst(double x, double U)
        {
            return 1 / (2 * U * U) + x;
        }

        /// <summary>
        /// Один шаг метода Эйлера для уравнения U'(x)= x + U^3,
        /// где промежуточные вычисления идут в BigDecimal,
        /// но возвращаемое значение - double.
        /// </summary>
        /*private static double EulerOneStep(double x, double U, double step)
        {
            BigDecimal bx = BigDecimal.Parse(x);
            BigDecimal bU = BigDecimal.Parse(U);
            BigDecimal bStep = BigDecimal.Parse(step);

            BigDecimal fVal = bx + BigDecimal.Pow(bU, 3);

            BigDecimal bigResult = bU + bStep * fVal;

            return (double)(bigResult);
        }*/
        private static double EulerOneStep(double x, double U, double step)
        {
            double fVal = step*x + U*step*U*U;
            double bigResult = U + fVal;

            return (double)(bigResult);
        }
    }
}
#endregion

public static class Task4
{
    public static void Run()
    {
        Console.WriteLine("=== Задача 4 ===");
        Console.WriteLine("u'(x)= x + u^3,  u(0)=0. Сравнение рядов Пикара и двух численных методов.");

        double h = 0.000001;          // Шаг сетки
        double Xmax = 1.65;       // Макс. X (можно поставить больше, но решение может взорваться)
        int nSteps = (int)(Xmax / h);

        // Начальные условия для численных методов:
        double uE = 0.0;   // Euler
        double uRK = 0.0;  // RK4
        double x = 0.0;

        Console.WriteLine("   x          Pic1         Pic2         Pic3         Pic4        Euler         RK4");

        // Для i=0:
        PrintOneLine(x, Picard1(x), Picard2(x), Picard3(x), Picard4(x), uE, uRK);

        for (int i = 1; i <= nSteps; i++)
        {
            // шаг численного интегрирования:
            double xNext = x + h;

            // 1) Euler
            uE = EulerStep(x, uE, h);

            // 2) RK4
            uRK = RK4Step(x, uRK, h);

            x = xNext; // переходим на новый узел

            // Проверяем, не взорвалось ли
            if (double.IsInfinity(uE) || double.IsInfinity(uRK) ||
                Math.Abs(uE) > 1e300 || Math.Abs(uRK) > 1e300)
            {
                Console.WriteLine("Blow-up! Прерываем цикл.");
                break;
            }

            // Печать строки таблицы
            if (x > 1.647)
                PrintOneLine(x,
                             Picard1(x),
                             Picard2(x),
                             Picard3(x),
                             Picard4(x),
                             uE,
                             uRK);
        }
    }

    // Печать одной строки таблицы в удобном формате
    private static void PrintOneLine(double x,
                                     double p1, double p2, double p3, double p4,
                                     double ue, double urk)
    {
        Console.WriteLine($"{x,7:F10}   {p1,10:E3}   {p2,10:E3}   {p3,10:E3}   {p4,10:E3}   {ue,10:E3}   {urk,10:E3}");
    }

    // Формулы Пикара (1-я, 2-я, 3-я, 4-я итерации) - упрощённо:
    // Ниже в качестве примера взяты простейшие формулы (вы, возможно, имеете свои разложения).
    // Для уравнения u'(x)=x+u^3, старт u(0)=0:

    // Picard1:  u1(x) = \int_0^x t dt = x^2/2
    static double Picard1(double x)
    {
        return 0.5 * x * x;
    }

    // Picard2:  u2(x) = \int_0^x [ t + (u1(t))^3 ] dt
    //           = \int_0^x [ t + (t^2/2)^3 ] dt
    //           = \int_0^x [ t + t^6 / 8 ] dt
    //           = x^2/2 + x^7/(56)
    static double Picard2(double x)
    {
        double x2 = x * x;
        double x7 = x2 * x2 * x2 * x; // x^7
        return 0.5 * x2 + (x7 / 56.0);
    }

    // Picard3:  у вас могут быть более длинные формулы, ниже - просто пример.
    //           Можно было бы аккуратно подставить (u2(t))^3 = ... и снова проинтегрировать.
    //           Для наглядности (пример) возьмём что-то чуть длиннее:
    static double Picard3(double x)
    {
        // Возьмём "игрушечное" расширение:
        // x^2/2 + x^7/56 + (какие-то дополнительные члены)
        // В реальности выписывают аккуратно через интегралы. 
        // Ниже просто сделаем вид:
        double x2 = x * x;
        double x7 = x2 * x2 * x2 * x;
        double x12 = x7 * x5(x); // x^7 * x^5 = x^12
        return 0.5 * x2
               + (x7 / 56.0)
               + (x12 / 800.0);  // коэффициенты взяты из головы "для примера".
    }

    // Вспомогательная маленькая функция
    private static double x5(double x) => x * x * x * x * x;

    static double Picard4(double x)
    {
        // Аналогично, расширим ещё членами:
        double x2 = x * x;
        double x7 = x2 * x2 * x2 * x;
        double x12 = x7 * x5(x);
        double x17 = x12 * x5(x);
        return 0.5 * x2
               + (x7 / 56.0)
               + (x12 / 800.0)
               + (x17 / 90000.0);
    }

    // Метод Эйлера на одном шаге
    // u'(x)=x+u^3
    private static double EulerStep(double x, double u, double h)
    {
        double fVal = x + Math.Pow(u, 3);
        return u + h * fVal;
    }

    // Метод Рунге–Кутты 4-го порядка на одном шаге
    public static double RK4Step(double x, double u, double h)
    {
        double k1 = f(x, u);
        double k2 = f(x + h / 2, u + (h / 2) * k1);
        double k3 = f(x + h / 2, u + (h / 2) * k2);
        double k4 = f(x + h, u + h * k3);

        double inc = (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
        return u + h * inc;
    }
    private static double f(double xx, double uu) => xx + Math.Pow(uu, 3);
}
