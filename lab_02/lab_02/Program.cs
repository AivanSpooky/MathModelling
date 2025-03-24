using System;
using System.Collections.Generic;
using System.Windows.Forms;
using ScottPlot;

namespace Lab02_SystemUF
{
    /*
     * ОБОЗНАЧЕНИЯ ПЕРЕМЕННЫХ И КОНСТАНТ:
     *
     * Ttable         - Массив температур (в Кельвинах), задающих точки таблицы.
     * KtableVariant1 - Массив коэффициентов поглощения для варианта 1.
     * KtableVariant2 - Массив коэффициентов поглощения для варианта 2.
     * Ktable         - Выбранный массив коэффициентов поглощения.
     *
     * R              - Радиус цилиндра (масштабная длина, используется для вычислений).
     * T0             - Температура в центре (при z=0).
     * Tw             - Температура на границе (при z=1).
     * pExp           - Показатель степени в формуле T(z)=T0+(Tw-T0)*z^(pExp).
     * c              - Скорость света.
     *
     * C1, C2         - Константы для функции Планка: uₚ(z)=C1/(exp(C2/T(z))-1).
     *
     * scaleF         - Дополнительный множитель для F'(z)=scaleF*3*k(T(z))*R*(uₚ(z)-U).
     *
     * U(z)           - Искомая функция (например, объёмная плотность излучения).
     * F(z)           - Поток излучения.
     * uₚ(z)          - Функция Планка для данной T(z).
     *
     * ksi            - Параметр граничного условия: U(0)=ksi*uₚ(0).
     *                  (В методе дихотомии подбирается такое ksi, чтобы на границе выполнялось
     *                   условие: F(1)=0.393*c*U(1), то есть Psi(ksi)=F(1)-0.393*c*U(1)=0.)
     *
     * Psi(ksi)       - Функция ошибки, по которой подбирается ksi.
     *
     * z              - Безразмерная переменная (изменяется от 0 до 1).
     *
     * slopes, intercepts - Коэффициенты для кусочно-линейной интерполяции в ln‑ln шкале для k(T).
     */

    public class Equation
    {
        private static double[] Ttable = { 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000 };
        private static double[] KtableVariant1 = {
            8.2e-03, 2.768e-02, 6.56e-02, 1.281e-01,
            2.214e-01, 3.516e-01, 5.248e-01, 7.472e-01, 1.025e+00
        };
        private static double[] KtableVariant2 = {
            1.6e+00, 5.4e+00, 1.28e+01, 2.5e+01,
            4.32e+01, 6.86e+01, 1.024e+02, 1.458e+02, 2.0e+02
        };
        // Здесь выбран вариант 1 (можно изменить на вариант 2 при необходимости)
        private static double[] Ktable = KtableVariant1;

        public double R = 0.35;
        public double T0 = 10000;
        public double Tw = 2000;
        public double pExp = 4.0;
        public double c = 3.0e10;
        public double C1 = 3.0084e-4;
        public double C2 = 4.799e4;
        public double scaleF = 3.3e+9;

        private double[] slopes;
        private double[] intercepts;

        public Equation()
        {
            // Здесь можно экспериментально корректировать scaleF
            scaleF = 2.5e+9;
            //scaleF = 1;
            ComputeInterpolationCoefficients();
        }

        private void ComputeInterpolationCoefficients()
        {
            int n = Ttable.Length;
            slopes = new double[n - 1];
            intercepts = new double[n - 1];
            for (int i = 0; i < n - 1; i++)
            {
                double T1 = Ttable[i];
                double T2 = Ttable[i + 1];
                double k1 = Ktable[i];
                double k2 = Ktable[i + 1];
                double x1 = Math.Log(T1);
                double x2 = Math.Log(T2);
                double y1 = Math.Log(k1);
                double y2 = Math.Log(k2);
                slopes[i] = (y2 - y1) / (x2 - x1);
                intercepts[i] = y1 - slopes[i] * x1;
            }
        }

        public double T_of_z(double z)
        {
            return T0 + (Tw - T0) * Math.Pow(z, pExp);
        }

        public double k_of_T(double T)
        {
            if (T <= Ttable[0])
                return Ktable[0];
            if (T >= Ttable[Ttable.Length - 1])
                return Ktable[Ktable.Length - 1];
            int i = 0;
            while (i < Ttable.Length - 1 && T > Ttable[i + 1])
                i++;
            double lx = Math.Log(T);
            double ly = intercepts[i] + slopes[i] * lx;
            return Math.Exp(ly);
        }

        public double u_p(double z)
        {
            double Tloc = T_of_z(z);
            return C1 / (Math.Exp(C2 / Tloc) - 1.0);
        }

        public double dU(double z, double F)
        {
            double kVal = k_of_T(T_of_z(z));
            return -(3.0 * kVal / c) * F;
        }

        public double dF(double z, double U, double F)
        {
            double kVal = k_of_T(T_of_z(z));
            double upVal = u_p(z);
            return scaleF * 3.0 * kVal * R * (upVal - U);
        }
    }

    public class Solver
    {
        private Equation eq;
        // Здесь интегрирование проводится адаптивно, поэтому равномерная сетка не задаётся заранее.
        private double z0, zMax;

        public Solver(Equation eq, double z0 = 0.0, double zMax = 1.0)
        {
            this.eq = eq;
            this.z0 = z0;
            this.zMax = zMax;
        }

        // Адаптивный метод Рунге–Кутта 4-го порядка с динамическим шагом.
        // Возвращает динамически сформированную сетку zOut, массивы U и F.
        public (double[] zOut, double[] U, double[] F) AdaptiveRungeKutta4(double ksi, double tolLocal = 1e-6)
        {
            List<double> zList = new List<double>();
            List<double> UList = new List<double>();
            List<double> FList = new List<double>();

            double currentZ = z0;
            double currentU = ksi * eq.u_p(z0);
            double currentF = 0.0;

            zList.Add(currentZ);
            UList.Add(currentU);
            FList.Add(currentF);

            double h = 0.01; // Начальный шаг
            double h_min = 1e-12;
            double h_max = 0.1;
            double safety = 0.9;

            while (currentZ < zMax)
            {
                if (currentZ + h > zMax)
                    h = zMax - currentZ;

                // Один шаг с шагом h
                var (U_big, F_big) = RK4Step(currentZ, currentU, currentF, h);
                // Два шага с шагом h/2
                var (U_mid, F_mid) = RK4Step(currentZ, currentU, currentF, h / 2.0);
                var (U_small, F_small) = RK4Step(currentZ + h / 2.0, U_mid, F_mid, h / 2.0);

                double errU = Math.Abs(U_small - U_big);
                double errF = Math.Abs(F_small - F_big);
                double normU = Math.Max(Math.Abs(U_small), 1e-10);
                double normF = Math.Max(Math.Abs(F_small), 1e-10);
                double relErr = Math.Max(errU / normU, errF / normF);

                if (relErr < tolLocal)
                {
                    // Шаг принят
                    currentZ += h;
                    currentU = U_small;
                    currentF = F_small;
                    zList.Add(currentZ);
                    UList.Add(currentU);
                    FList.Add(currentF);
                    // Увеличиваем шаг, если ошибка очень мала
                    double factor = safety * Math.Pow(tolLocal / relErr, 0.2);
                    h = Math.Min(h * factor, h_max);
                }
                else
                {
                    // Уменьшаем шаг
                    double factor = safety * Math.Pow(tolLocal / relErr, 0.2);
                    h = Math.Max(h * factor, h_min);
                    if (h <= h_min)
                    {
                        Console.WriteLine("Минимальный шаг достигнут. Прерывание интегрирования.");
                        break;
                    }
                }
            }
            Console.WriteLine($"h = {h}");
            return (zList.ToArray(), UList.ToArray(), FList.ToArray());
        }

        // Выполняет один шаг Рунге–Кутта 4-го порядка с шагом h.
        private (double U_next, double F_next) RK4Step(double zCurr, double U_curr, double F_curr, double h)
        {
            double k1 = eq.dU(zCurr, F_curr);
            double l1 = eq.dF(zCurr, U_curr, F_curr);

            double zHalf = zCurr + h / 2.0;
            double U_temp = U_curr + (h / 2.0) * k1;
            double F_temp = F_curr + (h / 2.0) * l1;
            double k2 = eq.dU(zHalf, F_temp);
            double l2 = eq.dF(zHalf, U_temp, F_temp);

            double U_temp2 = U_curr + (h / 2.0) * k2;
            double F_temp2 = F_curr + (h / 2.0) * l2;
            double k3 = eq.dU(zHalf, F_temp2);
            double l3 = eq.dF(zHalf, U_temp2, F_temp2);

            double zNext = zCurr + h;
            double U_temp3 = U_curr + h * k3;
            double F_temp3 = F_curr + h * l3;
            double k4 = eq.dU(zNext, F_temp3);
            double l4 = eq.dF(zNext, U_temp3, F_temp3);

            double U_next = U_curr + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
            double F_next = F_curr + (h / 6.0) * (l1 + 2 * l2 + 2 * l3 + l4);

            return (U_next, F_next);
        }

        // Функция ошибки: Psi(ksi)=F(zMax) - 0.393*c*U(zMax)
        public double Psi(double ksi)
        {
            var (zArray, UArray, FArray) = AdaptiveRungeKutta4(ksi);
            double U_last = UArray[UArray.Length - 1];
            double F_last = FArray[FArray.Length - 1];
            return F_last - 0.393 * eq.c * U_last;
        }

        /// <summary>
        /// Метод дихотомии для поиска ksi, при котором Psi(ksi)=0.
        /// Здесь используется адаптивный метод РК4 с динамическим шагом.
        /// </summary>
        public (double leftVal, double midVal, double rightVal, double[] zOut, double[] Usol, double[] Fsol)
            DichotomyMethod(double leftGuess, double rightGuess, double tol, int maxIter)
        {
            double psiLeft = Psi(leftGuess);
            double psiRight = Psi(rightGuess);
            if (psiLeft * psiRight > 0)
                Console.WriteLine("Warning: Psi(left)*Psi(right)>0 => нет гарантии корня!");

            double left = leftGuess;
            double right = rightGuess;
            double mid = 0, psiMid = 0;
            for (int i = 0; i < maxIter; i++)
            {
                mid = 0.5 * (left + right);
                psiMid = Psi(mid);
                Console.WriteLine($"i = {i}; left = {left}; mid = {mid}; right = {right}");
                Console.WriteLine($"psiLeft = {psiLeft}; psiMid = {psiMid}; psiRight = {psiRight}");
                if (Math.Abs(right - left) / Math.Abs(mid) < tol)
                    break;
                if (psiLeft * psiMid < 0)
                {
                    right = mid;
                    psiRight = psiMid;
                }
                else
                {
                    left = mid;
                    psiLeft = psiMid;
                }
            }
            var (zOut, Usol, Fsol) = AdaptiveRungeKutta4(mid);
            return (left, mid, right, zOut, Usol, Fsol);
        }
    }

    public static class Program
    {
        [STAThread]
        public static void Main()
        {
            Equation eq = new Equation();
            // Здесь решается задача с явным адаптивным методом РК4 с динамическим шагом.
            // Если бы мы решали неявным методом Эйлера или методом трапеций, то:
            //
            // – Неявный метод Эйлера обладает A‑устойчивостью и способен обрабатывать жёсткие задачи,
            //   но в каждом шаге требуется решение нелинейного уравнения (обычно методом Ньютона),
            //   что значительно увеличивает вычислительные затраты.
            //
            // – Метод трапеций (также неявный) имеет второй порядок точности и A‑устойчив,
            //   но он тоже требует решения системы нелинейных уравнений на каждом шаге,
            //   что может быть вычислительно затратным.
            //
            // Явный адаптивный метод Рунге–Кутта 4‑го порядка, как реализован здесь, 
            // часто оказывается эффективным для задач, не имеющих чрезмерной жёсткости.
            //
            // Решение неявными методами полезно, если система является жёсткой, и требуется 
            // высокая устойчивость, однако это приводит к необходимости использования итерационных
            // схем для решения алгебраических уравнений на каждом шаге.

            Solver solver = new Solver(eq, z0: 0.0, zMax: 1.0);

            // Задаём интервал поиска ksi (граничное условие: F(1)=0.393*c*U(1))
            double leftGuess = 0.0;
            double rightGuess = 1.0;
            var (leftVal, midVal, rightVal, zOut, Usol, Fsol) = solver.DichotomyMethod(leftGuess, rightGuess, tol: 1e-6, maxIter: 1000);

            Console.WriteLine($"ksi in [{leftVal}, {midVal}, {rightVal}]");

            // Для построения графиков используем динамически сформированную сетку zOut из интегратора.
            var pltU = new ScottPlot.Plot();
            pltU.AddScatter(zOut, Usol);
            pltU.Title("U(z)");
            pltU.XLabel("z");
            pltU.YLabel("U");
            pltU.YAxis.TickLabelFormat(x => x.ToString("E2"));
            var formU = new FormsPlotViewer(pltU) { Text = "U(z)" };
            formU.Show();

            var pltF = new ScottPlot.Plot();
            pltF.AddScatter(zOut, Fsol);
            pltF.Title("F(z)");
            pltF.XLabel("z");
            pltF.YLabel("F");
            pltF.YAxis.TickLabelFormat(x => x.ToString("E2"));
            var formF = new FormsPlotViewer(pltF) { Text = "F(z)" };
            formF.Show();

            Application.Run();
        }
    }
}
