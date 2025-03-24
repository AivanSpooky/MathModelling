using System;
using System.Collections.Generic;
using System.Windows.Forms;
using ScottPlot;

namespace Lab02_SystemUF_Boundary
{
    public class Equation
    {
        private static double[] Ttable = { 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000 };
        private static double[] KtableVariant2 = { 1.6e+00, 5.4e+00, 1.28e+01, 2.5e+01, 4.32e+01, 6.86e+01, 1.024e+02, 1.458e+02, 2.0e+02 };
        private static double[] Ktable = KtableVariant2;

        public double R = 0.35;
        public double T0 = 10000;
        public double Tw = 2000;
        public double pExp = 4.0;
        public double c = 3.0e10;
        public double C1 = 3.0084e-4;
        public double C2 = 4.799e4;

        private double[] slopes;
        private double[] intercepts;

        public Equation() => ComputeInterpolationCoefficients();

        private void ComputeInterpolationCoefficients()
        {
            int n = Ttable.Length;
            slopes = new double[n - 1];
            intercepts = new double[n - 1];
            for (int i = 0; i < n - 1; i++)
            {
                double T1 = Ttable[i], T2 = Ttable[i + 1];
                double k1 = Ktable[i], k2 = Ktable[i + 1];
                double x1 = Math.Log(T1), x2 = Math.Log(T2);
                double y1 = Math.Log(k1), y2 = Math.Log(k2);
                slopes[i] = (y2 - y1) / (x2 - x1);
                intercepts[i] = y1 - slopes[i] * x1;
            }
        }

        public double T_of_r(double r) => T0 + (Tw - T0) * Math.Pow(r / R, pExp);
        public double k_of_T(double T)
        {
            if (T <= Ttable[0]) return Ktable[0];
            if (T >= Ttable[Ttable.Length - 1]) return Ktable[Ttable.Length - 1];
            int i = 0;
            while (i < Ttable.Length - 1 && T > Ttable[i + 1]) i++;
            double lx = Math.Log(T);
            return Math.Exp(intercepts[i] + slopes[i] * lx);
        }
        public double u_p(double r) => C1 / (Math.Exp(C2 / T_of_r(r)) - 1.0);
    }

    public class BoundarySolver
    {
        private Equation eq;
        private int N; // Число узлов
        private double h; // Шаг сетки
        private double[] r; // Сетка по r
        private double[] U; // Решение U
        private double[] F; // Решение F

        public BoundarySolver(Equation eq, int N = 1000)
        {
            this.eq = eq;
            this.N = N;
            h = eq.R / (N - 1);
            r = new double[N];
            U = new double[N];
            F = new double[N];
            for (int i = 0; i < N; i++)
            {
                r[i] = i * h;
                U[i] = eq.u_p(r[i]); // Начальное приближение для U
                // Линейное начальное приближение для F
                F[i] = (r[i] / eq.R) * (0.39 * eq.c * eq.u_p(eq.R));
            }
        }

        public (double[] rOut, double[] U, double[] F) Solve(double tol = 1e-6, int maxIter = 100, double omega = 0.5)
        {
            for (int iter = 0; iter < maxIter; iter++)
            {
                double maxDelta = 0;
                double[] deltaU = new double[N];
                double[] deltaF = new double[N];

                // Формируем систему уравнений методом Ньютона
                for (int i = 1; i < N - 1; i++)
                {
                    double k = eq.k_of_T(eq.T_of_r(r[i]));
                    double up = eq.u_p(r[i]);

                    // Уравнение для U: (U_{i+1} - U_{i-1}) / (2h) = - (3 k / c) F_i
                    double G1_inner = (U[i + 1] - U[i - 1]) / (2 * h) + (3 * k / eq.c) * F[i];
                    // Уравнение для F: (F_{i+1} - F_{i-1}) / (2h) = 3 k R (u_p - U_i)
                    double G2_inner = (F[i + 1] - F[i - 1]) / (2 * h) - 3 * k * eq.R * (up - U[i]);

                    // Якобиан
                    double dG1_dUi_m1 = -1 / (2 * h);
                    double dG1_dUi_p1 = 1 / (2 * h);
                    double dG1_dFi = 3 * k / eq.c;

                    double dG2_dUi = 3 * k * eq.R;
                    double dG2_dFi_m1 = -1 / (2 * h);
                    double dG2_dFi_p1 = 1 / (2 * h);

                    // Решаем систему 2x2 для deltaU[i] и deltaF[i]
                    double det_inner = dG1_dUi_p1 * dG2_dFi_p1 - dG1_dFi * dG2_dUi;
                    deltaU[i] = (-G1_inner * dG2_dFi_p1 + G2_inner * dG1_dFi) / det_inner;
                    deltaF[i] = (-dG1_dUi_p1 * G2_inner + dG2_dUi * G1_inner) / det_inner;

                    maxDelta = Math.Max(maxDelta, Math.Max(Math.Abs(deltaU[i]), Math.Abs(deltaF[i])));
                }

                // Граничные условия
                F[0] = 0; // F(0) = 0
                // Корректируем U[N-1] и F[N-1] с учётом граничного условия
                double kR = eq.k_of_T(eq.T_of_r(r[N - 1]));
                double upR = eq.u_p(r[N - 1]);
                double G1 = (U[N - 1] - U[N - 2]) / h + (3 * kR / eq.c) * F[N - 1]; // Односторонняя разность для U
                double G2 = F[N - 1] - 0.39 * eq.c * U[N - 1]; // Граничное условие F(R)

                double dG1_dUN_m1 = -1 / h;
                double dG1_dUN = 1 / h;
                double dG1_dFN = 3 * kR / eq.c;

                double dG2_dUN = -0.39 * eq.c;
                double dG2_dFN = 1;

                double det = dG1_dUN * dG2_dFN - dG1_dFN * dG2_dUN;
                deltaU[N - 1] = (-G1 * dG2_dFN + G2 * dG1_dFN) / det;
                deltaF[N - 1] = (-dG1_dUN * G2 + dG2_dUN * G1) / det;

                maxDelta = Math.Max(maxDelta, Math.Max(Math.Abs(deltaU[N - 1]), Math.Abs(deltaF[N - 1])));

                //Обновляем U и F с релаксацией
                for (int i = 1; i < N; i++)
                {
                    U[i] = (1 - omega) * U[i] + omega * (U[i] + deltaU[i]);
                    F[i] = (1 - omega) * F[i] + omega * (F[i] + deltaF[i]);
                }

                // Проверяем сходимость
                if (maxDelta < tol)
                {
                    Console.WriteLine($"Сошлось на итерации {iter}, maxDelta = {maxDelta}");
                    break;
                }
            }

            return (r, U, F);
        }
    }

    public static class Program
    {
        [STAThread]
        public static void Main()
        {
            Equation eq = new Equation();
            BoundarySolver solver = new BoundarySolver(eq, N: 1000);

            var (r, U, F) = solver.Solve(tol: 1e-8, maxIter: 1000, omega: 0.1);

            // График U и F вместе
            var pltSol = new ScottPlot.Plot();
            pltSol.AddScatter(r, U, label: "U(r)");
            pltSol.AddScatter(r, F, label: "F(r)");
            pltSol.Title("Решение методом конечных разностей: U(r) и F(r)");
            pltSol.XLabel("r");
            pltSol.YLabel("Value");
            pltSol.Legend();
            var formSol = new FormsPlotViewer(pltSol) { Text = "U and F" };
            formSol.Show();

            // Отдельный график F
            var pltFOnly = new ScottPlot.Plot();
            pltFOnly.AddScatter(r, F, label: "F(r)");
            pltFOnly.Title("F(r) методом конечных разностей");
            pltFOnly.XLabel("r");
            pltFOnly.YLabel("F");
            pltFOnly.Legend();
            var formFOnly = new FormsPlotViewer(pltFOnly) { Text = "F(r) Only" };
            formFOnly.Show();

            // График u_p и u
            var pltUpU = new ScottPlot.Plot();
            pltUpU.AddScatter(r, U, label: "U(r)");
            double[] upValues = new double[r.Length];
            for (int i = 0; i < r.Length; i++) upValues[i] = eq.u_p(r[i]);
            pltUpU.AddScatter(r, upValues, label: "u_p(r)");
            pltUpU.Title("Сравнение U(r) и u_p(r)");
            pltUpU.XLabel("r");
            pltUpU.YLabel("Value");
            pltUpU.Legend();
            var formUpU = new FormsPlotViewer(pltUpU) { Text = "U vs u_p" };
            formUpU.Show();

            Application.Run();
        }
    }
}