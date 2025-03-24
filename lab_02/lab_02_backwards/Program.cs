using System;
using System.Collections.Generic;
//using static System.Net.Mime.MediaTypeNames;
using ScottPlot;

namespace Lab02_SystemUF
{
    public class Equation
    {
        private static double[] Ttable = { 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000 };
        private static double[] KtableVariant1 = {
            8.2e-03, 2.768e-02, 6.56e-02, 1.281e-01,
            2.214e-01, 3.516e-01, 5.248e-01, 7.472e-01, 1.025e+00
        };
        private static double[] KtableVariant2 = { 1.6e+00, 5.4e+00, 1.28e+01, 2.5e+01, 4.32e+01, 6.86e+01, 1.024e+02, 1.458e+02, 2.0e+02 };
        private static double[] Ktable = KtableVariant1;

        public double R = 0.35;
        public double T0 = 10000;
        public double Tw = 2000;
        public double pExp = 4.0;
        public double c = 3.0e10;
        public double C1 = 3.0084e-4;
        public double C2 = 4.799e4;
        public double scaleF = 1.0;

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
        public double dU(double r, double F) => -(3.0 * k_of_T(T_of_r(r)) / c) * F;
        public double dF(double r, double U, double F) => scaleF * 3.0 * k_of_T(T_of_r(r)) * R * (u_p(r) - U);

        // Производная u_p по r
        public double du_p_dr(double r)
        {
            double T = T_of_r(r);
            double dT_dr = pExp / R * (Tw - T0) * Math.Pow(r / R, pExp - 1);
            double expTerm = Math.Exp(C2 / T);
            double du_p_dT = C1 * (C2 / (T * T)) * expTerm / Math.Pow(expTerm - 1, 2);
            return du_p_dT * dT_dr;
        }
    }

    public class SolverBackward
    {
        private Equation eq;
        public SolverBackward(Equation eq) => this.eq = eq;

        public (double[] rOut, double[] U, double[] F) AdaptiveRungeKutta4_Backward(double xi, double tolLocal = 1e-6)
        {
            List<double> rList = new List<double>();
            List<double> UList = new List<double>();
            List<double> FList = new List<double>();

            double kAtR = eq.k_of_T(eq.T_of_r(eq.R));
            double U0 = -(3.0 / 0.39) * kAtR * xi;
            double F0 = xi;
            double currentR = eq.R, currentU = U0, currentF = F0;

            rList.Add(currentR); UList.Add(currentU); FList.Add(currentF);

            double h = -0.001, h_min = -1e-120, h_max = -0.1, safety = 0.9;

            while (currentR > 0)
            {
                if (currentR + h < 0) h = -currentR;

                var (U_big, F_big) = RK4Step(currentR, currentU, currentF, h);
                var (U_mid, F_mid) = RK4Step(currentR, currentU, currentF, h / 2.0);
                var (U_small, F_small) = RK4Step(currentR + h / 2.0, U_mid, F_mid, h / 2.0);

                double errU = Math.Abs(U_small - U_big), errF = Math.Abs(F_small - F_big);
                double relErr = Math.Max(errU / Math.Max(Math.Abs(U_small), 1e-10), errF / Math.Max(Math.Abs(F_small), 1e-10));

                if (relErr < tolLocal)
                {
                    currentR += h; currentU = U_small; currentF = F_small;
                    rList.Add(currentR); UList.Add(currentU); FList.Add(currentF);
                    h = Math.Max(h * safety * Math.Pow(tolLocal / relErr, 0.2), h_max);
                }
                else
                {
                    h = Math.Min(h * safety * Math.Pow(tolLocal / relErr, 0.2), h_min);
                    if (Math.Abs(h) < Math.Abs(h_min)) break;
                }
            }

            rList.Reverse(); UList.Reverse(); FList.Reverse();
            return (rList.ToArray(), UList.ToArray(), FList.ToArray());
        }

        private (double U_next, double F_next) RK4Step(double r, double U, double F, double h)
        {
            double k1 = eq.dU(r, F), l1 = eq.dF(r, U, F);
            double r2 = r + h / 2.0, U2 = U + h / 2.0 * k1, F2 = F + h / 2.0 * l1;
            double k2 = eq.dU(r2, F2), l2 = eq.dF(r2, U2, F2);
            double U3 = U + h / 2.0 * k2, F3 = F + h / 2.0 * l2;
            double k3 = eq.dU(r2, F3), l3 = eq.dF(r2, U3, F3);
            double r4 = r + h, U4 = U + h * k3, F4 = F + h * l3;
            double k4 = eq.dU(r4, F4), l4 = eq.dF(r4, U4, F4);

            return (U + h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4), F + h / 6.0 * (l1 + 2 * l2 + 2 * l3 + l4));
        }

        public double Psi_Backward(double xi)
        {
            var (_, _, FArray) = AdaptiveRungeKutta4_Backward(xi);
            return FArray[0];
        }

        public double FindXi(double lowerBound, double upperBound, double tol, int maxIter)
        {
            double psiLower = Psi_Backward(lowerBound), psiUpper = Psi_Backward(upperBound);
            Console.WriteLine($"Psi({lowerBound}) = {psiLower}, Psi({upperBound}) = {psiUpper}");
            if (psiLower * psiUpper > 0)
            {
                Console.WriteLine("Нет смены знака. Расширяем диапазон или используем другой подход.");
                return double.NaN;
            }
            for (int i = 0; i < maxIter; i++)
            {
                double mid = 0.5 * (lowerBound + upperBound);
                double psiMid = Psi_Backward(mid);
                if (Math.Abs(psiMid) < tol) return mid;
                if (psiLower * psiMid < 0) { upperBound = mid; psiUpper = psiMid; }
                else { lowerBound = mid; psiLower = psiMid; }
            }
            return 0.5 * (lowerBound + upperBound);
        }
    }

    public class Solver
    {
        private Equation eq;
        public Solver(Equation eq) => this.eq = eq;

        public (double[] rOut, double[] U, double[] F) AdaptiveRungeKutta4(double chi, double tolLocal = 1e-6, double alpha = 0.5, int maxIter = 10, double eps = 0.1)
        {
            // Первый прогон: вычисляем F(r), чтобы найти точку экстремума
            var (rOut, U_temp, F) = RungeKuttaStep(chi, tolLocal);

            // Находим точку экстремума F(r)
            int extremumIndex = 0;
            double maxF = double.MinValue;
            for (int i = 0; i < F.Length; i++)
            {
                if (F[i] > maxF)
                {
                    maxF = F[i];
                    extremumIndex = i;
                }
            }
            double rExtremum = rOut[extremumIndex];
            Console.WriteLine($"Точка экстремума F(r): r = {rExtremum}, F = {maxF}");

            // Полный прогон методом Рунге-Кутты
            var (rOutFinal, UFinal, FFinal) = RungeKuttaStep(chi, tolLocal);

            // Итеративная корректировка: до точки экстремума U(r) приближаем к u_p(r)
            for (int iter = 0; iter < maxIter; iter++)
            {
                bool converged = true;
                double[] U_new = new double[UFinal.Length];
                double[] F_new = new double[FFinal.Length];

                for (int i = 0; i < rOutFinal.Length; i++)
                {
                    double r = rOutFinal[i];
                    double up = eq.u_p(r);
                    double k = eq.k_of_T(eq.T_of_r(r));
                    double du_p_dr = eq.du_p_dr(r);

                    if (r <= rExtremum)
                    {
                        // До точки экстремума: корректируем U(r), чтобы оно было близко к u_p(r)
                        U_new[i] = alpha * UFinal[i] + (1 - alpha) * up;
                    }
                    else
                    {
                        // После точки экстремума: оставляем U(r) как есть (вычислено Рунге-Куттой)
                        U_new[i] = alpha * UFinal[i] + (1 - alpha) * up;
                    }

                    // Корректируем F(r)
                    F_new[i] = alpha * FFinal[i] + (1 - alpha) * (-eq.c / (3 * k)) * du_p_dr;

                    // Проверка сходимости
                    if (Math.Abs(U_new[i] - UFinal[i]) > eps || Math.Abs(F_new[i] - FFinal[i]) > eps)
                        converged = false;
                }

                UFinal = U_new;
                FFinal = F_new;
                if (converged) break;
            }

            return (rOutFinal, UFinal, FFinal);
        }

        private (double[] rOut, double[] U, double[] F) RungeKuttaStep(double chi, double tolLocal)
        {
            List<double> rList = new List<double>();
            List<double> UList = new List<double>();
            List<double> FList = new List<double>();

            double currentR = 0, currentU = chi * eq.u_p(0), currentF = 0;
            rList.Add(currentR); UList.Add(currentU); FList.Add(currentF);

            double h = 0.001, h_min = 1e-120, h_max = 0.1, safety = 0.9;

            while (currentR < eq.R)
            {
                if (currentR + h > eq.R) h = eq.R - currentR;

                var (U_big, F_big) = RK4Step(currentR, currentU, currentF, h);
                var (U_mid, F_mid) = RK4Step(currentR, currentU, currentF, h / 2.0);
                var (U_small, F_small) = RK4Step(currentR + h / 2.0, U_mid, F_mid, h / 2.0);

                double errU = Math.Abs(U_small - U_big), errF = Math.Abs(F_small - F_big);
                double relErr = Math.Max(errU / Math.Max(Math.Abs(U_small), 1e-10), errF / Math.Max(Math.Abs(F_small), 1e-10));

                if (relErr < tolLocal)
                {
                    currentR += h; currentU = U_small; currentF = F_small;
                    rList.Add(currentR); UList.Add(currentU); FList.Add(currentF);
                    h = Math.Min(h * safety * Math.Pow(tolLocal / relErr, 0.2), h_max);
                }
                else
                {
                    h = Math.Max(h * safety * Math.Pow(tolLocal / relErr, 0.2), h_min);
                    if (h < h_min) break;
                }
            }
            return (rList.ToArray(), UList.ToArray(), FList.ToArray());
        }

        private (double U_next, double F_next) RK4Step(double r, double U, double F, double h)
        {
            double k1 = eq.dU(r, F), l1 = eq.dF(r, U, F);
            double r2 = r + h / 2.0, U2 = U + h / 2.0 * k1, F2 = F + h / 2.0 * l1;
            double k2 = eq.dU(r2, F2), l2 = eq.dF(r2, U2, F2);
            double U3 = U + h / 2.0 * k2, F3 = F + h / 2.0 * l2;
            double k3 = eq.dU(r2, F3), l3 = eq.dF(r2, U3, F3);
            double r4 = r + h, U4 = U + h * k3, F4 = F + h * l3;
            double k4 = eq.dU(r4, F4), l4 = eq.dF(r4, U4, F4);

            return (U + h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4), F + h / 6.0 * (l1 + 2 * l2 + 2 * l3 + l4));
        }
    }

    public static class Program
    {
        [STAThread]
        public static void Main()
        {
            Equation eq = new Equation();
            SolverBackward solverBackward = new SolverBackward(eq);

            List<double> xiValues = new List<double>();
            Dictionary<double, (double[] r, double[] U, double[] F)> solutions = new Dictionary<double, (double[], double[], double[])>();

            for (double xi = 1e-12; xi <= 1e12; xi *= 5)
            {
                var (rOut, Usol, Fsol) = solverBackward.AdaptiveRungeKutta4_Backward(xi);
                xiValues.Add(xi);
                solutions[xi] = (rOut, Usol, Fsol);
                Console.WriteLine($"xi = {xi}, F(0) = {Fsol[0]}, U(0) = {Usol[0]}, u_p(0) = {eq.u_p(0)}");
            }

            // Задаем диапазон z (от 0 до 1)
            double zLeft = 0.0;  // Левая граница z (например, 0.0)
            double zRight = 1.0; // Правая граница z (например, 1.0, соответствует R = 0.35)

            // Преобразуем z в r с учетом R
            double rLeft = zLeft * eq.R;
            double rRight = zRight * eq.R;

            // Проверка корректности введенных значений
            if (zLeft < 0 || zRight > 1 || zLeft > zRight)
            {
                Console.WriteLine("Ошибка: zLeft и zRight должны быть в диапазоне [0, 1], и zLeft <= zRight.");
                return;
            }

            var pltU = new ScottPlot.Plot();
            foreach (var xi in xiValues)
            {
                var (r, U, _) = solutions[xi];
                // Фильтруем данные в диапазоне [rLeft, rRight]
                var rFiltered = r.Where(x => x >= rLeft && x <= rRight).ToArray();
                var uFiltered = U.SkipWhile((_, i) => r[i] < rLeft)
                                .Take(rFiltered.Length).ToArray();
                pltU.AddScatter(rFiltered, uFiltered, label: $"xi = {xi:E2}");
            }
            pltU.Title("U(r) для разных xi (справа-налево)");
            pltU.XLabel("r");
            pltU.YLabel("U");
            pltU.Legend();
            var formU = new FormsPlotViewer(pltU) { Text = "U(r) vs xi" };
            formU.Show();

            var pltF = new ScottPlot.Plot();
            foreach (var xi in xiValues)
            {
                var (r, _, F) = solutions[xi];
                // Фильтруем данные в диапазоне [rLeft, rRight]
                var rFiltered = r.Where(x => x >= rLeft && x <= rRight).ToArray();
                var fFiltered = F.SkipWhile((_, i) => r[i] < rLeft)
                                .Take(rFiltered.Length).ToArray();
                pltF.AddScatter(rFiltered, fFiltered, label: $"xi = {xi:E2}");
            }
            pltF.Title("F(r) для разных xi (справа-налево)");
            pltF.XLabel("r");
            pltF.YLabel("F");
            pltF.Legend();
            var formF = new FormsPlotViewer(pltF) { Text = "F(r) vs xi" };
            formF.Show();

            double xiOpt = solverBackward.FindXi(-1e-6, 1e+12, 1e-9, 10000);
            //xiOpt = 340;
            if (!double.IsNaN(xiOpt))
            {
                var (rOpt, UOpt, FOpt) = solverBackward.AdaptiveRungeKutta4_Backward(xiOpt);
                double chi = UOpt[0] / eq.u_p(0);
                Console.WriteLine($"Оптимальное xi = {xiOpt}, chi = {chi}");

                Solver solver = new Solver(eq);
                var (rForw, UForw, FForw) = solver.AdaptiveRungeKutta4(chi, alpha: 0.5, maxIter: 10, eps: 1e-4);
                double error = FForw.Last() - 0.39 * eq.c * UForw.Last();
                Console.WriteLine($"F(R) - 0.39*c*U(R) = {error}");

                // Фильтруем данные в диапазоне [rLeft, rRight] для всех последующих графиков
                var rForwFiltered = rForw.Where(x => x >= rLeft && x <= rRight).ToArray();
                var uForwFiltered = UForw.SkipWhile((_, i) => rForw[i] < rLeft)
                                        .Take(rForwFiltered.Length).ToArray();
                var fForwFiltered = FForw.SkipWhile((_, i) => rForw[i] < rLeft)
                                        .Take(rForwFiltered.Length).ToArray();

                Console.WriteLine($"r = {rForwFiltered[rForwFiltered.Length-1]}; F = {fForwFiltered[fForwFiltered.Length-1]}");
                // График U и F вместе
                var pltSol = new ScottPlot.Plot();
                pltSol.AddScatter(rForwFiltered, uForwFiltered, label: "U(r)");
                pltSol.AddScatter(rForwFiltered, fForwFiltered, label: "F(r)");
                pltSol.Title("Решение слева-направо: U(r) и F(r)");
                pltSol.XLabel("r");
                pltSol.YLabel("Value");
                pltSol.Legend();
                var formSol = new FormsPlotViewer(pltSol) { Text = "U and F" };
                formSol.Show();

                // Отдельный график F
                Console.WriteLine(fForwFiltered[fForwFiltered.Length - 1]);
                var pltFOnly = new ScottPlot.Plot();
                pltFOnly.AddScatter(rForwFiltered, fForwFiltered, label: "F(r)");
                pltFOnly.Title("F(r) слева-направо");
                pltFOnly.XLabel("r");
                pltFOnly.YLabel("F");
                pltFOnly.Legend();
                var formFOnly = new FormsPlotViewer(pltFOnly) { Text = "F(r) Only" };
                formFOnly.Show();

                // График u_p и u
                var pltUpU = new ScottPlot.Plot();
                pltUpU.AddScatter(rForwFiltered, uForwFiltered, label: "U(r)");
                double[] upValues = new double[rForwFiltered.Length];
                for (int i = 0; i < rForwFiltered.Length; i++) upValues[i] = eq.u_p(rForwFiltered[i]);
                pltUpU.AddScatter(rForwFiltered, upValues, label: "u_p(r)");
                pltUpU.Title("Сравнение U(r) и u_p(r)");
                pltUpU.XLabel("r");
                pltUpU.YLabel("Value");
                pltUpU.Legend();
                var formUpU = new FormsPlotViewer(pltUpU) { Text = "U vs u_p" };
                formUpU.Show();
            }
            else
            {
                Console.WriteLine("Не удалось найти xi, удовлетворяющее F(0) = 0.");
            }

            Application.Run();
        }
    }
}