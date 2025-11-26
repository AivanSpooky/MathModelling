using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Windows.Forms;
using ScottPlot;

namespace Lab3Modeling
{
    public static class Task2
    {
        private struct RadiativeTransfer
        {
            public double C;  // см/с
            public double R;  // см
            public double Tw; // К
            public double T0; // К
            public double K0; // 1/см
            public double P;
            public double M;
        }
        //private static readonly double[] TValues = { 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000 };
        //private static readonly double[] KValues = { 8.200E-03, 2.768E-02, 6.560E-02, 1.281E-01, 2.214E-01, 3.516E-01, 5.248E-01, 7.472E-01, 1.025E+00 };
        private static readonly double[] TValues = { 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000 };
        private static readonly double[] KValues = { 1.600E+00, 5.400E+00, 1.280E+01, 2.500E+01, 4.320E+01, 6.860E+01, 1.024E+02, 1.458E+02, 2.000E+02 };
        private const double EPS = 1e-4;

        private static double LogInterp(double TNew, double[] tValues, double[] kValues, int numPoints)
        {
            double[] logT = tValues.Select(t => Math.Log(t)).ToArray();
            double[] logK = kValues.Select(k => Math.Log(k)).ToArray();
            double logTNew = Math.Log(TNew);
            double logKNew = 0.0;

            for (int i = 0; i < numPoints - 1; i++)
            {
                if (logTNew >= logT[i] && logTNew <= logT[i + 1])
                {
                    logKNew = logK[i] + (logK[i + 1] - logK[i]) * (logTNew - logT[i]) / (logT[i + 1] - logT[i]);
                    break;
                }
            }
            return Math.Exp(logKNew);
        }

        private static double T(ref RadiativeTransfer rt, double z)
        {
            return (rt.Tw - rt.T0) * Math.Pow(z, rt.P) + rt.T0;
        }

        private static double K(ref RadiativeTransfer rt, double z)
        {
            double tNew = T(ref rt, z);
            return LogInterp(tNew, TValues, KValues, TValues.Length);
        }

        private static double DerivativeU(ref RadiativeTransfer rt, double z, double F)
        {
            return -3 * rt.R * K(ref rt, z) / rt.C * F;
        }

        private static double UP(ref RadiativeTransfer rt, double z)
        {
            double t = T(ref rt, z);
            if (t <= 0) throw new Exception("Температура T <= 0");
            return 3.084e-4 / (Math.Exp(4.799e4 / t) - 1);
        }

        private static double DerivativeF(ref RadiativeTransfer rt, double z, double F, double U)
        {
            double commonPart = rt.R * rt.C * K(ref rt, z) * (UP(ref rt, z) - U);
            return Math.Abs(z) > EPS ? (-F / z + commonPart) : (commonPart / 2);
        }

        private static double Lambda(ref RadiativeTransfer rt, double zn)
        {
            return rt.C / (3 * K(ref rt, zn));
        }

        private static double P(ref RadiativeTransfer rt, double zn)
        {
            return rt.C * K(ref rt, zn);
        }

        private static double F(ref RadiativeTransfer rt, double zn)
        {
            return rt.C * K(ref rt, zn) * UP(ref rt, zn);
        }

        private static double Kappa(ref RadiativeTransfer rt, double zn, double h)
        {
            return (Lambda(ref rt, zn) + Lambda(ref rt, zn + h)) / 2;
        }

        private static double V(ref RadiativeTransfer rt, double zn, double h)
        {
            return (Math.Pow(zn + h / 2.0, 2) - Math.Pow(zn - h / 2.0, 2)) / 2.0;
        }

        private static double A(ref RadiativeTransfer rt, double zn, double h)
        {
            return (zn - h / 2.0) * Kappa(ref rt, zn - h, h) / (rt.R * rt.R * h);
        }

        private static double D(ref RadiativeTransfer rt, double zn, double h)
        {
            return (zn + h / 2.0) * Kappa(ref rt, zn, h) / (rt.R * rt.R * h);
        }

        private static double B(ref RadiativeTransfer rt, double zn, double h)
        {
            return A(ref rt, zn, h) + D(ref rt, zn, h) + P(ref rt, zn) * V(ref rt, zn, h);
        }

        private static double FVal(ref RadiativeTransfer rt, double zn, double h)
        {
            return F(ref rt, zn) * V(ref rt, zn, h);
        }

        private static void ThomasAlgorithm(ref RadiativeTransfer rt, int N, double[] y)
        {
            double h = 1.0 / N;
            double[] ksi = new double[N + 1];
            double[] eta = new double[N + 1];

            double m0 = (Kappa(ref rt, 0.0, h) * h / 2.0) / (rt.R * rt.R * h) + P(ref rt, 0.0) * h * 0;
            double k0 = -(Kappa(ref rt, 0.0, h) * h / 2.0) / (rt.R * rt.R * h);
            double q0 = FVal(ref rt, 0, h) * h * 0;

            ksi[1] = -k0 / m0;
            eta[1] = q0 / m0;

            Console.WriteLine($"ksi[1] = {ksi[1]}, eta[1] = {eta[1]}");

            for (int n = 1; n <= N - 1; n++)
            {
                double xn = n * h;
                double denom = B(ref rt, xn, h) - A(ref rt, xn, h) * ksi[n];
                ksi[n + 1] = D(ref rt, xn, h) / denom;
                eta[n + 1] = (FVal(ref rt, xn, h) + A(ref rt, xn, h) * eta[n]) / denom;
            }

            double kR = K(ref rt, 1.0);
            double alpha = 0.39 * 3 * kR * h;
            y[N] = eta[N] / (1 + alpha - ksi[N]);
            Console.WriteLine($"y[N] = {y[N]}");

            for (int n = N - 1; n >= 0; n--)
                y[n] = ksi[n + 1] * y[n + 1] + eta[n + 1];

            Console.WriteLine($"y[0] = {y[0]}");
        }

        private static double IntegrateTrapezoid(ref RadiativeTransfer rt, double z, double[] y, int N, int currentIndex)
        {
            if (currentIndex == 0) return 0.0;
            double h = z / currentIndex;
            double sum = 0.0;

            for (int j = 0; j <= currentIndex; j++)
            {
                double xj = j * h;
                double kVal = K(ref rt, xj);
                double upVal = UP(ref rt, xj);
                double f = kVal * (upVal - y[j]) * xj;
                sum += (j == 0 || j == currentIndex) ? f / 2.0 : f;
            }

            return sum * h;
        }

        private static int HalfDivisionMethod(ref RadiativeTransfer rt, int initialN, double epsilon, out double[] yValues)
        {
            int N = initialN;
            int maxIterations = 100;
            double[] yPrev = new double[N + 1];
            double[] yCurr = new double[N + 1];
            ThomasAlgorithm(ref rt, N, yCurr);

            for (int iter = 0; iter < maxIterations; iter++)
            {
                N *= 2;
                yPrev = yCurr;
                yCurr = new double[N + 1];
                ThomasAlgorithm(ref rt, N, yCurr);

                double maxDiff = 0.0;
                for (int i = 0; i <= initialN; i++)
                {
                    int idx = i * (N / initialN);
                    double diff = Math.Abs(yCurr[idx] - yPrev[i]);
                    if (diff > maxDiff)
                        maxDiff = diff;
                }

                if (maxDiff < epsilon)
                {
                    yValues = yCurr;
                    return N;
                }
                initialN = N;
            }

            yValues = yCurr;
            return N;
        }

        private static Form CreatePlotForm(Plot plt, string title, string fileName)
        {
            // Сохранение графика
            plt.SaveFig(fileName);
            Console.WriteLine($"График {fileName} сохранён");

            // Создание формы для графика
            var form = new Form
            {
                Text = title,
                Size = new System.Drawing.Size(1200, 800)
            };
            var formsPlot = new ScottPlot.FormsPlot
            {
                Dock = DockStyle.Fill
            };
            foreach (var plottable in plt.GetPlottables())
            {
                formsPlot.Plot.Add(plottable);
            }
            formsPlot.Refresh();
            form.Controls.Add(formsPlot);
            return form;
        }

        public static void Execute()
        {
            try
            {
                Console.WriteLine("Начало Task2");
                RadiativeTransfer rt = new RadiativeTransfer
                {
                    C = 3e10,
                    R = 0.35,
                    Tw = 2000,
                    T0 = 10000,
                    K0 = 1e-4,
                    P = 4,
                    M = 0.39
                };
                Console.WriteLine("Инициализация rt завершена");

                double[] y;
                double eps = 1e-4;
                int N = HalfDivisionMethod(ref rt, 100, eps, out y);
                double h = 1.0 / N;
                double[] F = new double[N + 1];
                double[] uPArr = new double[N + 1];
                double[] kArr = new double[N + 1];
                double[] duDz = new double[N + 1];
                double[] dFDz = new double[N + 1];

                Console.WriteLine($"N = {N}");
                double z;
                F[0] = 0.0;
                for (int i = 1; i <= N; i++)
                {
                    z = i * h;
                    if (z < EPS)
                        F[i] = 0.0;
                    else
                    {
                        double integral = IntegrateTrapezoid(ref rt, z, y, N, i);
                        F[i] = (rt.R * rt.C / z) * integral;
                    }
                }

                for (int i = 0; i <= N; i++)
                {
                    z = i * h;
                    uPArr[i] = UP(ref rt, z);
                    kArr[i] = K(ref rt, z);
                    duDz[i] = DerivativeU(ref rt, z, F[i]);
                    dFDz[i] = i == 0 ? (F[i + 1] - F[i]) / h :
                              i == N ? (F[i] - F[i - 1]) / h :
                              (F[i + 1] - F[i - 1]) / (2 * h);
                }

                using (var writer = new StreamWriter("task2_data.dat"))
                {
                    for (int i = 0; i <= N; i++)
                    {
                        z = i * h;
                        writer.WriteLine($"{z} {y[i]} {uPArr[i]} {F[i]} {kArr[i]} {duDz[i]} {dFDz[i]}");
                    }
                }
                Console.WriteLine("Файл task2_data.dat сохранён");

                // Подготовка графиков
                var forms = new List<Form>();
                Application.EnableVisualStyles();
                Application.SetCompatibleTextRenderingDefault(false);

                var zz = Enumerable.Range(0, N + 1).Select(i => (double)i * h).ToArray();
                var plt = new Plot(1200, 800);
                plt.AddScatter(zz, y, label: "u(z)").LineWidth = 2;
                plt.Title("u(z)");
                plt.XLabel("z");
                plt.YLabel("u(z)");
                forms.Add(CreatePlotForm(plt, "u(z)", "Task2_u.png"));

                plt = new Plot(1200, 800);
                plt.AddScatter(zz, F, label: "F(z)").LineWidth = 2;
                plt.Title("F(z)");
                plt.XLabel("z");
                plt.YLabel("F(z)");
                forms.Add(CreatePlotForm(plt, "F(z)", "Task2_F.png"));

                plt = new Plot(1200, 800);
                plt.AddScatter(zz, uPArr, label: "u_p(z)").LineWidth = 2;
                plt.Title("u_p(z)");
                plt.XLabel("z");
                plt.YLabel("u_p(z)");
                forms.Add(CreatePlotForm(plt, "u_p(z)", "Task2_up.png"));

                plt = new Plot(1200, 800);
                plt.AddScatter(zz, kArr, label: "k(z)").LineWidth = 2;
                plt.Title("k(z)");
                plt.XLabel("z");
                plt.YLabel("k(z)");
                forms.Add(CreatePlotForm(plt, "k(z)", "Task2_k.png"));

                plt = new Plot(1200, 800);
                plt.AddScatter(zz, duDz, label: "du/dz(z)").LineWidth = 2;
                plt.Title("du/dz(z)");
                plt.XLabel("z");
                plt.YLabel("du/dz(z)");
                forms.Add(CreatePlotForm(plt, "du/dz(z)", "Task2_dudz.png"));

                plt = new Plot(1200, 800);
                plt.AddScatter(zz, dFDz, label: "dF/dz(z)").LineWidth = 2;
                plt.Title("dF/dz(z)");
                plt.XLabel("z");
                plt.YLabel("dF/dz(z)");
                forms.Add(CreatePlotForm(plt, "dF/dz(z)", "Task2_dFdz.png"));

                // Открытие всех форм одновременно
                foreach (var form in forms)
                {
                    form.Show();
                }
                Application.Run();
            }
            catch (Exception ex)
            {
                Console.WriteLine($"Ошибка в Task2: {ex.Message}");
            }
        }
    }
}