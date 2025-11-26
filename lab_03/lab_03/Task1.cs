using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Windows.Forms;
using ScottPlot;

namespace Lab3Modeling
{
    public static class Task1
    {
        private class Matrix
        {
            public int Size { get; }
            public double[][] Data { get; }
            public double[] Rhs { get; }

            public Matrix(int size)
            {
                Size = size;
                Data = new double[size][];
                for (int i = 0; i < size; i++)
                    Data[i] = new double[size];
                Rhs = new double[size];
            }
        }

        private static Matrix BuildSystemN3()
        {
            var m = new Matrix(3);
            //m.Data[0][0] = 38; m.Data[0][1] = 56; m.Data[0][2] = 98;
            //m.Rhs[0] = 35;
            //m.Data[1][0] = 37; m.Data[1][1] = 44; m.Data[1][2] = 42;
            //m.Rhs[1] = 21;
            //m.Data[2][0] = 31; m.Data[2][1] = 32; m.Data[2][2] = 22;
            //m.Rhs[2] = 14;

            //m.Data[0][0] = -11 / 15; m.Data[0][1] = -47 / 30; m.Data[0][2] = -86 / 35;
            //m.Rhs[0] = -5 / 12;
            //m.Data[1][0] = -560.0; m.Data[1][1] = -1248; m.Data[1][2] = -1995;
            //m.Rhs[1] = -366;
            //m.Data[2][0] = -66 / 35; m.Data[2][1] = -43 / 10; m.Data[2][2] = -146 / 21;
            //m.Rhs[2] = -7 / 6;

            m.Data[0][0] = -308.0; m.Data[0][1] = -658.0; m.Data[0][2] = -1032.0;
            m.Rhs[0] = -175.0;
            m.Data[1][0] = -4.0 / 3.0; m.Data[1][1] = -104.0 / 35.0; m.Data[1][2] = -19.0 / 4.0;
            m.Rhs[1] = -4.0 / 5.0;
            m.Data[2][0] = -396.0; m.Data[2][1] = -903.0; m.Data[2][2] = -1460.0;
            m.Rhs[2] = -245.0;
            return m;
        }

        private static void GaussianElimination(Matrix m, double[] solution)
        {
            int n = m.Size;
            for (int i = 0; i < n; i++)
            {
                int maxRow = i;
                for (int k = i + 1; k < n; k++)
                {
                    if (Math.Abs(m.Data[k][i]) > Math.Abs(m.Data[maxRow][i]))
                        maxRow = k;
                }
                if (maxRow != i)
                {
                    var tempRow = m.Data[i];
                    m.Data[i] = m.Data[maxRow];
                    m.Data[maxRow] = tempRow;
                    double tempRhs = m.Rhs[i];
                    m.Rhs[i] = m.Rhs[maxRow];
                    m.Rhs[maxRow] = tempRhs;
                }
                double pivot = m.Data[i][i];
                if (Math.Abs(pivot) < 1e-10)
                {
                    Console.WriteLine("Матрица сингулярна или почти сингулярна.");
                    return;
                }
                for (int j = i; j < n; j++)
                    m.Data[i][j] /= pivot;
                m.Rhs[i] /= pivot;
                for (int k = i + 1; k < n; k++)
                {
                    double factor = m.Data[k][i];
                    for (int j = i; j < n; j++)
                        m.Data[k][j] -= factor * m.Data[i][j];
                    m.Rhs[k] -= factor * m.Rhs[i];
                }
            }
            for (int i = n - 1; i >= 0; i--)
            {
                solution[i] = m.Rhs[i];
                for (int j = i + 1; j < n; j++)
                    solution[i] -= m.Data[i][j] * solution[j];
            }
        }

        private static double BasisF(double x, int i)
        {
            int integer = i + 2;
            return Math.Pow(x, integer) - x * integer;
        }

        private static double DerBasisF(double x, int i)
        {
            int integer = i + 2;
            return integer * Math.Pow(x, integer - 1) - integer;
        }

        private static double A(double xn, double h)
        {
            return (1.0 / (h * h)) + (xn / h);
        }

        private static double B(double xn, double h)
        {
            return (2.0 / (h * h)) - 2.0;
        }

        private static double D(double xn, double h)
        {
            return (1.0 / (h * h)) - (xn / h);
        }

        private static double F(double xn)
        {
            return -xn;
        }

        private static void ThomasAlgorithm(int N, double[] y)
        {
            double h = 1.0 / N;
            double[] ksi = new double[N + 1];
            double[] eta = new double[N + 1];

            ksi[1] = -1.0;
            eta[1] = 0.0;
            y[0] = 1.7520e-07;

            for (int n = 1; n <= N - 1; n++)
            {
                double xn = n * h;
                double denom = B(xn, h) - A(xn, h) * ksi[n];
                ksi[n + 1] = D(xn, h) / denom;
                eta[n + 1] = (F(xn) + A(xn, h) * eta[n]) / denom;
            }

            y[N] = (eta[N] + h) / (1 - ksi[N]);
            for (int n = N - 1; n >= 0; n--)
                y[n] = ksi[n + 1] * y[n + 1] + eta[n + 1];

            Console.WriteLine($"y[N] = {y[N]}");
        }

        private static int HalfDivisionMethod(int initialN, double epsilon, out double[] yValues)
        {
            int N = initialN;
            int maxIterations = 100;
            double[] yPrev = new double[N + 1];
            double[] yCurr = new double[N + 1];
            ThomasAlgorithm(N, yCurr);

            for (int iter = 0; iter < maxIterations; iter++)
            {
                N *= 2;
                yPrev = yCurr;
                yCurr = new double[N + 1];
                ThomasAlgorithm(N, yCurr);

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
                Size = new System.Drawing.Size(1000, 800)
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
                // Метод Галеркина (N=3)
                var m3 = BuildSystemN3();
                double[] solution3 = new double[m3.Size];
                GaussianElimination(m3, solution3);

                Console.WriteLine("Решения для N=3:");
                for (int i = 0; i < m3.Size; i++)
                    Console.WriteLine($"SolutionVector[{i}] = {solution3[i]}");

                // Подготовка данных для графиков
                var x = new double[101];
                var uN3 = new double[101];
                var duN3 = new double[101];
                var uN1 = new double[101];
                var duN1 = new double[101];

                for (int i = 0; i <= 100; i++)
                {
                    x[i] = i * 0.01;
                    uN3[i] = x[i];
                    duN3[i] = 1.0;
                    for (int j = 0; j < m3.Size; j++)
                    {
                        uN3[i] += solution3[j] * BasisF(x[i], j);
                        duN3[i] += solution3[j] * DerBasisF(x[i], j);
                    }
                    uN1[i] = (25.0 / 44.0) * (x[i] * x[i] - 2 * x[i]) + x[i];
                    duN1[i] = (25.0 / 44.0) * (2 * x[i] - 2) + 1;
                }

                // Разностный метод
                double[] y;
                int N = HalfDivisionMethod(10, 1e-4, out y);
                double h = 1.0 / N;
                double[] yX = new double[N + 1];
                double[] dyDx = new double[N];
                for (int i = 0; i <= N; i++)
                    yX[i] = i * h;
                for (int i = 0; i < N; i++)
                    dyDx[i] = (y[i + 1] - y[i]) / h;

                // Подготовка графиков
                var forms = new List<Form>();
                Application.EnableVisualStyles();
                Application.SetCompatibleTextRenderingDefault(false);

                var plt = new Plot(1000, 800);
                plt.AddScatter(x, uN3, label: "u(x) N=3").LineWidth = 2;
                plt.Title("Решение u(x) для N=3");
                plt.XLabel("x");
                plt.YLabel("u(x)");
                forms.Add(CreatePlotForm(plt, "Решение u(x) для N=3", "Task1_uN3.png"));

                plt = new Plot(1000, 800);
                plt.AddScatter(x, duN3, label: "du/dx N=3").LineWidth = 2;
                plt.Title("Производная du/dx для N=3");
                plt.XLabel("x");
                plt.YLabel("du/dx");
                forms.Add(CreatePlotForm(plt, "Производная du/dx для N=3", "Task1_duN3.png"));

                plt = new Plot(1000, 800);
                plt.AddScatter(x, uN1, label: "u(x) N=1").LineWidth = 2;
                plt.Title("Решение u(x) для N=1");
                plt.XLabel("x");
                plt.YLabel("u(x)");
                forms.Add(CreatePlotForm(plt, "Решение u(x) для N=1", "Task1_uN1.png"));

                plt = new Plot(1000, 800);
                plt.AddScatter(x, duN1, label: "du/dx N=1").LineWidth = 2;
                plt.Title("Производная du/dx для N=1");
                plt.XLabel("x");
                plt.YLabel("du/dx");
                forms.Add(CreatePlotForm(plt, "Производная du/dx для N=1", "Task1_duN1.png"));

                plt = new Plot(800, 600);
                plt.AddScatter(yX, y, label: "u(x)").LineWidth = 2;
                plt.Title($"Решение u(x) для N={N}");
                plt.XLabel("x");
                plt.YLabel("u(x)");
                forms.Add(CreatePlotForm(plt, $"Решение u(x) для N={N}", "Task1_fd_u.png"));

                plt = new Plot(800, 600);
                plt.AddScatter(yX.Take(N).ToArray(), dyDx, label: "du/dx").LineWidth = 2;
                plt.Title($"Производная du/dx для N={N}");
                plt.XLabel("x");
                plt.YLabel("du/dx");
                plt.SetAxisLimits(yMin: -1, yMax: 10);
                forms.Add(CreatePlotForm(plt, $"Производная du/dx для N={N}", "Task1_fd_du.png"));

                // Открытие всех форм одновременно
                foreach (var form in forms)
                {
                    form.Show();
                }
                Application.Run();
            }
            catch (Exception ex)
            {
                Console.WriteLine($"Ошибка в Task1: {ex.Message}");
            }
        }
    }
}