using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows.Forms;
using MathNet.Numerics;
using MathNet.Numerics.Integration;
using MathNet.Numerics.Interpolation;
using ScottPlot;
using ScottPlot.Plottable;

public class PlasmaSimulation
{
    private const double c = 3e10; // Speed of light
    private const double PlanckConstant = 3.084e-4;

    public class Configuration
    {
        public double R = 0.35;
        public double T0 = 8000;
        public double Tw = 1800;
        public double p = 2;
        public double Imax = 1000;
        public double Itmax = 80e-6;
        public double tmax = 80e-6;
        public double tau = 2e-6;
        public double zmax = 1;
        public double t0 = 0;
        public double z0 = 0;
        public double eps = 1e-8;

        public int zSteps = 100;
        public int tSteps = 200;
        public double Tau
        {
            get => tau;
            set
            {
                tau = value;
                CalculateSteps();
            }
        }

        public int TSteps
        {
            get => tSteps;
            set
            {
                tSteps = value;
                tau = (tmax - t0) / tSteps;
            }
        }
        public void CalculateSteps()
        {
            tSteps = (int)Math.Ceiling((tmax - t0) / tau) + 1;
        }



        public double[] Tarray = { 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000 };
        public double[] SigmaArray = { 0.309e-3, 0.309e-2, 0.309e-1, 0.270, 0.205e+1, 0.606e+1, 0.120e+2, 0.199e+2, 0.296e+2, 0.411e+2, 0.541e+2 };
        public double[] LambdaArray = { 0.381e-3, 0.381e-3, 0.381e-3, 0.448e-3, 0.577e-3, 0.733e-3, 0.131e-2, 0.218e-2, 0.358e-2, 0.562e-2, 0.832e-2 };
        public double[] cTArray = { 1.90e-3, 1.90e-3, 0.95e-3, 0.75e-3, 0.64e-3, 0.61e-3, 0.66e-3, 0.66e-3, 1.15e-3, 1.79e-3, 2.02e-3 };

        // For k interpolation
        public double[] TkArray = { 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000 };
        public int KVariant = 1;
        public double[] K1Array = { 8.200E-03, 2.768E-02, 6.560E-02, 1.281E-01, 2.214E-01, 3.516E-01, 5.248E-01, 7.472E-01, 1.025E+00 };
        public double[] K2Array = { 1.6, 5.4, 1.280E+01, 2.500E+01, 4.320E+01, 6.860E+01, 1.024E+02, 1.458E+02, 2.000E+02 };

        public double I(double t)
        {
            return Imax / Itmax * t * Math.Exp(-(t / Itmax - 1));
        }

        public double Sigma(double T)
        {
            var logT = Tarray.Select(x => Math.Log(x)).ToArray();
            var logSigma = SigmaArray.Select(x => Math.Log(x)).ToArray();
            return Math.Exp(Interpolate.Linear(logT, logSigma).Interpolate(Math.Log(T)));
        }

        public double Lambda(double T)
        {
            var logT = Tarray.Select(x => Math.Log(x)).ToArray();
            var logLambda = LambdaArray.Select(x => Math.Log(x)).ToArray();
            return Math.Exp(Interpolate.Linear(logT, logLambda).Interpolate(Math.Log(T)));
        }

        public double cT(double T)
        {
            var logT = Tarray.Select(x => Math.Log(x)).ToArray();
            var logcT = cTArray.Select(x => Math.Log(x)).ToArray();
            return Math.Exp(Interpolate.Linear(logT, logcT).Interpolate(Math.Log(T)));
        }

        public double k(double T)
        {
            var logT = TkArray.Select(x => Math.Log(x)).ToArray();
            double[] logK = KVariant == 1
                ? K1Array.Select(x => Math.Log(x)).ToArray()
                : K2Array.Select(x => Math.Log(x)).ToArray();
            return Math.Exp(Interpolate.Linear(logT, logK).Interpolate(Math.Log(T)));
        }

        public double E(double[] z, double[] T, double t)
        {
            if (z.Length != T.Length)
                throw new ArgumentException("z and T must have the same length");

            double integral = 0;
            for (int i = 0; i < z.Length - 1; i++)
            {
                double dz = z[i+1] - z[i];
                integral += (Sigma(T[i]) * z[i] + Sigma(T[i + 1]) * z[i + 1]) * dz / 2;
            }
            return I(t) / (2 * Math.PI * R * R * integral);
        }

        public double Tinit(double z)
        {
            return T0 + (Tw - T0) * Math.Pow(z, p);
        }

        public double Plank(double T)
        {
            return PlanckConstant / (Math.Exp(4.799e4 / T) - 1);
        }

        public double q(double[] z, double[] u, double[] T, int zindex)
        {
            if (z.Length != u.Length)
                throw new ArgumentException("z and u must have the same length");
            if (zindex >= z.Length)
                throw new ArgumentException("zindex must be less than z.Length");

            return c * k(T[zindex]) * (Plank(T[zindex]) - u[zindex]);
        }
    }

    public static double[] TridiagonalSolve(double[][] Mat, double[] F)
    {
        int n = Mat.Length;

        if (Mat.Any(row => row.Length != n))
            throw new ArgumentException("Matrix must be square");
        if (F.Length != n)
            throw new ArgumentException("F vector length must match matrix size");

        double[] A = new double[n];
        double[] B = new double[n];
        double[] C = new double[n];

        for (int i = 0; i < n; i++)
        {
            A[i] = i - 1 >= 0 ? Mat[i][i - 1] : 0;
            B[i] = Mat[i][i];
            C[i] = i + 1 < n ? Mat[i][i + 1] : 0;
        }

        double[] eps = new double[n - 1];
        double[] eta = new double[n - 1];
        double[] x = new double[n];

        eps[0] = -C[0] / B[0];
        eta[0] = F[0] / B[0];

        for (int i = 1; i < n - 1; i++)
        {
            double denominator = B[i] - A[i] * eps[i - 1];
            eps[i] = C[i] / denominator;
            eta[i] = (F[i] + A[i] * eta[i - 1]) / denominator;
        }

        x[n - 1] = (F[n - 1] - A[n - 1] * eta[n - 2]) / (A[n - 1] * eps[n - 2] + B[n - 1]);

        for (int i = n - 2; i >= 0; i--)
        {
            x[i] = eta[i] + eps[i] * x[i + 1];
        }

        return x;
    }

    public static double[] Solveu(Configuration config, int n, double[] T)
    {
        double[] z = MathNet.Numerics.Generate.LinearSpaced(n, config.z0, config.zmax);
        double h = z[1] - z[0];

        double[][] mat = new double[n][];
        for (int i = 0; i < n; i++) mat[i] = new double[n];

        Func<int, double> k = i => config.k(T[i]);
        Func<int, double> kappa = i => 2 / (k(i) + k(i + 1));

        mat[0][0] = 1 / ((k(0) + k(1)) * config.R);
        mat[0][1] = -mat[0][0];

        for (int i = 1; i < n - 1; i++)
        {
            double A = -(z[i] - h / 2) * kappa(i - 1) / (config.R * config.R * h);
            double C = -(z[i] + h / 2) * kappa(i) / (config.R * config.R * h);
            double B = (A + C - 3 * k(i) * z[i] * h);

            mat[i][i - 1] = A;
            mat[i][i] = B;
            mat[i][i + 1] = C;
        }

        mat[n - 1][n - 2] = -(1 / (k(n - 2) + k(n - 1)) / (config.R * config.R * h) * (2 - h));
        mat[n - 1][n - 1] = 3 * 0.39 / config.R - mat[n - 1][n - 2];

        double[] side = new double[n];
        for (int i = 1; i < n - 1; i++)
        {
            side[i] = -3 * k(i) * z[i] * h * config.Plank(T[i]);
        }

        return TridiagonalSolve(mat, side);
    }

    public static double[] Iteration(Configuration config, double[] iterT, double[] previousT,
                                   double[] z, double[] u, double tau, double t)
    {
        if (z.Length != u.Length)
            throw new ArgumentException("z and u must have the same length");
        if (z.Length != previousT.Length)
            throw new ArgumentException("z and previousT must have the same length");
        if (z.Length != iterT.Length)
            throw new ArgumentException("z and iterT must have the same length");

        double h = z[1] - z[0];

        Func<int, double> kappa = i => 2 * (config.Lambda(iterT[i - 1]) * config.Lambda(iterT[i])) /
                                     (config.Lambda(iterT[i - 1]) + config.Lambda(iterT[i]));

        double[][] mat = new double[z.Length][];
        for (int i = 0; i < z.Length; i++) mat[i] = new double[z.Length];

        mat[0][0] = tau / (config.R * config.R) *
                   (config.Lambda(iterT[0]) * config.Lambda(iterT[1])) /
                   (config.Lambda(iterT[0]) + config.Lambda(iterT[1]));
        mat[0][1] = -mat[0][0];

        double[] F = new double[z.Length];
        F[0] = 0;

        double E = config.E(z, iterT, t);

        for (int i = 1; i < z.Length - 1; i++)
        {
            double A = tau / (config.R * config.R) * (z[i] - h / 2) * kappa(i) / h;
            double C = tau / (config.R * config.R) * (z[i] + h / 2) * kappa(i + 1) / h;
            double B = A + C + z[i] * h * config.cT(iterT[i]);

            F[i] = z[i] * config.cT(iterT[i]) * previousT[i] * h +
                   tau * h * z[i] * config.Sigma(iterT[i]) * E * E -
                   tau * h * config.q(z, u, iterT, i) * z[i];

            mat[i][i - 1] = A;
            mat[i][i] = B;
            mat[i][i + 1] = C;
        }

        mat[z.Length - 1][z.Length - 2] = 0;
        mat[z.Length - 1][z.Length - 1] = 1;
        F[z.Length - 1] = config.Tw;

        return TridiagonalSolve(mat, F);
    }

    public static (double[] newT, double[] u) SolveTimeStep(Configuration config, double[] previousT,
                                                          double[] z, double t, double tau)
    {
        if (z.Length != previousT.Length)
            throw new ArgumentException("z and previousT must have the same length");

        double[] u = Solveu(config, z.Length, previousT);
        double[] iterT = (double[])previousT.Clone();

        int cnt = 0;
        while (true)
        {
            cnt++;
            double[] newT = Iteration(config, iterT, previousT, z, u, tau, t);

            if (newT.Any(double.IsNaN))
                break;

            bool converged = true;
            for (int i = 0; i < newT.Length; i++)
            {
                if (Math.Abs(newT[i] - iterT[i]) / newT[i] >= config.eps)
                {
                    converged = false;
                    break;
                }
            }

            if (converged) break;

            iterT = newT;
            u = Solveu(config, z.Length, iterT);
        }
        Console.WriteLine("AU");
        return (iterT, u);
    }

    public static (double[] zArr, double[] tArr, List<double[]> uSlices, List<double[]> timeSlices)
        Solve(Configuration config, int zSteps, int tSteps, bool verbose = true)
    {
        config.CalculateSteps();
        double[] zArr = MathNet.Numerics.Generate.LinearSpaced(zSteps, config.z0, config.zmax);
        double[] tArr = MathNet.Numerics.Generate.LinearSpaced(config.tSteps, config.t0, config.tmax);

        double tau = tArr[1] - tArr[0];
        double h = zArr[1] - zArr[0];

        var timeSlices = new List<double[]>();
        timeSlices.Add(zArr.Select(z => config.Tinit(z)).ToArray());

        var uSlices = new List<double[]>();
        uSlices.Add(Solveu(config, zSteps, timeSlices.Last()));

        foreach (double t in tArr.Skip(1))
        {
            var (temp, u) = SolveTimeStep(config, timeSlices.Last(), zArr, t, tau);

            if (verbose)
            {
                Console.WriteLine($"iter{timeSlices.Count} temp[0]={temp[0]}");
                Console.WriteLine($"iter{timeSlices.Count} u[0]={u[0]}");
            }

            timeSlices.Add(temp);
            uSlices.Add(u);
        }

        return (zArr, tArr, uSlices, timeSlices);
    }

    [STAThread]
    public static void Main()
    {
        var config = new Configuration
        {
            zSteps = 100,
            //tSteps = 200,
            //tmax = 80e-6,
            //tau = 80e-6 / 200,
            tau = 80e-6 / 500,
            tmax = 80e-6,

            Imax = 300
        };

        var (zArr, tArr, uSlices, timeSlices) = Solve(config, config.zSteps, config.tSteps, false);

        // Create plots using ScottPlot
        var plt1 = new Plot(1680, 800);
        for (int i = 0; i < timeSlices.Count; i += 50)
        {
            plt1.AddScatter(zArr, timeSlices[i], label: $"t = {tArr[i]:e2}");
        }
        plt1.Title("Temperature vs Position at Different Times");
        plt1.XLabel("Position z");
        plt1.YLabel("Temperature T");
        plt1.Legend();
        plt1.SaveFig("Temperature_vs_Position.png");

        Application.Run(new FormsPlotViewer(plt1, 1680, 800));

        var plt2 = new Plot(1680, 800);
        for (int i = 0; i < zArr.Length; i += 5)
        {
            double[] tempArr = timeSlices.Select(ts => ts[i]).ToArray();
            plt2.AddScatter(tArr, tempArr, label: $"z = {zArr[i]:.2f}");
        }
        plt2.Title("Temperature vs Time at Different Positions");
        plt2.XLabel("Time t");
        plt2.YLabel("Temperature T");
        plt2.Legend();
        plt2.SaveFig("Temperature_vs_Time.png");

        Application.Run(new FormsPlotViewer(plt2, 1680, 800));

        Console.WriteLine("Simulation completed successfully");
        Console.WriteLine("Plots saved as 'Temperature_vs_Position.png' and 'Temperature_vs_Time.png'");
    }
}