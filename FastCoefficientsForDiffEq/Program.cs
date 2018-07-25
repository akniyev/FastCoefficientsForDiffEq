﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using static System.Math;

namespace FastCoefficientsForDiffEq
{
    class Program
    {
        private static Complex[] _expArray;
        private static readonly Complex ImagNum = new Complex(0, 2);
        private static Complex[] InvFft(Complex[] array)
        {
            var n = array.Length;

            if (n % 2 != 0 && n != 1) throw new ArgumentException();

            if (n == 1)
            {
                return array;
            }

            var omegaN = new Complex(Cos(2 * PI / n), Sin(2 * PI / n));

            var omega = new Complex(1, 0);

            //Comment
            var even = new Complex[n / 2];
            var odd = new Complex[n / 2];

            for (var i = 0; i < n / 2; i++)
            {
                even[i] = array[2 * i];
                odd[i] = array[2 * i + 1];
            }

            var y0 = InvFft(even);
            var y1 = InvFft(odd);

            var y = new Complex[n];

            for (var k = 0; k < n / 2; k++)
            {
                y[k] = y0[k] + omega * y1[k];
                y[k + n / 2] = y0[k] - omega * y1[k];
                omega = omega * omegaN;
            }

            return y;
        }
        private static double[] Ettas(IReadOnlyList<double> a)
        {
            var n = a.Count;
            var b = new double[2 * n];
            b[0] = 0;
            for (var i = 1; i < n; i++)
            {
                b[i] = a[i] / i;
            }
            b[n] = 0;
            for (var i = 1; i < n; i++)
            {
                b[2 * n - i] = -a[i] / i;
            }
            var h = new Complex[n];
            for (int i = 0; i < n; i++)
            {
                h[i] = new Complex(b[2 * i], b[2 * i + 1]);
            }

            var a1 = InvFftForRealInput(h);
            var c = new Complex[n];

            for (var i = 0; i < c.Length; i++)
            {
                c[i] = Sqrt(2.0) / PI * a1[i] / ImagNum;
            }

            return c.Select(x => x.Real).ToArray();
        }
        private static Complex[] InvFftForRealInput(Complex[] a)
        {
            var n = a.Length;

            if (_expArray == null || _expArray.Length != n)
            {
                _expArray = new Complex[n];

                for (var k = 0; k < n; k++)
                {
                    _expArray[k] = new Complex(Cos(PI * k / n), Sin(PI * k / n));
                }
            }

            var imag = new Complex(0, 1);

            var b = new Complex[n];

            var a1 = InvFft(a);

            b[0] = new Complex(a1[0].Real, -a1[0].Imaginary);
            for (int i = 1; i < n; i++)
            {
                b[i] = new Complex(a1[n - i].Real, -a1[n - i].Imaginary);
            }
            var g1 = new Complex[n];
            var g2 = new Complex[n];
            var b1 = new Complex[2 * n];
            for (int i = 0; i < n; i++)
            {
                g1[i] = (a1[i] + b[i]) / 2;
                g2[i] = -imag * (a1[i] - b[i]) / 2;
            }
            for (int i = 0; i < n; i++)
            {
                b1[i] = g1[i] + _expArray[i] * g2[i];
                //b1[i] = g1[i] - _expArray[i] * g2[i];
            }
            return b1;
        }
        private static Complex[] Fft(Complex[] array)
        {
            var n = array.Length;

            if (n % 2 != 0 && n != 1) throw new ArgumentException();

            if (n == 1)
            {
                return array;
            }

            var omegaN = new Complex(Cos(2 * PI / n), -Sin(2 * PI / n));

            var omega = new Complex(1, 0);

            var even = new Complex[n / 2];
            var odd = new Complex[n / 2];
            for (var i = 0; i < n / 2; i++)
            {
                even[i] = array[2 * i];
                odd[i] = array[2 * i + 1];
            }

            var y0 = Fft(even);
            var y1 = Fft(odd);

            var y = new Complex[n];

            for (var k = 0; k < n / 2; k++)
            {
                y[k] = y0[k] + omega * y1[k];
                y[k + n / 2] = y0[k] - omega * y1[k];
                omega = omega * omegaN;
            }

            return y;
        }
        private static Complex[] FftForRealInput(Complex[] a)
        {
            var n = a.Length;
            // n = n / 2;

            var expArray = new Complex[n];

            for (var k = 0; k < n; k++)
            {
                expArray[k] = new Complex(Cos(PI * k / n), -Sin(PI * k / n));
            }
            var imag = new Complex(0, 1);

            var b = new Complex[n];

            var a1 = Fft(a);

            b[0] = new Complex(a1[0].Real, -a1[0].Imaginary);
            for (int i = 1; i < n; i++)
            {
                b[i] = new Complex(a1[n - i].Real, -a1[n - i].Imaginary);
            }
            var g1 = new Complex[n];
            var g2 = new Complex[n];
            var b1 = new Complex[n];
            for (int i = 0; i < n; i++)
            {
                g1[i] = (a1[i] + b[i]) / 2;
                g2[i] = -imag * (a1[i] - b[i]) / 2;
            }
            for (int i = 0; i < n; i = i + 1)
            {
                b1[i] = g1[i] + expArray[i] * g2[i];
                //b1[i + n] = g1[i] - expArray[i] * g2[i];
            }
            return b1;
        }
        private static double[] EttasWithInitialValues(IReadOnlyList<double> p, double y0)
        {
            return Ettas(p).Select(x => x + y0).ToArray();
        }

        // You should calculate ettas and hsin before calling this method
        //static Complex[] Cks(double[] f)
        //{
        //    var N = f.Length;

        //    var f1 = f.Append(0).Concat(f.Skip(1).Reverse()).ToArray();

        //    var f1_complex = Enumerable.Range(0, N).Select(i => new Complex(f1[2 * i], f1[2 * i + 1])).ToArray();

        //    var a1 = FftForRealInput(f1_complex);
        //    a1[0] = PI / N * a1[0] / 2;
        //    for (int i = 1; i < a1.Length; i++)
        //    {
        //        a1[i] = PI / N * (a1[i] + f[0]) / Sqrt(2);
        //    }
        //    //    .Select(x => (x + f[0]) / 2.0)
        //    //    .Select(x => Sqrt(2) * PI / N).ToArray();
        //    //a1[0] /= Sqrt(2); 
        //    return a1;
        //}
        static double[] Cks(double[] p, double h, double y0)
        {
            //y0 = 1;
            //h = 0.1;
            var N = p.Length;
            var ettas = EttasWithInitialValues(p, y0);
            var f = new Func<double, double, double>((x, y) => x + y);
            var indexes = Enumerable.Range(0, N).ToArray();

            var tau = indexes.Select(x => 1.0 * x / N).ToArray();

            var hsin = indexes.Select(j => h * Sin(PI * tau[j])).ToArray();
            var cos = indexes.Select(j => Cos(PI * tau[j])).ToArray();
            var dicrete_f = indexes.Select(j => f(hsin[j], ettas[j]) * cos[j]).ToArray();

            var f1 = dicrete_f.Append(0).Concat(dicrete_f.Skip(1).Reverse()).ToArray();

            var f1_complex = Enumerable.Range(0, N).Select(i => new Complex(f1[2 * i], f1[2 * i + 1])).ToArray();

            var a1 = FftForRealInput(f1_complex);
            a1[0] = PI / N * a1[0] / 2;
            for (int i = 1; i < a1.Length; i++)
            {
                a1[i] = PI / N * (a1[i] + f1[0]) / Sqrt(2);
            }
            //    .Select(x => (x + f[0]) / 2.0)
            //    .Select(x => Sqrt(2) * PI / N).ToArray();
            //a1[0] /= Sqrt(2);
            var a2 = new double[N];
            for (int i = 0; i < N; i++)
            {
                a2[i] = a1[i].Real;
            }
            return a2;
        }
        static double diff(double[] a, double[] b)
        {
            var result = new double [a.Length];
            
            for (int i = 0; i < a.Length; i++)
            {
                result[i] = (a[i] - b[i]);
            }
            var s =0.0;
            for (int i = 0; i < a.Length; i++)
            {
                s = s + Sqrt(result[i] * result[i]);
            }
            return s;
        }
        static void Main(string[] args)
        {
            //var N = 8;
            var p = new double[] { 1, 2, 3, 4, 5, 6, 7, 8};
            //var y0 = 1;
            //var h = 0.1;
            //var f = new Func<double, double, double>((x, y) => x + y);

            //var ettas = EttasWithInitialValues(p, y0);

            //var indexes = Enumerable.Range(0, N).ToArray();

            //var tau = indexes.Select(x => 1.0 * x / N).ToArray();

            //var hsin = indexes.Select(j => h * Sin(PI * tau[j])).ToArray();
            //var cos = indexes.Select(j => Cos(PI * tau[j])).ToArray();
            //var dicrete_f = indexes.Select(j => f(hsin[j], ettas[j])*cos[j]).ToArray();

            var eps = 0.01;
            var cks = Cks(p, 0.1, 1);
            var k = new double[p.Length];
            do
            {
                k = cks;
                cks = Cks(cks, 0.1, 1);
                Console.WriteLine(diff(cks, k));
            } while (diff(cks,k)>=eps);

            //Console.WriteLine(diff(cks, k));

            Console.WriteLine("Cks");

            foreach (var ck in cks)
            {
                Console.WriteLine(ck);
            }
            Console.ReadLine();
        }
        
        
    }
}