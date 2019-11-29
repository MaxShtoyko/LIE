using System;
using MathNet;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace LIE {
    class Program {
        static double PI = Math.PI;
        static double q = 7;
        static double n = 6;

        static void Main ( string[] args ) {

            var densities = QuadratureMethod ( );
            Console.WriteLine ( densities );

            Console.WriteLine ( GetApproximateResult ( densities, Vector.Build.DenseOfArray ( new double[] { 0, 0 } ) ) );
            for(int i =1; i<n;i++ ) {
               // Console.WriteLine ( $"Exact: {GetExactResult (  Z( S(i) )) } ----- Approximate {GetApproximateResult ( densities, Z ( S ( i ) ) )}" );
            }

            Console.ReadKey ( );
        }

        static Vector<double> QuadratureMethod ( ) {
            var matrix = Matrix<double>.Build.Dense ( ( int ) n - 1, ( int ) n - 1 );
            var b = Vector<double>.Build.Dense ( ( int ) n - 1 );

            for (int k = 1;k<n;k++ ) {
                for(int i = 1;i<n;i++ ) {
                    matrix[k - 1, i - 1] = GetQuadratureR ( Math.Abs ( i - k ) ) + GetQuadratureR ( i + k ) + 1 / n * GetL2 ( S ( k ), S ( i ) );
                }
                b[k - 1] = GetFunction ( );
            }

            Console.WriteLine ( matrix );

            return matrix.Solve ( b );
        }

        //Uex = Ф(x,y*)
        //static double GetExactResult (double[] x, double[] y ) {
        //    return 1 / ( 2 * PI ) * Math.Log ( 1 / GetR ( x, y ) );
        //}

        static double GetExactResult ( double[] x, double[] y ) {
            return 1;
        }

        static double GetFunction () {
            return 1;
        }

        static double GetApproximateResult ( Vector<double> densities, Vector<double> t ) {
            double result = 0;

            for(int k = 1; k<n;k++ ) {
                result += densities[k - 1] * Math.Log ( GetR ( t - Z ( S ( k ) ) ) );
            }

            return result * -1 / ( 2 * n );
        }

        static double S ( double j ) {
            return j * PI / n;
        }

        static double GetQuadratureR ( double j ) {

            double r = 0;

            r += 1;

            for ( double m = 1; m < n - 1; m++ ) {
                r += 2 / m * Math.Cos ( m * S ( j ) );
            }

            r += Math.Pow ( -1, j ) / Math.Pow ( n, 2 );

            return r * -1 / n;
        }


        static Vector<double> Z ( double s ) {
            return GetBoudaryFunction ( Y ( s ) );
        }

        static double Y ( double s ) {
            if ( s <= PI ) {
                return W ( s );
            }
            return PI + W ( s - PI );
        }

        static double V ( double s ) {
            return ( 1 / q - PI / 2 ) * Math.Pow ( ( PI - 2 * s ) / PI, 3 ) - 1 / q * ( PI - 2 * s ) / PI + PI / 2;
        }

        static double W ( double s ) {
            return PI * Math.Pow ( V ( s ), q ) / ( Math.Pow ( V ( s ), q ) + Math.Pow ( V ( PI - s ), q ) );
        }

        //EXAMPLE
        static double GetX1 ( double t ) {
            return Math.Sin ( t/2 );
        }
        
        //EXAMPLE
        static double GetX2 ( double t ) {
            return Math.Sin ( t );
        }

        //EXAMPLE
        static double GetDerivativeOfX1 ( double t ) {
            return 1 / 2 * Math.Cos ( t / 2 );
        }

        //EXAMPLE
        static double GetDerivativeOfX2 ( double t ) {
            return Math.Cos ( t / 2 );
        }

        // |x(s)-x(tau)|
        static double GetParametricR ( double s, double tau ) {
            return Math.Sqrt ( Math.Pow ( ( GetX1 ( s ) - GetX1 ( tau ) ), 2 ) + Math.Pow ( ( GetX2 ( s ) - GetX2 ( tau ) ), 2 ) );
        }

        static double GetR ( Vector<double> s ) {
            return Math.Sqrt ( Math.Pow ( s[0], 2 ) + Math.Pow ( s[1], 2 ) );
        }

        static double GetR ( double[] s, double[] tau ) {
            return Math.Sqrt ( Math.Pow ( s[0] - tau[0], 2 ) + Math.Pow ( s[1] - tau[1], 2 ) );
        }

        // |x'(s)|
        static double GetDerivativeOfR ( double s ) {
            return Math.Sqrt ( Math.Pow ( GetDerivativeOfX1 ( s ), 2 ) + Math.Pow ( GetDerivativeOfX2 ( s ), 2 ) );
        }

        static Vector<double> GetBoudaryFunction (double t) {
            return Vector.Build.DenseOfArray ( new[] { GetX1 ( t ), GetX2 ( t ) } );
        }

        static double GetL(double s, double tau ) {
            return GetL1 ( s, tau ) * Math.Log ( 4 / Math.E * Math.Pow ( Math.Sin ( ( s - tau ) / 2 ), 2 ) ) + GetL2 ( s, tau );
        }

        static double GetL2 (double s, double tau ) {
            if ( s - tau < 0.0001 ) {
                return 1 / 2 * Math.Log ( ( 4 * Math.Abs ( Math.Sin ( s ) ) ) / ( Math.E * GetDerivativeOfR ( s ) ) );
            }
            return 1 / 2 * Math.Log ( ( 4 * Math.Abs ( Math.Sin ( s ) ) ) / ( Math.E * Math.Pow ( GetParametricR ( s, tau ), 2 ) ) );
        }

        static double GetL1 ( double s, double tau ) {
            return -1 / 2;
        }
    }
}
