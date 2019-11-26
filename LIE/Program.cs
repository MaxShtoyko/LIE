using System;
using MathNet;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace LIE {
    class Program {
        static double PI = Math.PI;
        static double q = 5;
        static double n = 4;

        static void Main ( string[] args ) {

            var densities = QuadratureMethod ( );
            Console.WriteLine ( densities );

            Console.ReadKey ( );
        }


        static Vector<double> QuadratureMethod ( ) {
            var matrix = Matrix<double>.Build.Dense ( ( int ) n - 1, ( int ) n - 1 );
            var b = Vector<double>.Build.Dense ( ( int ) n - 1 );

            for (int k = 1;k<n;k++ ) {
                for(int i = 1;i<n;i++ ) {
                    matrix[k - 1, i - 1] = GetQuadratureR ( Math.Abs ( i - k ) ) + GetQuadratureR ( i + k ) + 1 / n * GetL2 ( GetSi ( k ), GetSi ( i ) );
                }
                b[k - 1] = GetFunction ( GetSi ( k ) );
            }

            Console.WriteLine ( matrix );

            return matrix.Solve ( b );
        }



        //Uex = x1^2 - x2^2
        static double GetExactResult (double t ) {
            return Math.Pow ( GetX1 ( t ), 2 ) - Math.Pow ( GetX2 ( t ), 2 );
        }

        //static double GetApproximateResult(double t ) {

        //}

        static double GetSi ( double j ) {
            return j * PI / n;
        }

        // Ri(s) = -1/2 { ... }
        static double GetQuadratureR ( double j ) {

            double r = 0;

            r += 1;

            for ( double m = 1; m < n - 1; m++ ) {
                r += 2 / m * Math.Cos ( m * GetSi ( j ) );
            }

            r += Math.Pow ( -1, j ) / Math.Pow ( n, 2 );

            return r * -1 / n;
        }

        //f 
        static double GetFunction (double s ) {
            return Math.Cos ( GetX1 ( s ) + GetX2 ( s ) );
        }

        static double V ( double s ) {
            return ( 1 / q - PI ) * Math.Pow ( ( PI - s ) / PI, 3 ) - 1 / q * ( PI - s ) / PI + PI;
        }

        static double W ( double s ) {
            return 2 * PI * Math.Pow ( V ( s ), q ) / ( Math.Pow ( V ( s ), q ) + Math.Pow ( V ( 2 * PI - s ), q ) );
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
        static double GetR(double s, double tau ) {
            return Math.Sqrt ( Math.Pow ( ( GetX1 ( s ) - GetX1 ( tau ) ), 2 ) + Math.Pow ( ( GetX2 ( s ) - GetX2 ( tau ) ), 2 ) );
        }

        // |x'(s)|
        static double GetDerivativeOfR ( double s ) {
            return Math.Sqrt ( Math.Pow ( GetDerivativeOfX1 ( s ), 2 ) + Math.Pow ( GetDerivativeOfX2 ( s ), 2 ) );
        }

        static Vector<double> GetParametricFunction (double t) {
            return Vector.Build.DenseOfArray ( new[] { GetX1 ( t ), GetX2 ( t ) } );
        }

        static double GetL(double s, double tau ) {
            return GetL1 ( s, tau ) * Math.Log ( 4 / Math.E * Math.Pow ( Math.Sin ( ( s - tau ) / 2 ), 2 ) ) + GetL2 ( s, tau );
        }

        static double GetL2 (double s, double tau ) {
            if ( s - tau < 0.0001 ) {
                return 1 / 2 * Math.Log ( ( 4 * Math.Abs ( Math.Sin ( s ) ) ) / ( Math.E * GetDerivativeOfR ( s ) ) );
            }
            return 1 / 2 * Math.Log ( ( 4 * Math.Abs ( Math.Sin ( s ) ) ) / ( Math.E * Math.Pow ( GetR ( s, tau ), 2 ) ) );
        }

        static double GetL1 ( double s, double tau ) {
            return -1 / 2;
        }
    }
}
