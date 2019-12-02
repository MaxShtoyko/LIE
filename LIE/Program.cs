using System;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Factorization;

namespace LIE {
    class Program {
        static double PI = Math.PI;
        static double n = 4;
        static Vector<double> y = Vector.Build.DenseOfArray ( new double[] { 5, 5 } );

        static void Main ( string[] args ) {

            double generalError = 0;

            for ( int j = 0; j < 6; j++ ) {
                n *= 2;

                var densities = QuadratureMethod ( );

                for ( int i = 1; i < n; i++ ) {
                    var x1 = i /(2*n)-0.25;
                    var x2 = i /(2*n)-0.25;
                    var point = Vector.Build.DenseOfArray ( new double[] { x1, x2 } );
                    var error = GetExactResult ( point ) - GetApproximateResult ( densities, point );
                    generalError += Math.Abs ( error );

                    //Console.WriteLine ( $"Exact: {GetExactResult ( point ) } ----- Approximate {GetApproximateResult ( densities, point )}" );
                    //Console.WriteLine ( $"(x1,x2) =({Math.Round ( x1, 4 )},{Math.Round ( x2, 4 )}) ---- Error : {Math.Abs ( error )}" );
                }

                Console.WriteLine ( $"n={n}, error={generalError / n}" );

                generalError = 0;
            }

            Console.ReadKey ( );
        }

        static Vector<double> QuadratureMethod ( ) {
            var matrix = Matrix<double>.Build.Dense ( ( int ) n - 1, ( int ) n - 1 );
            var b = Vector<double>.Build.Dense ( ( int ) n - 1 );

            for ( int k = 1; k < n; k++ ) {
                for ( int i = 1; i < n; i++ ) {
                    matrix[k-1, i-1] = GetL1 ( ) * R ( S ( k ), i ) + 1 / ( 2 * n ) * GetL2 ( S ( k ), S ( i ) );
                }
                b[k - 1] = GetParametricFunction ( S ( k ) );
            }

            return matrix.Solve ( b );
        }

        static double R ( double t, int j ) {

            double r = 0;

            r += 1;

            for ( double m = 1; m < n; m++ ) {
                r += 2.0 / m * Math.Cos ( m * ( t - S ( j ) ) );
            }

            r += Math.Cos ( n * ( t - S ( j ) ) ) / n;

            return ( r * -0.5 ) / n;
        }

        static double GetExactResult ( Vector<double> x ) {
            return 1 / ( 2 * PI ) * Math.Log ( 1 / ( GetR ( x - y ) ) );
        }

        static double GetFunction ( Vector<double> x ) {
            return 1 / ( 2 * PI ) * Math.Log ( 1 / ( GetR ( x - y ) ) );
        }


        static double GetParametricFunction ( double s ) {
            return GetFunction ( Z ( s ) );
        }

        static double GetApproximateResult ( Vector<double> densities, Vector<double> t ) {
            double result = 0;

            for(int k = 1; k<n; k++ ) {
                result += densities[k - 1] * Math.Log ( GetR ( t - Z ( S ( k ) ) ) );
            }

            return -0.5 * result / n;
        }

        static double S ( double j ) {
            return (j * PI) / n;
        }

        static Vector<double> Z ( double s ) {
            return GetBoudaryFunction ( ParametrizationHelper.Y ( s ) );
        }

        //EXAMPLE
        static double GetX1 ( double t ) {
            return 4 * Math.Cos ( t );
        }

        //EXAMPLE
        static double GetX2 ( double t ) {
            return 4 * Math.Sin ( t ) * Math.Cos ( t );
        }

        //EXAMPLE
        static double GetDerivativeOfX1 ( double t ) {
            return -4 * Math.Sin ( t );
        }

        //EXAMPLE
        static double GetDerivativeOfX2 ( double t ) {
            return 4 * Math.Cos ( 2*t );
        }

        // |x(s)|
        static double GetR ( Vector<double> s ) {
            return Math.Sqrt ( Math.Pow ( s[0], 2 ) + Math.Pow ( s[1], 2 ) );
        }

        // |x'(s)|
        static double GetDerivativeOfR ( Vector<double> s ) {
            return Math.Sqrt ( Math.Pow ( GetDerivativeOfX1 ( s[0] ), 2 ) + Math.Pow ( GetDerivativeOfX2 ( s[1] ), 2 ) );
        }

        static Vector<double> GetBoudaryFunction (double t) {
            return Vector.Build.DenseOfArray ( new[] { GetX1 ( t ), GetX2 ( t ) } );
        }

        static double GetL2 ( double s, double tau ) {
            if ( Math.Abs ( s - tau ) < 0.0000001 ) {
                return 0.5 * Math.Log ( Math.Abs ( Math.Sin ( s ) ) / ( Math.E * Math.Pow ( GetDerivativeOfR ( Z ( s ) ), 2 ) ) );
            }

            var v = 4 * Math.Pow ( Math.Sin ( ( tau - s ) / 2 ), 2 );
            var vv = Math.E * Math.Pow ( GetR ( Z ( s ) - Z ( tau ) ), 2 );

            return 0.5 * Math.Log ( v / vv );
        }

        static double B ( double s, double tau ) {
            if ( Math.Abs ( s - tau ) < 0.00001 ) {
                return 0.5* Math.Log ( ( 2 * Math.Abs ( Math.Sin ( s ) ) ) / ( Math.E * GetDerivativeOfR ( Z ( s ) ) ) );
            }

            var v = 2 * Math.Abs ( Math.Cos ( s ) - Math.Cos ( tau ) );
            var vv = ( Math.E * GetR ( Z ( s ) - Z ( tau ) ) );

            return Math.Log ( v / vv );
        }

        static double GetL1 (  ) {
            return -0.5;
        }

        //static Vector<double> QuadratureMethod ( ) {
        //    var matrix = Matrix<double>.Build.Dense ( ( int ) n - 1, ( int ) n - 1 );
        //    var b = Vector<double>.Build.Dense ( ( int ) n - 1 );

        //    for (int k = 1;k<n;k++ ) {
        //        for(int i = 1;i<n; i++ ) {
        //            matrix[k - 1, i - 1] = GetQuadratureR ( Math.Abs ( i - k ) ) + GetQuadratureR ( i + k ) + 1 / n * B ( S ( k ), S ( i ) );
        //        }
        //        b[k - 1] = 2*GetParametricFunction ( S ( k ) );
        //    }

        //    Console.WriteLine ( matrix );

        //    return matrix.Solve ( b );
        //}

        //static double GetQuadratureR ( double j ) {

        //    double r = 0;

        //    r += 1;

        //    for ( double m = 1; m < n; m++ ) {
        //        r += 2.0 / m * Math.Cos ( m * S ( j ) );
        //    }

        //    r += Math.Pow ( -1, j ) / Math.Pow ( n, 2 );

        //    return r * -1.0 / n;
        //}

        //static double GetX1 ( double t ) {
        //    return Math.Sin ( t / 2 );
        //}

        ////EXAMPLE
        //static double GetX2 ( double t ) {
        //    return Math.Sin ( t );
        //}

        ////EXAMPLE
        //static double GetDerivativeOfX1 ( double t ) {
        //    return 1 / 2 * Math.Cos ( t / 2 );
        //}

        ////EXAMPLE
        //static double GetDerivativeOfX2 ( double t ) {
        //    return Math.Cos ( t );
        //}

        //static double B(double s, double tau ) {
        //    if ( Math.Abs ( s - tau ) < 0.00001 ) {
        //        return Math.Log ( ( 2 * Math.Abs ( Math.Sin ( s ) ) ) / ( Math.E * GetDerivativeOfR ( Z ( s ) ) ) );
        //    }

        //    var v = 2 * Math.Abs ( Math.Cos ( s ) - Math.Cos ( tau ) );
        //    var vv = ( Math.E * GetR ( Z ( s ) - Z ( tau ) ) );

        //    return Math.Log ( v / vv );
        //}
    }
}
