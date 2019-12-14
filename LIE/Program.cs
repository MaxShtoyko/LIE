using System;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace LIE {
    class Program {
        static double PI = Math.PI;
        static double n = 2;
        static Vector<double> y = Vector.Build.DenseOfArray ( new double[] { 6, 6 } );

        static void Main ( string[] args ) {

            //Console.WriteLine ( "Для першої границі: " );

            //Console.WriteLine ( $"q = {ParametrizationHelper.q}" );
            //PrintResult ( 1.0, 0.0 );
            //PrintResult ( 0.75, 0.25 );

            //ParametrizationHelper.q = 5;
            //Console.WriteLine ( $"q = {ParametrizationHelper.q}" );
            //PrintResult ( 1.0, 0.0 );
            //PrintResult ( 0.75, 0.25 );

            Console.WriteLine ( "Для другої границі: " );

            Console.WriteLine ( $"q = {ParametrizationHelper.q}" );
            PrintResult ( 0.0, 0.1 );
            PrintResult ( -0.25, 0.5 );

            Console.WriteLine ( $"\n\n\nq = {ParametrizationHelper.q}" );
            PrintResult ( 0.0, 0.1 );
            PrintResult ( -0.25, 0.5 );

            Console.ReadKey ( );
        }


        static void PrintResult (double x1, double x2 ) {
            n = 2;

            Console.WriteLine ( $"\n(x1, x2) = ({ Math.Round ( x1, 4 )},{ Math.Round ( x2, 4 )})" );

            for ( int j = 0; j < 6; j++ ) {
                n *= 2;

                var densities = QuadratureMethod ( );

                var point = Vector.Build.DenseOfArray ( new double[] { x1, x2 } );
                var error = GetExactResult ( point ) - GetApproximateResult ( densities, point );

                Console.WriteLine ( $"n: {n}" );
                Console.WriteLine ( $"Error : {Math.Abs ( error )}" );
            }
        }
        //static Vector<double> QuadratureMethod ( ) {
        //    var matrix = Matrix<double>.Build.Dense ( 2 * ( int ) n - 1, 2 * ( int ) n - 1 );
        //    var b = Vector<double>.Build.Dense ( 2 * ( int ) n - 1 );

        //    for ( int k = 1; k < 2 * n; k++ ) {
        //        for ( int i = 1; i < 2 * n; i++ ) {
        //            matrix[k - 1, i - 1] = GetL1 ( ) * R ( S ( k ), i ) + 1.0 / ( 2 * n ) * B ( S ( i ), S ( k ) );
        //        }
        //        b[k - 1] = GetParametricFunction ( S ( k ) );
        //    }

        //    return matrix.Solve ( b );
        //}

        //static double GetApproximateResult ( Vector<double> densities, Vector<double> t ) {
        //    double result = 0;

        //    for ( int k = 1; k < 2 * n; k++ ) {
        //        result += densities[k - 1] * Math.Log ( GetR ( t - Z ( S ( k ) ) ) );
        //    }

        //    return result * -0.5 / n;
        //}

        static Vector<double> QuadratureMethod ( ) {
            var matrix = Matrix<double>.Build.Dense ( 2 * ( int ) n, 2 * ( int ) n );
            var b = Vector<double>.Build.Dense ( 2 * ( int ) n );

            for ( int k = 0; k < 2 * n; k++ ) {
                for ( int i = 0; i < 2 * n; i++ ) {
                    matrix[k, i] = GetL1 ( ) * R ( S ( k ), i ) + 1.0 / ( 2 * n ) * B ( S ( i ), S ( k ) );
                }
                b[k] = GetParametricFunction ( S ( k ) );
            }

            return matrix.Solve ( b );
        }

        static double GetApproximateResult ( Vector<double> densities, Vector<double> t ) {
            double result = 0;

            for ( int k = 0; k < 2 * n; k++ ) {
                result += densities[k] * Math.Log ( GetR ( t - Z ( S ( k ) ) ) );
            }

            return result * -0.5 / n;
        }

        static double B ( double s, double tau ) {
            if ( Math.Abs ( s - tau ) < 0.00001 ) {
                return 0.5 * Math.Log ( 1 / ( Math.E * Math.Pow ( GetR ( GetDerivativeOfZ ( s ) ), 2 ) ) );
            }

            var v = 4 * Math.Pow ( Math.Sin ( ( s - tau ) / 2 ), 2 );
            var vv = ( Math.E * Math.Pow ( GetR ( Z ( s ) - Z ( tau ) ), 2 ) );

            return 0.5 * Math.Log ( v / vv );
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

        //static double GetExactResult ( Vector<double> x ) {
        //    return 1;
        //}

        //static double GetFunction ( Vector<double> x ) {
        //    return 1;
        //}

        static double GetExactResult ( Vector<double> x ) {
            return 1.0 / ( 2.0 * PI ) * Math.Log ( 1.0 / ( GetR ( x - y ) ) );
        }

        static double GetFunction ( Vector<double> x ) {
            return 1.0 / ( 2.0 * PI ) * Math.Log ( 1.0 / ( GetR ( x - y ) ) );
        }

        static double GetParametricFunction ( double s ) {
            return GetFunction ( Z ( s ) );
        }

        static double S ( double j ) {
            return (j * PI) / n;
        }

        static Vector<double> Z ( double s ) {
            return GetBoudaryFunction ( ParametrizationHelper.Y ( s ) );
        }

        static Vector<double> GetDerivativeOfZ ( double s ) {
            return GetDerivativeOfBoudaryFunction ( ParametrizationHelper.Y ( s ) ) * ParametrizationHelper.Ydx ( s );
        }

        //static Vector<double> Z ( double s ) {
        //    return GetBoudaryFunction ( s );
        //}

        //static Vector<double> GetDerivativeOfZ ( double s ) {
        //    return GetDerivativeOfBoudaryFunction ( s );
        //}

        // |x(s)|
        static double GetR ( Vector<double> s ) {
            return Math.Sqrt ( Math.Pow ( s[0], 2 ) + Math.Pow ( s[1], 2 ) );
        }

        static Vector<double> GetBoudaryFunction ( double t ) {
            return Vector.Build.DenseOfArray ( new[] { GetX1 ( t ), GetX2 ( t ) } );
        }

        static Vector<double> GetDerivativeOfBoudaryFunction ( double t ) {
            return Vector.Build.DenseOfArray ( new[] { GetDerivativeOfX1 ( t ), GetDerivativeOfX2 ( t ) } );
        }

        static double GetL1 ( ) {
            return -0.5;
        }



        /////ПЕРША ГРАНИЦЯ
        //static double GetX1 ( double t ) {
        //    return 4.0 / Math.Sqrt ( 3 ) * Math.Sin ( t / 2 );
        //}

        ////EXAMPLE
        //static double GetX2 ( double t ) {
        //    return -2.0*Math.Sin ( t );
        //}

        ////EXAMPLE
        //static double GetDerivativeOfX1 ( double t ) {
        //    return 2.0 / Math.Sqrt ( 3.0 ) * Math.Cos ( t / 2 );
        //}

        ////EXAMPLE
        //static double GetDerivativeOfX2 ( double t ) {
        //    return -2.0*Math.Cos ( t );
        //}

        //ДРУГА ГРАНИЦЯ
        static double GetX1 ( double t ) {
            return -2.0 / 3.0 * Math.Sin ( 3 * t / 2 );
        }

        //EXAMPLE
        static double GetX2 ( double t ) {
            return -Math.Sin ( t );
        }

        //EXAMPLE
        static double GetDerivativeOfX1 ( double t ) {
            return -Math.Cos ( 3 * t / 2 );
        }

        //EXAMPLE
        static double GetDerivativeOfX2 ( double t ) {
            return -Math.Cos ( t );
        }


        ////EXAMPLE
        //static double GetX1 ( double t ) {
        //    return 2 * Math.Cos ( t );
        //}

        ////EXAMPLE
        //static double GetX2 ( double t ) {
        //    return 4 * Math.Sin ( t );
        //}

        ////EXAMPLE
        //static double GetDerivativeOfX1 ( double t ) {
        //    return -2 * Math.Sin ( t );
        //}

        ////EXAMPLE
        //static double GetDerivativeOfX2 ( double t ) {
        //    return 4 * Math.Cos ( t );
        //}
    }
}
