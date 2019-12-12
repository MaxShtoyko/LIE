using System;
namespace LIE {
    public class ParametrizationHelper {
        static double q = 3.0;
        static double PI = Math.PI;

        public static double Y ( double x ) {
            return W ( x );
        }

        private static double W ( double x ) {
            return 2.0 * PI * Math.Pow ( x, q ) / ( Math.Pow ( x, q ) + Math.Pow ( 2.0 * PI - x , q ) );
        }

        public static double Ydx ( double x ) {
            return Wdx ( x );
        }

        private static double Wdx ( double x ) {
            return 2.0 * PI * ( q * Math.Pow ( x, q - 1 ) * Math.Pow ( 2 * PI - x, q - 1 ) * ( ( 2.0 * PI - x ) + x ) )
                / Math.Pow ( Math.Pow ( x, q ) + Math.Pow ( 2.0 * PI - x , q ), 2.0 );
        }
    }
}
