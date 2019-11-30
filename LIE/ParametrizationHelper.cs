using System;
namespace LIE {
    public class ParametrizationHelper {
        static double q = 5;
        static double PI = Math.PI;

        public static double Y ( double s ) {
            if ( s <= PI ) {
                return W ( s );
            }
            return PI + W ( s - PI );
        }

        private static double V ( double s ) {
            return ( 1.0 / q - PI / 2.0 ) * Math.Pow ( ( PI - 2.0 * s ) / PI, 3 ) - ( 1.0 / q ) * ( PI - 2.0 * s ) / PI + PI / 2;
        }

        private static double W ( double s ) {
            return PI * Math.Pow ( V ( s ), q ) / ( Math.Pow ( V ( s ), q ) + Math.Pow ( V ( PI - s ), q ) );
        }
    }
}
