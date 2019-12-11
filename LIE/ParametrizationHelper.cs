using System;
namespace LIE {
    public class ParametrizationHelper {
        static double q = 7;
        static double PI = Math.PI;

        public static double Y ( double s ) {
            return W ( s );
        }

        private static double V ( double s ) {

            return ( 1.0 / q - PI / 2.0 ) * Math.Pow ( ( PI - s ) / PI, 3 ) - ( 1.0 / q ) * ( ( PI - s ) / PI ) + PI / 2;
        }

        private static double W ( double s ) {

            return 2 * PI * Math.Pow ( V ( s ), q ) / ( Math.Pow ( V ( s ), q ) + Math.Pow ( V ( 2 * PI - s ), q ) );
        }
    }
}
