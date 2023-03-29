namespace OS_convert_test
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        static double Marc(double bf0, double n, double PHI0, double PHI)
        {
            /*Compute meridional arc.
            'Input: -  ellipsoid semi major axis multiplied by central meridian scale factor (bf0) in meters; _
            'n (computed from a, b and f0); _
            'lat of false origin (PHI0) and initial or final latitude of point (PHI) IN RADIANS.

            'THIS FUNCTION IS CALLED BY THE - _
            '"Lat_Long_to_North" and "InitialLat" FUNCTIONS
            'THIS FUNCTION IS ALSO USED ON IT'S OWN IN THE "Projection and Transformation Calculations.xls" SPREADSHEET
            */

            double to_return;

            to_return = bf0 * (((1 + n + ((5 / 4) * (Math.Pow(n, 2))) + ((5 / 4) * (Math.Pow(n, 3)))) * (PHI - PHI0))
            - (((3 * n) + (3 * (Math.Pow(n, 2))) + ((21 / 8) * (Math.Pow(n, 3)))) * (Math.Sin(PHI - PHI0)) * (Math.Cos(PHI + PHI0)))
            + ((((15 / 8) * (Math.Pow(n, 2))) + ((15 / 8) * (Math.Pow(n, 3)))) * (Math.Sin(2 * (PHI - PHI0))) * (Math.Cos(2 * (PHI + PHI0))))
            - (((35 / 24) * (Math.Pow(n, 3))) * (Math.Sin(3 * (PHI - PHI0))) * (Math.Cos(3 * (PHI + PHI0)))));

            return to_return;
        }

        static double Lat_Long_to_North(double PHI, double LAM, double a, double b, double e0, double n0, double f0, double PHI0, double LAM0)
        {
            /*Project Latitude and longitude to Transverse Mercator northings
            'Input: - _
            ' Latitude (PHI) and Longitude (LAM) in decimal degrees; _
            ' ellipsoid axis dimensions (a & b) in meters; _
            ' eastings (e0) and northings (n0) of false origin in meters; _
            ' central meridian scale factor (f0); _
            'latitude (PHI0) and longitude (LAM0) of false origin in decimal degrees.

            'Convert angle measures to radians*/

            Double RadPHI, RadLAM, RadPHI0, RadLAM0;
            Double af0, bf0, e2, n, nu, rho, eta2, p, M, I, II, III, IIIA;
            double to_return;

            RadPHI = PHI * (Math.PI / 180);
            RadLAM = LAM * (Math.PI / 180);
            RadPHI0 = PHI0 * (Math.PI / 180);
            RadLAM0 = LAM0 * (Math.PI / 180);

            af0 = a * f0;
            bf0 = b * f0;

            e2 = ((Math.Pow(af0, 2)) - (Math.Pow(bf0, 2))) / (Math.Pow(af0, 2));
            n = (af0 - bf0) / (af0 + bf0);
            nu = af0 / (Math.Sqrt(1 - (e2 * (Math.Pow(Math.Sin(RadPHI), 2)))));

            rho = (nu * (1 - e2)) / (1 - (e2 * Math.Pow(Math.Sin(RadPHI), 2)));
            eta2 = (nu / rho) - 1;
            p = RadLAM - RadLAM0;
            M = Marc(bf0, n, RadPHI0, RadPHI);

            I = M + n0;
            II = (nu / 2) * (Math.Sin(RadPHI)) * (Math.Cos(RadPHI));
            III = ((nu / 24) * (Math.Sin(RadPHI)) * (Math.Pow(Math.Cos(RadPHI), 3))) * (5 - (Math.Pow(Math.Tan(RadPHI), 2)) + (9 * eta2));
            IIIA = ((nu / 720) * (Math.Sin(RadPHI)) * (Math.Pow(Math.Cos(RadPHI), 5))) * (61 - (58 * (Math.Pow(Math.Tan(RadPHI), 2))) + (Math.Pow(Math.Tan(RadPHI), 4)));

            to_return = I + ((Math.Pow(p, 2)) * II) + ((Math.Pow(p, 4)) * III) + ((Math.Pow(p, 6)) * IIIA);
            return to_return;
        }

        static double Lat_Long_to_East(double PHI,double  LAM, double a,double b,double e0,double  f0,double  PHI0, double LAM0)
        {
            /*Project Latitude and longitude to Transverse Mercator eastings.
            'Input: - _
            'Latitude (PHI) and Longitude (LAM) in decimal degrees; _
            'ellipsoid axis dimensions (a & b) in meters; _
            'eastings of false origin (e0) in meters; _
            'central meridian scale factor (f0); _
            'latitude (PHI0) and longitude (LAM0) of false origin in decimal degrees.

            'Convert angle measures to radians*/

            double RadPHI, RadLAM, RadPHI0, RadLAM0;
            double af0, bf0, e2, n, nu, rho, eta2, p, IV, V, VI;
            double to_return;

            RadPHI = PHI * (Math.PI / 180);
            RadLAM = LAM * (Math.PI / 180);
            RadPHI0 = PHI0 * (Math.PI / 180);
            RadLAM0 = LAM0 * (Math.PI / 180);

            af0 = a * f0;
            bf0 = b * f0;

            e2 = (Math.Pow((Math.Sin(RadPHI)) , 2));

            e2 = ((Math.Pow(af0, 2)) - (Math.Pow(bf0, 2))) / (Math.Pow(af0, 2));
            n = (af0 - bf0) / (af0 + bf0);
            nu = af0 / (Math.Sqrt(1 - (e2 * (Math.Pow((Math.Sin(RadPHI)), 2)))));
            rho = (nu * (1 - e2)) / (1 - (e2 * (Math.Pow((Math.Sin(RadPHI)), 2))));
            eta2 = (nu / rho) - 1;
            p = RadLAM - RadLAM0;

            IV = nu * (Math.Cos(RadPHI));
            V = (nu / 6) * (Math.Pow((Math.Cos(RadPHI)),3)) * ((nu / rho) - ((Math.Pow(Math.Tan(RadPHI), 2) )));
            VI = (nu / 120) * (Math.Pow(Math.Cos(RadPHI), 5)) * (5 - (18 * (Math.Pow(Math.Tan(RadPHI), 2))) + (Math.Pow(Math.Tan(RadPHI), 4)) + (14 * eta2) - (58 * (Math.Pow(Math.Tan(RadPHI), 2)) * eta2));
            to_return = e0 + (p * IV) + ((Math.Pow(p,3)) * V) + ((Math.Pow(p,5)) * VI);

            return to_return;
        }

        static void latlong_toOS (double DecimalLat, double DecimalLong, out double Eastings, out double northings)
        {
            double PHI, PHI0, LAM, LAM0, E0, N0, a, b, F0;
            //constants
            PHI0 = 49;
            LAM0 = -2;
            E0 = 400000;
            N0 = -100000.0;
            a = 6377563.396;
            b = 6356256.91;
            F0 = 0.9996012717;
            PHI = DecimalLat;
            LAM = DecimalLong;

            Eastings = Lat_Long_to_East(PHI, LAM, a, b, E0, F0, PHI0, LAM0);
            northings = Lat_Long_to_North(PHI, LAM, a, b, E0, N0, F0, PHI0, LAM0);
        }
        private void button1_Click(object sender, EventArgs e)
        {
            // see also https://github.com/dstl/osgb/tree/master/src/main/java/uk/gov/dstl/geo/osgb
            // and https://www.ordnancesurvey.co.uk/documents/resources/guide-coordinate-systems-great-britain.pdf

            double DecimalLat, DecimalLong, Eastings, Northings;
            string temp1;

            DecimalLat = 52; //Convert.ToDouble(TXT_LATITUDE.Text);
            DecimalLong = 0;  //Convert.ToDouble(TXT_LONGITUDE.Text);
            latlong_toOS(DecimalLat, DecimalLong, out Eastings, out Northings);
            TXT_EASTINGS.Text = Convert.ToString(Eastings);
            TXT_NORTHINGS.Text=Convert.ToString(Northings);
        }
    }
}