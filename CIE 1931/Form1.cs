using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace CIE_1931
{
    public partial class Form1 : Form
    {      
        public Form1()
        {
            InitializeComponent();
        }

        private void pictureBox1_Paint(object sender, PaintEventArgs e)
        {
            /*
            double step = 0.01;
            for (double X = 0; X < 1; X+=step)
            {
                for (double Y = 0; Y < 1; Y+=step)
                {
                    for (double Z = 0; Z < 1; Z+=step)
                    {
                        double[] xy = new double[2];
                        xy[0] = X / (X + Y + Z);
                        xy[1] = Y / (X + Y + Z);
                        double z = 1 - xy[0] - xy[1];
                        double[] XYZ = { X, Y, Z };
                        double[] lms = xyz2lms(XYZ);
                        double[] rgb = xyz2rgb(XYZ);
                        //rgb[0] = rgb[0] > 1 ? 1 : rgb[0];
                        //rgb[0] = rgb[0] < 0 ? 0 : rgb[0];
                        //rgb[1] = rgb[1] > 1 ? 1 : rgb[1];
                        //rgb[1] = rgb[1] < 0 ? 0 : rgb[1];
                        //rgb[2] = rgb[2] > 1 ? 1 : rgb[2];
                        //rgb[2] = rgb[2] < 0 ? 0 : rgb[2];
                        if (checkrgb(rgb))
                        {
                            Color color = Color.FromArgb((byte)(rgb[0] * 255), (byte)(rgb[1] * 255), (byte)(rgb[2] * 255));
                            e.Graphics.FillRectangle(new SolidBrush(color), (float)xy[0] * pictureBox1.Width, pictureBox1.Height - (float)xy[1] * pictureBox1.Height, 5, 5);
                        }
                        if (rgbisNaN(rgb))
                        {
                            //Console.WriteLine("Red: " + rgb[0] + ", Green: " + rgb[1] + ", Blue: " + rgb[2] + ", x: " + xyz[0] + ", y:" + xyz[1]);
                            Color color = Color.FromArgb(155, 155, 155);
                            e.Graphics.FillRectangle(new SolidBrush(color), (float)xy[0] * pictureBox1.Width, pictureBox1.Height - (float)xy[1] * pictureBox1.Height, 2, 2);
                        }
                    }
                }
            }*/
            /*
            double step = 0.005;
            double Y = 1;
            for (double i = 0; i < 1; i += step)
            {
                for (double j = 0; j < 1; j += step)
                {
                        double X = (Y / j) * i;
                        double Z = (Y / j) * (1 - j - i);
                        double[] XYZ = { X, Y, Z };
                        double[] rgb = xyz2rgb(XYZ);
                        if (checkrgb(rgb))
                        {
                            Color color = Color.FromArgb((int)(rgb[0] * 255), (int)(rgb[1] * 255), (int)(rgb[2] * 255));
                            e.Graphics.FillRectangle(new SolidBrush(color), (float)i * pictureBox1.Width, pictureBox1.Height - (float)j * pictureBox1.Height, 4, 4);
                        }
                        
                        else
                        {
                            Color color = Color.FromArgb(150, 150, 150);
                            e.Graphics.FillRectangle(new SolidBrush(color), (float)i * pictureBox1.Width, pictureBox1.Height - (float)j * pictureBox1.Height, 4, 4);
                        }
                }
            }
            */
            double step = 0.05;         
            for (double r = 0; r < 1; r += step)
            {
                progressBar1.Value = (int)(r * 100);
                for (double g = 0; g < 1; g += step)
                {
                    for (double b = 0; b < 1; b += step)
                    {
                        double X, Y, Z;
                        convolution(380, 780, out X, out Y, out Z, r, g, b);
                        double x = X / (X + Y + Z);
                        double y = Y / (X + Y + Z);
                        double z = 1 - x - y;
                        double[] XYZ = { X, Y, Z };
                        //double[] lms = xyz2lms(XYZ);
                        //double[] rgb = xyz2rgb(XYZ);
                        double[] rgb = { r, g, b };
                        if (checkrgb(rgb))
                        {
                            Color color = Color.FromArgb((byte)(rgb[0] * 255), (byte)(rgb[1] * 255), (byte)(rgb[2] * 255));
                            e.Graphics.FillRectangle(new SolidBrush(color), (float)x * pictureBox1.Width, pictureBox1.Height - (float)y * pictureBox1.Height, 5, 5);
                        }
                        if (rgbisNaN(rgb))
                        {
                            //Console.WriteLine("Red: " + rgb[0] + ", Green: " + rgb[1] + ", Blue: " + rgb[2] + ", x: " + xyz[0] + ", y:" + xyz[1]);
                            Color color = Color.FromArgb(155, 155, 155);
                            e.Graphics.FillRectangle(new SolidBrush(color), (float)x * pictureBox1.Width, pictureBox1.Height - (float)y * pictureBox1.Height, 2, 2);
                        }
                    }
                }
            }
            
            for (double lambda = 380; lambda < 640; lambda += step)
            {
                double X, Y, Z;
                convolution2(380, 780, out X, out Y, out Z, lambda);
                double x = X / (X + Y + Z);
                double y = Y / (X + Y + Z);
                double z = 1 - x - y;
                
                double hue = (650 - lambda) * 1.08;
                Color color = ColorFromHSV(hue, 1, 1);
                e.Graphics.FillRectangle(new SolidBrush(color), (float)x * pictureBox1.Width, pictureBox1.Height - (float)y * pictureBox1.Height, 5, 5);
            }
        }
        public static Color ColorFromHSV(double hue, double saturation, double value)
        {
            int hi = Convert.ToInt32(Math.Floor(hue / 60)) % 6;
            double f = hue / 60 - Math.Floor(hue / 60);

            value = value * 255;
            int v = Convert.ToInt32(value);
            int p = Convert.ToInt32(value * (1 - saturation));
            int q = Convert.ToInt32(value * (1 - f * saturation));
            int t = Convert.ToInt32(value * (1 - (1 - f) * saturation));

            if (hi == 0)
                return Color.FromArgb(255, v, t, p);
            else if (hi == 1)
                return Color.FromArgb(255, q, v, p);
            else if (hi == 2)
                return Color.FromArgb(255, p, v, t);
            else if (hi == 3)
                return Color.FromArgb(255, p, q, v);
            else if (hi == 4)
                return Color.FromArgb(255, t, p, v);
            else
                return Color.FromArgb(255, v, p, q);
        }
        private double Gauss(double x, double mu, double sigma1, double sigma2)
        {
            double g;
            if (x < mu)
            {
                g = Math.Exp(-0.5 * Math.Pow(x - mu, 2) / Math.Pow(sigma1, 2));
            }
            else
            {
                g = Math.Exp(-0.5 * Math.Pow(x - mu, 2) / Math.Pow(sigma2, 2));
            }
            return g;
        }        
        private void convolution(int la, int lb, out double X, out double Y, out double Z, double r, double g, double b)
        {
            X = 0; Y = 0; Z = 0;
            int N = 100;
            double delta = (lb - la) / N;
            double xBar, yBar, zBar;
            //double xMax = xBarMax(N, lb, la);
            //double yMax = yBarMax(N, lb, la);
            //double zMax = zBarMax(N, lb, la);
            for (double lambda = la; lambda < lb; lambda+=delta)
            {
                xBar = xBarF(lambda);
                yBar = yBarF(lambda);
                zBar = zBarF(lambda);

                X += xBar * radiance(lambda, r, g, b) * delta;
                Y += yBar * radiance(lambda, r, g, b) * delta;
                Z += zBar * radiance(lambda, r, g, b) * delta;
            }
            //X *= (lb - la) / N / xMax;
            //Y *= (lb - la) / N / yMax;
            //Z *= (lb - la) / N / zMax;
            //X: 541.199999999999, Y: 442.8, Z: 759.199999999999
            X /= 541.199999999999;
            Y /= 442.8;
            Z /= 759.199999999999;
        }
        private void convolution2(int la, int lb, out double X, out double Y, out double Z, double mu)
        {
            X = 0; Y = 0; Z = 0;
            int N = 100;
            double delta = (lb - la) / N;
            double xBar, yBar, zBar;
            //double xMax = xBarMax(N, lb, la);
            //double yMax = yBarMax(N, lb, la);
            //double zMax = zBarMax(N, lb, la);
            for (double lambda = la; lambda < lb; lambda += delta)
            {
                xBar = xBarF(lambda);
                yBar = yBarF(lambda);
                zBar = zBarF(lambda);
                
                X += xBar * laser(lambda, mu) * delta;
                Y += yBar * laser(lambda, mu) * delta;
                Z += zBar * laser(lambda, mu) * delta;
            }
            //X *= (lb - la) / N / xMax;
            //Y *= (lb - la) / N / yMax;
            //Z *= (lb - la) / N / zMax;
            //X: 541.199999999999, Y: 442.8, Z: 759.199999999999
            X /= 541.199999999999;
            Y /= 442.8;
            Z /= 759.199999999999;
        }        
        private double xBarMax(int N, double lb, double la)
        {
            double sum = 0;
            double delta = (lb - la) / N;
            for (int i = 0; i < N; i++)
            {
                sum += Math.Pow(xBarF(i * delta + la), 2);
            }
            return sum * (lb - la) / N;
        }
        private double yBarMax(int N, double lb, double la)
        {
            double sum = 0;
            double delta = (lb - la) / N;
            for (int i = 0; i < N; i++)
            {
                sum += Math.Pow(yBarF(i * delta + la), 2);
            }
            return sum * (lb - la) / N;
        }
        private double zBarMax(int N, double lb, double la)
        {
            double sum = 0;
            double delta = (lb - la) / N;
            for (int i = 0; i < N; i++)
            {
                sum += Math.Pow(zBarF(i * delta + la), 2);
            }
            return sum * (lb - la) / N;
        }
        private double xBarF(double lambda) 
        {
            double first = 1.056 * Gauss(lambda, 599.8, 37.9, 31.0);
            double second = 0.362 * Gauss(lambda, 442.0, 16.0, 26.7);
            double third = 0.065 * Gauss(lambda, 501.1, 20.4, 26.2);
            double result = first + second - third;
            return result;
        }
        private double yBarF(double lambda)
        {
            return 0.821 * Gauss(lambda, 568.8, 46.9, 40.5) + 0.286 * Gauss(lambda, 530.9, 16.3, 31.1);
        }
        private double zBarF(double lambda)
        {
            return 1.217 * Gauss(lambda, 437.0, 11.8, 36.0) + 0.681 * Gauss(lambda, 459.0, 26.0, 13.8);
        }
        private double laser(double lambda, double mu)
        {
            return Gauss(lambda, mu, 5, 5);
        }
        private double radiance(double lambda, double r, double g, double b)
        {
            double r0, g0, b0;
            r0 = 612; //SRGB
            g0 = 549; //SRGB
            b0 = 465; //SRGB
            double ds = 5;

            double result = r * Gauss(lambda, r0, ds, ds)
                + g * Gauss(lambda, g0, ds, ds)
                + b * Gauss(lambda, b0, ds, ds);
            //Random rnd = new Random();
            //result = rnd.NextDouble();
            //result = 1;
            return result;
        }
        private bool checkxyz(double[] xyz)
        {
            return xyz[0] < 0 || xyz[1] < 0 || xyz[2] < 0 || xyz[0] > 1 || xyz[1] > 1 || xyz[2] > 1;
        }

        private bool checkrgb(double[] rgb)
        {
            if (rgb[0] < 0 || rgb[0] > 1 || rgb[1] < 0 || rgb[1] > 1 || rgb[2] < 0 || rgb[2] > 1)
            {
                return false;
            }
            return true;
        }
        private bool rgbisNaN(double[] rgb)
        {
            if (Double.IsNaN(rgb[0]) || Double.IsNaN(rgb[1]) || Double.IsNaN(rgb[2]))
            {
                return true;
            }
            return false;
        }
        private double[] xyz2lms(double[] xyz)
        {
            double[] lms = new double[3];
            lms[0] = (xyz[0] * 0.38971) + (xyz[1] * -0.6889) + (xyz[2] * -0.07868);
            lms[1] = (xyz[0] * -0.22981) + (xyz[1] * 1.1834) + (xyz[2] * 0.04641);
            lms[2] = (xyz[0] * 0.0) + (xyz[1] * 0.0) + (xyz[2] * 1.0);
            return lms;
        }

        private double[] xyz2rgb(double[] xyz)
        {
            double[] rgb = new double[3];
            rgb[0] = (xyz[0] * 3.2406) + (xyz[1] * -1.5372) + (xyz[2] * -0.4986);
            rgb[1] = (xyz[0] * -0.9689) + (xyz[1] * 1.8758) + (xyz[2] * 0.0415);
            rgb[2] = (xyz[0] * 0.0557) + (xyz[1] * -0.2040) + (xyz[2] * 1.0570);
            return rgb;
        }

        private double[] rgbtoxyz(double[] rgb)
        {
            double[] xyz = new double[3];
            xyz[0] = (rgb[0] * 0.412453) + (rgb[1] * 0.357580) + (rgb[2] * 0.180423);
            xyz[1] = (rgb[0] * 0.212671) + (rgb[1] * 0.715160) + (rgb[2] * 0.072169);
            xyz[2] = (rgb[0] * 0.019334) + (rgb[1] * 0.119193) + (rgb[2] * 0.950227);
            return xyz;
        }

        private void pictureBox1_Resize(object sender, EventArgs e)
        {
            pictureBox1.Invalidate();
        }
    }
}
