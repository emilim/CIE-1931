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
            double step = 0.01;
            for (double r = 0; r < 1; r += step)
            {
                progressBar1.Value = (int)(r * 100);
                for (double g = 0; g < 1; g += step)
                {
                    for (double b = 0; b < 1; b += step)
                    {
                        double X, Y, Z;
                        convolution(380, 780, out X, out Y, out Z, r, g, b);
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
            }
        }

        private double Gauss(double x, double my, double sigma1, double sigma2)
        {
            double g = 0;
            if (x < my)
            {
                g = Math.Exp(-(1 / 2) * Math.Pow(x - my, 2) / Math.Pow(sigma1, 2));
            }
            else if (x >= my)
            {
                g = Math.Exp(-(1 / 2) * Math.Pow(x - my, 2) / Math.Pow(sigma2, 2));
            }
            return g;
        }        
        private void convolution(int la, int lb, out double X, out double Y, out double Z, double r, double g, double b)
        {
            X = 0; Y = 0; Z = 0;
            int N = 100;
            double deltalambda = (lb - la) / N;
            for (double i = la; i < lb; i+=deltalambda)
            {
                double xlambda = 1.056 * Gauss(i, 599.8, 37.9, 31.0) + 0.362 * Gauss(i, 442.0, 16.0, 26.7) - 0.065 * Gauss(i, 501.1, 20.4, 26.2);
                double ylambda = 0.821 * Gauss(i, 568.8, 46.9, 40.5) + 0.286 * Gauss(i, 530.9, 16.3, 31.1);
                double zlambda = 1.217 * Gauss(i, 437.0, 11.8, 36.0) + 0.681 * Gauss(i, 459.0, 26.0, 13.8);

                X += xlambda * light(i, r, g, b);
                Y += ylambda * light(i, r, g, b);
                Z += zlambda * light(i, r, g, b);
            }
            X /= (lb - la);
            Y /= (lb - la);
            Z /= (lb - la);
        }
        private double light(double lambda, double r, double g, double b)
        {
            double r0, g0, b0;
            r0 = 612; //SRGB
            g0 = 549; //SRGB
            b0 = 465; //SRGB

            double result = r * Gauss(lambda, r0, 30, 30)
                + g * Gauss(lambda, g0, 30, 30)
                + b * Gauss(lambda, b0, 30, 30);
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
