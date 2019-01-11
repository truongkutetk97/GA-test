using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.Runtime.InteropServices;

namespace Test_app_1
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {

        [DllImport("test_Jan_11_12_44.dll")]

        public static extern void DisplayHelloFromDLL();
        
            [DllImport("test_Jan_11_12_44.dll", CallingConvention =CallingConvention.Cdecl)]
        public static extern int Addd(int a, int b);
        [DllImport("test_Jan_11_12_44.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern int subb(int a, int b);
        [DllImport("test_Jan_11_12_44.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern int maiin(int b);

        public MainWindow()
        {
            InitializeComponent();
        }

        private void ss(object sender, RoutedEventArgs e)
        {

           DisplayHelloFromDLL();
            int c = Addd(3, 5);
            int j = 100;
            MessageBox.Show(maiin(j).ToString());

        }
    }
    
}
