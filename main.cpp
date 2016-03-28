#include <iostream>
#include <fstream>
#include <cmath>
#include <Eigen/Dense>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;


typedef std::vector<double> state_type;

class DoublePendulum
{
private:
    const double m1 = 0.5, m2 = 0.5, L1 = 0.1, L2 = 0.1, g = 9.81;

public:
    DoublePendulum() = default;
    void operator() (const state_type &y, state_type &dy, const double /* t */)
    {
        Eigen::MatrixXd M(2, 2), C(2, 1), B(2,1);

        double delta(y[0] - y[2]);

        // advanced initialization
        M <<  (m1 + m2)*L1, m2*L2*cos(delta), m2*L1*cos(delta),  m2*L2;

        C <<  m2*L1*L2*y[3]*y[3]*sin(delta) + g*(m1 + m2)*sin(y[0]),
                -m2*L1*y[1]*y[1]*sin(delta) + m2*g*sin(y[2]);

        //#####################( ODE Equations )################################
        dy[0] = y[1];
        dy[2] = y[3];

        B = M.inverse() * (-C);

        dy[1] = B(0,0);
        dy[3] = B(1,0);
    }
};



int main()
{
    const double dt = 0.01;
    runge_kutta_dopri5 < state_type > stepper;
    double pi = boost::math::constants::pi<double>();

    state_type y(4);
    // initial values
    y[0] = pi/3.0;  //  Theta 1
    y[1] = 0.0;     // dTheta 1
    y[2] = pi/4.0; //   Theta 2
    y[3] = 0.0;    //  dTheta 2

    DoublePendulum dp;

    string output_file;
    ofstream myfile;

    cout << "Please enter the name of the output file: " << endl;
    cin >> output_file;

    myfile.open(output_file + ".txt");

    for (double t(0.0); t <= 10; t += dt){

        myfile << fixed << setprecision(3);
        myfile << "t: " << t << "  Theta 1: " << y[0] << "   dTheta 1: " << y[1]
        << "  Theta 2: " << y[2] << "  dTheta 2: " << y[3] << endl;

        stepper.do_step(dp, y, t, dt);
    }

    myfile.close();

    cout << "The data has been saved to the output file." << endl;
    return 0;
}