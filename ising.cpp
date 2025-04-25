/* These is the new version
 *
 * by Eduardo Diaz Suarez
 *
*/
#include <iostream>
#include <sstream>
#include <iomanip>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <math.h>

#include "rng2.h"

long seed = (unsigned int)time(NULL);
//long seed2 = seed*(-1);


template<typename T>
T fromString(const std::string& s) {
  std::istringstream is(s);
  T t;
  is >> t;
  return t;
}

using namespace std;

class IsingSimul
{
    public:
        IsingSimul(int L)
        {
            srand((unsigned int)time(NULL));
            //int seed = (unsigned int)time(NULL);
            //int seed2 = seed*(-1);
            seed3= (-1)*time(NULL);
            mcs=0;
            lat_x  = L;
            lat_y  = L;
            this->spin_matrix_ = new double *[lat_y];
            for(int i = 0; i < lat_y; i++)
                spin_matrix_[i]= new double [lat_x];
        }
        ~IsingSimul()
        {
            for (int i = 0; i < lat_x; ++i)
              delete [] spin_matrix_[i];
            delete [] spin_matrix_;
        }

        void make_random_matrix_2d(void)
        {
            for(int i = 0; i < lat_y; i++)
            {
                for(int j = 0; j < lat_x; j++)
                {
                     double num = ran1(&seed3);
                     num >= 0.5 ? this->spin_matrix_[i][j] = 1.0 : this->spin_matrix_[i][j] = -1.0;
                }
            }
        }

        void show_()
        {
            for(int i = 0; i < lat_y; i++)
            {
                for(int j = 0; j < lat_x; j++)
                      cout<<"  "<<this->spin_matrix_[i][j];
                cout<<endl;
            }
        }

        double One_monte_carlo_step_per_spin( void );
        void precalc_w(double H, double T, double J);
        double spin_avg(void);
        double inter_energy(void);
        long MCS_(void){return mcs;}

    private:
        double ** spin_matrix_;
        double H;
        double J;
        long seed3;
        long mcs;
        int lat_x;
        int lat_y;

        int flip_spin(void);

        double wc[17][3];
};

using namespace std;

void run_metropolis_mc(double H, double Ti,double Tf, double dT ,double J, int L, unsigned int steps);

int main(int cant, char *args[])
{

    if(cant == 8 )
    {
        double H = fromString<double>(string(args[1]));
        double ti = fromString<double>(string(args[2]));
        double tf = fromString<double>(string(args[3]));
        double dt = fromString<double>(string(args[4]));

        int J = fromString<int>(string(args[5]));
        int L = fromString<int>(string(args[6]));
        unsigned int MC_steps_ = fromString<unsigned int>(string(args[7]));
        //run_monte_carlo_2d_dT(0,ti,tf,dt,1.0,Lx,Ly,MC_steps_);
        run_metropolis_mc(H,ti,tf,dt,J,L,MC_steps_);
        return 0;
    }
    cout<<"Ising 2D model, usage:   "<<endl;
    cout<<args[0]<<" H  Ti Tf dt J L MCsteps"<<endl;
    cout<<"Where H is the external magnetic field, J=1 is a ferromagetic sistem"<<endl;
    cout<<"Lx=Ly=L define the Heigth and Long of the 2D lattice, and "<<endl;
    cout<<"MCsteps are the number of Monte Carlos steps, we recomend 10^6"<<endl;
    cout<<"It is exclude the 10% of the initials values, to do the last"<<endl;
    cout<<"calculations parameters"<<endl<<endl;
    cout<<"sample: (where H=0,Ti=1.0, Tf=5.0, dt=0.1,  J=1, Lx=Ly=10, MCsteps=10^6)  "<<endl;
    cout<<args[0]<<" 0 1.0 5.0 0.1 1.0 10 1000000"<<endl;
    return -1;    
}

void run_metropolis_mc(double H, double Ti,double Tf, double dT ,double J, int L, unsigned int steps)
{
    //precalculate values

    cout<<"## Simulating Ising 2D  "<<L<<"x"<<L<<"  with J= "<<J<<endl;
    cout<<"##  H= "<<H<<",   From T= "<<Ti<<" to "<<Tf<<" K"<<endl;
    cout<<"#Tempt(K)   magnetization(m)    Susceptibility   Average_Energy_per_site    Specific_head"<<endl;
    for (double T=Ti;T<=Tf;T+=dT)
    {
        IsingSimul test_(L);
        test_.make_random_matrix_2d();
        test_.precalc_w(H,T,J);

        double M1_=0; double M2_ =0; double E_1=0; double E_2 = 0;
        unsigned int MCstep_10p =steps/10; // se desprecia el 10 % de los step MC
        unsigned int nsteps = steps-MCstep_10p;

        for(unsigned int MCstep_tmp=0;MCstep_tmp<MCstep_10p;MCstep_tmp++)
            test_.One_monte_carlo_step_per_spin();
        for (unsigned int MCstep=0;MCstep<nsteps;MCstep++)
        {
            test_.One_monte_carlo_step_per_spin();
            double magnet_total = test_.spin_avg();
            double E_t = test_.inter_energy();
            M1_ += fabs(magnet_total); M2_ += magnet_total*magnet_total; E_1 += E_t; E_2 += E_t*E_t;
        }
        //show the result!!!!!!!!!!!!!!!!!!!!!!
        double N = (double (L))*(double (L));
        M1_ = M1_/double (nsteps); M2_ =M2_/double (nsteps); E_1 = E_1/ double (nsteps);
        E_2 = E_2/double (nsteps);
        cout<<setw(4)<<T<<setw(15)<<M1_<<setw(15)<<(M2_  -  M1_*M1_)*N/T<<setw(15)<<E_1/N<<setw(15)<<(E_2 - E_1*E_1)/(T*T*N)<<endl;
    }
}


int IsingSimul::flip_spin(void)
{
    int Ly = this->lat_y;
    int Lx = this->lat_x;

    int x= int ( Lx * ran1(&(this->seed3)));
    int y= int ( Ly * ran1(&(this->seed3)));

    int xprev, xnext, yprev, ynext;
    if (x==0)    xprev =Lx-1; else xprev =x-1;
    if (x==Lx-1) xnext =0;    else xnext =x+1;
    if (y==0)    yprev =Ly-1; else yprev =y-1;
    if (y==Ly-1) ynext =0;    else ynext =y+1;

    //sums s values of the 4 neighbords
    int  s_neighbord = int (this->spin_matrix_[x][yprev] + this->spin_matrix_[x][ynext] +
                         this->spin_matrix_[xprev][y] + this->spin_matrix_[xnext][y] );
    //compare and change using the possible value calculate
    int delta_ss = 2*( int (this->spin_matrix_[x][y]))*( s_neighbord );
    //double delta_ss = 2.0*( spins[i][j])*(double ( s_neighbord) );

    double delta_E = (double(delta_ss))*J +  H*(this->spin_matrix_[x][y]);
    if (delta_E <= 0){this->spin_matrix_[x][y] = (-1.0)* ( this->spin_matrix_[x][y]); return 1;}

    double w_ = 0;
    w_ = this->wc[delta_ss+8][1 + ( int (this->spin_matrix_[x][y]))];

    double r =  ran1(&(this->seed3));

    if (r < w_)  // then change the configuration spin for s(i,j)
    {
        this->spin_matrix_[x][y] = (-1.0)* (this->spin_matrix_[x][y]);
        return 1;
    } else
    return 0;
}

double IsingSimul::One_monte_carlo_step_per_spin()
{
    int  N=( this->lat_x )* (this->lat_y);
    int accepts=0;
    for(int i=0;i<N;i++)
        if ( (this->flip_spin())==1 )
            ++accepts;
    double acceptanceRatio = double (accepts)/double (N);
    return acceptanceRatio;
    mcs++;
}


void IsingSimul::precalc_w(double H, double T, double J)
{
    this->H = H;
    this->J = J;
    for (int i = -8;i<=8;i+=4)
    {
        this->wc[i+8][0] = exp(-1.0*((double (i))*J + 2.0*H )/T);
        this->wc[i+8][1] = exp(-1.0*((double (i))*J )/T);
        this->wc[i+8][2] = exp(-1.0*((double (i))*J - 2.0*H )/T);
    }
}

double IsingSimul::spin_avg()
{
    double S=0;
    double N= double ( (this->lat_x ) * ( this->lat_y ) );
    for(int i = 0; i < this->lat_x; i++)
    {
        for(int j = 0; j < this->lat_y; j++)
        {
            S += this->spin_matrix_[i][j];
        }
    }
    return S/N;
}

double IsingSimul::inter_energy()
{
    double  sum =0;
    double s_sum=0;
    int Lx = this->lat_x;
    int Ly = this->lat_y;
    for(int i = 0; i < Lx; i++)
    {
        for(int j = 0; j < Ly; j++)
        {
            s_sum += this->spin_matrix_[i][j];
            //boundary conditions
            int i_next = (i == Lx - 1 ? 0 : i+1);
            int j_next = (j == Ly - 1 ? 0 : j+1);
            sum += (this->spin_matrix_[i][j])*(this->spin_matrix_[i_next][j] + this->spin_matrix_[i][j_next]);
        }
    }
    return -((this->J)*sum + (this->H)*s_sum);
}

//eof ising.cpp
