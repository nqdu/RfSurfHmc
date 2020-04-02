#include<random>
#include<iostream>
using std:: default_random_engine;
using std:: normal_distribution;
using std:: uniform_real_distribution;
using std:: uniform_int_distribution;

static default_random_engine e(11);
static normal_distribution<double> n(0,1);
static uniform_real_distribution<double> u(0,1);

/*
    generate random numbers with standard normal distribution
*/
double gaussian()
{
    return n(e);
}

double uniform(double a=0.0,double b=0.0)
{
    return a + (b-a) * u(e);
}

int uniformi(int a,int b){
    static uniform_int_distribution<int> ui(a,b);
    return ui(e);
}
/*
int main(){
    for(int i = 0; i <5;i++)
        std::cout << gaussian() << std::endl;
    for(int i = 0; i <5;i++)
        std::cout << gaussian() << std::endl;    
    return 0;
}
*/