#include <stdio.h>
#include <iostream>
#include <random>
#include <math.h>
#include <ctime>
#include <list>
#include <functional>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
using namespace std::placeholders;
using namespace std;
double k00=0.0;
double d=10;
double Pi = 3.14159265358979323846;

double gsl_sf_hyperg_2F1( double a, double b, double c, double x )
{
	const double TOLERANCE = 10e-5;
	double term = a * b * x / c;
	double value = 1.0 + term;
	int n = 1;
	while ( abs( term ) > TOLERANCE )
	{
		a++, b++, c++, n++;
		term *= a * b * x / c / n;
		value += term;
	}
	return value;
}


double normal_dist(double mean, double stdev)
{
	auto const seed = rand();
	std::default_random_engine engin {seed};
	std::normal_distribution<double> distribution(mean,stdev);
	double result = distribution(engin);
	return result;
}


std::vector<double> brownian(std::vector<double> x, double t)
{
	int siz= d;                     std::vector<double> y;
	for (int i=0;i<siz;i++)
	{
		y.push_back(x[i]+ sqrt(t) * normal_dist(0.0,1.0));
	}
	return y;
}


double uniform()
{
	static long s1 = 55555;
	static long s2 = 99999;
	static double factor = 1.0/2147483563.0;
	register long k,z;
	k = s1/53668;
	s1 = 40014*(s1%53668)-k*12211;
	if (s1<0)
		s1 += 2147483563;
	k = s2/52774;
	s2 = 40692*(s2%52774)-k*3791;
	if (s2<0)
		s2 += 2147483399;
	z = (s1-2147483563)+s2;
	if (z<1)
		z += 2147483562;
	return(((double)(z))*factor);
};


double uniform_dist(double a, double b)
{
	auto const seed = rand();
	std::default_random_engine engin {seed};
	std::uniform_real_distribution<double> distribution(a,b);
	double result = distribution(engin);
	return result;
}


double exponential_dist(double lambda)
{
	auto const seed = rand();
	std::default_random_engine engin(seed);
	std::exponential_distribution<double> distribution(lambda);
	double result = distribution(engin);
	return result;
}


double sum=0;
int count;
double subordinator(double x, double a, double e, double t)
{
	double Pi = 3.14159265358979323846;
	double lambda_0 = - Pi/2;
	double s=10;
	{
		double W = -log(uniform());
		double U = uniform_dist(-Pi/2, Pi/2);
		double A1 = sin((a/2)*(U-lambda_0)) / pow(cos(U),2/a);
		double A2 = pow(cos(U-(a/2)*(U-lambda_0))/W,(2-a)/a);
		s=x + pow(t,2/a) * A1 * A2;
	}
	if (s==0) {printf("%.80f\n",s);}
	return s;
}


double gamma_dist(double alpha, double beta)
{
	auto const seed = rand();
	std::default_random_engine engin {seed};
	std::gamma_distribution<double> distribution(alpha, beta);
	double result = distribution(engin);
	return result;
}


double gamma_pdf(double x, double a)
{
	return (pow(x,a-1)*exp(-x)/tgamma(a));
}


double min(double x, double y)
{
	if (x<y) return x; else return y;
}


double max(double x, double y)
{
	if (x<y) return y; else return x;
}


double norme_2 (std::vector<double> x)
{
	double r=0;
	int siz = d;                    for (int i=0;i<siz;i++)
	{
		r+=x[i]*x[i];
	}
	return sqrt(r);
}


double initial_condition(double a, std::vector<double> x)
{
	if (norme_2(x)>1)
	{
		return 0;
	}
	else
	{
		return pow(1-norme_2(x)*norme_2(x),k00+a/2.0);
	}
}


double constant(double a)
{
	return (- pow(2,a) * tgamma((d+a)/2.0)*tgamma(k00+1+a/2)/tgamma(d/2.0)/tgamma(k00+1));
}


double constanttilde(double a)
{
	return (- pow(2,a) * tgamma((d+a)/2.0)*tgamma(k00+1+a/2)/tgamma(-a/2)/tgamma(k00+1+(d+a)/2.0));
}


double sol_funct(double a, std::vector<double> x)
{
	if (norme_2(x)>1)
	{
		return -constanttilde(a)*gsl_sf_hyperg_2F1((d+a)/2.0,(2+a)/2.0,k00+1+(d+a)/2.0,1/norme_2(x)/norme_2(x))*pow(norme_2(x),-a-d);
	}
	else
	{
		return -constant(a)*gsl_sf_hyperg_2F1((d+a)/2.0,-k00,d/2.0,norme_2(x)*norme_2(x));
	}
}


double coefficient_0(double a, double t, std::vector<double> x)
{
	if (norme_2(x)>1)
	{
		return exp(-t)*sol_funct(a,x);
	}
	else
	{
		if (k00>0) return exp(-t)*sol_funct(a,x)-exp(-4*t)*pow(1-norme_2(x)*norme_2(x),4*k00+2*a);
		else return -exp(-t)*constant(a) - exp(-4*t) * pow(1-norme_2(x)*norme_2(x),4*k00+2*a);
	}
}


double integral(double(*f)(double alpha, double temps), double a, double x0, double x1, int n)
{
	double step = (x1 - x0)/n;
	double area = 0.0;
	for(int i = 0; i < n; ++i)
	{
		area += f(a,x0 + (i +0.5) * step) * step ;
	}
	return area;
}


double function_r0(double R, std::vector<double> x, std::vector<double> y)
{
	std::vector<double> result;
	result.reserve(x.size());
	std::transform(x.begin(),x.end(),y.begin(), std::back_inserter(result),std::minus<double>());
	return (pow(R,2)-pow(norme_2(x),2)) * (pow(R,2)-pow(norme_2(y),2) )/(pow(R,2)*pow(norme_2(result),2));
}


double integrand(double a, double t)
{
	return pow(t,a/2.0-1)/pow(t+1,d/2.0);
}


double green(double R, std::vector<double> x,std::vector<double> y, double a)
{
	double kappa = tgamma(d/2.0)/ (pow(2,a) * pow(Pi,d/2.0)* pow(tgamma(a/2.0),2.0));
	std::vector<double> result;
	result.reserve(x.size());
	std::transform(x.begin(),x.end(),y.begin(), std::back_inserter(result),std::minus<double>());
	double G = pow(norme_2(result),a-d) *integral(*integrand, a, 0.0 , function_r0(R,x,y), 5000);
	return kappa * G;
}


std::vector<double> nabla_G(double R, std::vector<double> x,std::vector<double> y, double a)
{
	double kappa = tgamma(d/2.0)/ (pow(2,a) * pow(Pi,d/2.0)* pow(tgamma(a/2.0),2.0));
	std::vector<double> result;
	result.reserve(x.size());
	std::transform(x.begin(),x.end(),y.begin(), std::back_inserter(result),std::minus<double>());
	std::vector<double> z(d);
	double N = norme_2(result);
	for (int i=0; i <d;i++)
	{
		z[i] = kappa * (2*(a-d)*(x[i]-y[i])*green(R,x,y,a)/N
			- (2*(x[i]-y[i])*function_r0(R,x,y)/pow(N,2) + 2*x[i] * (pow(R,2) -pow(norme_2(y),2))/(pow(R,2)*pow(N,2))) * pow(function_r0(R,x,y),a/2.0+1)/pow(function_r0(R,x,y)+1,d/2.0) );
	}
	return z;
}


std::vector<double> weight_green(double R, std::vector<double> x, std::vector<double> y, double a)
{
	std::vector<double> z = nabla_G(R,x,y,a);
	std::transform(z.begin(), z.end(), z.begin(),
		std::bind(std::multiplies<double>(), std::placeholders::_1, 1/green(R,x,y,a)));
	return z;
}


double poisson_kernel(double R, std::vector<double> x, std::vector<double> y, double a)
{
	double kappa= tgamma(d/2.0)*sin(Pi*a/2.0)/pow(Pi,d/2.0+1);
	std::vector<double> result;
	result.reserve(x.size());
	std::transform(x.begin(),x.end(),y.begin(), std::back_inserter(result),std::minus<double>());
	return kappa* pow( (pow(R,2)-pow(norme_2(x),2))/(pow(norme_2(y),2) - pow(R,2) ) , a/2.0) * pow(norme_2(result),-d);
}


double scalar_product(std::vector<double> x, std::vector<double> y)
{
	return std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
}


std::vector<double> coeff_b0(double R, std::vector<double> x)
{
	std::transform(x.begin(), x.end(), x.begin(),
		std::bind(std::multiplies<double>(), std::placeholders::_1, 1-pow(norme_2(x),2)));
	return x;
}


double solution_gradient_nonlin_elliptic_sample( double mu, double a, double e, std::vector<double> x,double S)
{
	int ifault;
	double W = gamma_dist(mu,1.0);
	double U = uniform_dist(0.0,1.0);
	int siz = d;                    std::vector<double> y = x;
	double dt=W/500/(1+0.5*sqrt(norme_2(x)));
	double t=0;
	while (t<W)
	{
		if (norme_2(y) > 1)
		{
			return (0.0);
		}
		y = brownian(y,subordinator(0,a,e,pow(2,a/2)*(dt)));
		t=t+dt;
	}
	if (norme_2(y) > 1)
	{
		return (0.0);
	}
	if (U < 0.5)
	{
		return 2*S*(sol_funct(a,y)+ pow((2*k00+a),2)*pow(norme_2(y),2*2)*pow(1-pow(norme_2(y),2), 2*k00+2*a/2.0)) /gamma_pdf(W,mu) ;
	}
	else
	{
		double u1 = solution_gradient_nonlin_elliptic_sample(mu, a,e, y, 1);
		return solution_gradient_nonlin_elliptic_sample(mu, a,e, y, 2*S*u1*scalar_product(coeff_b0(1,x),weight_green(1,x,y,a)) /gamma_pdf(W,mu));
	}
}


double solution_gradient_nonlin_elliptique(double mu, double a, double e, std::vector<double> x, int n)
{
	double result = 0;
	double v;
	int m=1;
	for (int i=0; i <n;i++)
	{
		double truc = max(-2000,min(2000,solution_gradient_nonlin_elliptic_sample(mu, a, e, x,1)));
		result += truc;
	}
	printf("result=%f\n",result/n);
	return (result/n);
}


int main()
{
	srand(0);
	double alpha = 1.75;
	double mu = 1.55;               double T = 1.0;
	double t = 0.9;
	double ee=0.15;
	double e=0;
	std::vector<double> x(d);
	FILE *datafile;
	datafile=fopen("datafile0.9","w");
	FILE *datacurve;
	datacurve=fopen("datacurve","w");
	int ifault;
	clock_t c1, c2;
	time_t rawtime;
	struct tm * timeinfo;
	int i,j,k;
	for (k=0;k<=20;k++)
	{
		x[0]=-1+k*0.1;for (j=1;j<d;j++){x[j]=0;}
		c1 = clock();
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		printf ( "Current time: %s", asctime (timeinfo) );
		printf("x=%f\n",x[0]);
		double y = solution_gradient_nonlin_elliptique(mu,alpha,0,x,5000000);
		c2 = clock();
		printf("x\t\t y\n");
		printf("%f\t %f \t Exact=%f\n",x[0],y,pow(1-min(norme_2(x)*norme_2(x),1),k00+alpha/2));
		fprintf(datafile,"%f\t%f\t%f\n",x[0],y,pow(1-min(norme_2(x)*norme_2(x),1),k00+alpha/2));
		printf("Time taken: %g minutes\n\n", (c2 - c1) / 60 / (double)CLOCKS_PER_SEC);
	}
	fclose(datafile);
	fclose(datacurve);
	return 0;
}
