//============================================================================
// Emanuel A. Fronhofer & Florian Altermatt
// Classical metapopulation dynamics and eco-evolutionary feedbacks in dendritic networks
// Ecography
// 2016
//
// procedures
//
//============================================================================

/*
	Copyright (C) 2016  Emanuel A. Fronhofer

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.
	
*/


const gsl_rng *gBaseRand;

//________________________________________________________________________________________
//------------------------------------------------------Initialize Random Number Generator

void specify_rng(unsigned long randSeed)
{
	gBaseRand = gsl_rng_alloc(gsl_rng_rand);

	srand(randSeed);
	unsigned long r = rand();
	gsl_rng_set(gBaseRand, r);
}

//________________________________________________________________________________________
//-------------------------------------------------------------------------Simplifications

//-------------------------------------------------Simplify Random Drawing between 0 and 1

double ran()
{
	return gsl_rng_uniform(gBaseRand);
}

//---------------------------------------------------------------Simplify Gaussian Randoms

double gauss(double sd)
{
	return gsl_ran_gaussian(gBaseRand,sd);
}

//-----------------------------------------------------------------Simplify Poisson Random

int poisson(double sd)
{
	return gsl_ran_poisson(gBaseRand,sd);
}

//-------------------------------------------------------------Simplify exponential Random
double expo(double sd)
{
	return gsl_ran_exponential(gBaseRand, 1/sd);
}

//---------------------------------------------------------------Simplify Lognormal Random

double lognorm(double zeta, double sigma)
{
	double var;										//variance of resulting randoms
	double mu,s_d;										//mean and sd of lognormal distr.
												//to be calculated by mean and
												//sigma of resulting randoms
	var = sigma*sigma;
	s_d = sqrt(log((var/(2*zeta))+1));
	mu = log(zeta)-0.5*(s_d*s_d);
	return gsl_ran_lognormal(gBaseRand,mu,s_d);
}

//---------------------------------------------------------------Simplify mean calculation
// for arrays
double mean(double data[], int n)
{
	return gsl_stats_mean(data,1,n);
}

double mean(float data[], int n)
{
	double data2[n];
	for (int i = 0; i < n; ++i) {
		data2[i] = data[i];
	}
	return gsl_stats_mean(data2,1,n);
}

//for vectors
double mean(vector<double> vec, int n)
{
	double data[n];
	for(int i = 0; i<n; i++) {
		data[i] = vec.at(i);
	}
	double m = gsl_stats_mean(data,1,n);
	return(m);
}

//-----------------------------------------------------Simplify integer
//median calculation

double median(vector<double> vec, int n)
{
	double data[n];
	for(int i = 0; i<n; i++) {
		data[i] = vec.at(i);
	}
	gsl_sort(data,1,n);
	double medi = gsl_stats_median_from_sorted_data(data,1,n);
	return(medi);
}

//-------------------------------------------------------------Simplify standard deviation
//for arrays
double sd(double data[], int n)
{
	return gsl_stats_sd(data,1,n);
}

//for vectors
double sd(vector<double> vec, int n)
{
	double data[n];
	for(int i = 0; i<n; i++) {
		data[i] = vec.at(i);
	}

	double s = gsl_stats_sd(data,1,n);
	return(s);
}

//------------------------------------------------------------- quartile

double uquant(vector<double> vec, int n)
{
	double data[n];
	for(int i = 0; i<n; i++) {
		data[i] = vec.at(i);
	}
	gsl_sort(data,1,n);
	double q = gsl_stats_quantile_from_sorted_data(data,1,n,0.75);
	return(q);
}

double lquant(vector<double> vec, int n)
{
	double data[n];
	for(int i = 0; i<n; i++) {
		data[i] = vec.at(i);
	}
	gsl_sort(data,1,n);
	double q = gsl_stats_quantile_from_sorted_data(data,1,n,0.25);
	return(q);
}

// ---------------------------------------------------- simplify maximum value calculation
double max(vector<double> data)
{
	return *max_element(data.begin(), data.end());
}

int max(vector<int> data)
{
	return *max_element(data.begin(), data.end());
}
