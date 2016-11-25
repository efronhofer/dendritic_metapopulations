//============================================================================
// Emanuel A. Fronhofer & Florian Altermatt
// Classical metapopulation dynamics and eco-evolutionary feedbacks in dendritic networks
// Ecography
// 2016
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

#include <iostream>
#include <cstdlib>								//standard C library
#include <ctime>								//access system time library
#include <fstream>								//file streaming library
#include <string>								//string library included
#include <sstream>								//string streaming for reading numbers from

#include <vector>
#include <cmath>								//standard math library

#include <gsl/gsl_rng.h>						//gsl random number generator
#include <gsl/gsl_randist.h>					//gsl random distributions
#include <gsl/gsl_statistics.h>					//gsl some statistical methods
#include <gsl/gsl_statistics_double.h> 			//gsl some double statistical methods
#include <gsl/gsl_sort_double.h> 				//gsl sort double arrays

#include <algorithm>

using namespace std;

#include "include/procedures.h"					//procedure simplifications
#include "include/classes.h"					//class definitions

//_____________________________________________________________________________
//------------------------------------------------------------ global variables
unsigned int sim_time;															// actual time in simulation
int max_runs;																	// no. of repeats
float mut_sd;																	// variance used for mutations
float mut_rate;																	// mutation rate
int capacity;																	// habitat capacity
float lambda_null;																// fertility
float mu0;																		// dispersal mortality
float alpha;																	// Allee effect strength
float sigma;																	// environmental stochasticity
float epsilon;																	// random patch extinctions

TPatch world[WORLDDIM];															// simulated world
float connectivity_matrix[WORLDDIM][WORLDDIM];									// connectivity matrix

float rel_metapopsize;															// relative metapopulation size
float occupancy;																// metapopulation occupancy
float rel_emigrants;															// relative number of emigrants

int no_neutral_alleles;															// pop gen: no. of alleles at these loci
float mut_rate_neutral;															// pop gen. mutation rate at these loci

// counter for turnover events
int TO_help_occ_cnt[WORLDDIM];
float turnover_rate;

//_____________________________________________________________________________
//------------------------------------------------------------------ procedures

//------------------------------------------------------------- read parameters
void readParameters(){
	ifstream parinfile("input/parameters.in");							//parameter input file
	string buffer;
	istringstream is;

	getline(parinfile,buffer); getline(parinfile,buffer);
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> sim_time;																		//simulation time
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> max_runs;																		//no. of repeats
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> mut_sd;																		//variance used for mutations
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> mut_rate;																		//mutation rate
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> no_neutral_alleles;															//pop gen: no. of alleles at these loci
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> mut_rate_neutral;																//pop gen. mutation rate at these loci
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> capacity;																		//habitat capacity
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> lambda_null;																	//fertility
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> sigma;																		//uncorrelated evironmental stochasticity
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> alpha;																		//Allee effect strength
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> epsilon;																		//random patch extinction probability
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> mu0;																			//dispersal mortality

	parinfile.close();
}

//------------------------------------------------------------- read connectivity matrix from file
void readConnectivityMatrix(){
	ifstream matinfile("input/connectivity_matrix_global.in");							//connectivity matrix input file

	for (int x = 0; x < WORLDDIM; ++x) {
		for (int y = 0; y < WORLDDIM; ++y) {
			matinfile >> connectivity_matrix[x][y]; // this is [row][column] or [to][from]
		}
	}

	matinfile.close();

	// now get probabilities for dispersing from column x to row y
	// help var for column sums
	float col_sums[WORLDDIM];

	// calculate column sums to normalize
	for (int x = 0; x < WORLDDIM; ++x) {
		// reset value
		col_sums[x] = 0;
		for (int y = 0; y < WORLDDIM; ++y) {
			// sum up to column
			col_sums[x] = col_sums[x] + connectivity_matrix[y][x];
		}
	}

	// now normalize
	for (int x = 0; x < WORLDDIM; ++x) {
		for (int y = 0; y < WORLDDIM; ++y) {
			connectivity_matrix[y][x] = connectivity_matrix[y][x] / col_sums[x];
		}
	}
}

//------------------------------------------------------- initialize simulation
void Initialize(){
	// initialize patches and individuals in patches
	for (int x = 0; x < WORLDDIM; ++x) {

		//first: clear females and males vector
		world[x].females.clear();
		world[x].newFemales.clear();
		world[x].males.clear();
		world[x].newMales.clear();

		// initialize individuals in this patch
		// females
		for (int f = 0; f < round(capacity/2); ++f) {
			TIndiv newfemale;
			for (int a = 0; a < N_ALLELES; ++a) {
				newfemale.dispRate[a] = ran();

				for (int n = 0; n < NO_NEUTRAL_LOCI; ++n) {
					newfemale.neutral_loci[a][n] = floor(ran()*no_neutral_alleles);

				}
			}

			world[x].females.push_back(newfemale);
		}

		// males
		for (int m = 0; m < round(capacity/2); ++m) {
			TIndiv newmale;
			for (int a = 0; a < N_ALLELES; ++a) {
				newmale.dispRate[a] = ran();

				for (int n = 0; n < NO_NEUTRAL_LOCI; ++n) {
					newmale.neutral_loci[a][n] = floor(ran()*no_neutral_alleles);

				}
			}
			world[x].males.push_back(newmale);
		}

		// save initial state for turnover rate calculation
		if (world[x].females.size() + world[x].males.size() > 0) {
			// occupied
			// save state for next comparision
			TO_help_occ_cnt[x] = 1;
		} else {
			// uncoccupied
			// save state for next comparison
			TO_help_occ_cnt[x] = 0;
		}

	}
}


// ------------------------------------------------ analyze population dynamics
void Analyze(unsigned int acttime, int actrun){
	//reset metapopulation size and occupancy
	rel_metapopsize = 0;
	occupancy = 0;

	unsigned int numberoccupied = 0;
	unsigned int metapopsize = 0;

	for (int x = 0; x < WORLDDIM; ++x) {
		unsigned int localpopsize = world[x].females.size() + world[x].males.size();
		metapopsize += localpopsize;
		if (localpopsize > 0) {
			++numberoccupied;
		}
	}
	// calculate occupancy
	occupancy = float(numberoccupied) / float(WORLDDIM);
	rel_metapopsize = float(metapopsize) / float(WORLDDIM*capacity);
}

// ---------------------------------------------------- save individual results
void saveResults(int actrun){
	// output file: individuals
	stringstream outputindiv_path_stream;
	outputindiv_path_stream << "output/output_individuals_run" << actrun << ".out";
	string outputindiv_path = outputindiv_path_stream.str();
	ofstream outputindiv(outputindiv_path.c_str());

	// headers
	outputindiv << "x" << "    " << "sex" << "    " << "dispRate_1" << "    " << "dispRate_2" << "    " ;

	for (int nl = 0; nl < NO_NEUTRAL_LOCI; ++nl) {
		for (int a = 0; a < N_ALLELES; ++a) {
			outputindiv << "neutral_locus_" << a+1 << "_" << nl+1 << "    ";
		}
	}

	outputindiv << endl;

	for (int x = 0; x < WORLDDIM; ++x) {
		for (unsigned int f = 0; f < world[x].females.size(); ++f) {
			// write metapop results to file
			outputindiv << x << "    " << "f" << "    " << world[x].females.at(f).dispRate[0] << "    " << world[x].females.at(f).dispRate[1] << "    ";

			for (int nl = 0; nl < NO_NEUTRAL_LOCI; ++nl) {
				for (int a = 0; a < N_ALLELES; ++a) {
					outputindiv << world[x].females.at(f).neutral_loci[a][nl] << "    ";
				}
			}

			outputindiv << endl;

		}
		for (unsigned int m = 0; m < world[x].males.size(); ++m) {
			// write metapop results to file
			outputindiv << x << "    " << "m" << "    " << world[x].males.at(m).dispRate[0] << "    " << world[x].males.at(m).dispRate[1] << "    ";

			for (int nl = 0; nl < NO_NEUTRAL_LOCI; ++nl) {
				for (int a = 0; a < N_ALLELES; ++a) {
					outputindiv << world[x].males.at(m).neutral_loci[a][nl] << "    ";
				}
			}

			outputindiv << endl;
		}
	}
	// close indiv output file
	outputindiv.close();
}

// --------------------------------------------- find new patch during dispersal
int findNewPatch(int x){

	// the individual is in patch x, let's choose a patch newx to disperse to
	// this uses a weighted lottery and the connectivity matrix
	float help = 1-ran(); // ran() generates uniform random numbers between 0 and 1 with 0 without 1; as I don't want zero to be included I take 1-ran()
	int newx = x;
	float connectivity_matrix_sum = 0;
	for (int y = 0; y < WORLDDIM; ++y) {
		connectivity_matrix_sum += connectivity_matrix[y][x]; // here y,x sothat the matrix read from file goes from col x to row newx
		if (connectivity_matrix_sum >= help) {
			newx = y;
			break;
		}
	}
	return(newx);
}

// -------------------------------------------------------- dispersal procedure
void Dispersal(){

	unsigned int no_emigrants = 0;
	unsigned int metapopsize = 0;
	rel_emigrants = 0;
	int no_TO_events = 0;

	for (int x = 0; x < WORLDDIM; ++x) {

		// counter for metapopsize
		metapopsize += world[x].females.size() + world[x].males.size();

		// start with females
		for (unsigned int f = 0; f < world[x].females.size(); ++f) {
			// should the individual disperse?
			if (ran() < mean(world[x].females.at(f).dispRate, N_ALLELES)){
				// this individual will disperse
				// increase counter
				++no_emigrants;
				// check whether this emigrant survives the dispersal process
				if (ran() > mu0){
					// find new patch (global dispersal)
					int newPatch = findNewPatch(x);

					// copy disperser into new patch
					TIndiv Disperser = world[x].females.at(f);
					world[newPatch].newFemales.push_back(Disperser);
				}
				// delete emigrant from natal patch
				world[x].females.at(f) = world[x].females.back();
				world[x].females.pop_back();
				// decrease female loop counter
				--f;
			}
		}
		// continue with males
		for (unsigned int m = 0; m < world[x].males.size(); ++m) {
			// should the individual disperse?
			if (ran() < mean(world[x].males.at(m).dispRate, N_ALLELES)){
				// this individual will disperse
				// increase counter
				++no_emigrants;
				// check whether this emigrant survives the dispersal process
				if (ran() > mu0){
					// find new patch (global dispersal)
					int newPatch = findNewPatch(x);
					// copy disperser into new patch
					TIndiv Disperser = world[x].males.at(m);
					world[newPatch].newMales.push_back(Disperser);
				}
				// delete emigrant from natal patch
				world[x].males.at(m) = world[x].males.back();
				world[x].males.pop_back();
				// decrease female loop counter
				--m;
			}
		}
	}

	// now that dispersal is over, merge philopatrics and residents
	for (int x = 0; x < WORLDDIM; ++x) {

		// first copy the females
		for (unsigned int f = 0; f < world[x].newFemales.size(); ++f) {
			world[x].females.push_back(world[x].newFemales.at(f));
		}
		// erase the "old" immigrants from newFemales
		world[x].newFemales.clear();
		// then copy the males
		for (unsigned int m = 0; m < world[x].newMales.size(); ++m) {
			world[x].males.push_back(world[x].newMales.at(m));
		}
		// erase the "old" immigrants from newFemales
		world[x].newMales.clear();

		//-----------------------------
		// quantify turnover rates

		if (world[x].females.size() + world[x].males.size() > 0) {
			// occupied
			if (TO_help_occ_cnt[x] == 0) {
				no_TO_events = no_TO_events + 1;
			}

			// save state for next comparision
			TO_help_occ_cnt[x] = 1;
		} else {
			// uncoccupied
			if (TO_help_occ_cnt[x] == 1) {
				no_TO_events = no_TO_events + 1;
			}

			// save state for next comparison
			TO_help_occ_cnt[x] = 0;
		}
	}

	turnover_rate = float(no_TO_events) / float(WORLDDIM);
	rel_emigrants = float(no_emigrants) / float(metapopsize);

}

// ------------------------------------------------------------------ mutations
float mutate(float allele){
	if(ran()< mut_rate){
		float newallele = allele + gauss(mut_sd);
		return(newallele);
	} else {
		return(allele);
	}
}

//------------------------------------------------------------------- mutation for neutral loci is completely random
int mutate_neutral(int allele){
	if(ran()< mut_rate_neutral){
		float newallele = floor(ran()*no_neutral_alleles);;
		return(newallele);
	} else {
		return(allele);
	}
}

// ------------------------------------------------------------ larval survival
float larvalSurvival(unsigned int x){
	// following the Beverton-Holt model
	// this version is the original Poethke & Hovestadt SISP
	float a = (float(lambda_null) - 1) / float(capacity);
	if (a < 0) {a = 0;}

	int n_adults = world[x].females.size() + world[x].males.size();

	float sqr_dens = float(n_adults)/float(capacity) * float(n_adults)/float(capacity) ;	//calculation of square-density for allee-effect
	float survival = (sqr_dens/(sqr_dens+(alpha*alpha)))/(1+a*float(n_adults));				//survival-probability of newborns based on theory of logistic growth


	return(survival);
}

// --------------------------------------------------------------- reproduction
void Reproduction(){
	for (int x = 0; x < WORLDDIM; ++x) {

		// just to be sure: resize new females and males vectors
		world[x].newFemales.clear();
		world[x].newMales.clear();

		// for each patch check whether there are females and males
		if (world[x].females.size() > 0 && world[x].males.size() > 0) {
			float survival = larvalSurvival(x);
			float lambda_local = lognorm(lambda_null, sigma);
			// females choose their mates
			for (unsigned int f = 0; f < world[x].females.size(); ++f) {
				// randomly choose male
				unsigned int m = floor(ran()*world[x].males.size());
				// calculate number of offspring from this mating (without density regulation)
				int no_offspring = poisson(2*lambda_local*survival);
				// loop over offspring
				for (int o = 0; o < no_offspring; ++o) {
					// sex ratio is 0.5
					if (ran() < 0.5) {
						// females
						// initialize new individual
						TIndiv newOffspring;
						// inherit emigration rate allele from mother
						int allele = floor(ran()*N_ALLELES);
						newOffspring.dispRate[0] = mutate(world[x].females.at(f).dispRate[allele]);

						// dispersal and neutral loci are not linked, but neural loci are
						allele = floor(ran()*N_ALLELES);
						for (int n = 0; n < NO_NEUTRAL_LOCI; ++n) {
							newOffspring.neutral_loci[0][n] = mutate_neutral(world[x].females.at(f).neutral_loci[allele][n]);
						}

						// inherit emigration rate allele from father
						allele = floor(ran()*N_ALLELES);
						newOffspring.dispRate[1] = mutate(world[x].males.at(m).dispRate[allele]);

						// dispersal and neutral loci are not linked, but neural loci are
						allele = floor(ran()*N_ALLELES);
						for (int n = 0; n < NO_NEUTRAL_LOCI; ++n) {
							newOffspring.neutral_loci[1][n] = mutate_neutral(world[x].males.at(m).neutral_loci[allele][n]);
						}

						// add new individual to new females vector
						world[x].newFemales.push_back(newOffspring);

					} else {
						//males
						// initialize new individual
						TIndiv newOffspring;
						// inherit emigration rate allele from mother
						int allele = floor(ran()*N_ALLELES);
						newOffspring.dispRate[0] = mutate(world[x].females.at(f).dispRate[allele]);
						// dispersal and neutral loci are not linked, but neural loci are
						allele = floor(ran()*N_ALLELES);
						for (int n = 0; n < NO_NEUTRAL_LOCI; ++n) {
							newOffspring.neutral_loci[0][n] = mutate_neutral(world[x].females.at(f).neutral_loci[allele][n]);
						}
						// inherit dispersal rate allele from father
						allele = floor(ran()*N_ALLELES);
						newOffspring.dispRate[1] = mutate(world[x].males.at(m).dispRate[allele]);
						// dispersal and neutral loci are not linked, but neural loci are
						allele = floor(ran()*N_ALLELES);
						for (int n = 0; n < NO_NEUTRAL_LOCI; ++n) {
							newOffspring.neutral_loci[1][n] = mutate_neutral(world[x].males.at(m).neutral_loci[allele][n]);
						}
						// add new individual to new males vector
						world[x].newMales.push_back(newOffspring);
					}
				}
			}
		}
	}
}



// -------------------------------------------------- death of annual organisms
void Death(){
	for (int x = 0; x < WORLDDIM; ++x) {

		int local_offspring_no = world[x].newFemales.size()+world[x].newMales.size();
		//cout << local_offspring_no << endl;
		if (local_offspring_no > 0) {
			// now clear adult vectors
			world[x].males.clear();
			world[x].females.clear();
			// include local patch extinctions
			if (ran() > epsilon){
				// now copy new females into females
				for (unsigned int nf = 0; nf < world[x].newFemales.size(); ++nf) {
					world[x].females.push_back(world[x].newFemales.at(nf));
				}
				// now copy surviving new males into males
				for (unsigned int nm = 0; nm < world[x].newMales.size(); ++nm) {
					world[x].males.push_back(world[x].newMales.at(nm));
				}
			}
			// clear new females vector
			world[x].newFemales.clear();
			// clear new males vector
			world[x].newMales.clear();
		} else {
			world[x].males.clear();
			world[x].females.clear();
		}
	}
}

//_____________________________________________________________________________
//------------------------------------------------------------------------ main

int main() {
	// random number generator
	specify_rng(time(NULL));
	//specify_rng(RS);

	//read parameters for all simulation runs
	readParameters();

	// read connectivity matrix
	readConnectivityMatrix();

	// repeat loop
	for (int actrun = 0; actrun < max_runs; ++actrun) {

		// output file: metapop
		stringstream outputmetapop_path_stream;
		outputmetapop_path_stream << "output/output_metapop_run" << actrun << ".out";
		string outputmetapop_path = outputmetapop_path_stream.str();
		ofstream outputmetapop(outputmetapop_path.c_str());

		// outputfile header
		outputmetapop << "time" << "    " << "rel_metapopsize" << "    " << "occupancy"<< "    " << "turnover"<< "    " << "emirate" << endl;

		Initialize();

		// time loop
		for (unsigned int acttime = 0; acttime < sim_time; ++acttime) {

			// natal dispersal
			Dispersal();

			// reproduction
			Reproduction();

			// density regulation and death of adults
			Death();

			// analyze metapopulation
			Analyze(acttime, actrun);

			// write metapop results to file
			outputmetapop << acttime << "    " << rel_metapopsize << "    " << occupancy << "    " << turnover_rate << "    " << rel_emigrants << endl;
			// write run and time to console
			//cout << actrun << "   " << acttime << "   " << rel_metapopsize << "    " << occupancy << "   " << turnover_rate << "   " << rel_emigrants << endl;

			// save individual results once at simulation end
			if (acttime == sim_time-1){
				saveResults(actrun);
			}
		}
		// close metapop output file
		outputmetapop.close();
	}
	return 0;
}
