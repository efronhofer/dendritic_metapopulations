//============================================================================
// Emanuel A. Fronhofer & Florian Altermatt
// Classical metapopulation dynamics and eco-evolutionary feedbacks in dendritic networks
// Ecography
// 2016
//
// classes & constants
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

//____________________________________________________________________________
//----------------------------------------------------------- define constants

const int WORLDDIM = 36;														// total world size (number of patches)
const int RS = 2;																// random seed
const int N_ALLELES = 2;														// at one locus; i.e. 1=haploid; 2=diploid

const int NO_NEUTRAL_LOCI = 10;													// no. of diploid neutral loci

//____________________________________________________________________________
//------------------------------------------------------------- define classes

// one individual ------------------------------------------------------------
class TIndiv {
public:
	TIndiv();
	float dispRate[N_ALLELES];
	int neutral_loci[N_ALLELES][NO_NEUTRAL_LOCI];
};

TIndiv::TIndiv() { //constructor for TIndiv
    for (int i = 0; i < N_ALLELES; ++i) {
		dispRate[i] = 0;
		for (int n = 0; n < NO_NEUTRAL_LOCI; ++n) {
			neutral_loci[i][n] = 0;
		}
	}
}

// one patch -----------------------------------------------------------------
class TPatch {
public:
	TPatch();
	vector<TIndiv> males;
	vector<TIndiv> newMales;
	vector<TIndiv> females;
	vector<TIndiv> newFemales;
};

TPatch::TPatch() {
	males.clear();
	newMales.clear();
	females.clear();
	newFemales.clear();
}
