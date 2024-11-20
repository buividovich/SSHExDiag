#ifndef _SSH_HPP_
#define _SSH_HPP_

using namespace std;

#include<random>
#include<iostream>
#include<iomanip>
#include<chrono>
#include<fstream>
#include<algorithm>
#include<vector>
#include<limits>
#include<arpackpp/arrscomp.h>
#include <boost/program_options.hpp>

#include "linalg.hpp"
#include "ansi_io.hpp"
#include "timing.hpp"

namespace po = boost::program_options;

static uint Lmax = 8; //Maximal possible value of L

//TODO: introduce an option for the interaction of the form (1 - \lambda*(x_i - x_{i+1})), this is the standard set-up in all 1D studies

//Some special values for indexing
static const uint64 NOT_A_STATE   = std::numeric_limits<uint64>::max();    //Code for any state that does not satisfy the P, Q constraints
static const uint64 BEYOND_CUTOFF = std::numeric_limits<uint64>::max()-1;  //Code for states which satisfy the P, Q constraints but go beyond cutoff

class SSHModel
{
	private:
		std::ranlux48 rng_engine;
		std::uniform_real_distribution<double> rng_uniform_dist;
		std::normal_distribution<double> rng_normal_dist{0.0, 1.0};
	public:
		SSHModel(int rng_engine_seed = 0);
		~SSHModel(){};
		//Option holder
		po::options_description ssh_options;
		void print_parameters();
		//Global parameters of the model
		uint         L      = 3;
		double       kappa  = 1.0;
		double       mu     = 0.0;
		double       lambda = 1.0;
		double       w      = 1.0;
		//Total momentum + total charge - conserved quantities to fix the sector
		uint         P      = 0; //P is defined in the range 0 ... L-1
		uint         Q      = 1; //Total charge, in the range 0 ... L
		//Approximation/truncation parameters
		uint         NM     = 5; //Max. number of bosonic eigenstates for each bosonic d.o.f.
		//Derived/calculated parameters
		uint         NS     = 0; //Total number of states in a given (P, Q) sector
		uint64       NS2    = 0; //NS*NS
		uint64       NPS    = 0; //Number of entries in the look-up table
		void InitBasis(); //Init the basis of states with fixed P, Q
		void PrintBasis(); //Print out the basis states
		bool CheckBasisConsistency();
		//Functions for converting indices
		uint64*      PIFNB  = NULL;
		uint64*      PIFNF  = NULL;
		//SI = sequential index of a basis state in a set of states with fixed P, Q
		//PI = "large" packed index: 2^L*(NM^(L-1)*nB[L-1] + NM^(L-2)*nB[L-2] + ... + nB[0] ) + 2^(L-1)*nF[L-1] + nF[0] 
		//ns = array of bosonic/fermionic occupation numbers {nB[0] ... nB[L-1], nF[0] .. nF[L-1]}
		void   PI2ns(uint64 pi, uint* nB, uint* nF);
		uint64 ns2PI(uint* nB, uint* nF);
		//Some intermediate variables for the current/hamiltonian
		double geff  = 0.0;
		t_complex* WH = NULL;
		t_complex* WJ = NULL;
		//Storage size estimators
		uint64       arpack_storage_size(uint nev, uint ncv){return (4*(uint64)NS + (uint64)NS*(uint64)ncv + (uint64)(ncv*ncv) + 8*(uint64)ncv + (uint64)nev)*sizeof(double);}; 
		uint64*      basis  = NULL;
		uint64*      lookup = NULL;
		/* Fermionic factors, needed for H and J */
		int          FermionicExchangeFactor(uint k1, uint k2, uint* nF);
		/* Hamiltonian definition */
		void         H(t_complex* in, t_complex* out);	 				//Hamiltonian
		t_complex*   HM(); //Returns the Hamiltonian matrix
		/* Electric current definition */
		void         J(t_complex* in, t_complex* out);	 				//Electric current, summed over all sites
		//Diagonalization routines
		void         diagonalize_H();
		int          LowestEigenstates(int anev, int ncv, double prec); //Saves nev lowest energy eigenstates to E and psi 
		void         CheckEigensystem(double* evec_err=NULL, double* ortho_err=NULL);
		uint         nev  = 0; 
		double*      E    = NULL;									//Eigenenergies
		t_complex*   psi  = NULL;									//The eigenvectors
		//Saving/reading eigensystems
		void       get_suffix0(char* cstr);
		uint       read_eigensystem(string datadir, bool read_evecs=true); //Returns the number of eigenvalues/vectors read in
		void       write_eigensystem(string datadir, bool write_evecs=true);
		//Tests
		double CheckHermiticity(uint ntrials); //Returns the hermiticity error
		void   HermiticityDetailedInvestigation(); //Prints out non-Hermitian matrix elements explicitly
		//Misc
		void        rand_vec(double*    out, uint n, bool normalize=true);	//Random real vector, filled with Gaussian numbers with unit dispersion
		void        rand_vec(t_complex* out, uint n, bool normalize=true);
		t_complex*  rand_vec(uint n, bool normalize=true);
		double      SPF_delta_contribution(uint m, uint n); //Contribution of a single delta-function for energy levels n = 0...nev-1, m = 0 .. nev -1
		double      Z(double beta);                         //Partition function for a given sector P, Q
};

#endif
