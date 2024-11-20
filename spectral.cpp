#include <iostream>
#include <omp.h>
#include <boost/program_options.hpp>
#include <vector>

#include "ssh.hpp"
#include "timing.hpp"

using namespace std;
namespace po = boost::program_options;

int main(int argc, char **argv)
{
	TIMING_INIT;

	int      nthreads        = 0; //If 0, the number of threads will be set automatically
	string   datadir         = "./data/";
	//Arnoldi params
	bool     use_arnoldi     = false;
	int      nev             = 0;
	int      ncv             = 0;
	double   arnoldi_prec    = 1.0E-8;
	//Physical temperature
	double   T               = 1.0;
	uint     NT              = 32; //Number of equidistant points where the current-current correlators are calculated
	//Histogramming of the spectral function params
	double   wmax            = 5.0; //Max. frequency for the histogram of the spectral functions
	double   wmin            = 1.0E-10; //Min. frequency for the histogram of the spectral functions
	double   dw              = 0.01;  //Frequency resolution for the histogram
	string   spf_hist_fname  = "";
	//Whether to save the full data, and how to cut away very small contributions
	bool     save_full_spf   = false;
	double   spf_cutoff      = 1.0E-8;
	string   full_spf_fname  = "";
	bool     many_body       = false;
	bool     check_eigsys    = false;
	
	po::options_description general_options("Algorithm options");
	general_options.add_options()
		("help,h", "produce this help message")
		("nthreads",	   po::value<int>(            &(nthreads))->default_value(        0), "Number of OpenMP threads to use, 0 = automatic choice"    		  	)
		("datadir", 	   po::value<string>(                       &datadir               ), "Directory for data output"                                		  	)
		("T",              po::value<double>(                &(T))->default_value(      1.0), "Temperature"                                              		  	)
		("NT",             po::value<uint>(                 &(NT))->default_value(       32), "Number of imaginary-time points for JJ correlator"           	  	)
		("nev",		 	   po::value<int>(                 &(nev))->default_value(        0), "Number of lowest eigenstates to consider"                	 	  	)
		("ncv", 		   po::value<int>(                 &(ncv))->default_value(        0), "NCV parameter of the Arnoldi algorithm"                  	 	  	)
		("arnoldi-prec",   po::value<double>(     &(arnoldi_prec))->default_value(  1.0E-10), "Target precision for the Arnoldi algorithm"               		  	)
		("wmax", 	   	   po::value<double>(             &(wmax))->default_value(      5.0), "Max. frequency for spectral function histogramming"      	 	  	)
		("wmin", 		   po::value<double>(             &(wmin))->default_value(   1.0E-5), "Max. frequency for spectral function histogramming"       		  	)
		("dw", 	           po::value<double>(               &(dw))->default_value(     0.01), "Bin width for spectral function histogramming"                   	)
		("spf-cutoff", 	   po::value<double>(       &(spf_cutoff))->default_value(  1.0E-10), "Min. weight of delta functions in spf"                            	)
		("many-body",      "Full many-body calculation (default: single-particle, Q=1 sector)" 	)
		("check-eigsys",   "Calculate the numerical eigensystem error" );
	
	//Diagonalizing the spin chain without any magnetic field ...
	SSHModel* SSH = new SSHModel();
	
	//Combine all the options
	po::options_description all_options;
	
	all_options.add(general_options).add(SSH->ssh_options);
	
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, all_options), vm);
	po::notify(vm);
	// Check if help is needed
	if(vm.count( "help" )){cout<<all_options<<endl; return 1;};
	if(vm.count( "check-eigsys" )){check_eigsys = true;};
	if(vm.count( "many-body" )){many_body = true;};

	use_arnoldi = (nev > 0); //If nev not set, use LAPACK
	
	ncv = (ncv==0? 2*nev + 1 : ncv); //Automatically set ncv if not specified
	
	double beta = 1.0/T;
	
	nthreads = (nthreads==0? omp_get_max_threads() : nthreads);
	//omp_set_num_threads(nthreads);
	openblas_set_num_threads(nthreads);
	
	cout << ansi::white << "Using " << ansi::magenta << nthreads << ansi::white << " OpenMP threads, openblas_num_threads = " << ansi::magenta << openblas_get_num_threads() << ansi::white << " omp_get_num_threads(): " << ansi::magenta << omp_get_max_threads() << ansi::reset << endl;
	
	cout << ansi::green << "Regime:                                             \t" << ansi::white << (many_body? "Full many-body Hamiltonian" : "Single-particle sector (Q=1)") << ansi::reset << endl;
	cout << ansi::green << "Directory for data output:                          \t" << ansi::white << datadir      << ansi::reset << endl;
	if(use_arnoldi)
	{
		cout << ansi::green << "Algorithm:                                          \t" << ansi::cyan    << "ARNOLDI"    << ansi::reset << endl;
		cout << ansi::green << "Number of lowest eigenstates to consider:           \t" << ansi::magenta << nev          << ansi::reset << endl;
		cout << ansi::green << "NCV parameter of the Arnoldi algorithm:             \t" << ansi::magenta << ncv          << ansi::reset << endl;
		cout << ansi::green << "Target precision for the Arnoldi algorithm:         \t" << ansi::magenta << arnoldi_prec << ansi::reset << endl;
	}
	else
	{
		cout << ansi::green << "Algorithm:                                          \t" << ansi::cyan    << "LAPACK" << ansi::yellow << " (finding all eigenstates)"    << ansi::reset << endl;
	};
	cout << ansi::green << "Temperature:                                        \t" << ansi::magenta << T                     ;
	cout << ansi::cyan  << " (beta = " << ansi::magenta << beta << ansi::cyan << ")" << ansi::reset << endl << flush; 
	cout << ansi::green << "Number of points for JJ correlators: \t"                << ansi::magenta << NT                    << ansi::reset << endl;
	cout << ansi::green << "Max. frequency for spectral function histogramming: \t" << ansi::magenta << wmax                  << ansi::reset << endl;
	cout << ansi::green << "Min. frequency for spectral function histogramming: \t" << ansi::magenta << wmin                  << ansi::reset << endl;
	cout << ansi::green << "Bin width for spectral function histogramming:      \t" << ansi::magenta << dw                    << ansi::reset << endl;
	cout << ansi::green << "Output file for the SPF histogram:                  \t" << ansi::white   << spf_hist_fname        << ansi::reset << endl;
	if(save_full_spf)
	{
		cout << ansi::green << "Output:                                             \t" << ansi::cyan    << "FULL SPF + HIST" << ansi::reset << endl;
		cout << ansi::green << "Min. weight of delta functions in spf:              \t" << ansi::magenta << spf_cutoff        << ansi::reset << endl;
		cout << ansi::green << "Output file for the full SPF data:                  \t" << ansi::white   << full_spf_fname    << ansi::reset << endl;
	};
	
	SSH->print_parameters();
	
	//File naming conventions
	char suffix0[512], suffix1[640], suffix2[780];
	SSH->get_suffix0(suffix0);
	if(use_arnoldi)
		sprintf(suffix1, "%s_%s_nev%u", suffix0, (many_body? "MB" : "SP"), nev);
	else
		sprintf(suffix1, "%s_%s", suffix0, (many_body? "MB" : "SP"));
	sprintf(suffix2, "%s_b%2.2lf", suffix1, beta);
	
	//Holder for the spectral function hist
	uint nbins = (uint)std::ceil(wmax/dw);
	double* SPF_hist = new double[nbins];
	std::fill(SPF_hist, SPF_hist + nbins, 0.0);
	
	//Holder for the Euclidean current-current correlator
	double* GE = new double[NT];
	std::fill(GE, GE + NT, 0.0);
	
	double Z0 = exp(-0.5*beta*SSH->w)/(1.0 - exp(-beta*SSH->w)); //Partition function of a free oscillator with frequency w
	
	//TODO: saving/loading of the eigenstates
	
	double Z  = (many_body? 2*SSH->L*Z0 : 0.0);
	
	uint Qmax = (many_body? SSH->L-1 : 1);
	
	for(SSH->Q=1; SSH->Q<=Qmax; SSH->Q++)
		for(SSH->P=0; SSH->P<=SSH->L/2; SSH->P++)
		{
			uint degeneracy_factor = ((SSH->P == (SSH->L - SSH->P)%SSH->L)? 1 : 2); //If P = -P, count contrib only once
			cout << ansi::cyan << "P = " << SSH->P << ", Q = " << SSH->Q << ", degeneracy factor: " << ansi::yellow << degeneracy_factor << ansi::reset << endl << flush;
			SSH->InitBasis();
			//Consistency check
			bool BasisConsistency = SSH->CheckBasisConsistency();
			cout << ansi::cyan << "Basis consistency check: " << (BasisConsistency? ansi::green : ansi::red);
			cout << (BasisConsistency? "PASS" : "FAIL") << ansi::reset << endl << flush;
			
			//Exact diagonalization
			if(use_arnoldi)
				SSH->LowestEigenstates(nev, ncv, arnoldi_prec);
			else
				SSH->diagonalize_H();
			
			cout << ansi::cyan << "  E[0] = " << ansi::magenta << SSH->E[0];
			cout << ansi::cyan << ", E[1] = " << ansi::magenta << SSH->E[1];
			cout << ansi::reset << endl << flush;

			if(check_eigsys)	
				SSH->CheckEigensystem();
			
			Z += degeneracy_factor*SSH->Z(beta);
			
			//Calculating the spectral function and the Euclidean current-current correlator
			t_complex* tmp   = new t_complex[SSH->NS];
			
			std::cout << ansi::white << "Calculating the spectral function and the Euclidean correlators, " << ansi::green << " ssh->nev = " << ansi::magenta << SSH->nev << ansi::reset << endl << flush;
			TIMING_START;
			for(uint m=0; m<SSH->nev; m++)
			{
				SSH->J(SSH->psi + m*SSH->NS, tmp);
				for(uint n=0; n<SSH->nev; n++) 
				{
					//Various components of the spectral function and the current-current correlator
					double Jmn2 = std::abs(scalar_prod(SSH->psi + n*SSH->NS, tmp, SSH->NS));
					Jmn2 = Jmn2*Jmn2;
					
					double SPF_contrib = Jmn2*(exp(-beta*SSH->E[m]) + exp(-beta*SSH->E[n]));
					double w = SSH->E[n] - SSH->E[m];
					//Current-current correlator
					for(uint it=0; it<NT; it++)
					{
						double tau = it*beta/(double)NT;
						GE[it] += Jmn2*exp(-tau*SSH->E[n] - (beta - tau)*SSH->E[m])*degeneracy_factor;
					};
					
					//Spectral function as a histogram
					int ibin = (int)std::floor(w/dw);
					if(ibin>=0 && ibin < nbins && w>wmin) 
						SPF_hist[ibin] += degeneracy_factor*SPF_contrib;
				};
			};
			TIMING_FINISH;
			std::cout << ansi::green << " ... Done in " << ansi::magenta << a_time << ansi::reset << " sec.\n" << endl;
			
			delete [] tmp;
			
			cout << endl;
		};
	
	cout << ansi::green << "Z = " << ansi::magenta << Z << ansi::reset << endl << flush;
	
	//Saving the spectral function
	char hist_fname[1100];
	sprintf(hist_fname, "%s/SPF_hist_%s.dat", datadir.c_str(), suffix2);
	FILE* hist_file = fopen(hist_fname, "w");
	if(hist_file==NULL) cout << ansi::red << "Could not open the file " << ansi::cyan << hist_fname << ansi::red << " for writing" << ansi::reset << endl << flush;
	
	for(uint ibin=0; ibin<nbins; ibin++)
	{
		SPF_hist[ibin] /= (Z*dw);
		if(hist_file!=NULL) fprintf(hist_file, "%2.4lf\t%2.4E\n", (ibin + 0.5)*dw, SPF_hist[ibin]); //TODO: change to more compact binary output
	};
	
	if(hist_file!=NULL)
	{
		fclose(hist_file);
		cout << ansi::green << "SPF histogram was saved to the file " << ansi::yellow << hist_fname << ansi::reset << endl << flush;
	};
	
	//Saving the current-current correlator
	char GE_fname[1100];
	sprintf(GE_fname, "%s/GE_%s.dat", datadir.c_str(), suffix2);
	FILE* GE_file = fopen(GE_fname, "w");
	if(GE_file==NULL) cout << ansi::red << "Could not open the file " << ansi::cyan << GE_fname << ansi::red << " for writing" << ansi::reset << endl << flush;
	
	for(uint it=0; it<NT; it++)
	{
		GE[it] /= Z;
		double tau = it*beta/(double)NT;
		if(GE_file!=NULL) fprintf(GE_file, "%2.4lf\t%2.4E\n", tau, GE[it]);
	};
	
	if(GE_file!=NULL)
	{
		fclose(GE_file);
		cout << ansi::green << "Current-current correlator was saved to the file " << ansi::yellow << GE_fname << ansi::reset << endl << flush;
	};
	
	//cout << ansi::cyan << "Storage space required by ARPACK with n = NS = " << NS << ", nev = " << nev << " : " << ansi::red << arpack_storage_size_GB << " Gb" << ansi::reset << endl << endl;
	cout << endl << endl;
	return EXIT_SUCCESS;
}

