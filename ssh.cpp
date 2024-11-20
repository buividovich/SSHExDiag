#include "ssh.hpp"

SSHModel::SSHModel(int rng_engine_seed)
{
	//Initialize the random number generator
	if(rng_engine_seed==0)
	{
		auto t1 = std::chrono::high_resolution_clock::now();
		int nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(t1.time_since_epoch()).count();
		rng_engine.seed(nanos);
	}
	else
		rng_engine.seed(rng_engine_seed);
	//Init options interface

	ssh_options.add_options()
		(      "L",	   po::value<uint>(  &(L)      )->default_value(    3), "Number of lattice sites")
		(  "kappa",	   po::value<double>(&(kappa)  )->default_value(  1.0), "Hopping amplitude")
		(     "mu",	   po::value<double>(&(mu)     )->default_value(  0.0), "Chemical potential")
		( "lambda",	   po::value<double>(&(lambda) )->default_value(  1.0), "Fermion-phonon coupling")
		(      "w",	   po::value<double>(&(w)      )->default_value(  1.0), "Phonon frequency")
		(      "P",	   po::value<uint>(  &(P)      )->default_value(    0), "Total momentum")
		(      "Q",	   po::value<uint>(  &(Q)      )->default_value(    2), "Total fermion number")
		(     "NM",	   po::value<uint>(  &(NM)     )->default_value(    5), "Cutoff for bosonic occupation numbers");
}

void SSHModel::print_parameters()
{
	cout << ansi::magenta << "SSH model parameters: " << ansi::reset << endl;
	cout << "\t" << ansi::green << "Number of lattice sites L: \t" << ansi::yellow << L       << ansi::reset << endl;
	cout << "\t" << ansi::green << "Hopping amplitude:         \t" << ansi::yellow << kappa   << ansi::reset << endl;
	cout << "\t" << ansi::green << "Chemical potential:        \t" << ansi::yellow << mu      << ansi::reset << endl;
	cout << "\t" << ansi::green << "Fermion-phonon coupling:   \t" << ansi::yellow << lambda  << ansi::reset << endl;
	cout << "\t" << ansi::green << "Phonon frequency:          \t" << ansi::yellow << w       << ansi::reset << endl;
	cout << "\t" << ansi::green << "Total momentum:            \t" << ansi::yellow << P       << ansi::reset << endl;
	cout << "\t" << ansi::green << "Total fermion number:      \t" << ansi::yellow << Q       << ansi::reset << endl;
	cout << "\t" << ansi::green << "Cutoff for nB[i]:          \t" << ansi::yellow << NM      << ansi::reset << endl;
	cout << endl;
}

void SSHModel::InitBasis()
{
	TIMING_INIT;
	if(L<3 || L>Lmax)
	{
		cerr << ansi::red << "Cannot initialize basis for L = " << ansi::cyan << L << ansi::red << "!!! Quitting..." << ansi::reset << endl << flush;
		exit(EXIT_FAILURE);
	};

	//Explicit formula for packing nB's and nF's into pi
	//pi = 2^L*(nB[L-1]*NM^{L-1} + nB[L-2]*NM^{L-2} + nB[0]) + 2^{L-1}*nF[L-1] + ... + nF[0];
	//pi = PIFNB[i]*nB[i] + PIFNF[i]*nF[i];
	//Conversion factors
	PIFNF = new uint64[L]; PIFNB = new uint64[L];
	PIFNF[0] = 1; PIFNB[0] = 2;
	for(uint i=1; i<L; i++){ PIFNF[i] = 2*PIFNF[i-1]; PIFNB[0]*=2; };
	//Initialize conversion factors
	for(uint i=1; i<L; i++)  PIFNB[i] = PIFNB[i-1]*NM;
	//Total number of possible pi values
	NPS = PIFNB[L-1]*NM;
	
	delete [] lookup; lookup = NULL;
	lookup = new uint64[NPS];
	std::fill(lookup, lookup + NPS, NOT_A_STATE);
	//First, construct the lookup table to convert pi to si
	cout << "Initializing the lookup table ..." << flush;
	NS = 0;
	uint64 nbi_count = 0;
	TIMING_START;
	for(uint64 pi=0; pi<NPS; pi++)
	{
		uint nB[Lmax], nF[Lmax];
		PI2ns(pi, nB, nF);
		//Calculate total charge, momentum and bosonic number
		uint aQ = 0, aP = 0, aNB = 0;
		for(uint k=0; k<L; k++)
		{
			aQ  += nF[k];
			aNB += nB[k];
			aP  += k*(nF[k] + nB[k]);
		};
		aP = aP % L; //Total momentum is only defined modulo L
		if(aQ==Q && aP==P && aNB<NM)
		{
			//Add this basis state to the list
			lookup[pi] = NS;
			//increase counter
			NS ++;
		};
		if(aNB<=NM) nbi_count ++;
	};
	
	TIMING_FINISH;
	cout << "... Done in " << ansi::magenta << a_time << ansi::reset << " sec., " << ansi::green << NS << ansi::reset << " basis vectors found! (would be " << ansi::magenta << nbi_count << ansi::reset << " without constraints)" << flush << endl;
	//Now construct the basis itself, to convert si to pi
	delete [] basis; basis = NULL;
	basis = new uint64[NS]; uint64 si = 0;

	for(uint64 pi=0; pi<NPS; pi++)
		if(lookup[pi]!=NOT_A_STATE)
		{
			basis[si] = pi;
			si ++;
		};
	//Consistency check
	if(si!=NS)
		cerr << "Something went wrong, si!=NS" << flush << endl;
	
	//Initializing some additional variables for the Hamiltonian/Current operators
	geff = -kappa*lambda/sqrt(2.0*w*(double)L);
	WH = new t_complex[L*L];
	for(uint k1=0; k1<L; k1++)
		for(uint k2=0; k2<L; k2++)
			WH[k1*L + k2] = exp(-2.0i*M_PI*(double)k1/(double)L) + exp(+2.0i*M_PI*(double)k2/(double)L);
	
	WJ = new t_complex[L*L];
	for(uint k1=0; k1<L; k1++)
		for(uint k2=0; k2<L; k2++)
			WJ[k1*L + k2] = 1.0i*(exp(-2.0i*M_PI*(double)k1/(double)L) - exp(+2.0i*M_PI*(double)k2/(double)L));
}

bool SSHModel::CheckBasisConsistency()
{
	bool res = true;
	for(uint64 si=0; si<NS; si++)
	{
		uint nB[Lmax], nF[Lmax];
		uint64 pi  = basis[si];
		uint64 si1 = lookup[pi];
		res = res && (si1 == si);
		PI2ns(pi, nB, nF);
		uint64 pi1 = ns2PI(nB, nF);
		res = res && (pi1 == pi);
		//Checking P and Q here
		uint aQ = 0, aP = 0, aNB = 0; 
		for(uint k=0; k<L; k++)
		{
			aQ  += nF[k];
			aNB += nB[k];
			aP  += k*(nF[k] + nB[k]);
		};
		aP = aP % L; //Total momentum is only defined modulo L
		res = res && aQ==Q && aP==P && aNB<NM;
	};
	return res;
}

void   SSHModel::PI2ns(uint64 pi, uint* nB, uint* nF)
{
	for(uint i=0; i<L; i++)
	{
		nB[i] = (pi/PIFNB[i])%NM;
		nF[i] = (pi/PIFNF[i])%2;
	};
}

uint64   SSHModel::ns2PI(uint* nB, uint* nF)
{
	uint64 pi = 0;
	for(uint i=0; i<L; i++)
	{
		pi += nB[i]*PIFNB[i];
		pi += nF[i]*PIFNF[i];
	};
	return pi;
}

void        SSHModel::get_suffix0(char* cstr)
{
	sprintf(cstr, "ssh_L%u_k%2.5lf_mu%2.5lf_l%2.5lf_w%2.5lf_N%u", L, kappa, mu, lambda, w, NM);
}

uint        SSHModel::read_eigensystem(string datadir, bool read_evecs)
{
	uint fnev = 0;
	char evals_fname[900], suffix[512], suffixPQ[640];
	get_suffix0(suffix);
	sprintf(suffixPQ, "%s_P%u_Q%u", suffix, P, Q);
	sprintf(evals_fname, "%s/eigensystems/evals_%s.dat", datadir.c_str(), suffixPQ);
	FILE* evals_file = fopen(evals_fname, "rb"); //Binary read
	if(evals_file!=NULL)
	{
		//Check the file size and determine how many eigenvalues it contains
		fseek(evals_file, 0L, SEEK_END);
		size_t fs = ftell(evals_file);
		rewind(evals_file);
		fnev = fs/sizeof(double);
		delete [] E;
		E = new double[fnev];
		if(fread(E, sizeof(double), fnev, evals_file)!=fnev)
		{
			cout << ansi::red << "Some error, Could not read " << fnev << " eigenvalues from the file " << evals_fname << "!!!" << ansi::reset << endl << flush;
			fclose(evals_file);
			return 0;
		};
		cout << ansi::white << "\t Loaded " << ansi::magenta << fnev << ansi::white << " eigenvalues from the file " << ansi::magenta << evals_fname << ansi::reset << endl << flush;
		fclose(evals_file);
	}
	else
	{
		cout << ansi::yellow << "The file " << ansi::magenta << evals_fname << ansi::yellow << " could not be opened!!!" << ansi::reset << endl << flush;
		return 0;
	};
	
	if(read_evecs)
	{
		char evecs_fname[900];
		sprintf(evecs_fname, "%s/eigensystems/evecs_%s.dat", datadir.c_str(), suffixPQ);
		FILE* evecs_file = fopen(evecs_fname, "rb"); //Binary write
		if(evecs_file!=NULL)
		{
			delete [] psi;
			psi = new t_complex[fnev*NS];
			cout << ansi::white << "\t Reading " << ansi::magenta << fnev << ansi::white << " eigenvectors from the file " << ansi::magenta << evecs_fname << ansi::white << " ..."<< ansi::reset << endl << flush;
			TIMING_INIT;
			TIMING_START;
			if(fread(psi, sizeof(double), NS*fnev, evecs_file)!=NS*fnev)
			{
				cerr << ansi::red << "Failed to read " << fnev << " eigenvectors from the file " << evecs_fname << "!!!" << ansi::reset << endl << flush;
				return 0;
			};
			fclose(evecs_file);
			TIMING_STOP;
			cout << ansi::white << "\t ... Done in " << ansi::magenta << a_time << ansi::white << " sec." << ansi::reset << endl << flush;
		}
		else
		{
			cerr << ansi::red << "Failed to open the file " << evecs_fname << " to read the eigenvectors!" << ansi::reset << endl << flush;
			return 0;
		};
	};
	
	this->nev = fnev;
	return fnev;
}

void       SSHModel::write_eigensystem(string datadir, bool write_evecs)
{
	char evals_fname[900], suffix[512], suffixPQ[640];
	get_suffix0(suffix);
	sprintf(suffixPQ, "%s_P%u_Q%u", suffix, P, Q);
	sprintf(evals_fname, "%s/eigensystems/evals_%s.dat", datadir.c_str(), suffixPQ);
	cout << ansi::white << "\t File to save eigenvalues: " << ansi::magenta << evals_fname << ansi::reset << flush << endl;
	FILE* evals_file = fopen(evals_fname, "wb"); //Binary write
	if(evals_file!=NULL)
	{
		size_t res = fwrite(E, sizeof(double), nev, evals_file);
		if(res!=nev)
			cerr << ansi::red << "Only " << res << " doubles out of " << nev << " could be written to the file " << evals_fname << "!!!" << ansi::reset << endl << flush;
		fclose(evals_file);
	}
	else
		cerr << ansi::red << "Failed to open the file " << evals_fname << " for binary writing!" << ansi::reset << endl << flush;
	
	if(!write_evecs) return;
		
	char evecs_fname[900];
	sprintf(evecs_fname, "%s/eigensystems/evecs_%s.dat", datadir.c_str(), suffixPQ);
	cout << ansi::white << "\t File to save eigenvectors: " << ansi::magenta << evecs_fname << ansi::reset << flush << endl;
	FILE* evecs_file = fopen(evecs_fname, "wb"); //Binary write
	if(evecs_file!=NULL)
	{
		if(fwrite(psi, sizeof(t_complex), nev*NS, evecs_file)!=nev*NS)
			cerr << ansi::red << "Failed to write " << nev*NS << " t_complex to the file " << evecs_fname << "!!!" << ansi::reset << endl << flush;
		fclose(evecs_file);
	}
	else
		cerr << ansi::red << "Failed to open the file " << evals_fname << " for binary writing!" << ansi::reset << endl << flush;	 
}

t_complex* SSHModel::HM()
{
	t_complex* res = new t_complex[NS*NS];
	t_complex* B   = new t_complex[NS];	t_complex* HB  = new t_complex[NS];
	
	for(uint si=0; si<NS; si++)
	{
		std::fill(B, B + NS, 0.0 + 0.0i);
		B[si] = 1.0 + 0.0i;
		H(B, HB);
		#pragma omp parallel for
		for(uint sj=0; sj<NS; sj++)
			res[sj*NS + si] = HB[sj];
	};
	
	delete [] HB;	delete [] B;
	return res;
}

int SSHModel::FermionicExchangeFactor(uint k1, uint k2, uint* nF)
{
	int res = 1; //Fermionic factor to reflect anti-commutativity
	uint ka = std::min(k1, k2);
	uint kb = std::max(k1, k2);
	for(uint k=ka+1; k<kb; k++)
		res *= (nF[k]==1? -1 : 1); 
	return res;
}

void SSHModel::H(t_complex* in, t_complex* out)
{
	
	#pragma omp parallel for
	for(uint64 si=0; si<NS; si++)
	{
		uint nB[Lmax], nF[Lmax];
		uint64 pi = basis[si];
		PI2ns(pi, nB, nF);
		//Calculating the diagonal part
		double H_diag = 0.0; uint nBt = 0;
		for(uint k=0; k<L; k++)
		{
			nBt += nB[k];
			H_diag += w*(nB[k] + 0.5) + 2.0*kappa*cos(2.0*M_PI*(double)k/(double)L)*nF[k] + mu*nF[k];
		};
		out[si] = H_diag*in[si];
		
		for(uint k1=0; k1<L; k1++)
			for(uint k2=0; k2<L; k2++)
				if( ((k1==k2) && nF[k1]>0) || ((k1!=k2) && nF[k1]>0 && nF[k2]<1) )
				{
					uint k1k2 = (2*L + k1 - k2)%L;
					uint k2k1 = (2*L + k2 - k1)%L;
					int FF = FermionicExchangeFactor(k1, k2, nF);
					
					nF[k1] = 1 - nF[k1];	nF[k2] = 1 - nF[k2];
					if(nB[k2k1]>0) //psi^+_k1 psi_k2 a^+_{k2-k1}
					{
						nB[k2k1] -= 1;
						uint64 pi1 = ns2PI(nB, nF);
						nB[k2k1]  += 1;
						uint64 si1 = lookup[pi1];
						if(si1==NOT_A_STATE) cout << ansi::red << "ALARM1!!!" << ansi::reset << endl << flush;
						out[si]   += geff*FF*sqrt((double)nB[k2k1])*WH[k1*L + k2]*in[si1];
					};
					if(nBt < NM-1) //psi^+_k1 psi_k2 a^_{k1-k2}
					{
						nB[k1k2]  += 1;
						uint64 pi1 = ns2PI(nB, nF);
						nB[k1k2]  -= 1;
						uint64 si1 = lookup[pi1];
						if(si1==NOT_A_STATE) cout << ansi::red << "ALARM2!!!" << ansi::reset << endl << flush;
						out[si]   += geff*FF*sqrt((double)(nB[k1k2]+1))*WH[k1*L + k2]*in[si1];
					};
					//Return to nF, nB
					nF[k1] = 1 - nF[k1];	nF[k2] = 1 - nF[k2];
				}; //Enf of if(nF ...) and loop over k1, k2 
	}; //End of loop over si 
}

void SSHModel::J(t_complex* in, t_complex* out)
{
	#pragma omp parallel for
	for(uint64 si=0; si<NS; si++)
	{
		uint nB[Lmax], nF[Lmax];
		uint64 pi = basis[si];
		PI2ns(pi, nB, nF);
		//Calculating the diagonal part
		double J_diag = 0.0; uint nBt = 0;
		for(uint k=0; k<L; k++)
		{
			nBt += nB[k];
			J_diag += 2.0*kappa*sin(2.0*M_PI*(double)k/(double)L)*nF[k];
		};
		out[si] = J_diag*in[si];
		
		for(uint k1=0; k1<L; k1++)
			for(uint k2=0; k2<L; k2++)
				if( ((k1==k2) && nF[k1]>0) || ((k1!=k2) && nF[k1]>0 && nF[k2]<1) )
				{
					uint k1k2 = (2*L + k1 - k2)%L;
					uint k2k1 = (2*L + k2 - k1)%L;
					int FF = FermionicExchangeFactor(k1, k2, nF);
					
					nF[k1] = 1 - nF[k1];	nF[k2] = 1 - nF[k2];
					if(nB[k2k1]>0) //psi^+_k1 psi_k2 a^+_{k2-k1}
					{
						nB[k2k1] -= 1;
						uint64 pi1 = ns2PI(nB, nF);
						nB[k2k1]  += 1;
						uint64 si1 = lookup[pi1];
						if(si1==NOT_A_STATE) cout << ansi::red << "ALARM1!!!" << ansi::reset << endl << flush;
						out[si]   += geff*FF*sqrt((double)nB[k2k1])*WJ[k1*L + k2]*in[si1];
					};
					if(nBt < NM-1)
					{
						nB[k1k2]  += 1;
						uint64 pi1 = ns2PI(nB, nF);
						nB[k1k2]  -= 1;
						uint64 si1 = lookup[pi1];
						if(si1==NOT_A_STATE) cout << ansi::red << "ALARM2!!!" << ansi::reset << endl << flush;
						out[si]   += geff*FF*sqrt((double)(nB[k1k2]+1))*WJ[k1*L + k2]*in[si1];
					}
					//Return to nF, nB
					nF[k1] = 1 - nF[k1];	nF[k2] = 1 - nF[k2];
				}; //Enf of if(nF ...) and loop over k1, k2 
	}; //End of loop over si 
}

void SSHModel::rand_vec(double* out, uint n, bool normalize)
{
	for(uint ips=0; ips<n; ips++)
		out[ips] = rng_normal_dist(rng_engine);
	if(normalize)
	{
		double nn = norm(out, n);
		rescale(out, 1.0/nn, n);
	};
}

void SSHModel::rand_vec(t_complex* out, uint n, bool normalize)
{
	for(uint ips=0; ips<n; ips++)
		out[ips] = rng_normal_dist(rng_engine) + 1.0i*rng_normal_dist(rng_engine);
	if(normalize)
	{
		double nn = norm(out, n);
		rescale(out, 1.0/nn + 0.0i, n);
	};
}

t_complex* SSHModel::rand_vec(uint n, bool normalize)
{
	t_complex* res = new t_complex[n];
	rand_vec(res, n, normalize);
	return res;
}

int    SSHModel::LowestEigenstates(int anev, int ncv, double prec)
{
	TIMING_INIT;
	ARrcCompStdEig<double> prob(NS, anev, "SR", ncv, prec, std::numeric_limits<int>::max());
	
	std::cout << "Running ARPACK to find " << anev << " lowest eigenstates, ncv = " << ncv << endl << flush;
	TIMING_START;
	uint64 itc = 0;
	while (!prob.ArnoldiBasisFound())
	{
		prob.TakeStep();
		if ((prob.GetIdo() == 1)||(prob.GetIdo() == -1)) 
		{
			H(prob.GetVector(), prob.PutVector());
			itc++;
       	};
    };
	TIMING_FINISH;
	std::cout << ansi::green << " ... Done in " << ansi::magenta << a_time << ansi::reset << " sec.\n" << endl;
    cout << "No. of iterations: " << prob.GetIter() << " (my own counter " << itc << ")" << endl;
    // Finding eigenvalues and eigenvectors.
    prob.FindEigenvectors();
    if(prob.EigenvaluesFound() && prob.ConvergedEigenvalues()==anev)
    {
    	//Re-packing the eigenvalues
		cout << ansi::green << "\t\t ARNOLDI finished successfully! " << ansi::reset << endl << flush;
		delete [] psi; psi = NULL;
		psi = new t_complex[NS*anev];
		t_complex* cE  = new t_complex[anev];
		prob.EigenValVectors(psi, cE);
		delete [] E; E = NULL;
		E = new double[anev];
		for(uint iev=0; iev<anev; iev++)
			E[iev] = cE[iev].real();
		delete [] cE;
		//Sorting the eigensystem
		sort_eigensystem(E, psi, anev, NS);
		this->nev = anev;
		return prob.GetIter();
	};
	cout << ansi::red << "\t\t ARNOLDI did not converge," << ansi::green << " prob.ConvergedEigenvalues = " << ansi::magenta << prob.ConvergedEigenvalues() << ansi::reset << endl << flush;
	return -1 - prob.ConvergedEigenvalues();
}

void    SSHModel::CheckEigensystem(double* evec_err, double* ortho_err)
{
		cout << "Hi from SSHModel::CheckEigensystem" << endl << flush;
		t_complex* tmp = new t_complex[NS];
		double my_evec_err = 0.0;
		for(uint iev=0; iev<nev; iev++)
		{
			H(psi + iev*NS, tmp);
			A_pluseq_bB(tmp, -E[iev] + 0.0i, psi + iev*NS, NS);
			my_evec_err = std::max(my_evec_err, norm(tmp, NS));
		};
		if(evec_err==NULL)
			cout << ansi::white << "Max. eigenvalue error is " << ansi::magenta << my_evec_err << ansi::reset << endl << flush;
		else
			*evec_err = my_evec_err;
		delete [] tmp;
		double my_ortho_err = orthonormality_norm(psi, nev, NS);
		if(ortho_err==NULL)
			cout << ansi::white << "Max. orthogonality error is " << ansi::magenta << my_ortho_err << ansi::reset << endl << flush;
		else
			*ortho_err = my_ortho_err;
}

void    SSHModel::diagonalize_H()					//find the eigenspectrum of the Hamiltonian
{
	TIMING_INIT;
	//Initialize storage for evecs
	delete [] E; delete [] psi; E = NULL; psi = NULL;
	
	E   = new double[NS];
	psi = HM(); //On input to LAPACKE_zheev, psi will contain the Hamiltonian matrix. Eigenvectors on output
	
	std::cout << "Running LAPACK_zheev to find all " << NS << " lowest eigenstates " << endl << flush;
	TIMING_START;
	int res = LAPACKE_zheev(LAPACK_ROW_MAJOR, 'V', 'U', NS, psi, NS, E);
	TIMING_FINISH;
	std::cout << ansi::green << " ... Done in " << ansi::magenta << a_time << ansi::reset << " sec.\n" << endl; 
	if(res!=0)
	{ 
		cerr << ansi::red << "Something went wrong, LAPACKE_zheev returned " << ansi::cyan << res << ansi::red << " !!!\n" << ansi::reset << endl << flush;
		std::exit(EXIT_FAILURE);
	};
	
	std::cout << "Transposing the eigensystem... " << endl << flush;
	TIMING_START;
	for(uint i=0; i<NS; i++)
		for(uint j=i+1; j<NS; j++)
		{
			t_complex tmp = psi[i*NS + j];
			psi[i*NS + j] = psi[j*NS + i];
			psi[j*NS + i] = tmp;
		};
	TIMING_FINISH;
	std::cout << ansi::green << " ... Done in " << ansi::magenta << a_time << ansi::reset << " sec.\n" << endl; 
	
	nev = NS;
	std::cout << "Sorting the eigensystem... " << endl << flush;
	TIMING_START;
	sort_eigensystem(E, psi, nev, NS);
	TIMING_FINISH;
	std::cout << ansi::green << " ... Done in " << ansi::magenta << a_time << ansi::reset << " sec.\n" << endl;
}

void SSHModel::PrintBasis() //Print out the basis states
{
	for(uint64 si=0; si<NS; si++)
	{
		uint nB[Lmax], nF[Lmax];
		PI2ns(basis[si], nB, nF);
		cout << "\t State " << ansi::cyan << si << ansi::green <<":\t " << "nB = (" << ansi::cyan;
		for(uint i=0; i<L-1; i++)	cout << nB[i] << ", ";
		cout << nB[L-1] << ansi::green << "),\t nF = (" << ansi::cyan;
		for(uint i=0; i<L-1; i++)	cout << nF[i] << ", ";
		cout << nF[L-1] << ansi::green << ")" << ansi::reset << endl << flush;
	};
}

double SSHModel::CheckHermiticity(uint ntrials)
{
	t_complex* psi = new t_complex[NS];
	t_complex* chi = new t_complex[NS];
	t_complex* tmp = new t_complex[NS];
	
	double err = 0.0;
	
	for(uint itrial=0; itrial<ntrials; itrial++)
	{
		rand_vec(psi, NS);
		rand_vec(chi, NS);

		H(psi, tmp); //tmp = H|psi>
		t_complex r1 = scalar_prod(chi, tmp, NS); //r1 = <chi|(H|psi>)
	
		H(chi, tmp); //tmp = H|psi>
		t_complex r2 = scalar_prod(tmp, psi, NS); //r1 = (H|chi>)^+ |psi>
		
		err = std::max(err, std::abs(r1-r2));
	};
	
	delete [] psi;
	delete [] chi;
	delete [] tmp;
	
	return err;
}

void SSHModel::HermiticityDetailedInvestigation()
{
	t_complex* aHM = HM();
	
	double max_err  = 0.0;
	for(uint si=0; si<NS; si++)
		for(uint sj=0; sj<NS; sj++)
		{
			double err = std::abs(aHM[si*NS + sj] - std::conj(aHM[sj*NS + si]));
			if(err>1.0E-10)
			{
				uint nB1[Lmax], nF1[Lmax], nB2[Lmax], nF2[Lmax];
				uint64 pi = basis[si]; uint64 pj=basis[sj];
				PI2ns(pi, nB1, nF1); PI2ns(pj, nB2, nF2); 
				cout << ansi::green << "State i: [" << ansi::cyan;
				for(uint k=0; k<L; k++) cout << nB1[k] << " ";
				for(uint k=0; k<L; k++) cout << nF1[k] << " ";
				cout << ansi::green << "]" << ansi::reset << endl << flush;
				cout << ansi::green << "State j: [" << ansi::cyan;
				for(uint k=0; k<L; k++) cout << nB2[k] << " ";
				for(uint k=0; k<L; k++) cout << nF2[k] << " ";
				cout << ansi::green << "]" << ansi::reset << endl << flush;
				cout << ansi::yellow << "\t" << aHM[si*NS + sj] << "\t" << aHM[sj*NS + si] << ansi::reset << endl << flush << endl;
			};
			max_err = std::max(max_err, err);
		};
	cout << ansi::red << max_err << ansi::reset << flush << endl;
}

double SSHModel::Z(double beta)
{
	double res  = 0.0;
	#pragma omp parallel for reduction (+:res)
	for(uint ie=0; ie<nev; ie++)
		res += exp(-beta*E[ie]);
	return res;
}
