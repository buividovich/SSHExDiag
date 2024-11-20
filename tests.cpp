/* This was a part of spectral.cpp */
	//Checking the matrix elements
	//1. Check that all matrix elements are zero in Q=0 and Q=L sectors
	double MaxErr = 0.0; uint Q = 0; double Err = 0.0; t_complex* JmnPQ = NULL;
	for(uint P=0; P<SSH->L; P++)
	{
		Q = 0;
		JmnPQ = Jmn + (Q*SSH->L + P)*SSH->nev*SSH->nev;
		Err = norm(JmnPQ, SSH->nev*SSH->nev);
		MaxErr = std::max(MaxErr, Err);
		
		Q = SSH->L;
		JmnPQ = Jmn + (Q*SSH->L + P)*SSH->nev*SSH->nev;
		Err = norm(JmnPQ, SSH->nev*SSH->nev);
		MaxErr = std::max(MaxErr, Err);
	};
	cout << ansi::green << "Max. norm of Jmn in Q = 0 and Q=L sectors: " << ansi::magenta << MaxErr << ansi::reset << endl << flush;
		
	//2. Check that Jmn coincide in P and -P sectors
	MaxErr = 0.0; double* Jmn1   = new double[SSH->nev*SSH->nev]; double* Jmn2   = new double[SSH->nev*SSH->nev];
	double MaxErr1 = 0.0; double MaxErr2 = 0.0;
	for(uint Q=1; Q<SSH->L; Q++)
		for(uint P=1; P<=SSH->L/2; P++)
		{
			uint P1 = P;
			uint P2 = (SSH->L - P)%SSH->L;
			cout << ansi::cyan << "Q = " << Q << ", P1 = " << P1 << ", P2 = " << P2 << ansi::reset << endl << flush;
			t_complex* JmnPQ1 = Jmn + (Q*SSH->L + P1)*SSH->nev*SSH->nev;
			t_complex* JmnPQ2 = Jmn + (Q*SSH->L + P2)*SSH->nev*SSH->nev;
			double*      EPQ1 = ES  + (Q*SSH->L + P1)*SSH->nev;
			double*      EPQ2 = ES  + (Q*SSH->L + P2)*SSH->nev;
			
			for(uint ie=0; ie<SSH->nev; ie++)
				for(uint je=0; je<SSH->nev; je++)
				{
					t_complex Jmn1c = JmnPQ1[ie*SSH->nev + je]*JmnPQ1[je*SSH->nev + ie];
					t_complex Jmn2c = JmnPQ2[ie*SSH->nev + je]*JmnPQ2[je*SSH->nev + ie];
					
					Jmn1[ie*SSH->nev + je] = Jmn1c.real();
					Jmn2[ie*SSH->nev + je] = Jmn2c.real();
					
					if(std::abs(Jmn1[ie*SSH->nev + je] - Jmn2[ie*SSH->nev + je])>1.0E-5)
					{
						cout << ansi::red << "(" << ie << "," << je << "): " << Jmn1[ie*SSH->nev + je] << " vs " << Jmn2[ie*SSH->nev + je];
					    cout << ansi::yellow << " E1[ie] = " << EPQ1[ie] << ", E1[je] = " << EPQ1[je]; 
						cout << ansi::yellow << " E2[ie] = " << EPQ2[ie] << ", E2[je] = " << EPQ2[je];
						cout << ansi::reset << endl << flush;
					};
					
					MaxErr1 = std::max(MaxErr1, Jmn1c.imag()*Jmn1c.imag());
					MaxErr1 = std::max(MaxErr1, Jmn2c.imag()*Jmn2c.imag());
				};
			
			for(uint ie=0; ie<SSH->nev; ie++)
				for(uint je=0; je<SSH->nev; je++)
				{	
					double e1 = Jmn1[ie*SSH->nev + je] - Jmn1[je*SSH->nev + ie];
					double e2 = Jmn2[ie*SSH->nev + je] - Jmn2[je*SSH->nev + ie];
					MaxErr2 = std::max(MaxErr2, e1*e1);
					MaxErr2 = std::max(MaxErr2, e2*e2);
				};
			
			Err = norm_diff(Jmn1, Jmn2, SSH->nev*SSH->nev);
			MaxErr = std::max(MaxErr, Err);
		};
	cout << ansi::green << "Max. difference of Jmn in P and -P sectors: " << ansi::magenta << MaxErr << ansi::reset << endl << flush;
	
	cout << ansi::white << MaxErr1 << ansi::reset << endl << flush;
	cout << ansi::white << MaxErr2 << ansi::reset << endl << flush;
	
	delete [] Jmn;
