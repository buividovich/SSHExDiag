#include "linalg.hpp"

void        rescale(double* A, double a, uint n)
{
	cblas_dscal(n, a, A, 1);
}

void		rescale(t_complex* A, t_complex a, uint n)
{
	cblas_zscal(n, &a, A, 1);
}

void  A_pluseq_bB(double* A, double b, double *B, uint n)
{
	cblas_daxpy(n,b,B,1,A,1);
}

void  A_pluseq_bB(t_complex* A, t_complex b, t_complex *B, uint n)
{
	cblas_zaxpy(n,&b,B,1,A,1);
}

void  A_eq_B_mult_C(double* A, double* B, double* C, uint n) //Matrix multiplication C = A*B
{
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, B, n, C, n, 0.0, A, n);
}

void  A_eq_B_mult_C(t_complex* A, t_complex* B, t_complex* C, uint n) //Matrix multiplication C = A*B
{
	t_complex alpha = 1.0 + 0.0i;
	t_complex beta  = 0.0 + 0.0i; 
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, &alpha, B, n, C, n, &beta, A, n);
}

void  commutator(t_complex* C, t_complex* A, t_complex* B, uint n) //C = [A, B] = A*B - B*A
{
	t_complex alpha = 1.0 + 0.0i; t_complex beta  = 0.0 + 0.0i; 
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, &alpha, A, n, B, n, &beta, C, n);
	alpha = -1.0 + 0.0i; beta = 1.0 + 0.0i;
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, &alpha, B, n, A, n, &beta, C, n);
}


void		A_pluseq_a_psi_dirprod_chi(double* A, double a, double* psi, double* chi, uint n)
{
	for(uint ij=0; ij<n*n; ij++)
	{
		uint i = ij/n; uint j = ij - i*n;
		A[ij] += a*psi[i]*chi[j];
	};
}

double  scalar_prod(double* psi1, double* psi2, uint n)
{
	return cblas_ddot(n, psi1, 1, psi2, 1);
}

t_complex   scalar_prod(t_complex* psi1, t_complex* psi2, uint n)
{
	t_complex res  = 0.0 + 0.0i;
	cblas_zdotc_sub(n, psi1, 1, psi2, 1, &res);
	return res;
}

double      norm(double*    psi, uint n)
{
	return cblas_dnrm2(n, psi, 1);
};

double   norm(t_complex* psi, uint n)
{
	return cblas_dznrm2(n, psi, 1);
};

double  norm_diff(double* psi1, double* psi2, uint n)
{
	double res = 0.0;
	for(uint i=0; i<n; i++)
	{
		double d = psi1[i] - psi2[i];
		res += d*d;
	};
	return sqrt(res);
}

double    norm_diff(t_complex* psi1, t_complex* psi2, uint n)
{
	double res = 0.0;
	for(uint i=0; i<n; i++)
	{
		t_complex d = (psi1[i] - psi2[i]);
		res += real(conj(d)*d);
	};
	return sqrt(res);
}

t_complex*  identity_matrix(uint n)
{
	t_complex* r = new t_complex[n*n];
	for(uint i=0; i<n; i++)
		for(uint j=0; j<n; j++)
			r[i*n + j] = (i==j? 1.0 + 0.0i : 0.0 + 0.0i);
	return r;
}

double      unitarity_norm(t_complex* U, uint n)
{
	double unorm = 0.0;
	for(uint i=0; i<n; i++)
		for(uint j=0; j<n; j++)
		{
			t_complex r = (i==j? -1.0 + 0.0i : 0.0 + 0.0i);
			for(uint k=0; k<n; k++)
				r += U[i*n + k]*conj(U[j*n + k]);
			unorm += real(r*conj(r));
		}
	return sqrt(unorm);
}

void check_eigensystem(double* A, double* evals, double* evecs, uint n, double* evec_err, double* ortho_err)
{
	double* tmp = new double[n];
	double max_evec_err = 0.0;
	//Eigenvector error
	for(uint ie=0; ie<n; ie++)
	{
		psi_eq_A_mult_chi(tmp, A, &(evecs[n*ie]), n);
		A_pluseq_bB(tmp, -evals[ie], &(evecs[n*ie]), n);
		max_evec_err = std::max(max_evec_err, norm(tmp, n));
	};
	delete [] tmp;
	if(evec_err!=NULL)
		*evec_err = max_evec_err;
	else
		cout << endl << "\t Max. eigenstate error: " << max_evec_err << endl; 
	//Orthonormality
	double max_ortho_err = 0.0;
	for(uint ie1=0; ie1<n; ie1++)
		for(uint ie2=0; ie2<n; ie2++)
		{
			double sp = scalar_prod(&(evecs[n*ie1]), &(evecs[n*ie2]), n) - (ie1==ie2? 1.0 : 0.0);
			max_ortho_err = std::max(max_ortho_err, std::abs(sp));
		};
	if(ortho_err!=NULL)
		*ortho_err = max_ortho_err;
	else
		cout << endl << "\t Max. orthogonality error: " << max_ortho_err << endl; 
}

void        double2complex(t_complex* out, double* in, uint n)
{
	for(uint i=0; i<n; i++)
		out[i] = in[i] + 0.0i;
}

t_complex*  double2complex(double* in, uint n)
{
	t_complex* res = new t_complex[n];
	double2complex(res, in, n);
	return res;
}

double      unitarity_err(t_complex* U, uint n)
{
	double max_err = 0.0;
	for(uint i=0; i<n; i++)
		for(uint j=0; j<n; j++)
		{
			t_complex UUd = 0.0 + 0.0i;
			for(uint k=0; k<n; k++)
				UUd += U[i*n + k]*conj(U[j*n + k]);
			UUd -= (i==j? 1.0 + 0.0i : 0.0i);
			max_err = std::max(max_err, std::abs(UUd));
		};
	return max_err;
}

void hermitian_conjugate(t_complex* out, t_complex* in, uint n)
{
	for(uint i=0; i<n; i++)
		for(uint j=0; j<n; j++)
		{
			out[i*n + j] = conj(in[j*n + i]);
		};
}

t_complex* hermitian_conjugate(t_complex* A, uint n)
{
	t_complex* r = new t_complex[n*n];
	hermitian_conjugate(r, A, n);
	return r;
}

