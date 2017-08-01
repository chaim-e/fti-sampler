#include "model.h"

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

long Petal::casson()
{
	long I = 0;

	for (int d = 0; d < 2 * n; ++d)
		for (int c = 0; c < d; ++c)
			for (int b = (d + 3) & 3; b < c; b += 4)
				if (P[d / 2] > P[b / 2])
					for (int a = (c + 3) & 3; a < b; a += 4)
						if (P[a / 2] > P[c / 2])
							I += 1 - 2 * ((a + d) & 1);
	
	return I;
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

long Petal::jones3()
{
	long I = 0;

	for (int g = 0; g < 2 * n; ++g)
		for (int f = 0; f < g; ++f)
			for (int d = 0; d < f; ++d)
				for (int c = (g + 3) & 3; c < d; c += 4)
					for (int b = (f + 3) & 3; b < c; b += 4)
						for (int a = (d + 3) & 3; a < b; a += 4)
						{
							if (P[g / 2] > P[c / 2] && P[d / 2] > P[a / 2] && P[b / 2] > P[f / 2]) I += 1 - 2 * ((g + d + b) & 1);
							if (P[g / 2] < P[c / 2] && P[d / 2] < P[a / 2] && P[b / 2] < P[f / 2]) I += 1 - 2 * ((a + c + f) & 1);
						}

	I *= 2;

	for (int g = 0; g < 2 * n; ++g)
		for (int f = 0; f < g; ++f)
			for (int d = (g + 3) & 3; d < f; d += 4)
				if (P[g / 2] > P[d / 2])
					for (int c = 0; c <= d; ++c)
						for (int b = (f + 3) & 3; b < c; b += 4)
							for (int a = (c + 3) & 3; a < b; a += 4)
								if (P[c / 2] > P[a / 2])
								{
									if (P[f / 2] > P[b / 2]) I += 1 - 2 * ((g + c + f) & 1);
									if (P[f / 2] < P[b / 2]) I += 1 - 2 * ((g + c + b) & 1);
								}
	
	for (int g = 0; g < 2 * n; ++g)
		for (int f = 0; f < g; ++f)
			for (int d = 0; d <= f; ++d)
				for (int c = (g + 3) & 3; c < d; c += 4)
					for (int b = (d + 3) & 3; b < c; b += 4)
						if (P[d / 2] > P[b / 2])
							for (int a = (f + 3) & 3; a <= b; a += 4)
								if (P[a / 2] > P[f / 2])
								{
									if (P[g / 2] > P[c / 2]) I += 1 - 2 * ((a + d + g) & 1);
									if (P[g / 2] < P[c / 2]) I += 1 - 2 * ((a + d + c) & 1);
								}
	
	for (int g = 0; g < 2 * n; ++g)
		for (int f = 0; f <= g; ++f)
			for (int d = 0; d < f; ++d)
				for (int c = (f + 3) & 3; c < d; c += 4)
					if (P[f / 2] > P[c / 2])
						for (int b = (g + 3) & 3; b <= c; b += 4)
							if (P[b / 2] > P[g / 2])
								for (int a = (d + 3) & 3; a < b; a += 4)
								{
									if (P[d / 2] > P[a / 2]) I += 1 - 2 * ((b + f + d) & 1);
									if (P[d / 2] < P[a / 2]) I += 1 - 2 * ((b + f + a) & 1);
								}
	
	return I/2;
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

long Petal::writhe()
{
	long I = 0;

	for (int b = 1; b < n; ++b)
		for (int a = 0; a < b; ++a)
			I += 1 - 2 * ((a + b + int(P[a] < P[b])) & 1);
	return I;
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

long Grid::casson()
{
	long I = 0;
	
	for (int a = 1; a <= n; ++a)
		for (int b = a; b < n; ++b)
			for (int c = b; c < n; ++c)
				if (((P[a]-P[c])*(P[a]-P[c+1]) < 0) && ((Q[c]-Q[a-1])*(Q[c]-Q[a]) < 0))
					for (int d = c+1; d <= n; ++d)
						if (((Q[b]-Q[d-1])*(Q[b]-Q[d]) < 0) && ((P[d]-P[b])*(P[d]-P[b+1]) < 0))
							if ((b!=c) || ((P[b+1]-P[b])*(P[a]-P[d]) > 0))
								I += 1 - 2*(1&(int(Q[a]>Q[a-1])+int(P[b+1]>P[b])+int(P[c+1]>P[c])+int(Q[d]>Q[d-1])));
	
	return I;
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

long Gridle::casson()
{
	long I = 0;
	int N = 4*n*n;
	
	for (int a = 1; a <= n; ++a)
		for (int c = 0; c < n; ++c)
			if (((P[a]-P[c])*(P[a]-P[c+1]) < 0) && ((Q[c]-Q[a-1])*(Q[c]-Q[a]) < 0))
				for (int b = a; b <= n; ++b)
				{
					int Ta = (4*a-1)*n + Q[c]*(1-2*int(Q[a-1] > Q[a]));
					int Tc = (4*c+1)*n + P[a]*(1-2*int(P[c] > P[c+1]));					
					int s1 = 1 - 2*(1&(int(Q[a-1]<Q[a])+int(P[c]<P[c+1])+int(a<c)));
					for (int d = int(a==b)*(c+1); d < n; ++d)
						if (((P[b]-P[d])*(P[b]-P[d+1]) < 0) && ((Q[d]-Q[b-1])*(Q[d]-Q[b]) < 0))
						{
							int Tb = (4*b-1)*n + Q[d]*(1-2*int(Q[b-1] > Q[b]));
							int Td = (4*d+1)*n + P[b]*(1-2*int(P[d] > P[d+1]));
							//cout << N << "  " << Ta << " " << Tc << "  " << Tb << " " << Td << " \t "; 
							//cout << a << " " << c << "  " << b << " " << d << " \t ";							
							if (
									((Ta < Tb) && (Tb < Tc) && (Tc < Td) && (S[P[a]*n+Q[c]] == 1) && (S[P[b]*n+Q[d]] == 0)) ||
									((Tc < Tb) && (Tb < Ta) && (Ta < Td) && (S[P[a]*n+Q[c]] == 0) && (S[P[b]*n+Q[d]] == 0)) ||
									((Ta < Td) && (Td < Tc) && (Tc < Tb) && (S[P[a]*n+Q[c]] == 1) && (S[P[b]*n+Q[d]] == 1)) ||
									((Tc < Td) && (Td < Ta) && (Ta < Tb) && (S[P[a]*n+Q[c]] == 0) && (S[P[b]*n+Q[d]] == 1)) ||
									((Tb < Ta) && (Ta < Td) && (Td < Tc) && (S[P[b]*n+Q[d]] == 1) && (S[P[a]*n+Q[c]] == 0)) ||
									((Td < Ta) && (Ta < Tb) && (Tb < Tc) && (S[P[b]*n+Q[d]] == 0) && (S[P[a]*n+Q[c]] == 0)) ||
									((Tb < Tc) && (Tc < Td) && (Td < Ta) && (S[P[b]*n+Q[d]] == 1) && (S[P[a]*n+Q[c]] == 1)) ||
									((Td < Tc) && (Tc < Tb) && (Tb < Ta) && (S[P[b]*n+Q[d]] == 0) && (S[P[a]*n+Q[c]] == 1)))
							{
								int s2 = 1 - 2*(1&(int(Q[b-1]<Q[b])+int(P[d]<P[d+1])+int(b<d)));
								I += s1*s2;								
								//cout << s1 << " " << s2;
							}
							//cout << endl;
						}					
				}
	
	return -I;
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

long Grid::defect()
{
	long I = 0;
	int N = 4*n*n;
	
	for (int a = 1; a <= n; ++a)
		for (int c = 0; c < n; ++c)
			if (((P[a]-P[c])*(P[a]-P[c+1]) < 0) && ((Q[c]-Q[a-1])*(Q[c]-Q[a]) < 0))
				for (int b = a; b <= n; ++b)
				{
					int Ta = (4*a-1)*n + Q[c]*(1-2*int(Q[a-1] > Q[a]));
					int Tc = (4*c+1)*n + P[a]*(1-2*int(P[c] > P[c+1]));					
					int s1 = 1 - 2*(1&(int(Q[a-1]<Q[a])+int(P[c]<P[c+1])+int(a<c)));
					for (int d = int(a==b)*(c+1); d < n; ++d)
						if (((P[b]-P[d])*(P[b]-P[d+1]) < 0) && ((Q[d]-Q[b-1])*(Q[d]-Q[b]) < 0))
						{
							int Tb = (4*b-1)*n + Q[d]*(1-2*int(Q[b-1] > Q[b]));
							int Td = (4*d+1)*n + P[b]*(1-2*int(P[d] > P[d+1]));
							if (((Tb-Ta+N)%N + (Tc-Tb+N)%N + (Td-Tc+N)%N + (Ta-Td+N)%N) != 2*N)
								I += s1*(1 - 2*(1&(int(Q[b-1]<Q[b])+int(P[d]<P[d+1])+int(b<d))));
						}
				}
	
	return -I;
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

long Grid::unnamed()
{
	long I = 0;
	int N = 4*n*n;
	
	for (int a = 1; a <= n; ++a)
		for (int c = 0; c < n; ++c)
			if (((P[a]-P[c])*(P[a]-P[c+1]) < 0) && ((Q[c]-Q[a-1])*(Q[c]-Q[a]) < 0))
				for (int b = a; b <= n; ++b)
				{
					int Ta = (4*a-1)*n + Q[c]*(1-2*int(Q[a-1] > Q[a]));
					int Tc = (4*c+1)*n + P[a]*(1-2*int(P[c] > P[c+1]));					
					int s1 = 1 - 2*(1&(int(Q[a-1]<Q[a])+int(P[c]<P[c+1])+int(a<c)));
					for (int d = int(a==b)*(c+1); d < n; ++d)
						if (((P[b]-P[d])*(P[b]-P[d+1]) < 0) && ((Q[d]-Q[b-1])*(Q[d]-Q[b]) < 0))
						{
							int Tb = (4*b-1)*n + Q[d]*(1-2*int(Q[b-1] > Q[b]));
							int Td = (4*d+1)*n + P[b]*(1-2*int(P[d] > P[d+1]));
							if (((Tb-Ta+N)%N + (Tc-Tb+N)%N + (Td-Tc+N)%N + (Ta-Td+N)%N) == 2*N)
								I += s1*(1 - 2*(1&(int(Q[b-1]<Q[b])+int(P[d]<P[d+1])+int(b<d)+
										int((max(Ta,Tc) > min(Tb,Td)) && (min(Ta,Tc) < max(Tb,Td))))));
						}
				}
	
	return I;
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

long Grid::cross()
{
	long I = 0;
	
	for (int a = 1; a <= n; ++a)
		for (int c = 0; c < n; ++c)
			if (((P[a]-P[c])*(P[a]-P[c+1]) < 0) && ((Q[c]-Q[a-1])*(Q[c]-Q[a]) < 0))
				I += 1;
	
	return I;
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

long Grid::whitney()
{
	long I = 0;
	
	for (int a = 1; a <= n; ++a)
		I += ((P[a] > P[a-1]) ^ (Q[a] > Q[a-1])) - ((P[(a+1)%n] > P[a]) ^ (Q[a] > Q[a-1]));
	
	return I/2;
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

long Star::casson()
{
	long I = 0;

	for (int d = 0; d < 2 * n; ++d)
		for (int c = 0; c < d; ++c)
			for (int b = (d + 3) & 3; b < c; b += 4)
				if (S[(b/2)*n+(d/2)] == 0)
					for (int a = (c + 3) & 3; a < b; a += 4)
						if (S[(a/2)*n+(c/2)] == 1)
							I += 1 - 2 * ((a + d) & 1);
	
	return I;
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

long Jump::casson()
{
	long I = 0;
	
	double* alpha = new double[n*n];
	double* beta = new double[n*n];	
	double* zdiff = new double[n*n];
	int* S = new int[n*n];
	for (int i = 0; i < n; ++i)
		for (int j = i+2; j < n; ++j) 
		{
			double det = (X[i+1] - X[i]) * (Y[j] - Y[j+1]) - (Y[i+1] - Y[i]) * (X[j] - X[j+1]);
			alpha[i*n+j] = ((X[j] - X[i]) * (Y[j] - Y[j+1]) - (Y[j] - Y[i]) * (X[j] - X[j+1])) / det;
			beta[i*n+j] = ((X[i+1] - X[i]) * (Y[j] - Y[i]) - (Y[i+1] - Y[i]) * (X[j] - X[i])) / det;
			zdiff[i*n+j] = (Z[j] + (Z[j+1]-Z[j]) * beta[i*n+j]) - (Z[i] + (Z[i+1]-Z[i]) * alpha[i*n+j]);
			if ((alpha[i*n+j] <= 0.0) || (alpha[i*n+j] >= 1.0) || (beta[i*n+j] <= 0.0) || (beta[i*n+j] >= 1.0))
				S[i*n+j] = 0;
			else
				S[i*n+j] = 1-2*(int(det < 0) ^ int(zdiff[i*n+j] < 0));
			// cout << i << " " << j << " " << det << " " << 
			//		alpha[i*n+j] << " " << beta[i*n+j] << " " << zdiff[i*n+j] << " " << S[i*n+j] << endl;
		}
	S[n-1] = 0;
	
	for (int a = 0; a < n; ++a)
		for (int b = a; b <= n; ++b)
			for (int c = b; c < n; ++c)
				for (int d = c; d < n; ++d)
					if ((c > a+1) && (d > b+1) &&
							(zdiff[a*n+c] > 0) && (zdiff[b*n+d] <= 0) && 
							((a!=b) || (alpha[a*n+c] < alpha[b*n+d])) && 
							((b!=c) || (beta[a*n+c] > alpha[b*n+d])) && 
							((c!=d) || (beta[a*n+c] < beta[b*n+d])))
					{
						I += S[n*a+c] * S[n*b+d];
						//cout << a << "~" << c << " " << b << "~" << d << " " << S[n*a+c] * S[n*b+d] << endl;
					}
	//cout << "I=" << I << endl;
	
	delete [] alpha;
	delete [] beta;
	delete [] zdiff;
	delete [] S;
	
	return I;
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

long Jump::whitney()
{
	double I = 0.0;
	double pi = acos(-1.0);	
	
	double* angles = new double[n+1];
	for (int i = 0; i < n; ++i) angles[i] = atan2(Y[i+1]-Y[i],X[i+1]-X[i]);
	angles[n] = angles[0];

	for (int i = 0; i < n; ++i) I += ((angles[i] > 0) - (angles[i+1] > 0)) * (abs(angles[i+1]-angles[i]) > pi); 
	
	delete [] angles;
	
	return int(round(I));
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

long Jump::cross()
{
	long I = 0;
	
	for (int i = 0; i < n; ++i)
		for (int j = i+2; j < n - int(i==0); ++j) 
		{
			double det = (X[i+1] - X[i]) * (Y[j] - Y[j+1]) - (Y[i+1] - Y[i]) * (X[j] - X[j+1]);
			double alpha = ((X[j] - X[i]) * (Y[j] - Y[j+1]) - (Y[j] - Y[i]) * (X[j] - X[j+1])) / det;
			double beta = ((X[i+1] - X[i]) * (Y[j] - Y[i]) - (Y[i+1] - Y[i]) * (X[j] - X[i])) / det;
			
			if ((alpha > 0.0) && (alpha < 1.0) && (beta > 0.0) && (beta < 1.0))
				I += 1;
		}
	
	return I;
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

long Jump::defect()
{
	long I = 0;
	
	double* alpha = new double[n*n];
	double* beta = new double[n*n];	
	int* S = new int[n*n];
	for (int i = 0; i < n; ++i)
		for (int j = i+2; j < n; ++j) 
		{
			double det = (X[i+1] - X[i]) * (Y[j] - Y[j+1]) - (Y[i+1] - Y[i]) * (X[j] - X[j+1]);
			alpha[i*n+j] = ((X[j] - X[i]) * (Y[j] - Y[j+1]) - (Y[j] - Y[i]) * (X[j] - X[j+1])) / det;
			beta[i*n+j] = ((X[i+1] - X[i]) * (Y[j] - Y[i]) - (Y[i+1] - Y[i]) * (X[j] - X[i])) / det;
			if ((alpha[i*n+j] <= 0.0) || (alpha[i*n+j] >= 1.0) || (beta[i*n+j] <= 0.0) || (beta[i*n+j] >= 1.0))
				S[i*n+j] = 0;
			else
				S[i*n+j] = 1-2*(int(det < 0));
		}
	S[n-1] = 0;
	
	for (int a = 0; a < n; ++a)
		for (int b = a; b < n; ++b)
			for (int c = b; c < n; ++c)
				for (int d = c; d < n; ++d)
					if ((c > a+1) && (d > b+1) &&							 
							((a!=b) || (alpha[a*n+c] < alpha[b*n+d])) && 
							((b!=c) || (beta[a*n+c] > alpha[b*n+d])) && 
							((c!=d) || (beta[a*n+c] < beta[b*n+d])))
						I -= S[n*a+c] * S[n*b+d];
	
	delete [] alpha;
	delete [] beta;
	delete [] S;
	
	return I;
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

long Jump::unnamed()
{
	long I = 0;
	
	double* alpha = new double[n*n];
	double* beta = new double[n*n];	
	int* S = new int[n*n];
	for (int i = 0; i < n; ++i)
		for (int j = i+2; j < n; ++j) 
		{
			double det = (X[i+1] - X[i]) * (Y[j] - Y[j+1]) - (Y[i+1] - Y[i]) * (X[j] - X[j+1]);
			if (det == 0.0) return SINGULARITY;
			alpha[i*n+j] = ((X[j] - X[i]) * (Y[j] - Y[j+1]) - (Y[j] - Y[i]) * (X[j] - X[j+1])) / det;
			beta[i*n+j] = ((X[i+1] - X[i]) * (Y[j] - Y[i]) - (Y[i+1] - Y[i]) * (X[j] - X[i])) / det;
			if ((alpha[i*n+j] <= 0.0) || (alpha[i*n+j] >= 1.0) || (beta[i*n+j] <= 0.0) || (beta[i*n+j] >= 1.0))
				S[i*n+j] = 0;
			else
				S[i*n+j] = 1-2*(int(det < 0));
		}
	S[n-1] = 0;
	
	for (int a = 0; a < n; ++a)
		for (int b = a; b < n; ++b)
			for (int c = b; c < n; ++c)
				for (int d = c; d < n; ++d)
				{
					if ((b > a+1) && (d > c+1) &&							 
							((b!=c) || (beta[a*n+b] < alpha[c*n+d]))) 
						I += S[n*a+b] * S[n*c+d];
					if ((c > b+1) && 							 
							((a!=b) || (alpha[a*n+d] < alpha[b*n+c])) && 
							((c!=d) || (beta[b*n+c] < beta[a*n+d]))) 
						I -= S[n*a+d] * S[n*b+c];					
				}
	
	delete [] alpha;
	delete [] beta;
	delete [] S;
	
	return I;
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

 


