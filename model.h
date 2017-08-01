#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <sstream>
#include <iomanip>
#include <random>
#include "time.h"

#define NOT_IMPLEMENTED 999999999
#define SINGULARITY 999999998

using namespace std;

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

class Model
{
public:
	
	int n;
	Model(int N) : n(N) { }
	
	// shuffling
	virtual void randomize(unsigned int seed) { srand(seed); }
	virtual void shuffle() = 0;
	
	// output
	friend ostream& operator<<(ostream& os, const Model& model) 
		{ return model.output(os); }
	virtual ostream& output(ostream& os) const 
		{ return os << "Output not implemented (n=" << n << ")"; }
	
	// invariants
	virtual long casson() { return NOT_IMPLEMENTED; }
	virtual long jones3() { return NOT_IMPLEMENTED; }
	virtual long writhe() { return NOT_IMPLEMENTED; }
	virtual long whitney() { return NOT_IMPLEMENTED; }
	virtual long cross() { return NOT_IMPLEMENTED; }
	virtual long defect() { return NOT_IMPLEMENTED; }
	virtual long unnamed() { return NOT_IMPLEMENTED; }

};

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

class Petal : public Model
{
public:
	
	vector<int> P;
	
	Petal(int N) : Model(N), P(N) 
	{ 
		for (int i = 0; i < n; ++i) P[i] = i;
	}
	
	void shuffle()
	{
		random_shuffle(P.begin(), P.end());
	}

	ostream& output(ostream& os) const
	{ 
		os << "Petaluma(["; 
		for (int i = 0; i < n; ++i) os << P[i] << ","; os << "])";
		return os; 
	}
	
	long casson();
	long jones3();
	long writhe();
	
};

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

class Grid : public Model
{
public:
	
	vector<int> P, Q;
	
	Grid(int N) : Model(N), P(N+1), Q(N+1) 
	{         
		for (int i = 0; i < n; ++i) P[i] = i;
		for (int i = 0; i < n; ++i) Q[i] = i;
	}
	
	void shuffle()
	{                
		random_shuffle(P.begin(), --P.end()); P[n]=P[0];
		random_shuffle(Q.begin(), --Q.end()); Q[n]=Q[0];
	}

	ostream& output(ostream& os) const
	{ 
		os << "Grid(["; 
		for (int i = 0; i < n; ++i) os << P[i] << ","; os << "],[";
		for (int i = 0; i < n; ++i) os << Q[i] << ","; os << "])";
		return os; 
	}	
	
	long casson(); 
	long whitney();
	long cross();
	long defect();
	long unnamed();
	
};

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

class Gridle : public Model
{
public:
	
	vector<int> P, Q, S;
	
	Gridle(int N) : Model(N), P(N+1), Q(N+1), S(N*N) 
	{         
		for (int i = 0; i < n; ++i) P[i] = i;
		for (int i = 0; i < n; ++i) Q[i] = i;
	}
	
	void shuffle()
	{                
		random_shuffle(P.begin(), --P.end()); P[n]=P[0];
		random_shuffle(Q.begin(), --Q.end()); Q[n]=Q[0];
		for (int x=0; x<n; ++x) for (int y=0; y<n; ++y) 
			S[y*n + x] = rand()>>30;
	}

	ostream& output(ostream& os) const
	{ 
		os << "Gridle(["; 
		for (int i = 0; i < n; ++i) os << P[i] << ","; os << "],[";
		for (int i = 0; i < n; ++i) os << Q[i] << ","; os << "],[";
		for (int i = 0; i < n*n; ++i) os << S[i] << ","; os << "])";
		return os; 
	}	
	
	long casson(); 
	
};

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

class Gridlock : public Grid
{
public:
	
	Gridlock(int N) : Grid(N) { } 

	void randomize(unsigned int seed) 
	{
		srand(seed);
		random_shuffle(Q.begin(), --Q.end()); Q[n]=Q[0];
	}
	
	void shuffle()
	{                
		random_shuffle(P.begin(), --P.end()); P[n]=P[0];
	}	
};

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

class Star : public Model
{
public:
	
	vector<int> S;
	
	Star(int N) : Model(N), S(N*N) { }
	
	void shuffle()
	{                
		for (int x=0; x<n; ++x) for (int y=0; y<x-1; ++y) 
			S[y*n + x] = rand()>>30;
	}
	
	long casson(); 
	
};

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

class Jump : public Model
{
public:
	
	vector<double> X,Y,Z;
	
	Jump(int N) : Model(N), X(N+1), Y(N+1), Z(N+1) { }

	ostream& output(ostream& os) const
	{ 
		os << "Jump([" << setprecision(10); 
		for (int i = 0; i < n; ++i) os << X[i] << ","; os << "],[";		
		for (int i = 0; i < n; ++i) os << Y[i] << ","; os << "],[";
		for (int i = 0; i < n; ++i) os << Z[i] << ","; os << "])";
		return os; 
	}		
	
	long casson();
	long whitney();	
	long cross();
	long defect();
	long unnamed();
	
};

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

class Cube : public Jump
{
public:

	Cube(int N) : Jump(N) { } 
	
	void shuffle()
	{                
		for (int i=0; i<n; ++i)
		{
            X[i] = rand(); Y[i] = rand(); Z[i] = rand();
		}
		X[n] = X[0]; Y[n] = Y[0]; Z[n] = Z[0];
	}	
};

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

class Sphere : public Jump
{
public:

	Sphere(int N) : Jump(N) { } 
	
	void shuffle()
	{                
		for (int i=0; i<n; ++i)
		{
			Z[i] = (rand() - (1<<30)) / double(1<<30);
			double r = sqrt(1 - Z[i]*Z[i]);
			double a = rand() / double(1UL<<31) * 6.283185307179586;
			X[i] = r * cos(a);
			Y[i] = r * sin(a);
		}
		X[n] = X[0]; Y[n] = Y[0]; Z[n] = Z[0];
	}	
};

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

class Disc : public Jump
{
public:

	Disc(int N) : Jump(N) { } 
	
	void shuffle()
	{                
		for (int i=0; i<n; ++i)
		{
			Z[i] = (rand() - (1<<30)) / double(1<<30);
			double r = sqrt(rand() / double(1UL<<31));
			double a = rand() / double(1UL<<31) * 6.283185307179586;
			X[i] = r * cos(a);
			Y[i] = r * sin(a);
		}
		X[n] = X[0]; Y[n] = Y[0]; Z[n] = Z[0];
	}	
};

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

class Gaussian : public Jump
{
public:
	
	mt19937 rnd;
	normal_distribution<double> normal;		

	Gaussian(int N) : Jump(N) { }
	
	void randomize(unsigned int seed) { rnd.seed(seed); }
	
	void shuffle()
	{                
		for (int i=0; i<n; ++i)
		{			
			X[i] = normal(rnd);
			Y[i] = normal(rnd);
			Z[i] = normal(rnd);
		}
		X[n] = X[0]; Y[n] = Y[0]; Z[n] = Z[0];
	}	
};

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

#endif
