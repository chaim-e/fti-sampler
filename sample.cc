#include "model.h"

#include <unordered_map>

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

int main(int argc, char** argv)
{
	// Input
	if (argc<4) 
	{ 
		cerr << "usage: " << argv[0] << " model inv1,inv2,...,invD size [iters] [prints] [seed] [outfile]\n";
		return 1;
	}

	// Initialize parameters
	int n(atoi(argv[3]));
	int samples(argc>4 ? atoi(argv[4]) : 1000);
	int prints(argc>5 ? atoi(argv[5]) : samples / 100);
	unsigned int seed(argc>6 ? atoi(argv[6]) : time(0));
	
	// Initialize model
	string model_name(argv[1]);
	Model* model;
	if      (model_name == "petal") 	model = new Petal(n);
	else if (model_name == "grid")  	model = new Grid(n);
	else if (model_name == "gridle")  	model = new Gridle(n);
	else if (model_name == "gridlock") 	model = new Gridlock(n);
	else if (model_name == "star")  	model = new Star(n); 
	else if (model_name == "cube")  	model = new Cube(n); 
	else if (model_name == "sphere")  	model = new Sphere(n); 
	else if (model_name == "disc")  	model = new Disc(n);
	else if (model_name == "gaussian")  model = new Gaussian(n); 
	else 
	{
		cerr << "no random model named " << model_name << endl;
		return 1;
	}
	model->randomize(seed);	
	
	// Initialize invariant
	string invariant_names(argv[2]);
	int dim = count(invariant_names.begin(), invariant_names.end(), ',') + 1;
	stringstream invariant_names_stream(invariant_names);
	long (Model::*invariants[dim])();
	string name;
	for (int i = 0; (i < dim) && getline(invariant_names_stream, name, ','); ++i)
	{
		if      (name == "casson") invariants[i] = &Model::casson;
		else if (name == "jones3") invariants[i] = &Model::jones3;
		else if (name == "writhe") invariants[i] = &Model::writhe;
		else if (name == "whitney")invariants[i] = &Model::whitney; 
		else if (name == "cross")  invariants[i] = &Model::cross;
		else if (name == "defect") invariants[i] = &Model::defect;
		else if (name == "unnamed")invariants[i] = &Model::unnamed;
		else 
		{
			cerr << "no invariant named " << name << endl;
			return 1;
		}
	}

	// Summarize input 
	stringstream summary;
	summary << "model = " << model_name << ", invariants = " << invariant_names << " (" << dim << ")" <<  
			", n = " << n << ", samples = " << samples <<  ", prints = " << prints << ", seed = " << seed;
	cerr << summary.str() << endl;
	
	// Main loop
	unordered_map<string, long> count;
	time_t rawtime;
	for (int counter = 0; counter < samples; ++counter)
	{
		model->shuffle();
		string key;
		for (int i = 0; i < dim; ++i) key += to_string((model->*invariants[i])()) + ",";
		count[key] += 1;
		
		if (prints == 0) 
			cerr << *model << " : " << key << endl;
		else if ((counter+1) % prints == 0) 				
			cerr << setw(8) << count.size() << setw(12) << (counter+1) << " " << (time(&rawtime)*0 + ctime(&rawtime));
	}
	
	// Output
	ofstream fout; 
	if (argc>7) 
	{ 
		fout.open(argv[7]); 
		if (fout.fail()) { cerr << "couldn't open file" << argv[7] << endl; return 1; }
	}
	ostream* out(argc>7 ? &fout : &cout);
	*out << summary.str() << endl;
	*out << "{";
	for (auto& x : count) *out << "(" << x.first << "):" << x.second << ",";
	*out << "}" << endl;
	
	delete model;
	
	return 0;
}

