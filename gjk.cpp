#include <iostream>
#include <math.h>
#include <vector>
#include <random>
#include <fstream>

#define DIMENSION 3
#define EPSILON 0.00001
using namespace std;

uniform_real_distribution<double> unif(-5,5);
default_random_engine re;


double dotProduct(vector<double> a, vector<double> b);
vector<double> negateVec(vector<double> d);
vector<double> addVec(vector<double> a, vector<double> b, double scale);
bool checkPoint(vector<double> oldVec, vector<double> newVec);
vector<double> crossProduct(vector<double>a,vector<double> b);
void printResult(bool value);

class Polytype{
public:
		int numpoints; // Num of points.
		vector<vector<double>> vertices; // vertices
		vector<double> support_pt; // support point
};

class Simplex{
public:
		int npoints=0;
		vector<vector<double>> vertices;

		bool addPoint(vector<double> a, vector<double> b, vector<double> d){
			vector<double> Minkow_pt = addVec(a, b, -1.0);
			vertices.push_back(Minkow_pt);
			npoints++;
			if(npoints-1){	
				if(checkPoint(d, Minkow_pt)){
					return false;
				}
				for(int i=0;i<npoints-1;i++){
					int c=0;
					for(int j=0;j<DIMENSION;j++){
						if(fabs(vertices[i][j]-Minkow_pt[j]) < EPSILON){
				                c++;
						}
					}
					if(c==3){
						printResult(0);
						exit(0)	;
					}
				}
			}
			return true;
		}
		vector<double> LineNormal(){
		    vector<double> p3,p5;

		    p3 = addVec(vertices[0], vertices[1], -1.0); // O(d) ~ O(1).
		    p5 = negateVec(vertices[1]); // O(d) ~ O(1).
		    
		    vector<double> p4;
		   	p4=crossProduct(p3,p5); // O(d) ~ O(1).
		   	vector<double> ans;
		   	ans=crossProduct(p4,p3); // O(d) ~ O(1).
		   	return ans;
		}

		vector<double> PlaneNormal(int a, int b, int c)
		{
		    vector<double> vec1 = addVec(vertices[c], vertices[a], -1.0);
		    vector<double> vec2 = addVec(vertices[b], vertices[a], -1.0);

		    vector<double> vec3 = crossProduct(vec1, vec2);
			vector<double> vec4 = negateVec(vertices[a]); 
		    
		    if(dotProduct(vec3, vec4) > 0)
		        return vec3;
		    else
		        return negateVec(vec3);
		}

		bool InOrigin()
		{
		   	vector<double> vec2 = PlaneNormal(0, 1, 3);
		   	vector<double> opp = addVec(vertices[2], vertices[3], -1.0);

		    if(dotProduct(vec2, opp) < 0)
		    {
		        vertices.erase(vertices.begin()+2);
		        return false;
		    }

		    vec2 = PlaneNormal(0, 2, 3);
		    opp = addVec(vertices[1], vertices[3], -1.0);

		    if(dotProduct(vec2, opp) < 0)
		    {
		        vertices.erase(vertices.begin()+1);
		        return false;
		    }

		    vec2 = PlaneNormal(1, 2, 3);
		    opp = addVec(vertices[0], vertices[3], -1.0);

		    if(dotProduct(vec2, opp) < 0)
		    {
		        vertices.erase(vertices.begin());
		        return false;
		    }
		    return true;
		}
};

void support(Polytype &body, vector<double> d);
void readinp(Polytype &b1, Polytype &b2, string filename);

int main(void){
	Polytype b1;
	Polytype b2;
	Simplex simp;
	
	string filename;
	cout <<"Enter file: ";
	cin >> filename;

	readinp(b1, b2, filename);
	
	vector<double>Direction;
	vector<double>NDirection;

	for(int i =0; i<DIMENSION; i++){
		Direction.push_back(unif(re));
		NDirection.push_back(-Direction[i]);
	}

	// First point.
	support(b1,Direction); // O(d*N1) ~ O(N1)
	support(b2,NDirection); // O(d*N2) ~ O(N2)


	// O(d) ~ O(1).
	if(!simp.addPoint(b1.support_pt, b2.support_pt, Direction)){ // First point in MKD.
		printResult(0);
		exit(0);
	}

	Direction = negateVec(simp.vertices[0]); // O(d) ~ O(1).
	NDirection = simp.vertices[0];


	// Second Point.
	support(b1, Direction);  // O(d*N1) ~ O(N1).
	support(b2, NDirection);  // O(d*N2) ~ O(N2).

	// O(d) ~ O(1).
	if(!simp.addPoint(b1.support_pt, b2.support_pt, Direction)){ // Second point in MKD.
		printResult(0);
		exit(0);
	}

	Direction = simp.LineNormal(); // O(d) ~ O(1).
	NDirection = negateVec(Direction);  // O(d) ~ O(1).


	// Third Point.
	support(b1, Direction); // O(d*N1) ~ O(N1).
	support(b2, NDirection); // O(d*N2) ~ O(N2).

	// O(d) ~ O(1).
	if(!simp.addPoint(b1.support_pt, b2.support_pt, Direction)){ // Third point in MKD.
		printResult(0);
		exit(0);
	}

	do{
		Direction = simp.PlaneNormal(0, 1, 2); // O(d) ~ O(1).
		NDirection = negateVec(Direction); // O(d) ~ O(1).

		// Fourth point
		support(b1, Direction); // O(d*N1) ~ O(N1).
		support(b2, NDirection); // O(d*N2) ~ O(N2).

		// O(d) ~ O(1).
		if(!simp.addPoint(b1.support_pt, b2.support_pt, Direction)){
			printResult(0);
			exit(0);
		}	
		
		// O(d) ~ O(1).
		if(simp.InOrigin()){
			printResult(1);
			exit(0);
		}
	}while(1);
}


void readinp(Polytype &b1, Polytype &b2, string filename){
	ifstream f(("inputs/" + filename).c_str());
	if(!f.good()){
		cout << "No file named " << filename << " found\n";
		cout << "File should be at : inputs/" << filename << endl;
		exit(0);
	}
	f.clear();
	streambuf *cinbuf = cin.rdbuf();
	ifstream in("inputs/" + filename);
    cin.rdbuf(in.rdbuf());

	//cout << "Number of points for body1\n";
	cin >> b1.numpoints;
	//cout << "\nEnter Vertices of body1"<<endl;
	double temp;
	for(int i=0; i<b1.numpoints; i++){
		vector<double> tmp;
		for(int j=0; j<DIMENSION; j++){
			cin >> temp;
			tmp.push_back(temp);
		}
		b1.vertices.push_back(tmp);
	}

	//cout << "Number of points for body2\n";
	cin >> b2.numpoints;
	//cout << "\nEnter Vertices of body2"<<endl;
	for(int i=0; i<b2.numpoints; i++){
		vector<double> tmp;
		for(int j=0; j<DIMENSION; j++){
			cin >> temp;
			tmp.push_back(temp);
		}
		b2.vertices.push_back(tmp);
	}
	cin.rdbuf(cinbuf);
	return;
}

double dotProduct(vector<double> a, vector<double> b){
	double ans=0.0;
	for(int i=0; i<a.size(); i++){
		ans += a[i] * b[i];
	}
	return ans;
}


vector<double> crossProduct(vector<double>a,vector<double> b){
    vector<double >result;
    result.push_back(a[1] * b[2] - a[2] * b[1]);
    result.push_back(a[2] * b[0] - a[0] * b[2]);
    result.push_back(a[0] * b[1] - a[1] * b[0]);
    return result;
}

vector<double> negateVec(vector<double> d){
	vector<double>nd;
		for(int i=0;i<d.size();i++){
			nd.push_back(-d[i]);
		}
	return nd;
}

vector<double> addVec(vector<double> a, vector<double> b, double scale){
	vector<double> ans;
	for(int i=0; i<a.size(); i++){
		ans.push_back(a[i]+b[i]*scale);
	}
	return ans;
}

void support(Polytype &body, vector<double> d){
	double maxi=-INFINITY;
	int index=0;
	for(int i=0; i<body.numpoints; i++){
		double dot = dotProduct(body.vertices[i], d);
		if(dot > maxi){
			maxi = dot;
			index = i;
		}
	}
	body.support_pt.assign(body.vertices[index].begin(), body.vertices[index].end());
}


bool checkPoint(vector<double> oldVec, vector<double> newVec){
	if(dotProduct(oldVec, newVec) < 0.0){
		return true;
	}
	return false;
}

void printResult(bool value){
	if(value){
		cout << "Collision Detected" << endl;
		return;
	}
	cout << "No Collision" <<endl;
}