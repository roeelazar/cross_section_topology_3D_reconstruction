#include "mex.h"

#include <limits>
#include <cmath>



#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<vector>
#include <unordered_map>

#define TRUE 1
#define FALSE 0
#define INFTY std::numeric_limits<double>::max()



using namespace std;


class UnionFind {
	
public:
	vector<int> _parent;
	vector<int> _elements;
	unordered_map<int, int> _ele2index;
	
	
	UnionFind() {}
	~UnionFind() {}
	// Add the initial elements for union-find. This function is supposed to be
	// called only once. If called twice, everything will be reset.
	void SetElements(const vector<int>& elements){
		_elements = elements;
		_ele2index.clear();
		_parent.resize(_elements.size(), 0);
		for (int i = 0; i < _elements.size(); ++i) {
			_ele2index[_elements[i]] = i;
			_parent[i] = i;
		}
	}	
	// Add additional element for union-find. It can be useful if initial elements
	// are unkown in the beginning.
	void AddElement(int ele) {
		_ele2index[ele] = _elements.size();
		_elements.push_back(ele);
		_parent.push_back(_parent.size());
	}
	int Find(int ele) { return Find_index(_ele2index[ele]); }
	int Find_index(int ind) {
		if (ind != _parent[ind]) {
			_parent[ind] = Find_index(_parent[ind]);
		}
		return _parent[ind];
	}
	
	// Union with elements.
	void Union(int ele0, int ele1) { Union_roots(Find(ele0), Find(ele1)); }
	// Union with index instead of elements.
	void Union_indexes(int ind0, int ind1) {
		_parent[Find_index(ind0)] = Find_index(ind1);
	}
	// Union with root instead of elements or index.
	void Union_roots(int ind0, int ind1) { _parent[ind0] = ind1; }
	
	int getIndex(int ele) {
		// if (_ele2index.find(ele) == _ele2index.end()) return -1;
		return _ele2index[ele];
	}
	int getElement(int ind) {
		// if (ind >= _elements.size()) return -1;
		return _elements[ind];
	}
	void ExtractComponents(vector<vector<int> >& comps){
		comps.clear();
		unordered_map<int, int> root2comp;
		for (int i = 0; i < _elements.size(); ++i) {
			int r = Find_index(i);
	  // A new comp is found.
			if (root2comp.find(r) == root2comp.end()) {
				root2comp[r] = comps.size();
				comps.push_back(vector<int>(1,_elements[i]));
			} else {
				comps[root2comp[r]].push_back(_elements[i]);
			}
		}
	}
	
	void ExtractSingleComponents(int ele, vector<int>&comp){
		comp.clear();	
		int ind = getIndex(ele);
		int r_be = Find_index(ind);
		for (int i = 0; i < _elements.size(); ++i) {
			if(r_be == Find_index(i)){
				comp.push_back(i);
			}
		}
	}
	int GetNumOfComponents(){
		int nComp = 0;
		for (int i=0; i<_parent.size(); ++i) {
			if (i == _parent[i]) {
				nComp++;
			}
		}
		return nComp;
	}
};




class minCut{
	
	
private:
	vector<int>tmp_sets;
	

	
	
	vector<vector<double> >W;	// the weight matrix
	vector<int>Del;
	int n;				// the number of veritces
	int m;				// the number of edges
	
	UnionFind unifind;
	
	int recordSource;
	int recordTarget;

public:
	
	minCut(){};
	
private:	
	double minCutPhase(int V){
		int i = 0, j = 0;
		int s[2];
		if(V  == 2) {
			for( i = 0; i < n; i++){
				if(Del[i] == FALSE){
					s[j] = i; j++;
				}
			}
			recordSource = s[0];
			recordTarget = s[1];
			return W[s[0]][s[1]];
		}
		double L[n];
		int T[n];
		memset(L, 0, n*sizeof(double));
		memset(T, FALSE, n*sizeof(int));
		i = 1;	
		j = 0;
		int v,u;
		while( i <= V){
			v = maxStickiness(T,L);
			T[v] = TRUE;
			for(u = 0; u < n; u++){
				if(W[v][u] != 0 && Del[u] == FALSE && T[u] == FALSE){
					L[u] = L[u] + W[u][v];
				}
			}
			if( i >= V-1){
				s[j] = v; j++;
			}
			i++;
		}
		recordSource = s[0];
		recordTarget = s[1];
		return L[s[1]];	
	}
	
	void merge(int s, int t){
		int v = 0;
		for(v = 0; v < n; v++){
			if(Del[v] == FALSE && v != s && v!= t){
				W[s][v] = W[s][v] + W[v][t];
				W[v][s] = W[s][v];
			}
		}
		Del[t] = TRUE;
		unifind.Union(s,t);
	}
	
	int maxStickiness(int *T, double *L){
		int i = 0;
		int v = 0;
		double max = 0;
		//for(int k=0;k<n;++k)if(T[k] == FALSE){v=k;max = L[k];break;}	
		for(i = 0; i < n; i++){
			if(Del[i] == FALSE && T[i] == FALSE && max <= L[i]){
				v = i;
				max = L[i];
			} 
		}
		return v;
	}
	
	void ExtractCurrentComponents(double *sets){
		
		
		unifind.ExtractSingleComponents(recordTarget,tmp_sets);
		

		for(int i=0;i<n;++i)sets[i] = 0.0;
		for(auto a:tmp_sets)sets[a] = 1.0;
	
	
	
	}

	
public:
	
	void StoerWagner(double *in_Wmatrix, mwSize n_node, double *mincut, double *sets){
		int V  = n = n_node;
		double C = INFTY;
		Del.resize(n_node);
		memset(Del.data(), FALSE, n*sizeof(int));
		W.resize(n_node,vector<double>(n_node));
		for(int i=0;i<n_node;++i){
			for(int j=0;j<n_node;++j){
				W[i][j] = in_Wmatrix[i*n_node+j];
			}
		}
		
		vector<int>ele(n_node);for(int i=0;i<n_node;++i)ele[i] = i;
		unifind.SetElements(ele);
		
		
		
		for(V = n; V > 1; V--){
			double cutValue = minCutPhase(V);
			if(C>cutValue){
				C = cutValue;
				if(sets){					
					ExtractCurrentComponents(sets);
				}
			}
			merge(recordSource,recordTarget);
		}
		
		if(mincut)*mincut = C;

		return;
	}
	
	
	
};



void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])
{
	if(nrhs > 1) {
		
		mexErrMsgIdAndTxt("StoerWagner:nrhs","At most one inputs.");
		
	}
	
	
	if(nlhs > 2) {
		mexErrMsgIdAndTxt("StoerWagner", "Too many output arguments, expected 1 - 2");
	}
	
	
	
	
	
	double *inMatrix;       /* 1xN input matrix */
	mwSize ncols;           /* size of matrix */
	
	
	double *outMatrix;      /* output matrix */
	double *outWeight;      /* output matrix */
	
	
	
	/* get the value of the scalar input  */
	
	
	
	/* create a pointer to the real data in the input matrix  */
	inMatrix = (double *)mxGetData(prhs[0]);
	
	
	/* get dimensions of the input matrix */
	ncols = mxGetN(prhs[0]);
	
	
	/* create the output matrix */
	if(nlhs > 0)plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
	if(nlhs > 1)plhs[1] = mxCreateDoubleMatrix(ncols,1,mxREAL);
	
	
	/* get a pointer to the real data in the output matrix */
	outWeight = nlhs > 0?(double *)mxGetData(plhs[0]):NULL;
	outMatrix = nlhs > 1?(double *)mxGetData(plhs[1]):NULL;
	
	
	minCut mc;
	mc.StoerWagner(inMatrix,ncols,outWeight,outMatrix);
	
	
	/* call the computational routine */
	//arrayProduct(multiplier,inMatrix,outMatrix,ncols);
	
	
}
