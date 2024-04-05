#include "base_cell_complex.h"
#include <map>
#include <queue>

template<class Element, class FType>
class UnionFind {


	struct UFpair {
		Element x;
		FType value;
		Element parent;
		//unsigned char rank;
	};

	map<Element, UFpair> e2ufp;

public:
		
	UnionFind() {
		e2ufp.clear();
	}

	void MakeSet(Element x, FType value) {
		UFpair p; p.x=x; p.value=value; p.parent=x;
		e2ufp[x] = p;
	}

	void Union(Element x, Element y) {
		Element xroot = Find(x);
		Element yroot = Find(y);
		if (xroot == yroot) return;

		if (e2ufp[xroot].value < e2ufp[yroot].value) {
			e2ufp[xroot].parent = yroot;
		} else {
			e2ufp[yroot].parent = xroot;
		}
	}

	Element Find(Element x) {
		Element xp = e2ufp[x].parent;	
		if (xp == x) return x;
		Element xroot = Find(xp);
		e2ufp[x].parent = xroot;
		return xroot;
	}

	void print() {
		typename map<Element, UFpair>::iterator it = e2ufp.begin();
		while (it != e2ufp.end()) {
			printf("%d(%f)->%d\n", (int) (*it).second.x, (float) (*it).second.value,
				Find((*it).second.x));
			it++;
		}
	}

};


//template<class Element, class FType, unsigned char Dim>
//Branch {
//
//	Branch<Element, FType, Dim>* parent;
//	Branch<Element, 
//
//public:
//
//	Branch(Branch* parent, vector<index_type>* geom) {
//
//
//
//	}
//
//
//
//};


template<class Element, class FType, unsigned char Dim>
class MergeTree {
	BaseCellComplex<Dim, FType>* bcc;
	UnionFind<index_type, FType> uf;


	struct sc_element {
		index_type cellid;
		unsigned int insertion_time;
		FType value;
	};

	struct sc_comparator {
		bool operator() (sc_element a, sc_element b) {
			// next value
			if (a.value > b.value) return true;
			if (a.value < b.value) return false;
			// next insertion time
			if (a.insertion_time < b.insertion_time) return true; 
			if (a.insertion_time > b.insertion_time) return false;	
			// finally cell index
			return a.cellid > b.cellid;
		}
	};


	index_type otherVertex(index_type vertex, index_type edge) {
		IndexIterator facetiter = bcc->getFacetIterator(edge);
		index_type other = *facetiter.loc;
		if (other != vertex) return other;
		facetiter++;
		return *facetiter.loc;
	};

	inline bool isLower(index_type a_id, FType a_val, index_type b_id, FType b_val) {
		if (a_val < b_val) return true;
		if (a_val > b_val) return false;
		return a_id < b_id;
	};

	priority_queue<sc_element, vector<sc_element>, sc_comparator > cell_queue;

	int seedQueueWithMinima() {

		int ct = 0;
		IndexIterator vertiter = bcc->getCellIterator(0);
		while(vertiter.isValid()){
			index_type vert = *vertiter.loc;
			FType vert_val = bcc->getValue(vert);
			unsigned char bval = bcc->getBoundary(vert);
			//printf("is %d critical?\n", vert);

			bool is_lowest = true;
			IndexIterator edgeiter = bcc->getCofacetIterator(vert);
			while (is_lowest && edgeiter.isValid()) {
				index_type other_vert = otherVertex(vert, *edgeiter.loc);
				//printf("%d's other is %d\n", vert, other_vert);
				if (bval == bcc->getBoundary(other_vert) &&
					isLower(other_vert, bcc->getValue(other_vert), vert, vert_val))
					is_lowest = false;
				edgeiter++;
			}
			if (is_lowest) {
				sc_element si;
				si.cellid = vert;
				si.value = vert_val;
				si.insertion_time = 0;
				cell_queue.push(si);
				//pqlist.push_back(vert);
				ct++;
				//printf("it IS! %d %f\n", vert, (float) vert_val);
			}
			vertiter++;
		}
		return ct;

	};




	struct MTelement {
		Element x;
		index_type parent;
		//FType value;
		unsigned char num_children;
		index_type branch;
		bool assigned;
	};


	struct Branch {
		float color[3];
		index_type root;
	};

	map<index_type, MTelement> mt;
	vector<Branch> branches;

public:
	MergeTree(BaseCellComplex<Dim, FType>* bcc) : bcc(bcc) {
		mt.clear();
	}


	void computeMT() {

		// initialize bcc with unassigneds.
		IndexIterator celliter = bcc->getCellIterator(0);	
		while (celliter.isValid()) {
			index_type cellid = *celliter.loc;
			bcc->setCritical(cellid, false);
			bcc->setAssigned(cellid, false);
			celliter++;
		}


		// initialize with potential minima
		seedQueueWithMinima();
		unsigned int itime = 1;

		// do assignment!
		while(! cell_queue.empty()) {
			sc_element si = cell_queue.top();
			cell_queue.pop();

			if (bcc->getAssigned(si.cellid)) continue;

			index_type x = si.cellid;
			bcc->setAssigned(x, true);
			uf.MakeSet(x, si.value);
			MTelement mtx; mtx.x=x; mtx.parent=x; mtx.num_children=0;
			mt[x] = mtx;

			IndexIterator fiter = bcc->getCofacetIterator(x);
			while (fiter.isValid()) {
				index_type medge = *fiter.loc;
				index_type ov = otherVertex(x, medge);

				if (! bcc->getAssigned(ov)) {
					// add to queue
					sc_element ni;
					ni.cellid = ov;
					ni.value = bcc->getValue(ov);
					ni.insertion_time = itime++;
					cell_queue.push(ni);
					fiter++; continue;
				}

				index_type ovroot = uf.Find(ov);
				if (ovroot != x) {
					// add mt edge
					mt[ovroot].parent = x;
					mt[x].num_children++;
				}
				uf.Union(x, ovroot);
				fiter++;
			}

		}
	}

	void decompose() {

		typename map<index_type, MTelement>::iterator it = mt.begin();
		while (it != mt.end()) {
			MTelement& mte = (*it).second;
			mte.assigned = false;
			it++;
		}

		it = mt.begin();
		while (it != mt.end()) {
			MTelement& mte = (*it).second;
			if (mte.assigned) { it++; continue; }

			if (mte.num_children != 1) {

				mte.assigned = true;
				mte.branch = mte.x;

				MTelement mtn = mt[mte.parent];
				while (mtn.num_children == 1 && mtn.parent != mtn.x) {
					mt[mtn.x].assigned = true;
					mt[mtn.x].branch = mte.branch;
					mtn = mt[mtn.parent];
				}
			}
			it++;
		}
	}

	void gl_fill_vertices2(vector<index_type>& verts) {
		typename map<index_type, MTelement>::iterator it = mt.begin();
		while (it != mt.end()) {
			MTelement& mte = (*it).second;
			//if (mte.num_children == 1 && mte.x != mte.parent) {
			//	it++; continue;
			//}
			verts.push_back(mte.x);	
			verts.push_back(mte.branch);
			it++;
		}
	}
	void gl_fill_vertices(vector<index_type>& verts) {
		typename map<index_type, MTelement>::iterator it = mt.begin();
		while (it != mt.end()) {
			MTelement& mte = (*it).second;
			if (mte.num_children == 1 && mte.x != mte.parent) {
				it++; continue;
			}
			verts.push_back(mte.x);	
			it++;
		}
	}
	void gl_fill_arcs(vector<index_type>& arcs) {
		typename map<index_type, MTelement>::iterator it = mt.begin();
		while (it != mt.end()) {
			MTelement& mte = (*it).second;
			if (mte.num_children == 1) {
				it++; continue;
			}
			
			MTelement mtn = mt[mte.parent];
			while (mtn.num_children == 1) {
				if (mtn.parent == mtn.x) break;
				mtn = mt[mtn.parent];
			}

			arcs.push_back(mte.x);
			arcs.push_back(mtn.x);

			it++;
		}
	}

	void dumpDot(char* filename) {
		FILE* fout = fopen(filename, "w");
		fprintf(fout, "digraph mergetree {\n");
		fprintf(fout, "\trankdir=BT;\n");

		//output nodes
		typename map<index_type, MTelement>::iterator it = mt.begin();
		while (it != mt.end()) {
			MTelement& mte = (*it).second;
			if (mte.num_children == 1 && mte.x != mte.parent) {
				it++; continue;
			}
			
			fprintf(fout, "\t%d [label=\"%.2f\"];\n", 
				mte.x, (float) bcc->getValue(mte.x));
			it++;
		}


		it = mt.begin();
		while (it != mt.end()) {
			MTelement& mte = (*it).second;
			if (mte.num_children == 1) {
				it++; continue;
			}
			
			MTelement mtn = mt[mte.parent];
			while (mtn.num_children == 1) {
				if (mtn.parent == mtn.x) break;
				mtn = mt[mtn.parent];
			}

			fprintf(fout, "\t%d -> %d [len=%f];\n", mte.x, mtn.x,
				(float) mtn.parent - (float) mte.parent);
			it++;
		}
		fprintf(fout, "}\n");
		fclose(fout);
	}
};
