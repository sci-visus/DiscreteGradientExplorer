#ifndef MSCOMPLEX
#define MSCOMPLEX
#include "base_cell_complex.h"
#include <map>

//#include <ext/hash_map>
#include <queue>

// this class models an edge that has an upper and lower linked list
//template <class type>
//class d_linked_node {
//
//	d_linked_node* u_prev;
//	d_linked_node* u_next;
//	d_linked_node* l_prev;
//	d_linked_node* l_next;
//
//public:
//	type element;
//	
//	d_linked_node(type element) : element(element) {
//		u_prev = u_next = l_prev = l_next = NULL;
//	}
//
//	void remove() {
//		if (prev != NULL) {
//			prev->next = next;
//		} 
//		if (next != NULL) {
//			next->prev = prev;
//		}
//	}
//	void insert(d_linked_node* p, d_linked_node* n) {
//		if (p != NULL) {
//			p->next = this;
//		} 
//		if (n != NULL) {
//			n->prev = this;
//		}
//		next = n; prev = p;
//	}
//	void insert_before(d_linked_node* n) {
//		if (n == NULL) return insert(NULL, n);
//		insert(n->prev, n);
//	}
//	d_linked_node* next() {
//		return next;
//	}
//	d_linked_node* prev() {
//		return prev;
//	}
//};


// a simple list node structure.  
template <class type>
class list_node {
	list_node* ln_next;
	type ln_element;
public:
	list_node(type element) : ln_element(element), ln_next(NULL) {}

	void insert_before(list_node* e) {
		if (e == NULL) return;
		ln_next = e;
	}

	type element() {
		return ln_element;
	}

	list_node* next() {
		return ln_next;
	}
};

struct geometry {
	geometry* merged[3];
	vector<index_type> base;
};

template <class FType>
class node;

template <class FType>
class arc {
public:
	node<FType>* lower;
	node<FType>* upper;
	geometry* geom;
	FType persistence;
	int created;
	int destroyed;

	arc(node<FType>* upper, node<FType>* lower, vector<index_type>* v) {
		this->created = -1;
		this->destroyed = -1;
		this->upper = upper;
		if (this->upper != NULL) {
			list_node< arc<FType>* >* na = new list_node< arc<FType>* >(this);
			na->insert_before(this->upper->arcs);
			this->upper->arcs = na;
		}

		this->lower = lower;
		if (this->lower != NULL) {
			list_node< arc<FType>* >* na = new list_node< arc<FType>* >(this);
			na->insert_before(this->lower->arcs);
			this->lower->arcs = na;
		}		
		this->persistence = upper->value - lower->value;
		if(this->persistence< 0) printf("ERROR, PERSISTence = %f\n", this->persistence);
		this->geom = new geometry();
		this->geom->base.assign(v->begin(), v->end());
	}

	// create from a1.upper->a1.lower=a2.lower->a2.upper==a3.upper->a3.lower
	arc(arc<FType>* a1, arc<FType>* a2, arc<FType>* a3, int created ) {
		printf("%f -> %f=%f -> %f=%f -> %f\n", a1->upper->value, a1->lower->value, 
			a2->lower->value, a2->upper->value, 
			a3->upper->value, a3->lower->value);
		this->created = created;
		this->destroyed = -1;
		this->upper = a1->upper;
		if (this->upper != NULL) {
			list_node< arc<FType>* >* na = new list_node< arc<FType>* >(this);
			na->insert_before(this->upper->arcs);
			this->upper->arcs = na;
		}

		this->lower = a3->lower;
		if (this->lower != NULL) {
			list_node< arc<FType>* >* na = new list_node< arc<FType>* >(this);
			na->insert_before(this->lower->arcs);
			this->lower->arcs = na;
		}		
		this->persistence = this->upper->value - this->lower->value;
		this->geom = new geometry();
		this->geom->merged[0] = a1->geom;
		this->geom->merged[1] = a2->geom;
		this->geom->merged[2] = a3->geom;
	}


};

template <class FType>
class node {
public:
	index_type id;
	list_node< arc<FType>* >* arcs;	
	unsigned char dim;
	FType value;
	int destroyed;
	int boundary;

	node(index_type id, unsigned char dim, FType value) :
	id(id), dim(dim), value(value), destroyed(-1) {
		arcs = NULL;
	}

};

using namespace std;

template <unsigned char Dim, class FType>
class MSComplex {
    int msc_num_cancelled;
public:
	// store minimum and maximum values of critical points
	FType mymin, mymax;

	// the number canceled and a structure to translate
	// a number canceled to persistence canceled
	unsigned int num_canceled;
	vector< FType > num_to_value_canceled;

	// the nodes and arcs
	map< index_type, node<FType>* > nodes;
	vector< arc<FType>* > arcs;

	// underlying mesh storing discrete gradient
	BaseCellComplex<Dim, FType>* bcc;


	void msc_rec_trace_down(node<FType>* n, index_type id, vector<index_type>* path) {
		//printf("did=%d\n", bcc->getDim(id));
		path->push_back(id);
		//printf("s= %d\n", path->size());
		// end of recursion
		if (bcc->getCritical(id) || nodes.count(id) != 0) {
			// create a new arc
			printf("s= %d\n", path->size());
			arc<FType>* a = new arc<FType>(n, nodes[id], path);
			arcs.push_back(a);
		} else {

			// recurse on facets of end of arrow
			index_type pair = bcc->getPair(id);
			if (bcc->getDim(pair) != bcc->getDim(id) + 1) {
				// then this is the "head" of an arrow, so path ends
				path->pop_back();
				return;
				//&&
				//bcc->getDim(pair) != bcc->getDim(id) - 1) {
				//printf("ERROR pair=%d id=%d\n", bcc->getDim(pair), bcc->getDim(id));
				//path->pop_back();		return;
			}
			IndexIterator fit = bcc->getFacetIterator(pair);
			while (fit.isValid() ) {
				if (bcc->getDim(pair) != bcc->getDim(*fit.loc) +1)
					printf("WHOA2 ndim=%d iddim=%d cdim=%d fdim=%d\n", 
					n->dim, bcc->getDim(pair), bcc->getDim(*fit.loc));
/*				if (n->dim == (Dim - (bcc->getDimAscMan(*fit.loc) - 1))) {
					fit++; continue;
				}*/		
				index_type nid = *fit.loc;
				if (nid != id) {
					msc_rec_trace_down(n, nid, path);
				}
				fit++;
			}
		}
		path->pop_back();
	}

	void msc_trace_down(node<FType>* n) {
		if (n->dim == 0) {
			return;	
		}

		vector<index_type> v;
		v.push_back(n->id);

		// seed the downward traversal in each direction
		IndexIterator fit = bcc->getFacetIterator(n->id);
		while (fit.isValid()) {
			if (n->dim != bcc->getDim(*fit.loc) +1)
				printf("WHOA ndim=%d cdim=%d fdim=%d\n", 
				n->dim, bcc->getDim(n->id), bcc->getDim(*fit.loc));
			//if (n->dim == (Dim - (bcc->getDimAscMan(*fit.loc) - 1))) {
			//	fit++; continue;
			//}
			printf("n=%d s=%d\n", n->dim, bcc->getDim(*fit.loc));
			 msc_rec_trace_down(n,*fit.loc, &v);
			fit++;
		}
	}



	struct sortedArc{
		arc<FType>* a;
		
		bool operator()(const sortedArc& _Left, const sortedArc& _Right) {
			if(_Left.a->persistence < _Right.a->persistence) return false;
			//if(_Left.persistence > _Right.persistence) return true;
			return true;
			//return _Left.a > _Right.a;
		}
	};
	
	
	int msc_arc_count_multiplicity(arc<FType>* a) {
		int counter = 0;
		list_node< arc<FType>* >* larcs = a->lower->arcs;
		while (larcs != NULL) {
			if (larcs->element()->lower == a->lower &&
				larcs->element()->upper == a->upper) counter++;
			larcs = larcs->next();
		}
		return counter;
	}

	bool msc_arc_valid_to_cancel(arc<FType>* a) {
		if (a->destroyed != -1) return false;
		if (a->lower->boundary != a->upper->boundary) return false;
		int counter = this->msc_arc_count_multiplicity(a);
		//printf("cc = %d\n", counter);
		if(counter != 1) return false;

		return true;
	}

	void msc_cancel(arc<FType>* a, 
		priority_queue<sortedArc, vector< sortedArc>, sortedArc > &tocancel){
		
			bool res = this->msc_arc_valid_to_cancel(a);
			if (! res) return;
		// is "a" valid to cancel?
		
			msc_num_cancelled++;
			num_to_value_canceled.push_back(a->persistence);
			//cancel
			node<FType>* n;

			list_node< arc<FType>* >* larcs = a->lower->arcs;
			while (larcs != NULL) {

				if (larcs->element()->destroyed != -1 ||
					larcs->element()->upper->dim != a->upper->dim ||
					larcs->element() == a) {
						larcs = larcs->next(); 
						continue;
				}
				//larcs->element()->destroyed = msc_num_cancelled;

				list_node< arc<FType>* >* uarcs = a->upper->arcs;
			
				while (uarcs != NULL) {

					if(uarcs->element()->destroyed != -1 ||
						uarcs->element()->lower->dim != a->lower->dim ||
						uarcs->element() == a) {
							uarcs = uarcs->next();
							continue;
					}

					//new arc
					arc<FType>* ua = uarcs->element();
					arc<FType>* la = larcs->element();
					arc<FType>* na = new arc<FType>(la, a, ua, msc_num_cancelled);
					arcs.push_back(na);
					sortedArc nna = { na };
					if (nna.a->persistence < 0){
						printf("ERROR: p = %f\n", nna.a->persistence);
					}
					tocancel.push(nna);
					uarcs = uarcs->next();
				}
				larcs = larcs->next();
			}
			larcs = a->lower->arcs;
			while (larcs != NULL) {
				if(larcs->element()->destroyed == -1)
					larcs->element()->destroyed = msc_num_cancelled;
				larcs = larcs->next(); 
			}
			list_node< arc<FType>* >* uarcs = a->upper->arcs;
			while (uarcs != NULL) {
					if (uarcs->element()->destroyed == -1)
						uarcs->element()->destroyed = msc_num_cancelled;
					uarcs = uarcs->next();
			}
			a->lower->destroyed = msc_num_cancelled;
			a->upper->destroyed = msc_num_cancelled;
	}




	MSComplex(BaseCellComplex<Dim, FType>* bcc) : bcc(bcc) {
		this->num_canceled = 0;
	}

	void computeComplex() {
		printf("computing complex!\n");
		// add all nodes to the complex, record min,max

		bool firsttime = true;
		int ncounts[4]; for (int i=0;i<=Dim;i++) ncounts[i]=0;
		for (int i = 0; i <= Dim; i++) {
			IndexIterator mit = bcc->getCellIterator(i);
			if (firsttime) {
				this->mymin = this->mymax = bcc->getValue(*mit.loc);
				firsttime = false;
			}
			while (mit.isValid()) {

				// current id
				index_type m = *mit.loc;

				// current value, set global min max
				FType nv = bcc->getValue(m);
				if (nv > this->mymax) this->mymax = nv;
				if (nv < this->mymin) this->mymin = nv;

				// if its critical, make new critical point
				if (bcc->getCritical(m)) {
					ncounts[bcc->getDim(m)]++;
					node<FType>* nn = 
						new node<FType>(m, bcc->getDim(m), nv);
					nn->boundary = bcc->getBoundary(m);
					this->nodes[m] = nn;	
				}
				mit++; 
			}
		}
		printf("Found Nodes!:\n");
		for (int i=0; i<=Dim; i++)
			printf("\t# Index-%d = %d\n", i, ncounts[i]);

		// add all arcs
		// iterate over all nodes and trace down
		typename map< index_type, node<FType>* >::iterator it;
		for(it= nodes.begin(); it!= nodes.end(); it++) {
			// current node
			node<FType>* n = ((node<FType>*) (*it).second);
			printf("tracing from %d\n", n->dim);
			int barf = this->arcs.size();
			msc_trace_down(n);
			printf("found %d arcs\n",  this->arcs.size() - barf);
		}
		
		printf("MS Complex: %d nodes, %d arcs\n", nodes.size(), arcs.size());
	}

	void computeHeirarchy(){
		msc_num_cancelled = 0;
		printf("Computing MSC hierarchy...\n");
		priority_queue<sortedArc, vector< sortedArc>, sortedArc > tocancel;
		for (unsigned int i = 0; i < this->arcs.size(); i++) {
			sortedArc na = { arcs[i] };
			if (na.a->persistence < 0){
				printf("ERROR: p = %f\n", na.a->persistence);
			}
			tocancel.push(na);
		}

		while(! tocancel.empty()) {
			printf("tocance.size=%d\n", tocancel.size());
			arc<FType>* a = tocancel.top().a;
			tocancel.pop();

			this->msc_cancel(a, tocancel);
			//return;
		}
	}

	void msc_fill_geometry(geometry* geom, bool direction, vector<index_type>& v) {
		if(geom->base.size() == 0) {
			this->msc_fill_geometry(geom->merged[0], direction, v);
			this->msc_fill_geometry(geom->merged[1], ! direction, v);
			this->msc_fill_geometry(geom->merged[2], direction, v);
		} else {
			if (direction) {
				for (int i = 0; i < geom->base.size(); i++) {
					if (v.size() > 0 && v[v.size()-1] == geom->base[i]) {
						v.pop_back();
					} else {
						v.push_back(geom->base[i]);
					}
				}
			} else {
				for (int i = geom->base.size() - 1; i >= 0; i--) {
					if (v.size() > 0 && v[v.size()-1] == geom->base[i]) {
						v.pop_back();
					} else {
						v.push_back(geom->base[i]);
					}
				}
			}
		}

	}
	void fill_geometry(geometry* geom, vector<index_type>& v){
		msc_fill_geometry(geom, true, v);
	}
	void dump_obj(char* fname) {
		FILE* f = fopen(fname, "w");
		
		int counter = 0;
		map<index_type, int> index_map;


	}


};

#endif