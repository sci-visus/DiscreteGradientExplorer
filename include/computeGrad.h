#include "base_cell_complex.h"
#include <queue>
#define S_EPSILON 0.00001


// creates a data value for every cell in the dataset
template <unsigned char Dim, class FType>
  void setAugmentedFunction(BaseCellComplex<Dim, FType>* bcc){
  
  // set the value of each cell to be the max of its facets
  // vertices are already set, so set 1, 2, 3... cells
  for (unsigned int i = 1; i <= Dim; i++) {
    //printf("gothere\n");
    // iterate through all cells of dimension i
    IndexIterator dimiter = bcc->getCellIterator(i);
    //printf("gothere2\n");
    while (dimiter.isValid()) {
      index_type cellid = *dimiter.loc;
      //printf("gothere3\n");
      // for each cell, get max value of facets
      IndexIterator facetiter = bcc->getFacetIterator(cellid);
      //printf("gothere 4:%d %d valid=%d\n",cellid, bcc->getDim(cellid), facetiter.isValid());
      FType value = bcc->getValue(*facetiter.loc);  // value of first vertex
      FType min_value = value;
      //printf("gothere 5\n");
      while (facetiter.isValid()) {
	FType nvalue = bcc->getValue(*facetiter.loc);
	if (nvalue > value) value = nvalue;
	if (nvalue < min_value) min_value = nvalue;
	facetiter++;
      }
      bcc->setValue(cellid, value);
      bcc->setMinValue(cellid, min_value);
      //printf("set dim%d %d to %f\n", i, cellid, (float)value);
      dimiter++;
    }
  }
}; 

template <unsigned char Dim, class FType> 
  void initializeCells(BaseCellComplex<Dim, FType>* bcc){
  
  // iterate through all cells and set number of cofacets
  for (unsigned int i = 0; i <= Dim; i++) {
    IndexIterator celliter = bcc->getCellIterator(i);	
    while (celliter.isValid()) {
      index_type cellid = *celliter.loc;

      
      index_type size = bcc->getFacetIterator(cellid).size;
      bcc->setNumUFacets(cellid, size);
	  bcc->setCritical(cellid, false);
	  bcc->setAssigned(cellid, false);
	  
      //printf("setting num U cofacets: dim=%d, stuff=%d, size=%d\n", bcc->getDim(cellid), bcc->getNumUFacets(cellid), size);
      

      celliter++;
    }
  }

  
  
};

// sorted cells elements are structs to keep track of region growing sorted
// order simulation of simplicity.
template <class FType> 
  struct sc_element {
	  
    index_type cellid;
    unsigned int insertion_time;
    FType value;
    FType slope;
    unsigned char dim;	
    unsigned char boundary; // boundary ordering
    unsigned char numleft;
  };

template <class FType>
struct sc_comparator {
	  bool operator() (sc_element<FType> a, sc_element<FType> b) {
			// first do boundaries
			  //if (a.boundary < b.boundary) return true;
			  //if (a.boundary > b.boundary) return false;
			


 			  // next value
			  if (a.value > b.value) return true;
			  if (a.value < b.value) return false;
			  

			  // next do dimension - resolve higher dim stuff if possible
			  if (a.dim < b.dim) return true;
			  if (a.dim > b.dim) return false;

			  // next high slope
			  //if (a.slope < b.slope) return true;
			  //if (a.slope > b.slope) return false;


			  // next insertion time
			  if (a.insertion_time > b.insertion_time) return true; 
			  if (a.insertion_time < b.insertion_time) return false;	
			  
			  // finally cell index
			  return a.cellid > b.cellid;

		  }
	  };

template <unsigned char Dim, class FType>
inline index_type otherVertex(index_type vertex, index_type edge, BaseCellComplex<Dim, FType>* bcc) {

	 IndexIterator facetiter = bcc->getFacetIterator(edge);
	 index_type other = *facetiter.loc;
	 if (other != vertex) return other;
	 facetiter++;
	 return *facetiter.loc;
};

template<class FType>
inline bool isLower(index_type a_id, FType a_val, index_type b_id, FType b_val) {
	//printf("islower: %d %f %d %d\n", a_id, (float) a_val, b_id, (float) b_val);
	if (a_val < b_val) return true;
	if (a_val > b_val) return false;
	return a_id < b_id;

};

template <unsigned char Dim, class FType>
int seedQueueWithMinima(priority_queue<sc_element<FType>, vector<sc_element<FType> >, 
						 sc_comparator<FType> >& cell_queue, BaseCellComplex<Dim, FType>* bcc, vector<index_type>& pqlist) {

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
			index_type other_vert = otherVertex<Dim, FType>(vert, *edgeiter.loc, bcc);
			//printf("%d's other is %d\n", vert, other_vert);
			if (bval == bcc->getBoundary(other_vert) &&
				isLower<FType>(other_vert, bcc->getValue(other_vert), vert, vert_val))
				is_lowest = false;
			edgeiter++;
		}
		if (is_lowest) {
			sc_element<FType> si;
			si.cellid = vert;
			si.value = vert_val;
			si.slope = 0;
			si.dim = 0;
			si.insertion_time = 0;
			si.boundary = bcc->getBoundary(vert);
			cell_queue.push(si);
			pqlist.push_back(vert);
			ct++;
			//printf("it IS! %d %f\n", vert, (float) vert_val);
		}


		vertiter++;
	}
	return ct;

};

template <unsigned char Dim, class FType>
int decNeighUFacetCount(index_type cellid, index_type& counter, BaseCellComplex<Dim, FType>* bcc, 
						 priority_queue<sc_element<FType>, vector<sc_element<FType> >,
						 sc_comparator<FType> >& cell_queue , vector<index_type>& pqlist) {
	
    IndexIterator cfiter = bcc->getCofacetIterator(cellid);
	int ct = 0;
	while (cfiter.isValid()) {
		index_type cf = *cfiter.loc;

		if (bcc->getAssigned(cf)) {
			cfiter++; continue;
		}

		index_type numUF = bcc->getNumUFacets(cf);

		if (numUF == 0) printf("LOGICAL ERROR. should never get here\n");
		
		numUF--;
		bcc->setNumUFacets(cf, numUF);
		// add stuff to list if only one ununpaired facet
		if (numUF == 1) {

		  //get potential neighbor
		  IndexIterator asdf = bcc->getFacetIterator(cf);
		  index_type freeone;
		  while(asdf.isValid()) {
		    if (*asdf.loc != cellid && ! bcc->getAssigned(*asdf.loc))
		      freeone = *asdf.loc;
		    asdf++;
		  }
		  sc_element<FType> si;
		  si.cellid = cf;
		  si.value = bcc->getValue(cf);
		  si.slope = (si.value - bcc->getMinValue(cf)) /
		    max(S_EPSILON, bcc->distance(cf, freeone));
		  si.dim = bcc->getDim(cf);
		  si.insertion_time = counter++;
		  si.boundary = bcc->getBoundary(cf);
		  cell_queue.push(si);
		  pqlist.push_back(cf);
		  ct++;
		}

		cfiter++;
	}
	return ct;
};

template <unsigned char Dim, class FType>
void setCritical(index_type cellid, BaseCellComplex<Dim, FType>* bcc) {
	
	bcc->setAssigned(cellid, true);

	bcc->setCritical(cellid, true);
	bcc->setDimAscMan(cellid, Dim - bcc->getDim(cellid));
};

// returns cellID of pair
// sets pair, assigned, and dimA
template <unsigned char Dim, class FType>
index_type findPairAndPair(index_type cellid, BaseCellComplex<Dim, FType>* bcc) {

	IndexIterator fiter = bcc->getFacetIterator(cellid);
	index_type pair;
	unsigned char dimA = Dim;
	int found_pair = 0;
	while (fiter.isValid()) {
		index_type potential_pair = *fiter.loc;
		if (! bcc->getAssigned(potential_pair)) {
			found_pair++;
			pair = potential_pair;
		} else {
			//it IS assigned
			dimA = min(dimA, bcc->getDimAscMan(potential_pair));
		}
		fiter++;
	}
	if (found_pair != 1) printf("ERROR, found %d pairs!\n", found_pair);

	// now pair them
	bcc->setAssigned(pair, true);
	bcc->setPair(pair, cellid);
	bcc->setDimAscMan(pair, dimA);

	bcc->setAssigned(cellid, true);
	bcc->setPair(cellid, pair);
	bcc->setDimAscMan(cellid, dimA);

	return pair;
};

template <unsigned char Dim, class FType>
  void buildGradient(BaseCellComplex<Dim, FType>* bcc, vector<index_type>& order, vector<int>& counts, vector<index_type>& pqlist){

	priority_queue<sc_element<FType>, vector<sc_element<FType> >, sc_comparator<FType> > cell_queue;

	int nummin = seedQueueWithMinima<Dim, FType>(cell_queue, bcc, pqlist);
	counts.push_back(nummin);

	index_type insert_count = 1;

	while(! cell_queue.empty()) {
		sc_element<FType> si = cell_queue.top();
		//printf("%f\n", (float) bcc->getValue(si.cellid));
		cell_queue.pop();

		// if it has been assigned, skip it
		if (bcc->getAssigned(si.cellid)) continue;

		// if there are no unassigned facets, then mark it critical
		if (bcc->getNumUFacets(si.cellid) == 0) {

		  //printf("crit:%d %f\n", si.cellid, bcc->getValue(si.cellid));
			setCritical<Dim, FType>(si.cellid, bcc);
			int ct = decNeighUFacetCount<Dim, FType>(si.cellid, insert_count, bcc, cell_queue, pqlist);
			order.push_back(si.cellid);
			counts.push_back(ct);
		} else if (bcc->getNumUFacets(si.cellid) == 1) {
			// find its pair and pair it
			index_type pair = findPairAndPair(si.cellid, bcc);
			//printf("pair:%d %f - %d %f\n", si.cellid, bcc->getValue(si.cellid),
			//	 pair, bcc->getValue(pair));
			int ct1 = decNeighUFacetCount<Dim, FType>(pair, insert_count, bcc, cell_queue, pqlist);
			int ct2 = decNeighUFacetCount<Dim, FType>(si.cellid, insert_count, bcc, cell_queue, pqlist);
			order.push_back(si.cellid);
			counts.push_back(ct1);
			order.push_back(pair);
			counts.push_back(ct2);
		} else {
			printf("error, should never have nonzero or one num unassigned facets\n");
		}
	}

	//validation
	IndexIterator it = bcc->getCellIterator();
	while (it.isValid()) {
		index_type cel = *it.loc;
		if (! bcc->getCritical(cel)) {
			index_type pair = bcc->getPair(cel);
			index_type me = bcc->getPair(pair);
			if (bcc->getDim(pair) == bcc->getDim(cel)) printf("PPPLDLDLLDLL\n");
			if (me != cel) printf ("B:DSKFJ:DSLKJSDKLJFDS:L\n");
		}
		it++;
	}

};


template <unsigned char Dim, class FType> 
  void computeGradient(BaseCellComplex<Dim, FType>* bcc, vector<index_type>& order, vector<int>& counts, vector<index_type>& pqlist){
  
  printf("setting function\n");
  setAugmentedFunction<Dim, FType>(bcc);
  
  printf("initializing cells\n");
  initializeCells<Dim, FType>(bcc);

  printf("Beginning grdient computation\n");
  buildGradient<Dim, FType>(bcc, order, counts, pqlist);  

  printf("done!\n");
  
};
