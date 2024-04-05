#ifndef UNSTRUCTURED_CELL_COMPLEX
#define UNSTRUCTURED_CELL_COMPLEX


#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <set>
#include <algorithm>
#include <map>
#include <cmath>

using namespace std;

#include "base_cell_complex.h"











template <class FType> class Simplex {
public:
   Simplex() :
      verts(0), numberOfVerts(0),
      iFacets(0,0), iCofacets(0,0), 
      fValue(0), dim(0), is_bdry(0) {};

   Simplex(FType v, unsigned char d, bool b) :
      verts(0), numberOfVerts(0),
      iFacets(0,0), iCofacets(0,0), 
      fValue(v), dim(d), is_bdry(b) {};

   Simplex(const Simplex &s) :
      iFacets(s.iFacets), iCofacets(s.iCofacets),
      fValue(s.fValue), dim(s.dim), is_bdry(s.is_bdry) 
   {
      numberOfVerts = s.numberOfVerts;
      if (numberOfVerts != 0) {
         verts = new index_type[numberOfVerts];
         for (int i=0; i<numberOfVerts; i++) {
            verts[i] = s.verts[i];
         }
      } else {
         verts = 0;
      }
   };



   FType function_value() {
      return fValue;
   }

   bool is_boundary() {
      return is_bdry;
   }

   unsigned char dimension() {
      return dim;
   };

   inline index_type* getFacet(index_type f) {
      return iFacets.begin+f;
   };
   inline index_type* getCofacet(index_type cf) {
      return iCofacets.begin+cf;
   };
   inline IndexIterator& getFacetIterator() {
      return iFacets;
   };
   inline IndexIterator& getCofacetIterator() {
      return iCofacets; 
   };


   index_type* verts;
   index_type numberOfVerts;

   IndexIterator iFacets;
   IndexIterator iCofacets;

   FType fValue;
   FType fMinValue;
   unsigned char dim;
   bool is_bdry;
   bool is_assigned;
   index_type pair;
   index_type numUFacets;
   bool is_critical;

   inline bool isCritical() {
     return is_critical;
     //max(iFacets.size, iCofacets.size) == pair;
   };

   inline void setCritical(bool value) {
     is_critical = value;
     //if (value) pair = max(iFacets.size, iCofacets.size);
   };

   inline bool operator<(const Simplex<FType> &s) const;
   inline Simplex<FType>& operator=(const Simplex<FType> &s);
};


template <class FType>
inline bool Simplex<FType>::operator<(
      const Simplex<FType> &s) const 
{
   if (numberOfVerts != s.numberOfVerts) {
      return numberOfVerts < s.numberOfVerts;
   } else {
      for (int i=0; i<numberOfVerts; i++) {
         if (verts[i] == s.verts[i]) {
            continue;
         } else {
            return verts[i] < s.verts[i];
         }
      }
   }

   return false;
}

template <class FType>
inline Simplex<FType>& Simplex<FType>::operator=(
      const Simplex<FType> &s) 
{
   if (this != &s) {
      this->Simplex<FType>::~Simplex<FType>();
      new (this) Simplex<FType>(s);
   }
   return *this;
};










template <unsigned char Dim, class FType> class UnstructuredSimplicialComplex : 
   public BaseCellComplex<Dim, FType> {
   
public:
   //constructors
   UnstructuredSimplicialComplex() :
      verts(0), numberOfVerts(0), cells(0), numberOfCells(0),
      facets(0), numberOfFacets(0), cofacets(0), numberOfCofacets(0) 
      {};

   //loaders
   void loadFromObj(const char* filename, float scale = 1.0);
   void loadFromOff(const char* filename, float scale = 1.0);
   void addLowDimSimplices(set<Simplex<FType> > setOfCell);
   void loadFromBrick(const char* filename, 
      int size_x, int size_y, int size_z, 
      int shift = 1, float scale = 1.0);



   //data members
   Simplex<FType>* cells;
   index_type numberOfCells;
   //[0],...,[Dim-1] are only off a particular dimensions
   //[Dim] is for all the facets
   index_type* cell_ids;
   IndexIterator cellOffsets[Dim+2];

   index_type* facets;
   index_type numberOfFacets;
 
   index_type* cofacets;
   index_type numberOfCofacets;



   //vertex array, stores positions
   BaseVertex<Dim>* verts;
   index_type numberOfVerts;





	// access functions to cells
   //all cells is the default
   //all cells of a particular dimension otherwise
   inline IndexIterator getCellIterator(unsigned char dim = Dim) {
      return cellOffsets[dim];
   }
   
   //facets of a cell
   inline IndexIterator getFacetIterator(index_type cellid) {
      return cells[cellid].iFacets;
   }
   //cofacets of a cell
   inline IndexIterator getCofacetIterator(index_type cellid) {
      return cells[cellid].iCofacets;
   }


   // function value of a cell
   inline void setValue(index_type cellid, FType value) {
      cells[cellid].fValue = value;
   };
	inline FType getValue(index_type cellid) { 
      return cells[cellid].fValue; 
   };

	
   // function value of a cell
   inline void setMinValue(index_type cellid, FType value) {
      cells[cellid].fMinValue = value;
   };
	inline FType getMinValue(index_type cellid) { 
      return cells[cellid].fMinValue; 
   };

	
	
	// dimension of a cell
	inline unsigned char getDim(index_type cellid) { 
      return cells[cellid].dimension();
   };




	
	// set and get the boundary value - this is a byte value, however,
	// the possible values are 0 through MAX_DIM
	inline void setBoundary(index_type cellid, unsigned char value) { 
		cells[cellid].is_bdry = value;
	}; 
	inline unsigned char getBoundary(index_type cellid) { 
		return cells[cellid].is_bdry; 
	};

	// set and get the assigned bit - used to mark if something has
	// been processed by the algorithm. 
	inline void setAssigned(index_type cellid, bool value) {
		cells[cellid].is_assigned = value;
	};
	inline bool getAssigned(index_type cellid) { 
		return cells[cellid].is_assigned;
	};

	// set and get the number of unpaired facets: could be 0 through
	// the maximum number of facets of any cell. NOTE, this is
	// only valid while the cell is NOT assigned 
	inline void setNumUFacets(index_type cellid, index_type value) {
		cells[cellid].numUFacets = value;
	};
	inline index_type getNumUFacets(index_type cellid) { 
		return cells[cellid].numUFacets; 
	};

	// set and get the pair (offset) of a cell: could be 0 through
	// the maximum number of cofacets of any cell. NOTE this is only
	// valid AFTER the cell is assigned.
	inline void setPair(index_type cellid, index_type value) {		
		cells[cellid].pair = value;
	};
	inline index_type getPair(index_type cellid) { 
		return cells[cellid].pair;
	};

	// set and get whether the cell is critical: Note this is only 
	// valid after the cell is assigned. 
	inline void setCritical(index_type cellid, bool value) {
		cells[cellid].setCritical(value);
	};
	inline bool getCritical(index_type cellid) { 
		return cells[cellid].isCritical();
	};

	// set and get the dimension of the ascending manifold this
	// cell is part of: range is 0 to MAX_DIM, NOTE only valid
	// after a cell as been assigned
	inline void setDimAscMan(index_type cellid, unsigned char value) {
		cells[cellid].numUFacets = value;
	};
	inline unsigned char getDimAscMan(index_type cellid) { 
		return cells[cellid].numUFacets;
	};




   inline double distance(index_type cellid0, index_type cellid1) {
      BaseVertex<Dim> c0, c1;

      if (cellid0 < 0 || cellid1 < 0 ||
            cellid0 > numberOfCells || cellid1 > numberOfCells) {
         cerr << "distance(): invalid cellid passed in: ";
         cerr << cellid0 << " " << cellid1 << endl;
         return 9e99;
      }

      //compute centroids of each cell
      for (int i=0; i<Dim; i++) {
         for (int j=0; j<cells[cellid0].numberOfVerts; j++) {
            c0.position[i] += verts[cells[cellid0].verts[j]].position[i];
         }
         c0.position[i] /= cells[cellid0].numberOfVerts;
      }

      for (int i=0; i<Dim; i++) {
         for (int j=0; j<cells[cellid1].numberOfVerts; j++) {
            c1.position[i] += verts[cells[cellid1].verts[j]].position[i];
         }
         c1.position[i] /= cells[cellid1].numberOfVerts;
      }

      //compute distance

      double d = 0;

      for (int i=0; i<Dim; i++) {
         d += (c0.position[i] - c1.position[i])*(c0.position[i] - c1.position[i]);
      }

      return sqrt(d);
   }

};


template <unsigned char Dim, class FType>
void UnstructuredSimplicialComplex<Dim, FType>::loadFromObj(const char* filename, float scale /*= 1.0*/) {

	char buffer[128];
	//make sure we're clean
	numberOfVerts = 0;
	if (verts != 0) {
		delete[] verts;
	}

	numberOfCells = 0;
	if (cells != 0) {
		delete[] cells;
	}

	printf("about to read...\n");
	// get number of verts and number of faces
	numberOfVerts = 0;
	int numberOfMaxDimCells = 0;
	set<Simplex<FType> > setOfCell;

	char firstchar;
	FILE* f = fopen(filename, "r");
	while (! feof(f)) {
		fscanf(f, "%c %s\n", &firstchar, buffer);
		//printf("%c %s\n", firstchar, buffer);
		if (firstchar == 'v') numberOfVerts++;
		if (firstchar == 'f') numberOfMaxDimCells++;
	}
	fclose(f);
	printf("counted %d verts, %d faces\n", numberOfVerts,
		numberOfMaxDimCells);

	verts = new BaseVertex<Dim>[numberOfVerts];

	f = fopen(filename, "r");
	for (int counter = 0; counter < numberOfVerts; counter++) {   
		fscanf(f, "%c ", &firstchar);
		for (int i = 0; i < Dim; i++)
			fscanf(f, "%f ", &(verts[counter].position[i]));

		float value;
		fscanf(f, "%f\n", &value);

		Simplex<FType> uc;
		uc.numberOfVerts = 1;
		uc.verts = new index_type[uc.numberOfVerts];
		uc.verts[0] = counter;
		uc.fValue = scale*value;
		uc.dim = 0;
		setOfCell.insert(uc);




	}
	for (int counter = 0; counter < numberOfMaxDimCells; counter++) {   
		fscanf(f, "%c ", &firstchar);

		int coords[Dim + 1];
		for (int i = 0; i < Dim + 1; i++)
			fscanf(f, "%d ", &(coords[i]));

		fscanf(f, "\n");

		//reading cells

		Simplex<FType> uc;
		uc.numberOfVerts =  Dim + 1;
		uc.verts = new index_type[uc.numberOfVerts];
		uc.dim = Dim;

		for (int i=0; i<uc.numberOfVerts; i++) {
			uc.verts[i] = coords[i] - 1;
		}

		sort(uc.verts, uc.verts+uc.numberOfVerts);

		setOfCell.insert(uc);

	}
	fclose(f);

	addLowDimSimplices(setOfCell);

}



template <unsigned char Dim, class FType>
void UnstructuredSimplicialComplex<Dim, FType>::loadFromOff(const char* filename, float scale /*= 1.0*/) {
   ifstream input(filename);


   //make sure we're clean
   numberOfVerts = 0;
   if (verts != 0) {
      delete[] verts;
   }

   numberOfCells = 0;
   if (cells != 0) {
      delete[] cells;
   }



   //this is temporary storage we will allocate and delete
   int numberOfMaxDimCells = 0;
   set<Simplex<FType> > setOfCell;



   if (input.is_open()) {
      //tokenize, line by line
      //stage 0 = header, 1 = verts, 2 = cells;
      int stage = 0;
      int v_index = 0;
      int c_index = 0;
      
      
      while(input.good() && stage < 3) {
         string line;
         getline(input, line);
         //clean dos endlines too
         int cr_idx = line.find('\r', 0);
         if (cr_idx != string::npos) {
            line = line.substr(0, cr_idx);
         }

         vector<string> tokens;
         string token;
         istringstream iss(line);
         while (getline(iss, token, ' ')) {
            tokens.push_back(token);
         }

         if (tokens.size() > 0) {
            //skip comments
            if (tokens[0][0] != '#') {
               if (stage == 0) {
                  //try to read the header
                  if (tokens[0] == "OFF") {
                     if (tokens.size() == 1) {
                        //read the next line
                        tokens.clear();
                        getline(input, line);
                        istringstream iss2(line);
                        while (getline(iss2, token, ' ')) {
                           tokens.push_back(token);
                        }

                        //read the counts
                        numberOfVerts = atoi(tokens[0].c_str());
						if (numberOfVerts == 0) exit(1);
                        numberOfMaxDimCells = atoi(tokens[1].c_str());

                     } else {
                        //read the counts
                        numberOfVerts = atoi(tokens[0].c_str());
                        numberOfMaxDimCells = atoi(tokens[1].c_str());
                     }

                     //cout << numberOfVerts << " " << numberOfMaxDimCells << endl;
                     stage = 1;
                  }
               } else if (stage == 1) {
                  //reading verts

                  //initialize storage
                  if (v_index == 0) {
                     verts = new BaseVertex<Dim>[numberOfVerts];
                  }

                  //read until we've hit the cap
                  for (int i=0; i<Dim; i++) {
                     float v = atof(tokens[i].c_str());
                     verts[v_index].position[i] = v;
                  }

                  Simplex<FType> uc;
                  uc.numberOfVerts = 1;
                  uc.verts = new index_type[uc.numberOfVerts];
                  uc.verts[0] = v_index;
                  uc.fValue = scale*atof(tokens[tokens.size()-1].c_str());
                  uc.dim = 0;
                  setOfCell.insert(uc);

                  v_index++;

                  if (v_index == numberOfVerts) {
                     stage = 2;
                  }
               } else if (stage == 2) {
                  //reading cells

                  Simplex<FType> uc;
                  uc.numberOfVerts =  atoi(tokens[0].c_str());
                  uc.verts = new index_type[uc.numberOfVerts];
                  uc.dim = uc.numberOfVerts-1;

                  for (int i=0; i<uc.numberOfVerts; i++) {
                     uc.verts[i] = atoi(tokens[i+1].c_str());
                  }

                  sort(uc.verts, uc.verts+uc.numberOfVerts);

                  setOfCell.insert(uc);

                  c_index++;

                  if (c_index == numberOfMaxDimCells) {
                     stage = 3;
                  }
               }
            }
         }
      }

   } else {
      cerr << "File " << filename << " could not be opened" << endl;
   }


   input.close();

	addLowDimSimplices(setOfCell);

}
   
template <unsigned char Dim, class FType>
void UnstructuredSimplicialComplex<Dim, FType>::addLowDimSimplices(
	set<Simplex<FType> > setOfCell) {


   //at this point we have facets varied and Vertex positions
   //let's build the boundary faces to make it a complex!
   //In non-simplicial settings, this step needs to be replaced with 
   //all boundaries explicit


   for (int i=Dim; i>1; i--) {
      typename set< Simplex<FType> >::iterator iter = setOfCell.begin();
      while (iter != setOfCell.end()) {
         if (iter->dim == i) {
            //add in the boundary faces;

            for (int j=0; j<iter->numberOfVerts; j++) {
               //simplicial
               //add in the face that takes out vertex j

               Simplex<FType> uc;
               uc.numberOfVerts = iter->numberOfVerts-1;
               uc.verts = new index_type[uc.numberOfVerts];
               uc.dim = i-1;

               int v_index = 0;
               for (int k=0; k<iter->numberOfVerts; k++) {
                  if (k != j) {
                     uc.verts[v_index] = iter->verts[k];
                     v_index++;
                  }
               }

               setOfCell.insert(uc);
            }
         }


         iter++;
      }
   }


   //setOfCell now has all cells, in order to boot
   //the sort on these cells is increasing dimension, 
   //and then by vertex orderings
   //cout << setOfCell.size() << endl;

   numberOfCells = setOfCell.size();
   cells = new Simplex<FType>[numberOfCells];
   cell_ids = new index_type[numberOfCells];

   map<Simplex<FType>, int> mapOfCells;

   typename set< Simplex<FType> >::iterator iter = setOfCell.begin();
   int cell_index = 0;

   index_type* begin = cell_ids;
   index_type iter_size = 0;
   char current_dim = 0;

   while (iter != setOfCell.end()) {
      cells[cell_index] = *iter;
      mapOfCells[*iter] = cell_index;
      cell_ids[cell_index] = cell_index;


      if (cells[cell_index].dim != current_dim) {
         //initialize the iterator
         cellOffsets[current_dim] = IndexIterator(begin,iter_size);

         current_dim = cells[cell_index].dim;
         begin = begin + iter_size;
         iter_size = 0;
      }

      iter_size++;
      cell_index++;
      iter++;
   }

   //set the last one for simplices in setOfCells
   cellOffsets[current_dim] = IndexIterator(begin, iter_size);
   current_dim++;

   //if there happen to be no higher dimensional cells left, make
   //blank iterators to the end of the array
   while (current_dim < Dim+1) {
      cellOffsets[current_dim] = IndexIterator(begin+iter_size, 0);
      current_dim++;
   }

   //and then finally set the one for all simplices
   cellOffsets[current_dim] = IndexIterator(cell_ids, numberOfCells);


   cout << "Simplicial complex read with the following counts" << endl;
   cout << "for elements of each dimension: " << endl;
   for (int i=0; i<Dim+1; i++) {
      cout << i << "  " << cellOffsets[i].size << endl;
   }
   cout << "Total number of elements: " << cellOffsets[Dim+1].size << endl;


   /* Debug code

   for (int i=0; i<Dim+2; i++) {
      cout << (cellOffsets[i].begin - cell_ids) ;
      cout << " " << cellOffsets[i].size << endl;
   }

   */




   //let's build the adjacencies

   vector<vector<index_type> > fct_vec;
   vector<vector<index_type> > cofct_vec;
   fct_vec.resize(numberOfCells);
   cofct_vec.resize(numberOfCells);
   
   //increase counts
   numberOfFacets = 0;
   numberOfCofacets = 0;

   iter = setOfCell.begin();

   while (iter != setOfCell.end()) {
      if (iter->dim > 0) {
         //add each facet, and to the cofacet list too
         for (int i=iter->numberOfVerts-1; i>=0; i--) {
            Simplex<FType> s;
            s.numberOfVerts = iter->numberOfVerts-1;
            s.verts = new index_type[s.numberOfVerts];
            int v_count = 0;
            for (int j=0; j<iter->numberOfVerts; j++) {
               if (i != j) {
                  s.verts[v_count] = iter->verts[j];
                  v_count++;
               }
            }
            fct_vec[mapOfCells[*iter]].push_back(mapOfCells[s]);
            cofct_vec[mapOfCells[s]].push_back(mapOfCells[*iter]);

            numberOfFacets++;
            numberOfCofacets++;

            delete[] s.verts;
         }
      }

      iter++;
   }


   ///allocate some space for these, mark in the vertices points
   facets = new index_type[numberOfFacets];
   cofacets = new index_type[numberOfCofacets];


   int f_index = 0;

   for (int i=0; i<fct_vec.size(); i++) {
      if (!fct_vec[i].empty()) {
         index_type* begin = facets+f_index;
         index_type size = fct_vec[i].size();
         for (int j=0; j<size; j++) {
            facets[f_index] = fct_vec[i][j];
            f_index++;
         }
         
         //set the iterator
         cells[i].iFacets = IndexIterator(begin, size);
      }
   }

   f_index = 0;

   for (int i=0; i<cofct_vec.size(); i++) {
      if (!cofct_vec[i].empty()) {
         index_type* begin = cofacets+f_index;
         index_type size = cofct_vec[i].size();
         for (int j=0; j<size; j++) {
            cofacets[f_index] = cofct_vec[i][j];
            f_index++;
         }
         
         //set the iterator
         cells[i].iCofacets = IndexIterator(begin, size);
      }
   }


   //mark boundary!

   //start with the cells of dimenion d-1
   IndexIterator ii = getCellIterator(Dim-1);

   while (ii.isValid()) {
      //check its cofacets
      IndexIterator ii2 = getCofacetIterator(*ii.loc);

      if (ii2.size != 2) {
         //it is a boundary element!
         cells[*ii.loc].is_bdry = true;
      }

      ii++;
   }


   //now mark all cells based on adjacency
   for (int i=Dim-2; i>=0; i--) {
      IndexIterator ii = getCellIterator(i);

      while (ii.isValid()) {
         //check its cofacets
         IndexIterator ii2 = getCofacetIterator(*ii.loc);

         while(ii2.isValid()) {
            if (cells[*ii2.loc].is_bdry) {
               cells[*ii.loc].is_bdry = true;
               break;
            }
            ii2++;
         }

         ii++;
      }
   }






   /* Debug Codea


   for (int i=0; i<numberOfCells; i++) {
//      if (cells[i].is_bdry) {
         cout << i << " ";
         cout << int(cells[i].dim) << " ";

         for (int j=0; j<cells[i].numberOfVerts; j++) {
            cout << cells[i].verts[j] << ",";
         }
         cout << endl;
//      }
   }


   for (int i=0; i<numberOfCells; i++) {
      IndexIterator ii = cells[i].iFacets;
      if (ii.isValid()) {
         cout << i << " ";
         while (ii.isValid()) {
            cout << *ii.loc << ",";
            ii++;
         }
         cout << endl;
      }
   }

   for (int i=0; i<numberOfCells; i++) {
      IndexIterator ii = cells[i].iCofacets;
      if (ii.isValid()) {
         cout << i << " ";
         while (ii.isValid()) {
            cout << *ii.loc << ",";
            ii++;
         }
         cout << endl;
      }
   }

   */

}





template <unsigned char Dim, class FType>
void UnstructuredSimplicialComplex<Dim, FType>::loadFromBrick(const char* filename,
      int size_x, int size_y, int size_z, 
      int shift /* = 1*/, float scale /*= 1.0*/) {
   ifstream input(filename);


   //make sure we're clean
   numberOfVerts = 0;
   if (verts != 0) {
      delete[] verts;
   }

   numberOfCells = 0;
   if (cells != 0) {
      delete[] cells;
   }



   //this is temporary storage we will allocate and delete
   set<Simplex<FType> > setOfCell;

   FILE* file = fopen(filename, "rb");
   if (file == 0) 
   {
      fprintf(stderr, "loadFromBrick: Error loading file %s\n", filename);
      exit(1);
   }

   numberOfVerts = size_x*size_y*size_z;

   float* data = new float[numberOfVerts];
   fread(data, sizeof(float), numberOfVerts, file);
   fclose(file);
   
   size_x /= shift;
   size_y /= shift;
   size_z /= shift;
   
   numberOfVerts = size_x*size_y*size_z;

   //this is the number of edges
   int numberOfMaxDimCells = (size_x-1)*size_y*size_z;
   numberOfMaxDimCells += size_x*(size_y-1)*size_z;
   numberOfMaxDimCells += size_x*size_y*(size_z-1);

   numberOfCells = numberOfVerts + numberOfMaxDimCells;

   cout << numberOfVerts << " " << numberOfCells << endl;

   //first build the verts
   verts = new BaseVertex<Dim>[numberOfVerts];

   //next add the cells
   cells = new Simplex<FType>[numberOfCells];
   cell_ids = new index_type[numberOfCells];


   int c_index = 0;

   for (int k=0; k<size_z; k++) {
      for (int j=0; j<size_y; j++) {
         for (int i=0; i<size_x; i++) {
            int i0 = shift*i;
            int j0 = shift*j;
            int k0 = shift*k;

            int v_index = k0*size_y*size_x + j0*size_x + i0;
            //int v_index = k*size_y*size_x + j*size_x + i;
            
            verts[c_index].position[0] = i0;
            verts[c_index].position[1] = j0;
            verts[c_index].position[2] = k0;

                  
            Simplex<FType> uc;
            uc.numberOfVerts = 1;
            uc.verts = new index_type[uc.numberOfVerts];
            uc.verts[0] = c_index;
            uc.fValue = scale*data[v_index];
            uc.dim = 0;
            
            cells[c_index] = uc;
            c_index++;
         }
      }
   }

   delete[] data;




   numberOfFacets = 0;
   numberOfCofacets = 0;
   vector<vector<index_type> > fct_vec;
   vector<vector<index_type> > cofct_vec;
   fct_vec.resize(numberOfCells);
   cofct_vec.resize(numberOfCells);


   //add left cells

   for (int k=0; k<size_z; k++) {
      for (int j=0; j<size_y; j++) {
         for (int i=0; i<size_x-1; i++) {
            int v_index = k*size_y*size_x + j*size_x + i;
            int v_index2 = k*size_y*size_x + j*size_x + i+1;

            Simplex<FType> uc;
            uc.numberOfVerts = 2;
            uc.verts = new index_type[uc.numberOfVerts];
            uc.dim = uc.numberOfVerts-1;

            uc.verts[0] = v_index;
            uc.verts[1] = v_index2;

            cells[c_index] = uc;

            fct_vec[c_index].push_back(uc.verts[0]);
            fct_vec[c_index].push_back(uc.verts[1]);
            cofct_vec[uc.verts[0]].push_back(c_index);
            cofct_vec[uc.verts[1]].push_back(c_index);

            numberOfFacets += 2;
            numberOfCofacets += 2;

            c_index++;
         }
      }
   }

   //add up cells

   for (int k=0; k<size_z; k++) {
      for (int j=0; j<size_y-1; j++) {
         for (int i=0; i<size_x; i++) {
            int v_index = k*size_y*size_x + j*size_x + i;
            int v_index2 = k*size_y*size_x + (j+1)*size_x + i;

            Simplex<FType> uc;
            uc.numberOfVerts = 2;
            uc.verts = new index_type[uc.numberOfVerts];
            uc.dim = uc.numberOfVerts-1;

            uc.verts[0] = v_index;
            uc.verts[1] = v_index2;

            cells[c_index] = uc;

            fct_vec[c_index].push_back(uc.verts[0]);
            fct_vec[c_index].push_back(uc.verts[1]);
            cofct_vec[uc.verts[0]].push_back(c_index);
            cofct_vec[uc.verts[1]].push_back(c_index);

            numberOfFacets += 2;
            numberOfCofacets += 2;

            c_index++;
         }
      }
   }

   //add back cells

   for (int k=0; k<size_z-1; k++) {
      for (int j=0; j<size_y; j++) {
         for (int i=0; i<size_x; i++) {
            int v_index = k*size_y*size_x + j*size_x + i;
            int v_index2 = (k+1)*size_y*size_x + j*size_x + i;

            Simplex<FType> uc;
            uc.numberOfVerts = 2;
            uc.verts = new index_type[uc.numberOfVerts];
            uc.dim = uc.numberOfVerts-1;

            uc.verts[0] = v_index;
            uc.verts[1] = v_index2;

            cells[c_index] = uc;

            fct_vec[c_index].push_back(uc.verts[0]);
            fct_vec[c_index].push_back(uc.verts[1]);
            cofct_vec[uc.verts[0]].push_back(c_index);
            cofct_vec[uc.verts[1]].push_back(c_index);

            numberOfFacets += 2;
            numberOfCofacets += 2;

            c_index++;
         }
      }
   }



   //at this point we have facets varied and Vertex positions
   //let's build the boundary faces to make it a complex!
   //In non-simplicial settings, this step needs to be replaced with 
   //all boundaries explicit
/*

   for (int i=Dim; i>1; i--) {
      typename set< Simplex<FType> >::iterator iter = setOfCell.begin();
      while (iter != setOfCell.end()) {
         if (iter->dim == i) {
            //add in the boundary faces;

            for (int j=0; j<iter->numberOfVerts; j++) {
               //simplicial
               //add in the face that takes out vertex j

               Simplex<FType> uc;
               uc.numberOfVerts = iter->numberOfVerts-1;
               uc.verts = new index_type[uc.numberOfVerts];
               uc.dim = i-1;

               int v_index = 0;
               for (int k=0; k<iter->numberOfVerts; k++) {
                  if (k != j) {
                     uc.verts[v_index] = iter->verts[k];
                     v_index++;
                  }
               }

               setOfCell.insert(uc);
            }
         }


         iter++;
      }
   }*/


   //setOfCell now has all cells, in order to boot
   //the sort on these cells is increasing dimension, 
   //and then by vertex orderings
   //cout << setOfCell.size() << endl;
   
   //first set up the cell offsets and ids
   char current_dim = 0;
   index_type* begin = cell_ids;
   index_type iter_size = 0;

   for (int i=0; i<numberOfCells; i++) {
      cell_ids[i] = i;
      if (cells[i].dim != current_dim) {
         cellOffsets[current_dim] = IndexIterator(begin,iter_size);

         current_dim = cells[i].dim;
         begin = begin + iter_size;
         iter_size = 0;
      }

      iter_size++;
   }
   //set the last one for simplices in setOfCells
   cellOffsets[current_dim] = IndexIterator(begin, iter_size);
   current_dim++;

   //if there happen to be no higher dimensional cells left, make
   //blank iterators to the end of the array
   while (current_dim < Dim+1) {
      cellOffsets[current_dim] = IndexIterator(begin+iter_size, 0);
      current_dim++;
   }

   //and then finally set the one for all simplices
   cellOffsets[current_dim] = IndexIterator(cell_ids, numberOfCells);

/*

   numberOfCells = setOfCell.size();
   cells = new Simplex<FType>[numberOfCells];
   cell_ids = new index_type[numberOfCells];

   map<Simplex<FType>, int> mapOfCells;

   typename set< Simplex<FType> >::iterator iter = setOfCell.begin();
   int cell_index = 0;

   index_type* begin = cell_ids;
   index_type iter_size = 0;
   char current_dim = 0;

   while (iter != setOfCell.end()) {
      cells[cell_index] = *iter;
      mapOfCells[*iter] = cell_index;
      cell_ids[cell_index] = cell_index;


      if (cells[cell_index].dim != current_dim) {
         //initialize the iterator
         cellOffsets[current_dim] = IndexIterator(begin,iter_size);

         current_dim = cells[cell_index].dim;
         begin = begin + iter_size;
         iter_size = 0;
      }

      iter_size++;
      cell_index++;
      iter++;
   }

   //set the last one for simplices in setOfCells
   cellOffsets[current_dim] = IndexIterator(begin, iter_size);
   current_dim++;

   //if there happen to be no higher dimensional cells left, make
   //blank iterators to the end of the array
   while (current_dim < Dim+1) {
      cellOffsets[current_dim] = IndexIterator(begin+iter_size, 0);
      current_dim++;
   }

   //and then finally set the one for all simplices
   cellOffsets[current_dim] = IndexIterator(cell_ids, numberOfCells);

*/
   cout << "Simplicial complex read with the following counts" << endl;
   cout << "for elements of each dimension: " << endl;
   for (int i=0; i<Dim+1; i++) {
      cout << i << "  " << cellOffsets[i].size << endl;
   }
   cout << "Total number of elements: " << cellOffsets[Dim+1].size << endl;


   /* Debug code

   for (int i=0; i<Dim+2; i++) {
      cout << (cellOffsets[i].begin - cell_ids) ;
      cout << " " << cellOffsets[i].size << endl;
   }

   */




   //let's build the adjacencies
   ///allocate some space for these, mark in the vertices points
   facets = new index_type[numberOfFacets];
   cofacets = new index_type[numberOfCofacets];


   int f_index = 0;

   for (int i=0; i<fct_vec.size(); i++) {
      if (!fct_vec[i].empty()) {
         index_type* begin = facets+f_index;
         index_type size = fct_vec[i].size();
         for (int j=0; j<size; j++) {
            facets[f_index] = fct_vec[i][j];
            f_index++;
         }
         
         //set the iterator
         cells[i].iFacets = IndexIterator(begin, size);
      }
   }

   f_index = 0;

   for (int i=0; i<cofct_vec.size(); i++) {
      if (!cofct_vec[i].empty()) {
         index_type* begin = cofacets+f_index;
         index_type size = cofct_vec[i].size();
         for (int j=0; j<size; j++) {
            cofacets[f_index] = cofct_vec[i][j];
            f_index++;
         }
         
         //set the iterator
         cells[i].iCofacets = IndexIterator(begin, size);
      }
   }


   //mark boundary!

   //start with the cells of dimenion d-1
   IndexIterator ii = getCellIterator(Dim-1);

   while (ii.isValid()) {
      //check its cofacets
      IndexIterator ii2 = getCofacetIterator(*ii.loc);

      if (ii2.size != 2) {
         //it is a boundary element!
         cells[*ii.loc].is_bdry = true;
      }

      ii++;
   }


   //now mark all cells based on adjacency
   for (int i=Dim-2; i>=0; i--) {
      IndexIterator ii = getCellIterator(i);

      while (ii.isValid()) {
         //check its cofacets
         IndexIterator ii2 = getCofacetIterator(*ii.loc);

         while(ii2.isValid()) {
            if (cells[*ii2.loc].is_bdry) {
               cells[*ii.loc].is_bdry = true;
               break;
            }
            ii2++;
         }

         ii++;
      }
   }

}


#endif
