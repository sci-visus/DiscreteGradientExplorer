#ifndef BASE_CELL_COMPLEX
#define BASE_CELL_COMPLEX




#include "index_type.h"
#include "index_iterator.h"


template <unsigned char Dim> class BaseVertex {
public:
   BaseVertex() {
      for (int i=0; i<Dim; i++) {
         position[i] = 0;
      }
   };

   float position[Dim];
};






template <unsigned char Dim, class FType> class BaseCellComplex {
   
public:
	// access functions to cells
   //all cells is the default
   //all cells of a particular dimension otherwise
   virtual IndexIterator getCellIterator(unsigned char dim = Dim) = 0;
   
   //facets of a cell
   virtual IndexIterator getFacetIterator(index_type cellid) = 0;
   //cofacets of a cell
   virtual IndexIterator getCofacetIterator(index_type cellid) = 0;



   //vertex array, stores positions
   BaseVertex<Dim>* verts;
   index_type numberOfVerts;





    // function value of a cell
	virtual void setValue(index_type cellid, FType value) = 0;
	virtual FType getValue(index_type cellid) = 0;

	// minimum function value of verts of a cell
	virtual void setMinValue(index_type cellid, FType value) = 0;
	virtual FType getMinValue(index_type cellid) = 0;
	
	// dimension of a cell
	virtual unsigned char getDim(index_type cellid) = 0;
	
	// set and get the boundary value - this is a byte value, however,
	// the possible values are 0 through MAX_DIM
	virtual void setBoundary(index_type cellid, unsigned char value) = 0;
	virtual unsigned char getBoundary(index_type cellid) = 0;

	// set and get the assigned bit - used to mark if something has
	// been processed by the algorithm. 
	virtual void setAssigned(index_type cellid, bool value) = 0;
	virtual bool getAssigned(index_type cellid) = 0;

	// set and get the number of unpaired cofacets: could be 0 through
	// the maximum number of cofacets of any cell. NOTE, this is
	// only valid while the cell is NOT assigned 
	virtual void setNumUFacets(index_type cellid, index_type value) = 0;
	virtual index_type getNumUFacets(index_type cellid) = 0;

	// set and get the pair (offset) of a cell: could be 0 through
	// the maximum number of cofacets of any cell. NOTE this is only
	// valid AFTER the cell is assigned.
	virtual void setPair(index_type cellid, index_type value) = 0;
	virtual index_type getPair(index_type cellid) = 0;

	// set and get whether the cell is critical: Note this is only 
	// valid after the cell is assigned. 
	virtual void setCritical(index_type cellid, bool value) = 0;
	virtual bool getCritical(index_type cellid) = 0;

	// set and get the dimension of the ascending manifold this
	// cell is part of: range is 0 to MAX_DIM, NOTE only valid
	// after a cell as been assigned
	virtual void setDimAscMan(index_type cellid, unsigned char value) = 0;
	virtual unsigned char getDimAscMan(index_type cellid) = 0;


   virtual double distance(index_type cellid0, index_type cellid1) = 0;

};



#endif
