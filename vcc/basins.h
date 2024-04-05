#include "base_cell_complex.h"
#include <map>

//#include <ext/hash_map>
#include <queue>

#define INT_INFTY 4294967295


template <class FType>
struct basin {
  index_type minimum;
  float color[4];
  unsigned int destroyed;
  basin<FType>* parent;
  FType value;
  unsigned int vcount;

};

template <class FType>
struct basin_edge {
	basin<FType>* b1;
	basin<FType>* b2;
	FType persistence;
	FType value;



};


using namespace std;

template <unsigned char Dim, class FType>
class basins {
	
	bool basin_LT(const basin<FType>& left, const basin<FType>& right) {
		if(left.value < right.value) return false;
		if(left.value > right.value) return true;
		return left.minimum > right.minimum;
	}	


	template<class FType2>
	struct edge_LT {
		bool operator()(const basin_edge<FType2>& left, const basin_edge<FType2>& right) {
			if(left.persistence < right.persistence) return false;
			if(left.persistence > right.persistence) return true;
			if(left.b1 < right.b1) return false;
			return true;
		}
	};

	BaseCellComplex<Dim, FType>* bcc;
	priority_queue<basin_edge<FType> ,vector<basin_edge<FType> >, edge_LT<FType> > be_sorter; 	
	map<index_type, basin<FType>*> min_to_basin_map;
	map<index_type, index_type> cell_to_min_map;

	void rec_asc_man_assign(index_type start_cell, index_type curr_cell) {

	}
	unsigned int g_dest_count;
	vector<FType> g_int_to_f_pers;

	basin<FType>* createBasin(index_type m) {
		basin<FType>* be = new basin<FType>();	
		be->minimum = m;
		be->parent = NULL;
		be->destroyed = INT_INFTY;
		be->value = bcc->getValue(m);
		for (int i = 0; i < 3; i++) 
			be->color[i] = 0.01f * (rand()%100) + 0.01f;
		be->color[3] = 1.0f;
		return be;
	}
	
	void fill_edge(index_type m1, index_type m2, 
		index_type eid, basin_edge<FType> &e) {
		basin<FType>* b1 = min_to_basin_map[m1];
		basin<FType>* b2 = min_to_basin_map[m2];
		basin<FType>* bu; basin<FType>* bl;
		if (basin_LT(*b1, *b2)) {
			bl = b2; bu = b1;
		} else {
			bu = b2; bl = b1;
		}
		e.value = bcc->getValue(eid);
		e.persistence = e.value - bu->value;
		e.b1 = b1; e.b2 = b2;
	}

	index_type rec_trace_down(index_type vert) {
		if (cell_to_min_map.count(vert) != 0) {
			return cell_to_min_map[vert];	
		}


		vector<index_type> v;
		
		while(cell_to_min_map.count(vert) == 0) {
		  if (bcc->getCritical(vert)) {
		    cell_to_min_map[vert] = vert;
		    continue;
		  }

		  v.push_back(vert);

		  index_type e = bcc->getPair(vert);
		  IndexIterator fit = bcc->getFacetIterator(e);
		  while (fit.isValid()) {
		    index_type nv = *fit.loc;
		    if (nv != vert) {
		      vert = nv;
		      break;
		    }
		    fit++;
		  }
		}
		index_type res = cell_to_min_map[vert];
		for (int i = 0; i < v.size(); i++) {
		  cell_to_min_map[v[i]]=res;
		}
		return res;
		
	}

	index_type rec_trace_down2(index_type vert) {
		if (cell_to_min_map.count(vert) != 0) {
			return cell_to_min_map[vert];	
		}


		if (bcc->getCritical(vert)) {
			cell_to_min_map[vert]=vert;
			return vert;
		}
		index_type e = bcc->getPair(vert);
		IndexIterator fit = bcc->getFacetIterator(e);
		while (fit.isValid()) {
			index_type nv = *fit.loc;
			if (nv != vert) {
				index_type badsf = rec_trace_down(nv);
				cell_to_min_map[vert] = badsf;
				cell_to_min_map[e] = badsf;
				return badsf;
			}
			fit++;
		}
		printf("WARNING: trace down failed!\n");
		return 0;
	}
	
	index_type rec_assign_down(index_type cel) {
			//printf("%d\n", cel);
//		printf("gg\n");
			if (bcc->getDim(cel) == 0 && bcc->getCritical(cel)) {
				cell_to_min_map[cel] = cel;
				return cel;
			}

		if (cell_to_min_map.count(cel) != 0) {
			//printf("ee\n");
			return cell_to_min_map[cel];
		}
				//printf ("a\n");

		int ind = bcc->getDim(cel);
		//printf("ind=%d\n", ind);
		if (bcc->getCritical(cel)) {
					//		printf ("b\n");

			IndexIterator fit = bcc->getFacetIterator(cel);
			index_type res = rec_assign_down(*fit.loc);
			cell_to_min_map[cel] = res;
							//printf ("b\n");
			return res;
		} else if (ind > bcc->getDim(bcc->getPair(cel))) {
					//		printf ("c\n");

			// then return the res of first facet of this head
			FType val = bcc->getValue(cel);
			//printf("cc\n");
			IndexIterator fit = bcc->getFacetIterator(cel);
			//printf("ccc\n");
			while (fit.isValid()) {
				//printf("ccccc\n");
				index_type ncel = *fit.loc;
					//printf("ii\n");
				if (ncel != bcc->getPair(cel) && bcc->getValue(ncel) <= val) {
					//printf("hh\n");
					index_type res = rec_assign_down(ncel);
					cell_to_min_map[cel] = res;
								//printf ("c\n");
				return res;
				}
						//printf("jj\n");
			fit++;
				//		printf("cccc\n");

			}
			//printf("barf\n");
		} else {
				//	printf ("d\n");
			//printf("%d\n", cel);
			//this is tail of arrow
			index_type res = rec_assign_down(bcc->getPair(cel));
			cell_to_min_map[cel] = res;
			//printf("d\n");
			return res;
		}		
	}

	FType mymin;
	FType mymax;

	inline index_type otherVertex(index_type vertex, index_type edge) {

		IndexIterator facetiter = bcc->getFacetIterator(edge);
		index_type other = *facetiter.loc;
		if (other != vertex) return other;
		facetiter++;
		return *facetiter.loc;
	};



public: 



	void dumpBasinCountPerPersistence(char* filename) {
	  
	  FILE* fout = fopen(filename, "w");
	  unsigned int s = this->g_int_to_f_pers.size();
	  for (int i = 0; i < s; i++)
	    fprintf(fout, "%u %f\n", 1 + s - i, this->g_int_to_f_pers[i]);
	  fclose(fout);

	}

	void dumpVertBCount(char* filename, unsigned int cancelcount) {

		FILE* fout = fopen(filename, "wb");

		IndexIterator it = bcc->getCellIterator(0);
		while (it.isValid()) {
			index_type vert = *it.loc;

			index_type my_min; // = rec_assign_down(vert);
			if (cell_to_min_map.count(vert) == 0) {
			  rec_assign_down(vert);
			}
			basin<FType>* b = min_to_basin_map[cell_to_min_map[vert]];
			while (b->destroyed < cancelcount) 
			  b = b->parent;
			my_min = b->minimum;

			index_type marray[8];
			unsigned char numm = 1;
			marray[0] = my_min;

			unsigned char num_neg = 0;
			IndexIterator fit = bcc->getCofacetIterator(vert);
			
			// count number of neighbors with different vals
			while (fit.isValid()) {
				index_type edge = *fit.loc;
				index_type ov = this->otherVertex(vert, edge);


				index_type ov_min; // = rec_assign_down(ov);
				
				if (cell_to_min_map.count(ov) == 0) {
				  rec_assign_down(ov);
				}
				basin<FType>* b2 = min_to_basin_map[cell_to_min_map[ov]];
				while (b2->destroyed < cancelcount) 
				  b2 = b2->parent;
				ov_min = b2->minimum;

				bool isin = false;
				for (unsigned char aa = 0; aa < numm; aa++)
				  isin = isin || (ov_min == marray[aa]);

				if (! isin) {
				  marray[numm] = ov_min;
				  numm++;
				}

				//if (ov_min != my_min) num_neg++;
				
				fit++;
			}
			num_neg = numm-1;

			fwrite(&num_neg, sizeof(unsigned char), 1, fout);

			it++;
		}
		fclose(fout);
	}

	void dumpBasinSizes(char* filename, unsigned int cancelcount) {
	  int barf;
	
	  typename map<index_type, basin<FType>* >::iterator it;

	  for(it= min_to_basin_map.begin(); it!= min_to_basin_map.end(); it++) {
	    ((basin<FType>*) (*it).second)->vcount = 0;
	  }

	  IndexIterator vit = bcc->getCellIterator(0);
	  while(vit.isValid()) {
	    index_type cel = *vit.loc;
	    
	    if (cell_to_min_map.count(cel) == 0) {
	      rec_assign_down(cel);
	    }
	    basin<FType>* b = min_to_basin_map[cell_to_min_map[cel]];
	    while (b->destroyed < cancelcount) 
	      b = b->parent;

	    b->vcount++;
	    vit++;
	  }

	  FILE* fout = fopen(filename, "w");	    
	  for(it= min_to_basin_map.begin(); it!= min_to_basin_map.end(); it++) {
	    unsigned int vcount = ((basin<FType>*) (*it).second)->vcount;
	    if (vcount != 0)
	      fprintf(fout, "%u\n", vcount);
	  }
	  
	  fclose(fout);
	  
  	}
	

	//void dumpBasinFDist(char* filename, unsigned int cancelcount) {
	//  FILE* fout = fopen(filename, "w");
	//  unsigned int s = this->g_int_to_f_pers.size();
	//  for (int i = cancelcount; i < s; i++)
	//    fprintf(fout, "%u %f\n", 1 + s - i, this->g_int_to_f_pers[i]);
	//  fclose(fout);
	//  
 // 	}

	
	unsigned int getDestrCount(float percent){
	  FType value = ((float) mymax - (float) mymin) * percent;
	  for (int i = 0; i < this->g_int_to_f_pers.size(); i++)
	    if (g_int_to_f_pers[i] > value) return i;
	  return g_int_to_f_pers.size();
	}
	
	

	int num_destroyed() {
	  return g_dest_count;
	}
	
	float* cellColor(index_type cel) {
	  if (cell_to_min_map.count(cel) == 0) {
			rec_assign_down(cel);
			//printf("ERROR DID NOT FIND MAP\n");
		}
		basin<FType>* b = min_to_basin_map[cell_to_min_map[cel]];
		//printf("bmin=%d\n", b->minimum);
		return b->color;
	}
	float livingCellValue(index_type cel, unsigned int destroyed) {
	  if (cell_to_min_map.count(cel) == 0) {
	    rec_assign_down(cel);
	    //printf("ERROR DID NOT FIND MAP\n");
	  }
	  basin<FType>* b = min_to_basin_map[cell_to_min_map[cel]];
	  while (b->destroyed < destroyed) 
	    b = b->parent;
	  //printf("bmin=%d\n", b->minimum);
	  return b->value;
	}


	float* livingCellColor(index_type cel, unsigned int destroyed) {
		if (cell_to_min_map.count(cel) == 0) {
			rec_assign_down(cel);
			//printf("ERROR DID NOT FIND MAP\n");
		}
		basin<FType>* b = min_to_basin_map[cell_to_min_map[cel]];
		while (b->destroyed < destroyed) 
			b = b->parent;
		//printf("bmin=%d\n", b->minimum);
		return b->color;
	}
	basins(BaseCellComplex<Dim, FType>* bcc) {
		this->bcc = bcc;
		this->g_dest_count = 0;

	}

	void computeBasins() {
		printf("creating basins...\n");
		// for each minimum create a basin
		int basin_count = 0;
		cell_to_min_map.clear();
		IndexIterator mit = bcc->getCellIterator(0);
		mymin = mymax = bcc->getValue(*mit.loc);

		while (mit.isValid()) {
			index_type m = *mit.loc;
			
			// set global min max
			FType nv = bcc->getValue(m);
			if (nv > mymax) mymax = nv;
			if (nv < mymin) mymin = nv;

			if (bcc->getCritical(m)) {
				basin_count++;
				min_to_basin_map[m] = createBasin(m);
				cell_to_min_map[m] = m;
			}
			mit++; 
		}
		printf("found %d basins!\n", basin_count);

		printf("creating connections...\n");
		// for each 1-saddle
		IndexIterator it = bcc->getCellIterator(1);
		while (it.isValid()) {
			//printf("doing %u\n", *it.loc);
			index_type sad = *it.loc;
			if (! bcc->getCritical(sad)) {
				it++; continue;
			}

			//	printf ("1\n");
			// trace down each lower path
			index_type mins[2];
			int mincount = 0;
			IndexIterator fit = bcc->getFacetIterator(sad);
			while (fit.isValid() && mincount < 2) {
			//	printf ("2\n");
				index_type vert = *fit.loc;
				cell_to_min_map[sad] = this->rec_trace_down(vert); //rec_assign_down(vert);
			//	printf ("3\n");
				if (this->cell_to_min_map.count(vert) == 0) {
					printf("EROROORROORO\n");
				}
				mins[mincount] = this->cell_to_min_map[vert];
//printf ("4\n");
				//rec_trace_down(vert);
				mincount++;	fit++;
			}
			// now mins contains mins at end of downpaths from this saddle
				//	printf ("5\n");
			basin_edge<FType> ed;
			fill_edge(mins[0], mins[1], sad, ed);
			be_sorter.push(ed);
			//	printf ("6\n");

			it++;
		}	
		printf("Done!\n");
		
		printf("Assigning cells to basins...\n");
		IndexIterator cit = bcc->getCellIterator();
		while (cit.isValid()) {
			index_type cel = *cit.loc;
			rec_assign_down(cel);
			cit++;
		}
		printf("Done!\n");
	
		printf("Creating hierarchy...\n");
		while (cancel_lowest());
		printf("Done!\n");

	}

	bool cancel_lowest() {
		if (be_sorter.empty()) return false;

		basin_edge<FType> be = be_sorter.top();
		if (be.b1 == be.b2) {
			be_sorter.pop();
			return true;
		}
		if (be.b1->parent != NULL) {
			basin_edge<FType> ne;
			ne.value = be.value;
			ne.b1 = be.b1->parent;
			ne.persistence = 
				ne.value - max(be.b1->parent->value, be.b2->value);
			ne.b2 = be.b2;
			be_sorter.pop();
			be_sorter.push(ne);
			return true;
		}
		if (be.b2->parent != NULL) {
			basin_edge<FType> ne;
			ne.value = be.value;
			ne.b2 = be.b2->parent;
			ne.persistence = 
				ne.value - max(be.b2->parent->value, be.b1->value);
			ne.b1 = be.b1;
			be_sorter.pop();
			be_sorter.push(ne);
			return true;
		}
		// otherwise be is now the lowest and a valid cancellation
		// set the parent
		basin<FType>* bupper;
		basin<FType>* blower;
		if (basin_LT(*be.b1, *be.b2)) {
			bupper = be.b1;
			blower = be.b2;
		} else {
			bupper = be.b2;
			blower = be.b1;
		}
		g_int_to_f_pers.push_back(be.persistence);
		bupper->destroyed = g_dest_count++;
		bupper->parent = blower;
		be_sorter.pop();
		return true;
	}


   void print_persistence() {
/*      typename map<index_type, basin<FType>*>::iterator b_iter;
      b_iter = min_to_basin_map.begin();

      while (b_iter != min_to_basin_map.end()) {
         cout << b_iter->first << " " << b_iter->second->value << endl;
         b_iter++;
      }
*/
      IndexIterator  it = bcc->getCellIterator(0);
      FType min = bcc->getValue(*it.loc);
      FType max = bcc->getValue(*it.loc);
      while (it.isValid()) {
         index_type cellid = *it.loc;
         FType val = bcc->getValue(cellid);

         if (val < min) {
            min = val;
         }
         if (val > max) {
            max = val;
         }
         it++;
      }

      double max_f = max - min;
      double num_verts = bcc->getCellIterator(0).size;
      double num_edges = bcc->getCellIterator(1).size;
      int i;
      for (i=1; i<g_int_to_f_pers.size(); i++) {
         cerr << g_int_to_f_pers[i-1]/max_f << " ";
         cerr << g_int_to_f_pers[i-1] << " ";
         cerr << float(min_to_basin_map.size()-i)/num_edges << " ";
         cerr << float(min_to_basin_map.size()-i)/num_verts << " ";
         cerr << min_to_basin_map.size()-i << endl;

         cerr << g_int_to_f_pers[i]/max_f << " ";
         cerr << g_int_to_f_pers[i] << " ";
         cerr << float(min_to_basin_map.size()-i)/num_edges << " ";
         cerr << float(min_to_basin_map.size()-i)/num_verts << " ";
         cerr << min_to_basin_map.size()-i << endl;
      }
      cerr << g_int_to_f_pers[i-1]/max_f << " ";
      cerr << g_int_to_f_pers[i-1] << " 0 0" << endl;
   }

};
