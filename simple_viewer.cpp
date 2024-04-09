


#include <iostream>
#include "unstructured_simplicial_complex.h"
#include "index_iterator.h"
#include "computeGrad.h"
#include "basins.h"
//#include "dump_image.h"
#include "mergeTree.h"


#ifdef __APPLE__
    #include "TargetConditionals.h"
    #ifdef TARGET_OS_MAC
        #include <GLUT/glut.h>
        #include <OpenGL/OpenGL.h>
    #endif
#elif defined _WIN32 || defined _WIN64
    #include <GL\glut.h>
#endif 

#include <stdio.h>
#include "Matrix4.h"
//#include "CImg.h"

#include <cmath>
#include <algorithm>

#define EPSILON 0.0 //0.0001

UnstructuredSimplicialComplex<2, float>* usc;
basins<2, float>* bsn;
bool drawbasincolor = true;
int glivingcounter = 0;
MergeTree<index_type, float, 2>* mt;


float red_scale(float s) {
	Vector4 v;
	v.vals[0] = s; //max(s, 1.0f-s);
	v.vals[1] = 1.0f - s;
	v.vals[2] = max(s, 1.0f - s);
	Normalize3(&v);
	return v.vals[0];
	//return min(max(4.0*(0.5-s), 0.0), 1.0);
}

float green_scale(float s) {
	Vector4 v;
	v.vals[0] = s;//max(s, 1.0f-s);
	v.vals[1] = 1.0f - s;
	v.vals[2] = max(s, 1.0f - s);
	Normalize3(&v);
	return v.vals[1];
	//return 0.2f;
	//return min(max(4.0*fabs(s-0.25)-1.0, 0.0), 1.0);
}

float blue_scale(float s) {
	Vector4 v;
	v.vals[0] = s;//max(s, 1.0f-s);
	v.vals[1] = 1.0f - s;
	v.vals[2] = max(s, 1.0f - s);
	Normalize3(&v);
	return v.vals[2];
	// return  1.0-s;
   //return min(max(4.0*(0.75-s), 0.0), 1.0);
}


using namespace std;

vector<index_type> total_order;
vector<int> total_counts;
vector<index_type> total_pqlist;

void Initialize_GLUT(int argc, char** argv);


void draw_earth(float x, float y, float z, float size)
{
	glEnable(GL_LIGHTING);
	glPushMatrix();
	glTranslatef(x, y, z);
	//glScalef(size, size, size);
	GLUquadricObj* q = gluNewQuadric();
	gluQuadricDrawStyle(q, GLU_FILL);
	gluQuadricNormals(q, GLU_SMOOTH);
	gluSphere(q, size, 10, 10);
	gluDeleteQuadric(q);
	glPopMatrix();
	glDisable(GL_LIGHTING);
}



int XMIN = 2;
int YMIN = 2;
int ZMIN = 2;


int XMAX = 2;
int YMAX = 2;
int ZMAX = 2;

float ballsize = 0.1f;
float percp = 0.0f;

template<class FType> struct gminmax {
	FType minval;
	FType maxval;
	float xmin;
	float xmax;
	float ymin;
	float ymax;

};

index_type cutoff = 0;

bool gusecutoffhack = false;
bool draw_flat = false;
bool draw_edges = true;
gminmax<float> fminmax;

template <unsigned char Dim, class FType>
void setGlobalMinMax(UnstructuredSimplicialComplex<Dim, FType>* bcc) {
	IndexIterator  it = bcc->getCellIterator(0);
	fminmax.minval = bcc->getValue(*it.loc);
	fminmax.maxval = bcc->getValue(*it.loc);
	XMIN = XMAX = bcc->verts[*it.loc].position[0];
	ZMIN = ZMAX = bcc->verts[*it.loc].position[1];


	while (it.isValid()) {
		index_type cellid = *it.loc;
		FType val = bcc->getValue(cellid);
		if (val > fminmax.maxval) fminmax.maxval = val;
		if (val < fminmax.minval) fminmax.minval = val;

		if (XMIN > bcc->verts[*it.loc].position[0])
			XMIN = bcc->verts[*it.loc].position[0];
		if (XMAX < bcc->verts[*it.loc].position[0])
			XMAX = bcc->verts[*it.loc].position[0];

		if (ZMIN > bcc->verts[*it.loc].position[1])
			ZMIN = bcc->verts[*it.loc].position[1];
		if (ZMAX < bcc->verts[*it.loc].position[1])
			ZMAX = bcc->verts[*it.loc].position[1];

		it++;
	}

	YMIN = fminmax.minval;
	YMAX = fminmax.maxval;

};

bool draw_gradient = false;
float arrow_width = 0.1;
float line_width = 2.0;

template <unsigned char Dim, class FType>
void centroid(UnstructuredSimplicialComplex<Dim, FType>* bcc, index_type cellid, float* verts) {

	for (int i = 0; i < 3; i++) verts[i] = 0.0f;
	Simplex<FType>& s = bcc->cells[cellid];
	for (int i = 0; i < s.numberOfVerts; i++) {
		BaseVertex<Dim>& v = bcc->verts[s.verts[i]];
		verts[0] += v.position[0];
		if (!draw_flat)
			verts[1] += bcc->getValue(s.verts[i]);
		else
			verts[1] += bcc->getValue(cellid);
		verts[2] += v.position[1];
	}
	for (int i = 0; i < 3; i++) verts[i] /= s.numberOfVerts;

};

template <unsigned char Dim, class FType>
void drawArrow(UnstructuredSimplicialComplex<Dim, FType>* bcc, index_type cellid) {
	if (bcc->getCritical(cellid)) return;
	if (!bcc->getAssigned(cellid)) return;
	index_type pair = bcc->getPair(cellid);
	if (bcc->getDim(pair) > bcc->getDim(cellid)) return;
	float start[3];
	float end[3];
	centroid<Dim, FType>(bcc, cellid, end);
	centroid<Dim, FType>(bcc, pair, start);

	//instead of a line
	//  glVertex3f(start[0], start[1]+ 2*EPSILON, start[2]);
	//  glVertex3f(end[0], end[1]+ 2*EPSILON, end[2]);


	//let's use code for computing the worlds most expensive arrow.
	double diff[3];
	diff[0] = end[0] - start[0];
	diff[1] = end[1] - start[1];
	diff[2] = end[2] - start[2];
	double len = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
	len = sqrt(len);
	diff[0] /= len;
	diff[1] /= len;
	diff[2] /= len;

	const float RADDEG = 57.29578;

	float Q = atan2(diff[1], diff[0]) * RADDEG;
	float P = acos(diff[2]) * RADDEG;

	glColor3f(0.1, 0.1, 0.1);

	glPushMatrix();
	GLUquadricObj* quadObj = gluNewQuadric();
	GLUquadricObj* quadObj2 = gluNewQuadric();
	GLUquadricObj* quadObj3 = gluNewQuadric();
	GLUquadricObj* quadObj4 = gluNewQuadric();

	gluQuadricDrawStyle(quadObj, GLU_FILL);
	gluQuadricNormals(quadObj, GLU_SMOOTH);
	gluQuadricDrawStyle(quadObj2, GLU_FILL);
	gluQuadricNormals(quadObj2, GLU_SMOOTH);
	gluQuadricDrawStyle(quadObj3, GLU_FILL);
	gluQuadricNormals(quadObj3, GLU_SMOOTH);
	gluQuadricDrawStyle(quadObj4, GLU_FILL);
	gluQuadricNormals(quadObj4, GLU_SMOOTH);


	glTranslatef(start[0], start[1], start[2]);
	glRotatef(Q, 0, 0, 1);
	glRotatef(P, 0, 1, 0);
	gluCylinder(quadObj, 0.4 * arrow_width, 0.4 * arrow_width, len - 2 * arrow_width, 10, 10);
	gluDisk(quadObj2, 0, 0.4 * arrow_width, 10, 10);

	glTranslatef(0, 0, len - 2.5 * arrow_width);
	gluCylinder(quadObj3, arrow_width, 0, 2.5 * arrow_width, 10, 10);
	gluDisk(quadObj4, 0, arrow_width, 10, 10);

	gluDeleteQuadric(quadObj);
	gluDeleteQuadric(quadObj2);
	gluDeleteQuadric(quadObj3);
	gluDeleteQuadric(quadObj4);
	glPopMatrix();



};

template <unsigned char Dim, class FType>
void drawTube(UnstructuredSimplicialComplex<Dim, FType>* bcc, index_type cellid) {
	
	IndexIterator fit = bcc->getFacetIterator(cellid);
	//fit.begin();
	std::vector<index_type> verts;
	while (fit.isValid()) {
		verts.push_back(*fit.loc);
		fit++;
	}
	if (verts.size() != 2) {
		printf("WHOA verts size is %d\n", verts.size());
		return;
	}

	float start[3];
	float end[3];
	centroid<Dim, FType>(bcc, verts[0], end);
	centroid<Dim, FType>(bcc, verts[1], start);
	
	// need flat func.
	
	//instead of a line
	//  glVertex3f(start[0], start[1]+ 2*EPSILON, start[2]);
	//  glVertex3f(end[0], end[1]+ 2*EPSILON, end[2]);


	//let's use code for computing the worlds most expensive arrow.
	double diff[3];
	diff[0] = end[0] - start[0];
	diff[1] = end[1] - start[1];
	diff[2] = end[2] - start[2];
	double len = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
	len = sqrt(len);
	diff[0] /= len;
	diff[1] /= len;
	diff[2] /= len;

	const float RADDEG = 57.29578;

	float Q = atan2(diff[1], diff[0]) * RADDEG;
	float P = acos(diff[2]) * RADDEG;

	glColor3f(1, 0.1, 0.1);

	glPushMatrix();
	GLUquadricObj* quadObj = gluNewQuadric();
	GLUquadricObj* quadObj2 = gluNewQuadric();
	GLUquadricObj* quadObj4 = gluNewQuadric();

	gluQuadricDrawStyle(quadObj, GLU_FILL);
	gluQuadricNormals(quadObj, GLU_SMOOTH);
	gluQuadricDrawStyle(quadObj2, GLU_FILL);
	gluQuadricNormals(quadObj2, GLU_SMOOTH);
	gluQuadricDrawStyle(quadObj4, GLU_FILL);
	gluQuadricNormals(quadObj4, GLU_SMOOTH);


	glTranslatef(start[0], start[1], start[2]);
	glRotatef(Q, 0, 0, 1);
	glRotatef(P, 0, 1, 0);
	gluCylinder(quadObj, 0.7 * arrow_width, 0.7 * arrow_width, len, 10, 10);
	gluDisk(quadObj2, 0, 0.7 * arrow_width, 10, 10);

	glTranslatef(0, 0, len);
	gluDisk(quadObj4, 0, .7* arrow_width, 10, 10);

	gluDeleteQuadric(quadObj);
	gluDeleteQuadric(quadObj2);
	gluDeleteQuadric(quadObj4);
	glPopMatrix();



};


bool draw_mt = false;
bool flat_funct = false;
bool redrawstuff = true;
bool draw_unassigned = true;
GLuint drawlist = -1;

template <unsigned char Dim, class FType>
void drawStuff(UnstructuredSimplicialComplex<Dim, FType>* bcc) {




	if (!redrawstuff) {
		glCallList(drawlist);
		return;
	}

	if (drawlist != -1)
		glDeleteLists(drawlist, 1);

	drawlist = glGenLists(1);
	glNewList(drawlist, GL_COMPILE_AND_EXECUTE);

	redrawstuff = false;

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.0, 1.0);


	// TRIANGLES	

	   //  glColor3f(.2, .2, .2); 
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
	if (!draw_flat)  glEnable(GL_LIGHTING);
	glBegin(GL_TRIANGLES);
	IndexIterator  it = bcc->getCellIterator(Dim);
	while (it.isValid()) {
		index_type cellid = *it.loc;
		//total hack!
		Simplex<FType>& s = bcc->cells[cellid];
		if (s.numberOfVerts != 3) {
			printf("SDJKFHDS:LJFSKL:D\n");
		}

		BaseVertex<Dim>& v1 = bcc->verts[s.verts[0]];
		BaseVertex<Dim>& v2 = bcc->verts[s.verts[1]];
		BaseVertex<Dim>& v3 = bcc->verts[s.verts[2]];

		bool is_assigned = bcc->getAssigned(cellid);

		if (!draw_unassigned && !is_assigned) {
			it++; continue;
		}
		//if (bcc->getNumUFacets(cellid) == 1) {
	  //	  glColor3f(0,1,0);
		//}
		if (is_assigned) {


			unsigned char asif = bcc->getDimAscMan(cellid);
			if (asif == 2) {
				glColor3f(0.8, 0.8, 0.0);
			}
			else if (asif == 1) {
				glColor3f(0.3, 0.3, 0.9);
			}
			else {
				glColor3f(1.0, 0, 0);
			}
		}

		if (drawbasincolor) {
			float* c = bsn->livingCellColor(cellid, glivingcounter);
			glColor3f(c[0], c[1], c[2]);
		}

		if (!draw_flat) {
			// normal for lighting
			Vector4 vv1;
			vv1.vals[0] = v1.position[0];
			vv1.vals[1] = bcc->getValue(s.verts[0]);
			vv1.vals[2] = v1.position[1];
			vv1.vals[3] = 1.0f;
			Vector4 vv2;
			vv2.vals[0] = v2.position[0];
			vv2.vals[1] = bcc->getValue(s.verts[1]);
			vv2.vals[2] = v2.position[1];
			vv2.vals[3] = 1.0f;
			Vector4 vv3;
			vv3.vals[0] = v3.position[0];
			vv3.vals[1] = bcc->getValue(s.verts[2]);
			vv3.vals[2] = v3.position[1];
			vv3.vals[3] = 1.0f;
			for (int i = 0; i < 4; i++) {
				vv1.vals[i] = vv3.vals[i] - vv1.vals[i];
				vv2.vals[i] = vv3.vals[i] - vv2.vals[i];
			}
			Normalize4(&vv1);
			Normalize4(&vv2);
			Vector4 n = Cross(vv1, vv2);
			Normalize4(&n);
			glNormal3f(n.vals[0], n.vals[1], n.vals[2]);
		}

		unsigned char asif = bcc->getDimAscMan(cellid);
		if (!drawbasincolor && !is_assigned /*&& bcc->getNumUFacets(cellid) == 0*/ && flat_funct) {
			float nval = ((float)bcc->getValue(cellid) - (float)fminmax.minval)
				/ ((float)fminmax.maxval - (float)fminmax.minval);
			float r = red_scale(nval);
			float g = green_scale(nval);
			float b = blue_scale(nval);
			glColor3f(r, g, b);


		}
		if (!drawbasincolor && !is_assigned  /*&& bcc->getNumUFacets(cellid) == 0*/ && !flat_funct) {
			float nval = ((float)bcc->getValue(s.verts[0]) - (float)fminmax.minval)
				/ ((float)fminmax.maxval - (float)fminmax.minval);
			float r = red_scale(nval);
			float g = green_scale(nval);
			float b = blue_scale(nval);
			glColor3f(r, g, b);
		}

		if (!draw_flat)
			glVertex3f(v1.position[0], bcc->getValue(s.verts[0]), v1.position[1]);
		else
			glVertex3f(v1.position[0], bcc->getValue(cellid), v1.position[1]);

		if (!drawbasincolor && !is_assigned  /*&& bcc->getNumUFacets(cellid) == 0 */ && !flat_funct) {
			float nval = ((float)bcc->getValue(s.verts[1]) - (float)fminmax.minval)
				/ ((float)fminmax.maxval - (float)fminmax.minval);
			float r = red_scale(nval);
			float g = green_scale(nval);
			float b = blue_scale(nval);
			glColor3f(r, g, b);
		}
		if (!draw_flat)
			glVertex3f(v2.position[0], bcc->getValue(s.verts[1]), v2.position[1]);
		else
			glVertex3f(v2.position[0], bcc->getValue(cellid), v2.position[1]);

		if (!drawbasincolor && !is_assigned /*&& bcc->getNumUFacets(cellid) == 0*/ && !flat_funct) {
			float nval = ((float)bcc->getValue(s.verts[2]) - (float)fminmax.minval)
				/ ((float)fminmax.maxval - (float)fminmax.minval);

			float r = red_scale(nval);
			float g = green_scale(nval);
			float b = blue_scale(nval);
			glColor3f(r, g, b);
		}

		if (!draw_flat)
			glVertex3f(v3.position[0], bcc->getValue(s.verts[2]), v3.position[1]);
		else
			glVertex3f(v3.position[0], bcc->getValue(cellid), v3.position[1]);



		it++;
	}
	glEnd();
	if (!draw_flat) glDisable(GL_LIGHTING);





	it = bcc->getCellIterator(0);
	glPointSize(line_width * 2 + 2);
	glBegin(GL_POINTS);
	while (it.isValid()) {
		index_type cellid = *it.loc;
		//total hack!
		Simplex<FType>& s = bcc->cells[cellid];

		bool is_assigned = bcc->getAssigned(cellid);
		if (!draw_unassigned && !is_assigned) {
			it++; continue;
		}

		BaseVertex<Dim>& v1 = bcc->verts[s.verts[0]];

		if (s.verts[0] != cellid) printf("ERORORORORORORO\n");


		if (gusecutoffhack && bsn->livingCellValue(cellid, glivingcounter) < 0.01f) {
			it++;
			continue;
		}

		if (drawbasincolor && !bcc->getCritical(cellid)) {
			float* c = bsn->livingCellColor(cellid, glivingcounter);
			glColor3f(c[0], c[1], c[2]);
			glVertex3f(v1.position[0], bcc->getValue(cellid) + EPSILON, v1.position[1]);

		}
		else if (bcc->getAssigned(cellid)) {

			if (!bcc->getCritical(cellid)) {
				glColor3f(.6, .6, 0);
				glVertex3f(v1.position[0], bcc->getValue(cellid) + EPSILON, v1.position[1]);
			}
			//printf("reder min%d = %f %f %f\n", cellid,v1.position[0],bcc->getValue(s.verts[0]),v1.position[1] );

		}
		else {
			float nval = ((float)bcc->getValue(s.verts[0]) - (float)fminmax.minval)
				/ ((float)fminmax.maxval - (float)fminmax.minval);
			float r = red_scale(nval) - .1;
			float g = green_scale(nval) - .1;
			float b = blue_scale(nval) - .1;
			glColor3f(r, g, b);
			glVertex3f(v1.position[0], bcc->getValue(cellid) + EPSILON, v1.position[1]);//draw_earth(v1.position[0],bcc->getValue(s.verts[0]),v1.position[1], ballsize/2.0);
		}

		it++;
	}
	glEnd();




	// LINES
	glDisable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(0.0, 0.0);

	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_BLEND);
	glLineWidth(line_width);



	if (draw_mt) {
		vector<index_type> vts;
		mt->gl_fill_vertices2(vts);

		for (int i = 0; i < vts.size(); i += 2) {
			BaseVertex<Dim>& v1 = bcc->verts[vts[i]];
			index_type b = vts[i + 1];
			glColor3f(
				(((float)((b * 321 + 93) % 31))) / 31.0f,
				(((float)((b * 421 + 91) % 51))) / 51.0f,
				(((float)((b * 221 + 92) % 71))) / 71.0f);
			draw_earth(v1.position[0], bcc->getValue(vts[i]), v1.position[1], ballsize * 3);
		}
		vts.clear();
		mt->gl_fill_arcs(vts);
		glBegin(GL_LINES);
		for (int i = 0; i < vts.size(); i++) {
			BaseVertex<Dim>& v1 = bcc->verts[vts[i]];
			glVertex3f(v1.position[0], bcc->getValue(vts[i]) + EPSILON, v1.position[1]);
		}
		glEnd();
		vts.clear();
	}


	std::vector<index_type> saddle_ids;
	if (draw_edges) {

		glBegin(GL_LINES);
		glColor4f(0.6, 0.5, 0.5, 1.0);
		it = bcc->getCellIterator(1);
		while (it.isValid()) {
			index_type cellid = *it.loc;
			//total hack!
			Simplex<FType>& s = bcc->cells[cellid];
			if (s.numberOfVerts != 2) {
				printf("SDJKFHDS:LJFSKL:D\n");
			}
			bool is_assigned = bcc->getAssigned(cellid);
			if (!draw_unassigned && !is_assigned) {
				it++; continue;
			}

			if (gusecutoffhack && bsn->livingCellValue(cellid, glivingcounter) < 0.01f) {
				it++;
				continue;
			}

			BaseVertex<Dim>& v1 = bcc->verts[s.verts[0]];
			BaseVertex<Dim>& v2 = bcc->verts[s.verts[1]];
			if (drawbasincolor) {
				float* c = bsn->livingCellColor(cellid, glivingcounter);
				glColor3f(c[0], c[1], c[2]);

				glVertex3f(v1.position[0], bcc->getValue(s.verts[0]) + EPSILON, v1.position[1]);
				glVertex3f(v2.position[0], bcc->getValue(s.verts[1]) + EPSILON, v2.position[1]);

			}
			else if (!draw_flat) {
				if (bcc->getNumUFacets(cellid) == 1) {

					//glColor4f(0,1, 0, 1.0);


					if (bcc->getAssigned(cellid)) {


						if (bcc->getCritical(cellid)) {
							glColor4f(1, 0, 0, 1.0);
							saddle_ids.push_back(cellid);
						}
						else {


							unsigned char asif = bcc->getDimAscMan(cellid);
							if (asif == 2) {
								glColor4f(.6, .6, 0, 1.0);
							}
							else if (asif == 1) {
								glColor4f(0.3, 0.3, 0.9, 1.0);
							}
							else {
								glColor4f(1.0, 0, 0, 1.0);
							}
						}
					}


					// float c[3];
					// centroid<Dim,FType>(bcc, cellid, c);
					// draw_earth(c[0],c[1],c[2],ballsize);
					glVertex3f(v1.position[0], bcc->getValue(s.verts[0]) + EPSILON, v1.position[1]);


					glVertex3f(v2.position[0], bcc->getValue(s.verts[1]) + EPSILON, v2.position[1]);

				}
				else {
					if (!flat_funct) {
						float nval = ((float)bcc->getValue(s.verts[0]) - (float)fminmax.minval)
							/ ((float)fminmax.maxval - (float)fminmax.minval);
						float r = red_scale(nval);
						float g = green_scale(nval);
						float b = blue_scale(nval);
						float nval2 = ((float)bcc->getValue(s.verts[1]) - (float)fminmax.minval)
							/ ((float)fminmax.maxval - (float)fminmax.minval);
						float r2 = red_scale(nval2);
						float g2 = green_scale(nval2);
						float b2 = blue_scale(nval2);

						float a = 0.3f;
						glColor4f(r + a, g + a, b + a, 0.2);
						glVertex3f(v1.position[0], bcc->getValue(s.verts[0]) + EPSILON, v1.position[1]);

						glColor4f(r2 + a, g2 + a, b2 + a, 0.2);
						glVertex3f(v2.position[0], bcc->getValue(s.verts[1]) + EPSILON, v2.position[1]);
					}
					else {
						float nval = ((float)bcc->getValue(cellid) - (float)fminmax.minval)
							/ ((float)fminmax.maxval - (float)fminmax.minval);
						float r = red_scale(nval);
						float g = green_scale(nval);
						float b = blue_scale(nval);
						float a = 0.3f;
						glColor4f(r + a, g + a, b + a, 0.2);
						glVertex3f(v1.position[0], bcc->getValue(s.verts[0]) + EPSILON, v1.position[1]);
						glVertex3f(v2.position[0], bcc->getValue(s.verts[1]) + EPSILON, v2.position[1]);
					}
				}


			}
			else {
				float nval = ((float)bcc->getValue(cellid) - (float)fminmax.minval)
					/ ((float)fminmax.maxval - (float)fminmax.minval);
				float r = red_scale(nval);
				float g = green_scale(nval);
				float b = blue_scale(nval);

				float a = 0.3f;
				glColor4f(r + a, g + a, b + a, 0.4);
				glVertex3f(v1.position[0], bcc->getValue(cellid) + EPSILON, v1.position[1]);
				glVertex3f(v2.position[0], bcc->getValue(cellid) + EPSILON, v2.position[1]);
			}


			it++;
		}
		glEnd();

	}
	glDisable(GL_BLEND);


	// VERTICES
	it = bcc->getCellIterator(0);
	index_type sizeshit = it.size;
	while (it.isValid()) {
		index_type cellid = *it.loc;
		//total hack!
		Simplex<FType>& s = bcc->cells[cellid];

		bool is_assigned = bcc->getAssigned(cellid);
		if (!draw_unassigned && !is_assigned) {
			it++; continue;
		}

		BaseVertex<Dim>& v1 = bcc->verts[s.verts[0]];


		if (drawbasincolor && bcc->getCritical(cellid)) {
			float* c = bsn->livingCellColor(cellid, glivingcounter);
			glColor3f(c[0], c[1], c[2]);
			draw_earth(v1.position[0], bcc->getValue(s.verts[0]), v1.position[1], ballsize);

		}
		else if (bcc->getAssigned(cellid)) {

			if (bcc->getCritical(cellid)) {
				glColor3f(1, 0, 0);
				draw_earth(v1.position[0], bcc->getValue(s.verts[0]), v1.position[1], ballsize);
			}
			else {
				glColor3f(.6, .6, 0);
				//draw_earth(v1.position[0],bcc->getValue(s.verts[0]),v1.position[1], ballsize/2.0);
			}
			//printf("reder min%d = %f %f %f\n", cellid,v1.position[0],bcc->getValue(s.verts[0]),v1.position[1] );

		}
		else {
			float nval = ((float)bcc->getValue(s.verts[0]) - (float)fminmax.minval)
				/ ((float)fminmax.maxval - (float)fminmax.minval);
			float r = red_scale(nval);
			float g = green_scale(nval);
			float b = blue_scale(nval);
			glColor3f(r, g, b);
			//draw_earth(v1.position[0],bcc->getValue(s.verts[0]),v1.position[1], ballsize/2.0);
		}

		it++;
	}

	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
	glEnable(GL_LIGHTING);


	it = bcc->getCellIterator(1);
	while (it.isValid()) {
		index_type cellid = *it.loc;
		if (bcc->getCritical(cellid)) {
			drawTube(bcc, cellid);
		}//drawArrow<Dim, FType>(bcc, cellid);
		it++;
	}

	glDisable(GL_LIGHTING);
	//      glEnd();

	if (draw_gradient) {
		//      glLineWidth(arrow_width);
		//      glColor3f(0, 0, 0);
		//      glBegin(GL_LINES);

		glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
		glEnable(GL_LIGHTING);

		for (int dd = 1; dd <= 2; dd++) {
			it = bcc->getCellIterator(dd);
			while (it.isValid()) {
				index_type cellid = *it.loc;
				drawArrow<Dim, FType>(bcc, cellid);
				it++;
			}
		}
		glDisable(GL_LIGHTING);
		//      glEnd();
	}
	glEndList();

};

int main(int argc, char** argv) {


	// const unsigned char purple[] = {255, 0, 255};
	// cimg_library::CImg<unsigned char>(640, 400, 1, 3, 0).draw_text(100, 100, "Hello World", purple).display("my first casasfd");


	/*
	   int array[10];

	   for (int i=0; i<10; i++) {
	   array[i] = 2*i+5;
	   }

	   IndexIterator ii(array, 5);

	   while (ii.isValid()) {
	   cout << *ii.loc << endl;
	   ii++;
	   }*/
	printf("Welcome to the unstructures MS Complex computer!\n");

	usc = new UnstructuredSimplicialComplex<2, float>();
	if (argc > 2) {
		//treat the second parameter as a scale to function value
		usc->loadFromOff(argv[1], atof(argv[2]));
	}
	else {
		usc->loadFromOff(argv[1]);
	}

	printf("number of verts read %d\n", usc->numberOfVerts);

	printf("testing uf\n");

	UnionFind<index_type, int> myuf;
	myuf.MakeSet(0, 10);
	myuf.MakeSet(1, 2);
	myuf.MakeSet(2, 15);
	myuf.MakeSet(3, 21);
	myuf.Union(0, 2);
	myuf.Union(1, 3);
	myuf.Union(2, 3);
	myuf.print();

	printf("testing mt\n");

	mt = new MergeTree<index_type, float, 2>(usc);
	mt->computeMT();
	mt->decompose();
	mt->dumpDot("mt.dot");

	printf("done\n");


	computeGradient<2, float>(usc, total_order, total_counts, total_pqlist);

	bsn = new basins<2, float>(usc);
	bsn->computeBasins();
	glivingcounter = bsn->num_destroyed();

	for (int i = 0; i < total_order.size(); i++) {
		usc->setAssigned(total_order[i], false);
	}

	int numcps[3];
	numcps[0] = 0;  numcps[1] = 0; numcps[2] = 0;
	for (int kk = 0; kk <= 2; kk++) {
		IndexIterator it = usc->getCellIterator(kk);
		while (it.isValid()) {
			index_type itv = *it.loc;
			if (usc->getCritical(itv)) numcps[usc->getDim(itv)]++;
			//if (! usc->getAssigned(itv)) printf("ERROR NOT ASSIGNED\n");
			it++;
		}
	}
	printf("mins=%d sad=%d maxs=%d\n", numcps[0], numcps[1], numcps[2]);

	/*
	   UnstructuredSimplicialComplex<3, float> usc;
	   usc.loadFromOff("../test_data/test3d.off");
	 */

	int win;

	Initialize_GLUT(argc, argv);

	return 0;
}


/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//////////                                       ////////////////
//////////    GLUT STUFF DOWN HERE               ////////////////
//////////                                       ////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

Matrix4 mat;
float translateoff[3];
bool clearcolor = 1;
bool g_mouseActive;
int mainWindow;

// Initial size of graphics window.
const int WIDTH = 640;
const int HEIGHT = 480;

// Current size of window.
int width = WIDTH;
int height = HEIGHT;

// Mouse positions, normalized to [0,1].
double xMouse = 0.0;
double yMouse = 0.0;

// Bounds of viewing frustum.
double nearPlane = 0.1;
double farPlane = 30000;

// Viewing angle.
double fovy = 40.0;

// Variables.
double alpha = 0;                                  // Set by idle function.
double beta = 0;                                   // Set by mouse X.
double mdistance = -(farPlane - nearPlane) / 2;    // Set by mouse Y.
float dist;



GLfloat light_diffuse[] = { 0.8, 0.8, 0.8, 1.0 };  /* Red diffuse light. */
GLfloat light_position[] = { 1.0, 1.0, 1.0, 0.0 };  /* Infinite light location. */

float diffuseLight[] = { 0.8f, 0.8f, 0.8f, 1.0f };
float specularLight[] = { 1.0f, 1.0f, 1.0f, 1.0f };
float LightPosition[] = { 1.1f, 0.0f, 8.0f, 1.0f };

struct Quat {
	float x, y, z, w;
};

Quat times(Quat a, Quat b) {
	Quat res;
	res.w = a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z;
	res.x = a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y;
	res.y = a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x;
	res.z = a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w;
	return res;
}

Quat rotation(Vector4 v, float angle) {
	Quat res;
	Normalize3(&v);
	res.w = cosf(angle / 2.0);
	float tmp = sinf(angle / 2.0);
	res.x = v.vals[0] * tmp;
	res.y = v.vals[1] * tmp;
	res.z = v.vals[2] * tmp;
	return res;
}

Matrix4 rotmatrix;
Quat rotationtotal;
void resetQuat(Quat& q) {
	q.x = 0;
	q.y = 0;
	q.z = 0;
	q.w = 1;
}

void setRotQuat(Matrix4& m, Quat& q) {
	{
		float x2 = 2.0f * q.x, y2 = 2.0f * q.y, z2 = 2.0f * q.z;

		float xy = x2 * q.y, xz = x2 * q.z;
		float yy = y2 * q.y, yw = y2 * q.w;
		float zw = z2 * q.w, zz = z2 * q.z;

		m.vals[0] = 1.0f - (yy + zz);
		m.vals[1] = (xy - zw);
		m.vals[2] = (xz + yw);
		m.vals[3] = 0.0f;

		float xx = x2 * q.x, xw = x2 * q.w, yz = y2 * q.z;

		m.vals[4] = (xy + zw);
		m.vals[5] = 1.0f - (xx + zz);
		m.vals[6] = (yz - xw);
		m.vals[7] = 0.0f;

		m.vals[8] = (xz - yw);
		m.vals[9] = (yz + xw);
		m.vals[10] = 1.0f - (xx + yy);
		m.vals[11] = 0.0f;

		m.vals[12] = 0.0f;
		m.vals[13] = 0.0f;
		m.vals[14] = 0.0f;
		m.vals[15] = 1.0f;
	}

}

void reset_view() {
	resetQuat(rotationtotal);
	mat = MIdentity();
	dist = 300 * (ZMAX - ZMIN);
}

void render_view() {
	resetQuat(rotationtotal);
	mat = MIdentity();
	dist = 240 * (ZMAX - ZMIN);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	//glMultMatrixf(mat.vals);
	//glRotatef(-50, 0,1,0);
	Vector4 v; v.vals[0] = 0; v.vals[1] = 1; v.vals[2] = 0;
	rotationtotal = rotation(v, -50 * PI / 180.0);
	setRotQuat(rotmatrix, rotationtotal);
	glMultMatrixf(rotmatrix.vals);

	//glMultMatrixf(mat.vals);

	//glGetFloatv(GL_MODELVIEW_MATRIX, mat.vals);
}

void myidle() {
	//cimg_library::cimg::sleep(10);

}

bool dbbox = true;

void display()
{
	glPointSize(4.0);
	if (clearcolor) {
		glClearColor(1, 1, 1, 1);
	}
	else {
		glClearColor(0, 0, 0, 1);
	}
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// enable crystal ball transform here

	glTranslatef(0, 0, -dist / 100);
	glTranslatef(-translateoff[0], -translateoff[1], -translateoff[2]);

	setRotQuat(rotmatrix, rotationtotal);
	glMultMatrixf(rotmatrix.vals);


	glTranslatef(-(XMAX - XMIN) / 2.0 - XMIN,
		-(YMAX - YMIN) / 2.0 - YMIN,
		-(ZMAX - ZMIN) / 2.0 - ZMIN);

	glEnable(GL_DEPTH_TEST);


	// draw the complex

	// if (dbbox) {
	// 	   if (clearcolor)
	// 		glColor3f(.05, .05, .05);
	// 	   else 
	// 		   glColor3f(.9,.9,.9);
	// 	 glBegin(GL_LINES);
	// 	 glVertex3f(0, 0, 0);
	// 	 glVertex3f(XDIM-1 , 0, 0);
	// 	 glVertex3f(0, YDIM-1, 0);
	// 	 glVertex3f(XDIM-1 , YDIM-1, 0);
	// 	 glVertex3f(0, 0, ZDIM-1);
	// 	 glVertex3f(XDIM-1 , 0, ZDIM-1);
	// 	 glVertex3f(0, YDIM-1, ZDIM-1);
	// 	 glVertex3f(XDIM-1 , YDIM-1, ZDIM-1);

	// 	glVertex3f(0, 0, 0);
	// 	glVertex3f(0, YDIM-1, 0);
	// 	glVertex3f(XDIM-1, 0, 0);
	// 	glVertex3f(XDIM-1, YDIM-1, 0);
	// 	glVertex3f(0, 0, ZDIM-1);
	// 	glVertex3f(0, YDIM-1, ZDIM-1);
	// 	glVertex3f(XDIM-1, 0, ZDIM-1);
	// 	glVertex3f(XDIM-1, YDIM-1, ZDIM-1);

	// 	 glVertex3f(0, 0, 0);
	// 	glVertex3f(0, 0, ZDIM-1);
	// 		 glVertex3f(XDIM-1, 0, 0);
	// 	glVertex3f(XDIM-1, 0, ZDIM-1);
	// 		 glVertex3f(0, YDIM-1, 0);
	// 	glVertex3f(0, YDIM-1, ZDIM-1);
	// 		 glVertex3f(XDIM-1, YDIM-1, 0);
	// 	glVertex3f(XDIM-1, YDIM-1, ZDIM-1);

	// 	glEnd();
	// }

	drawStuff<2, float>(usc);
	glutSwapBuffers();


}



void* gl_copy_pixels = NULL;



void reshapeMainWindow(int newWidth, int newHeight)
{
	width = newWidth;
	height = newHeight;
	if (gl_copy_pixels != NULL) {
		//delete(gl_copy_pixels);
		gl_copy_pixels = NULL;
	}
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(fovy, GLfloat(width) / GLfloat(height), nearPlane, farPlane);
	glutPostRedisplay();
}

void outputimage(char* name) {
	//unsigned char* rpixels = new unsigned char[width * height * 3];
	//unsigned char* mpixels = new unsigned char[width * height * 3];
	//glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, rpixels);
	//for (int x = 0; x < width; x++) {
	//	for (int y = 0; y < height; y++) {
	//		int id1 = x + y * width;
	//		int id2 = x + (height - 1 - y) * width;
	//		mpixels[id1] = rpixels[id2 * 3];
	//		mpixels[id1 + width * height] = rpixels[id2 * 3 + 1];
	//		mpixels[id1 + 2 * width * height] = rpixels[id2 * 3 + 2];
	//	}
	//}
	//cimg_library::CImg< unsigned char > img(mpixels, width, height, 1, 3, false);

	//img.save_bmp(name);
	//delete(rpixels);
	//delete(mpixels);
	//return;

}

void makeSequence() {

	for (int i = 0; i < total_order.size(); i++) {
		usc->setAssigned(total_order[i], false);
	}
	for (int i = 0; i < total_pqlist.size(); i++) {
		usc->setNumUFacets(total_pqlist[i], 0);
	}
	int nc = total_counts[0];
	for (int i = 0; i < nc; i++) {
		usc->setNumUFacets(total_pqlist[i], 1);
	}

	char fname[1024];

	int counter = 0;
	for (int i = 0; i < 700/*total_order.size()*/; i++) {
		sprintf(fname, "out/full_%05d.png", counter++);
		printf("outputtting: %s\n", fname);


		if (usc->getCritical(total_order[i])) {
			for (int j = nc; j < nc + total_counts[i + 1]; j++) {
				usc->setNumUFacets(total_pqlist[j], 1);
			}
			nc += total_counts[i + 1];

			usc->setAssigned(total_order[i], true);
			usc->setNumUFacets(total_order[i], 0);
			glutPostRedisplay();
			display();
			glutSwapBuffers();
			// outputimage(fname);

		}
		else {
			for (int j = nc; j < nc + total_counts[i + 1]; j++) {
				usc->setNumUFacets(total_pqlist[j], 1);
			}
			nc += total_counts[i + 1];

			usc->setAssigned(total_order[i++], true);

			for (int j = nc - 1; j < nc + total_counts[i]; j++) {
				usc->setNumUFacets(total_pqlist[j], 1);
			}
			nc += total_counts[i];
			usc->setAssigned(total_order[i], true);

			display();
			glutSwapBuffers();
			//glutPostRedisplay();  
		   // outputimage(fname);
		}
	}
}


void graphicKeys(unsigned char key, int x, int y)
{

	switch (key)
	{

	case 'T':
		draw_unassigned = !draw_unassigned;
		redrawstuff = true;
		break;
	case 't':
		draw_mt = !draw_mt;
		redrawstuff = true;
		//bsn->dumpBasinCountPerPersistence("basins.txt");
		break;

	case 'y':
		char tmpname[1024];
		sprintf(tmpname, "bsizes_%.6u.txt", glivingcounter);
		bsn->dumpBasinSizes(tmpname, glivingcounter);
		break;

	case 'h':
		gusecutoffhack = !gusecutoffhack;
		redrawstuff = true;
		break;

	case '1':
		glivingcounter = bsn->getDestrCount(0.0f);
		redrawstuff = true;
		break;
	case '2':
		glivingcounter = bsn->getDestrCount(0.001f);
		redrawstuff = true;
		break;
	case '3':
		glivingcounter = bsn->getDestrCount(0.005f);
		redrawstuff = true;
		break;
	case '4':
		glivingcounter = bsn->getDestrCount(0.01f);
		redrawstuff = true;
		break;
	case '5':
		glivingcounter = bsn->getDestrCount(0.05f);
		redrawstuff = true;
		break;
	case '6':
		glivingcounter = bsn->getDestrCount(0.1f);
		redrawstuff = true;
		break;
	case '7':
		glivingcounter = bsn->getDestrCount(0.2f);
		redrawstuff = true;
		break;
	case '8':
		glivingcounter = bsn->getDestrCount(0.5f);
		redrawstuff = true;
		break;
	case '9':
		glivingcounter = bsn->getDestrCount(1.0f);
		redrawstuff = true;
		break;

	case 'c':
		drawbasincolor = !drawbasincolor;
		redrawstuff = true;
		break;
	case 'j':
		if (glivingcounter > 0)
			glivingcounter--;
		redrawstuff = true;
		break;
	case 'k':
		if (glivingcounter < bsn->num_destroyed())
			glivingcounter++;
		redrawstuff = true;
		break;

	case 'b':
		clearcolor = !clearcolor;
		printf("Background color is (%d, %d, %d)\n", clearcolor, clearcolor, clearcolor);
		break;

	case 'm':
	{
		for (int i = 0; i < 100; i++) {
			if (cutoff < total_order.size()) {
				if (usc->getCritical(total_order[cutoff])) {
					usc->setAssigned(total_order[cutoff], true);
					if (cutoff < total_order.size() - 1) cutoff++;
				}
				else {
					usc->setAssigned(total_order[cutoff], true);
					cutoff++;
					usc->setAssigned(total_order[cutoff], true);
					cutoff++;
				}
			}
		}
	}
	redrawstuff = true;
	break;
	case 'n':
	{
		for (int i = 0; i < 100; i++) {
			if (cutoff > 0) {
				if (usc->getCritical(total_order[cutoff])) {
					usc->setAssigned(total_order[cutoff], false);
					cutoff--;
				}
				else {
					usc->setAssigned(total_order[cutoff], false);
					cutoff--;
					usc->setAssigned(total_order[cutoff], false);
					cutoff--;
				}
			}
		}
	}
	redrawstuff = true;
	break;
	case 'M':
	{
		for (int i = 0; i < 1; i++) {
			if (cutoff < total_order.size()) {
				if (usc->getCritical(total_order[cutoff])) {
					usc->setAssigned(total_order[cutoff], true);
					if (cutoff < total_order.size() - 1) cutoff++;
				}
				else {
					usc->setAssigned(total_order[cutoff], true);
					cutoff++;
					usc->setAssigned(total_order[cutoff], true);
					cutoff++;
				}
			}
		}
	}
	redrawstuff = true;
	break;
	case 'N':
	{
		for (int i = 0; i < 1; i++) {
			if (cutoff > 0) {
				if (usc->getCritical(total_order[cutoff])) {
					usc->setAssigned(total_order[cutoff], false);
					cutoff--;
				}
				else {
					usc->setAssigned(total_order[cutoff], false);
					cutoff--;
					usc->setAssigned(total_order[cutoff], false);
					cutoff--;
				}
			}
		}
	}
	redrawstuff = true;
	break;


	case 's':
		makeSequence();
		break;



	case '[':
		line_width *= 1.2f;
		redrawstuff = true;
		break;
	case ']':
		line_width /= 1.2f;
		redrawstuff = true;
		break;
	case 'e':
		draw_edges = !draw_edges;
		redrawstuff = true;
		break;
	case 'F':
		flat_funct = !flat_funct;
		redrawstuff = true;
		break;
	case 'f':
		draw_flat = !draw_flat;
		redrawstuff = true;
		break;
	case 'g':
		draw_gradient = !draw_gradient;
		redrawstuff = true;
		break;

	case '(':
		arrow_width *= 1.2f;
		redrawstuff = true;
		break;

	case ')':
		arrow_width /= 1.2f;
		redrawstuff = true;
		break;

	case 'x':
		ballsize *= 1.2f;
		redrawstuff = true;
		break;

	case 'z':
		ballsize /= 1.2f;
		redrawstuff = true;
		break;
	case 'w':
	{
		char filename[1024];
		sprintf(filename, "out.bmp");
		outputimage(filename);

		break;
	}
	case 'd':
		render_view();
		break;
	case 'r':
		reset_view();
		break;
	case 'q':
		exit(1);
		break;
	}
	glutPostRedisplay();
}

int xold;
int yold;
Vector4 vstart;
int buttonstate;
void mouseClicked(int button, int state, int x, int y)
{
	xold = x;
	yold = y;

	//printf("Mouse Clicked! button = %d, %d  x = %d   y = %d\n", state, button, x, y);
	buttonstate = button;
	if (state == GLUT_DOWN) {
		g_mouseActive = true;
	}
	else {
		g_mouseActive = false;
	}

	glutPostRedisplay();
}


void mouseMovement(int mx, int my)
{
	//printf("%d, %d\n", mx, my);
	if (!g_mouseActive) return;
	if (buttonstate == 2) {
		int dd = XMAX - XMIN;
		dist += dd * (my - yold);
		yold = my;

	}

	if (buttonstate == 1) {

		float ldim = (width > height ? width : height);


		vstart.vals[0] = xold - .5 * width;
		vstart.vals[1] = .5 * height - yold;

		if ((.5 * ldim) * (.5 * ldim) - vstart.vals[0] * vstart.vals[0] - vstart.vals[1] * vstart.vals[1] <= 0) return;
		vstart.vals[2] = sqrt((.5 * ldim) * (.5 * ldim) - vstart.vals[0] * vstart.vals[0] - vstart.vals[1] * vstart.vals[1]);
		Normalize3(&vstart);

		xold = mx;
		yold = my;

		// calculate the vectors;
		Vector4 vend;


		vend.vals[0] = mx - .5 * width;
		vend.vals[1] = .5 * height - my;

		if ((.5 * ldim) * (.5 * ldim) - vend.vals[0] * vend.vals[0] - vend.vals[1] * vend.vals[1] <= 0) return;
		vend.vals[2] = sqrt((.5 * ldim) * (.5 * ldim) - vend.vals[0] * vend.vals[0] - vend.vals[1] * vend.vals[1]);
		Normalize3(&vend);

		translateoff[0] += -(XMAX - XMIN) * (vend.vals[0] - vstart.vals[0]);
		translateoff[1] += -(XMAX - XMIN) * (vend.vals[1] - vstart.vals[1]);
		translateoff[2] += -0 * (vend.vals[2] - vstart.vals[2]);
		glutPostRedisplay();
	}

	if (buttonstate == 0) {
		float ldim = (width > height ? width : height);

		if (xold == mx && yold == my) return;

		vstart.vals[0] = xold - .5 * width;
		vstart.vals[1] = .5 * height - yold;

		if ((.5 * ldim) * (.5 * ldim) - vstart.vals[0] * vstart.vals[0] - vstart.vals[1] * vstart.vals[1] < 0) return;
		vstart.vals[2] = sqrt((.5 * ldim) * (.5 * ldim) - vstart.vals[0] * vstart.vals[0] - vstart.vals[1] * vstart.vals[1]);
		if (Normalize3(&vstart) == 0) return;
		xold = mx;
		yold = my;

		// calculate the vectors;
		Vector4 vend;


		vend.vals[0] = mx - .5 * width;
		vend.vals[1] = .5 * height - my;

		if ((.5 * ldim) * (.5 * ldim) - vend.vals[0] * vend.vals[0] - vend.vals[1] * vend.vals[1] <= 0) return;
		vend.vals[2] = sqrt((.5 * ldim) * (.5 * ldim) - vend.vals[0] * vend.vals[0] - vend.vals[1] * vend.vals[1]);
		if (Normalize3(&vend) == 0) return;



		float alpha = InteriorAngle(vstart, vend);
		if (alpha < 0.01) return;

		Vector4 cp = Cross(vstart, vend);
		if (Normalize3(&cp) == 0) return;

		//update the crystal ball matrix
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

		rotationtotal = times(rotationtotal, rotation(cp, alpha * PI / -180.0));

		setRotQuat(rotmatrix, rotationtotal);
		glMultMatrixf(rotmatrix.vals);

	}


	glutPostRedisplay();

}

void initScene() {
	glDepthMask(true);
	if (clearcolor) {
		glClearColor(1, 1, 1, 1);
	}
	else {
		glClearColor(0, 0, 0, 1);
	}

	glDisable(GL_CULL_FACE);
	glDepthFunc(GL_LESS);
	glEnable(GL_DEPTH_TEST);

	glShadeModel(GL_SMOOTH);
	glLineWidth(1);
	glEnable(GL_LINE_SMOOTH);


	// Enable our light.

	glEnable(GL_LIGHT0);

	// Set up the material information for our objects.  Again this is just for show.
	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specularLight);
	glMateriali(GL_FRONT_AND_BACK, GL_SHININESS, 128);


	// initialize crystal ball
	mat = MIdentity();
	dist = 300 * (ZMAX - ZMIN);

}



void Initialize_GLUT(int argc, char** argv) {


	translateoff[0] = 0.0;
	translateoff[1] = 0.0;
	translateoff[2] = 0.0;
	resetQuat(rotationtotal);
	g_mouseActive = false;

	// GLUT initialization.
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(width, height);
	glutInitWindowPosition(570, 80);
	mainWindow = glutCreateWindow("Morse3d Efficient Computation");

	//		glewInit();

	glutDisplayFunc(display);
	glutReshapeFunc(reshapeMainWindow);
	glutKeyboardFunc(graphicKeys);
	glutMouseFunc(mouseClicked);
	//glutSpecialFunc(functionKeys);
	glutMotionFunc(mouseMovement);
	glutIdleFunc(myidle);

	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);

	initScene();
	setGlobalMinMax(usc);



	glutMainLoop();
}
