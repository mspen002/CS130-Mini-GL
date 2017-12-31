/**
 * minigl.cpp
 * -------------------------------
 * Implement miniGL here.
 *
 * You may include minigl.h and any of the standard C++ libraries.
 * No other includes are permitted.  Other preprocessing directives
 * are also not permitted.  These requirements are strictly
 * enforced.  Be sure to run a test grading to make sure your file
 * passes the sanity tests.
 *
 * The behavior of the routines your are implenting is documented here:
 * https://www.opengl.org/sdk/docs/man2/
 * Note that you will only be implementing a subset of this.  In particular,
 * you only need to implement enough to pass the tests in the suite.
 */

#include "minigl.h"
#include "vec.h"
#include "mat.h"
#include <algorithm>
#include <cassert>
#include <stack>
#include <cmath>
#include <vector>
#include <cstdio>

using namespace std;

/**
 * Useful data types
 */
typedef mat<MGLfloat,4,4> mat4; //data structure storing a 4x4 matrix, see mat.h
typedef mat<MGLfloat,3,3> mat3; //data structure storing a 3x3 matrix, see mat.h
typedef vec<MGLfloat,4> vec4;   //data structure storing a 4 dimensional vector, see vec.h
typedef vec<MGLfloat,3> vec3;   //data structure storing a 3 dimensional vector, see vec.h
typedef vec<MGLfloat,2> vec2;   //data structure storing a 2 dimensional vector, see vec.h

//MGLfloat zero = 0;

MGLpoly_mode draw_mode;
MGLmatrix_mode mat_mode;

vec3 current_color;

struct vertex{
vec4 loc;
vec3 color;
};

vector<vertex> list;
vector<vec3> col;

struct triangle{
vertex a;
vertex b;
vertex c;
};

vector<triangle> trilist;
vector<triangle> screenList;

mat4 projmatrix;
mat4 modelmatrix;
stack<mat4> proj;
stack<mat4> model;

MGLfloat xmin;
MGLfloat xmax;
MGLfloat ymin;
MGLfloat ymax;

triangle t;
triangle q;
vertex v;


vertex matrixMult(vertex v, mat4 m)
{
	vertex store = v;
	MGLfloat total;
	MGLfloat f;

	//cout << m(0,0);// << proj.at(0)(0,0);
	
	for(int x = 0; x <= 3; x++)
	{
	total = 0;
	for(int y = 0; y <= 3; y++)
	{
		f = v.loc[y];
		total += m(x,y) * f;
		//cout << m(x,y)<< ", ";
	}
	store.loc[x] = total;
	//cout << store.loc[x] << ", ";
	}
	return store;
}

void norm(MGLfloat &x, MGLfloat &y, MGLfloat &z)
{
	MGLfloat check = sqrt(x*x + y*y + z*z);
	if(check == 0)
		return;
	x = x/check;
	y = y/check;
	z = z/check;
}
vertex screenScale(MGLsize width, MGLsize height, vertex v)
{
	if(v.loc[3] != 0)
	{
		v.loc[0] = v.loc[0]/v.loc[3];
		v.loc[1] = v.loc[1]/v.loc[3];
		v.loc[2] = v.loc[2]/v.loc[3];
		//v.loc[3] = v.loc[3]/v.loc[3];
	}
	v.color = v.color;

	v.loc[0] = ((v.loc[0]+1)*(width/2));
	v.loc[1] = ((v.loc[1]+1)*(height/2));
	//cout << (int)v.loc[0] << ", " << (int)v.loc[1] << ", ";
	//v.loc[0] = floor((v.loc[0] + 1)/(width/2));	

	//cout << (int)v.loc[0] << ", " << (int)v.loc[1] << ", ";
	if(v.loc[0] > width)
		v.loc[0] = width;
	else if(v.loc[0] < 0)
		v.loc[0] = 0;
	if(v.loc[1] > height)
		v.loc[1] = height;
	else if(v.loc[1] < 0)
		v.loc[1] = 0;
	return v;
	
}


MGLfloat bary(vertex a, vertex b, vertex c)
{
	return (b.loc[0] - a.loc[0]) * (c.loc[1] -  a.loc[1]) - (b.loc[1] - a.loc[1]) * (c.loc[0] - a.loc[0]);
}

MGLpixel colorpicker(MGLfloat alpha, MGLfloat beta, MGLfloat gamma, vertex v1, vertex v2, vertex v3)
{
	MGLfloat red = alpha*v1.color[0] + beta*v2.color[0] + gamma*v3.color[0];
	MGLfloat green = alpha*v1.color[1] + beta*v2.color[1] + gamma*v3.color[1];
	MGLfloat blue = alpha*v1.color[2] + beta*v2.color[2] + gamma*v3.color[2];
	
	//cout << red << " ";
	return Make_Pixel(red, green, blue);
}

/**
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description) {
    printf("%s\n", description);
    exit(1);
}


/**
 * Read pixel data starting with the pixel at coordinates
 * (0, 0), up to (width,  height), into the array
 * pointed to by data.  The boundaries are lower-inclusive,
 * that is, a call with width = height = 1 would just read
 * the pixel at (0, 0).
 *
 * Rasterization and z-buffering should be performed when
 * this function is called, so that the data array is filled
 * with the actual pixel values that should be displayed on
 * the two-dimensional screen.
 */
void mglReadPixels(MGLsize width,
                   MGLsize height,
                   MGLpixel *data)
{

	for(int i = 0; i < (int)trilist.size(); i++)
	{
		trilist.at(i).a = screenScale(width, height, trilist.at(i).a);
		trilist.at(i).b = screenScale(width, height, trilist.at(i).b);
		trilist.at(i).c = screenScale(width, height, trilist.at(i).c);
	}
	
	//cout << xmax<< endl;
	int pixels = 0;
	MGLfloat zip = 0;
	MGLfloat w = width-1;
	MGLfloat h = height-1;
	vector<int> zbuff(width*height);
	fill(zbuff.begin(),zbuff.end(), 999999);
	for(int i = 0; i < (int)trilist.size(); i++)
	{
	//MGLfloat xmin = max(zip, min(min(trilist.at(i).a.loc[0],trilist.at(i).b.loc[0]),trilist.at(i).c.loc[0]));
	//MGLfloat ymin = max(zip, min(min(trilist.at(i).a.loc[1],trilist.at(i).b.loc[1]),trilist.at(i).c.loc[1]));
	//MGLfloat xmax = min(w, max(max(trilist.at(i).a.loc[0],trilist.at(i).b.loc[0]),trilist.at(i).c.loc[0]));
	//MGLfloat ymax = min(h, max(max(trilist.at(i).a.loc[1],trilist.at(i).b.loc[1]),trilist.at(i).c.loc[1]));
	for(MGLfloat y = 0; y < height; y++)
	{
	for(MGLfloat x = 0; x < width; x++) 
	{
		pixels = x + y * width;
		//cout << "Running! \n";
		vertex n;
		n.loc[0] = x;
		n.loc[1] = y;
		n.loc[2] = 1;
		n.loc[3] = 1;
		MGLfloat area = trilist.at(i).a.loc[0]*(trilist.at(i).b.loc[1]-trilist.at(i).c.loc[1])
				 + trilist.at(i).a.loc[1]*(trilist.at(i).c.loc[0]-trilist.at(i).b.loc[0])
				 + (trilist.at(i).b.loc[0]*trilist.at(i).c.loc[1]-trilist.at(i).b.loc[1]
				 * trilist.at(i).c.loc[0]);

		MGLfloat a1 = (n.loc[0]*(trilist.at(i).b.loc[1]-trilist.at(i).c.loc[1])
				 + n.loc[1]*(trilist.at(i).c.loc[0]-trilist.at(i).b.loc[0])
				 + (trilist.at(i).b.loc[0]*trilist.at(i).c.loc[1] - trilist.at(i).b.loc[1]
				 * trilist.at(i).c.loc[0]))/area; 

		MGLfloat b1 = (trilist.at(i).a.loc[0]*(n.loc[1]-trilist.at(i).c.loc[1])
				 + trilist.at(i).a.loc[1]*(trilist.at(i).c.loc[0]-n.loc[0])
				 + (n.loc[0]*trilist.at(i).c.loc[1] - n.loc[1]*trilist.at(i).c.loc[0]))/area;

		MGLfloat c1 = (trilist.at(i).a.loc[0]*(trilist.at(i).b.loc[1]-n.loc[1])
				 + trilist.at(i).a.loc[1]*(n.loc[0]-trilist.at(i).b.loc[0])
				 + (trilist.at(i).b.loc[0]*n.loc[1] - trilist.at(i).b.loc[1]*n.loc[0]))/area;

		MGLfloat za = trilist.at(i).a.loc[2]/trilist.at(i).a.loc[3];
		MGLfloat zb = trilist.at(i).b.loc[2]/trilist.at(i).b.loc[3];
		MGLfloat zc = trilist.at(i).c.loc[2]/trilist.at(i).c.loc[3];
		//norm(a1,b1,c1);
		MGLfloat z = a1*za + b1*zb + c1*zc;
		//cout << a1 << b1 << c1;
		if(a1 >= 0 && b1 >= 0 && c1 >= 0 && z <= zbuff.at(pixels) && z <= 1 && z >= -1)
		{
			MGLfloat pa = a1;
//(a1*trilist.at(i).a.loc[3]/(a1*trilist.at(i).a.loc[3] + b1*trilist.at(i).b.loc[3] + c1*trilist.at(i).c.loc[3]));
			MGLfloat pb = b1;
//(b1*trilist.at(i).b.loc[3]/(a1*trilist.at(i).a.loc[3] + b1*trilist.at(i).b.loc[3] + c1*trilist.at(i).c.loc[3]));
			MGLfloat pc = c1;
//(c1*trilist.at(i).c.loc[3]/(a1*trilist.at(i).a.loc[3] + b1*trilist.at(i).b.loc[3] + c1*trilist.at(i).c.loc[3]));
			a1 = (pa/trilist.at(i).a.loc[3]/(pa/trilist.at(i).a.loc[3] + pb/trilist.at(i).b.loc[3] + pc/trilist.at(i).c.loc[3]));
			b1 = (pb/trilist.at(i).b.loc[3]/(pa/trilist.at(i).a.loc[3] + pb/trilist.at(i).b.loc[3] + pc/trilist.at(i).c.loc[3]));
			c1 = (pc/trilist.at(i).c.loc[3]/(pa/trilist.at(i).a.loc[3] + pb/trilist.at(i).b.loc[3] + pc/trilist.at(i).c.loc[3]));
			zbuff.at(pixels) = z;
			MGLpixel pix = colorpicker(a1, b1, c1, trilist.at(i).a, trilist.at(i).b, trilist.at(i).c);
			data[pixels] = pix;
		}
	}
	}
	pixels = 0;
	}
trilist.clear();
}

/**
 * Start specifying the vertices for a group of primitives,
 * whose type is specified by the given mode.
 */
void mglBegin(MGLpoly_mode mode)
{
	draw_mode = mode;
}


/**
 * Stop specifying the vertices for a group of primitives.
 */
void mglEnd()
{
if(draw_mode == MGL_TRIANGLES)
{
	for(unsigned int i = 0; i < list.size(); i+=3)
	{
		t.a = list.at(i);
		t.b = list.at(i+1);
		t.c = list.at(i+2);
		trilist.push_back(t);
	}
list.clear();
}
if(draw_mode == MGL_QUADS)
{
	for(unsigned int i = 0; i < list.size(); i+=4)
	{
		t.a = list.at(i);
		t.b = list.at(i+1);
		t.c = list.at(i+2);
		trilist.push_back(t);
		q.a = list.at(i);
		q.b = list.at(i+2);
		q.c = list.at(i+3);
		trilist.push_back(q);
	}
list.clear();
}
}

/**
 * Specify a two-dimensional vertex; the x- and y-coordinates
 * are explicitly specified, while the z-coordinate is assumed
 * to be zero.  Must appear between calls to mglBegin() and
 * mglEnd().
 */
void mglVertex2(MGLfloat x,
                MGLfloat y)
{
	mglVertex3(x,y,0);
}

/**
 * Specify a three-dimensional vertex.  Must appear between
 * calls to mglBegin() and mglEnd().
 */
void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{
	v.loc = vec4(x,y,z,1);
	v = matrixMult(v, modelmatrix);
	v = matrixMult(v, projmatrix);
	//cout << v.loc[0];
	//cout << proj.at(0) << endl;
	v.color = current_color;
	list.push_back(v);
}

/**
 * Set the current matrix mode (modelview or projection).
 */
void mglMatrixMode(MGLmatrix_mode mode)
{
	mat_mode = mode;
}

/**
 * Push a copy of the current matrix onto the stack for the
 * current matrix mode.
 */
void mglPushMatrix()
{
	if(mat_mode == MGL_MODELVIEW)
	{
		model.push(modelmatrix);
	}
	if(mat_mode == MGL_PROJECTION)
	{
		proj.push(projmatrix);
	}
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
	if(mat_mode == MGL_MODELVIEW)
	{
		modelmatrix = model.top();
		model.pop();
	}
	if(mat_mode == MGL_PROJECTION)
	{
		projmatrix = proj.top();
		proj.pop();
	}
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
	mat4 identity;
	identity.make_zero();
	identity(0,0) = 1;
	identity(1,1) = 1;
	identity(2,2) = 1;
	identity(3,3) = 1;
	//current_matrix()=identity;
	//matrix = identity;
	if(mat_mode == MGL_MODELVIEW)
	{
		modelmatrix.make_zero();
		modelmatrix = identity;
	}
	if(mat_mode == MGL_PROJECTION)
	{
		projmatrix.make_zero();
		projmatrix = identity;
	}
}
/*
mat4& current_matrix()
{
if(mat_mode==MGL_MODELVIEW)
   return model.end();
else 
   return proj.end();
}*/

/**
 * Replace the current matrix with an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglLoadMatrix(const MGLfloat *matrix)
{
	mat4 load;
	load.make_zero();
	
	load(0,0) = matrix[0];
	load(1,0) = matrix[1];
	load(2,0) = matrix[2];
	load(3,0) = matrix[3];	
	load(0,1) = matrix[4];
	load(1,1) = matrix[5];
	load(2,1) = matrix[6];
	load(3,1) = matrix[7];
	load(0,2) = matrix[8];
	load(1,2) = matrix[9];
	load(2,2) = matrix[10];
	load(3,2) = matrix[11];
	load(0,3) = matrix[12];
	load(1,3) = matrix[13];
	load(2,3) = matrix[14];
	load(3,3) = matrix[15];
	
	if(mat_mode == MGL_MODELVIEW)
	{
		modelmatrix = load;
	}
	if(mat_mode == MGL_PROJECTION)
	{
		projmatrix = load;
	}
}

/**
 * Multiply the current matrix by an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglMultMatrix(const MGLfloat *matrix)
{
	mat4 mult;
	mult.make_zero();
	
	mult(0,0) = matrix[0];
	mult(1,0) = matrix[1];
	mult(2,0) = matrix[2];
	mult(3,0) = matrix[3];	

	mult(0,1) = matrix[4];
	mult(1,1) = matrix[5];
	mult(2,1) = matrix[6];
	mult(3,1) = matrix[7];

	mult(0,2) = matrix[8];
	mult(1,2) = matrix[9];
	mult(2,2) = matrix[10];
	mult(3,2) = matrix[11];

	mult(0,3) = matrix[12];
	mult(1,3) = matrix[13];
	mult(2,3) = matrix[14];
	mult(3,3) = matrix[15];
	
	if(mat_mode == MGL_MODELVIEW)
	{
		modelmatrix = modelmatrix*mult;
	}
	if(mat_mode == MGL_PROJECTION)
	{
		projmatrix = projmatrix*mult;
	}
}

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
	mat4 tran;
	tran.make_zero();
	
	tran(0,0) = 1;
	tran(1,1) = 1;
	tran(2,2) = 1;
	tran(3,3) = 1;	
	tran(0,3) = x;
	tran(1,3) = y;
	tran(2,3) = z;
	
	if(mat_mode == MGL_MODELVIEW)
	{
		modelmatrix = modelmatrix*tran;
	}
	if(mat_mode == MGL_PROJECTION)
	{
		projmatrix = projmatrix*tran;
	}
}

/**
 * Multiply the current matrix by the rotation matrix
 * for a rotation of (angle) degrees about the vector
 * from the origin to the point (x, y, z).
 */
void mglRotate(MGLfloat angle,
               MGLfloat x,
               MGLfloat y,
               MGLfloat z)
{
	MGLfloat s = sin(angle*M_PI/180);
	MGLfloat c = cos(angle*M_PI/180);
	mat4 rotate;
	rotate.make_zero();

	norm(x,y,z);

	rotate(0,0) = x*x*(1-c) + c;
	rotate(0,1) = x*y*(1-c) - z*s;
	rotate(0,2) = x*z*(1-c) + y*s;

	rotate(1,0) = y*x*(1-c) + z*s;
	rotate(1,1) = y*y*(1-c) + c;
	rotate(1,2) = y*z*(1-c) - z*s;

	rotate(2,0) = x*z*(1-c) - y*s;
	rotate(2,1) = y*z*(1-c) + x*s;
	rotate(2,2) = z*z*(1-c) + c;

	rotate(3,3) = 1;
	
	if(mat_mode == MGL_MODELVIEW)
	{
		modelmatrix = modelmatrix*rotate;
	}
	if(mat_mode == MGL_PROJECTION)
	{
		projmatrix = projmatrix*rotate;
	}
}

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
	mat4 scale;
	scale.make_zero();
	
	scale(0,0) = x;
	scale(1,1) = y;
	scale(2,2) = z;
	scale(3,3) = 1;	
	
	if(mat_mode == MGL_MODELVIEW)
	{
		modelmatrix = modelmatrix*scale;
	}
	if(mat_mode == MGL_PROJECTION)
	{
		projmatrix = projmatrix*scale;
	}
}

/**
 * Multiply the current matrix by the perspective matrix
 * with the given clipping plane coordinates.
 */
void mglFrustum(MGLfloat left,
                MGLfloat right,
                MGLfloat bottom,
                MGLfloat top,
                MGLfloat near,
                MGLfloat far)
{
	MGLfloat x1,y1,z1;
	
	x1 = (right + left)/(right - left);
	y1 = (top + bottom)/(top - bottom);
	z1 = -2*(far*near)/(far - near);
	
	mat4 frust;
	frust.make_zero();
	frust(0,0) = (2*near)/(right - left);
	frust(1,1) = (2*near)/(top - bottom);
	frust(2,2) = -(far + near)/(far - near);
	frust(3,2) = -1;
	
	
	frust(0,2) = x1;
	frust(1,2) = y1;
	frust(2,3) = z1;	
	
	if(mat_mode == MGL_MODELVIEW)
	{
		modelmatrix = modelmatrix*frust;
	}
	if(mat_mode == MGL_PROJECTION)
	{
		projmatrix = projmatrix*frust;
	}
}

/**
 * Multiply the current matrix by the orthographic matrix
 * with the given clipping plane coordinates.
 */
void mglOrtho(MGLfloat left,
              MGLfloat right,
              MGLfloat bottom,
              MGLfloat top,
              MGLfloat near,
              MGLfloat far)
{
	MGLfloat x1,y1,z1;
	
	x1 = -(right + left)/(right - left);
	y1 = -(top + bottom)/(top - bottom);
	z1 = -(far + near)/(far - near);
	
	mat4 orth;
	orth.make_zero();
	orth(0,0) = (MGLfloat)2/(right - left);
	orth(1,1) = (MGLfloat)2.0/(top - bottom);
	orth(2,2) = (MGLfloat)-2.0/(far - near);
	orth(3,3) = 1;
	
	orth(0,3) = x1;
	orth(1,3) = y1;
	orth(2,3) = z1;	
	
	if(mat_mode == MGL_MODELVIEW)
	{
		modelmatrix = modelmatrix*orth;
	}
	if(mat_mode == MGL_PROJECTION)
	{
		projmatrix = projmatrix*orth;
		//cout << proj.at(0)(0,0);
	}
}

/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLfloat red,
              MGLfloat green,
              MGLfloat blue)
{
	current_color = vec3(red*255, green*255, blue*255);
}
