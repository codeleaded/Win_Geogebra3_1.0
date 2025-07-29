#ifndef MATH3D_H
#define MATH3D_H

#include "/home/codeleaded/System/Static/Container/Vector.h"

typedef struct vec3d{
	float x;
	float y;
	float z;
	float w;
} vec3d;

vec3d vec3d_new(float x,float y,float z){
	return (vec3d){ x,y,z,1.0f };
}
vec3d vec3d_Add(vec3d v1, vec3d v2)
{
	return vec3d_new( v1.x + v2.x, v1.y + v2.y, v1.z + v2.z );
}
vec3d vec3d_Sub(vec3d v1, vec3d v2)
{
	return vec3d_new( v1.x - v2.x, v1.y - v2.y, v1.z - v2.z );
}
vec3d vec3d_Mul(vec3d v1, float k)
{
	return vec3d_new( v1.x * k, v1.y * k, v1.z * k );
}
vec3d vec3d_Div(vec3d v1, float k)
{
	if(k==0.0f) return vec3d_new( 0.0f,0.0f,0.0f );
	return vec3d_new( v1.x / k, v1.y / k, v1.z / k );
}
float vec3d_DotProduct(vec3d v1, vec3d v2)
{
	return v1.x*v2.x + v1.y*v2.y + v1.z * v2.z;
}
float vec3d_Length(vec3d v)
{
	return sqrtf(vec3d_DotProduct(v, v));
}
vec3d vec3d_Normalise(vec3d v){
	float l = vec3d_Length(v);
	return vec3d_new( v.x / l, v.y / l, v.z / l );
}
vec3d vec3d_Perp(vec3d v){
	return vec3d_new( -v.z, v.y, v.x );
}
vec3d vec3d_Round(vec3d v){
	return vec3d_new( roundf(v.x),roundf(v.y),roundf(v.z) );
}

vec3d vec3d_CrossProduct(vec3d v1, vec3d v2){
	vec3d v = vec3d_new(0.0f,0.0f,0.0f);
	v.x = v1.y * v2.z - v1.z * v2.y;
	v.y = v1.z * v2.x - v1.x * v2.z;
	v.z = v1.x * v2.y - v1.y * v2.x;
	return v;
}
vec3d vec3d_IntersectPlane(vec3d plane_p, vec3d* plane_n, vec3d lineStart, vec3d lineEnd){
	*plane_n = vec3d_Normalise(*plane_n);
	float plane_d = -vec3d_DotProduct(*plane_n, plane_p);
	float ad = vec3d_DotProduct(lineStart, *plane_n);
	float bd = vec3d_DotProduct(lineEnd, *plane_n);
	float t = (-plane_d - ad) / (bd - ad);
	vec3d lineStartToEnd = vec3d_Sub(lineEnd, lineStart);
	vec3d lineToIntersect = vec3d_Mul(lineStartToEnd, t);
	return vec3d_Add(lineStart, lineToIntersect);
}

typedef struct triangle{
	vec3d p[3];
	vec3d n;
	unsigned int c;
} triangle;

void Triangle_CalcNorm(triangle* t){
	vec3d line1 = vec3d_Sub(t->p[1],t->p[0]);
	vec3d line2 = vec3d_Sub(t->p[2],t->p[0]);
	vec3d normal = vec3d_CrossProduct(line1,line2);

	t->n = vec3d_Normalise(normal);
}
void Triangle_ShadeNorm(triangle* t,vec3d dirLight){
	float l = F32_Max(0.2f,vec3d_DotProduct(t->n,dirLight));
	t->c = Pixel_Mulf(t->c,l);
}

#define CUBE_SIDE_SOUTH		0
#define CUBE_SIDE_EAST		1
#define CUBE_SIDE_NORTH		2
#define CUBE_SIDE_WEST		3
#define CUBE_SIDE_TOP		4
#define CUBE_SIDE_BOTTOM	5

#define CUBE_SIDE_000		0
#define CUBE_SIDE_001		1
#define CUBE_SIDE_100		2
#define CUBE_SIDE_101		3
#define CUBE_SIDE_002		4
#define CUBE_SIDE_102		5
#define CUBE_SIDE_200		6
#define CUBE_SIDE_201		7
#define CUBE_SIDE_202		8
#define CUBE_SIDE_010		9
#define CUBE_SIDE_011		10
#define CUBE_SIDE_110		11
#define CUBE_SIDE_111		12
#define CUBE_SIDE_012		13
#define CUBE_SIDE_112		14
#define CUBE_SIDE_210		15
#define CUBE_SIDE_211		16
#define CUBE_SIDE_212		17
#define CUBE_SIDE_020		18
#define CUBE_SIDE_021		19
#define CUBE_SIDE_120		20
#define CUBE_SIDE_121		21
#define CUBE_SIDE_022		22
#define CUBE_SIDE_122		23
#define CUBE_SIDE_220		24
#define CUBE_SIDE_221		25
#define CUBE_SIDE_222		26

vec3d Neighbour_Side(int s){
	switch (s){
	case CUBE_SIDE_SOUTH: 	return (vec3d){ 0.0f, 0.0f,-1.0f,1.0f };
	case CUBE_SIDE_EAST: 	return (vec3d){ 1.0f, 0.0f, 0.0f,1.0f };
	case CUBE_SIDE_NORTH: 	return (vec3d){ 0.0f, 0.0f, 1.0f,1.0f };
	case CUBE_SIDE_WEST: 	return (vec3d){-1.0f, 0.0f, 0.0f,1.0f };
	case CUBE_SIDE_TOP: 		return (vec3d){ 0.0f, 1.0f, 0.0f,1.0f };
	case CUBE_SIDE_BOTTOM: 	return (vec3d){ 0.0f,-1.0f, 0.0f,1.0f };
	}
	return (vec3d){ 0.0f,0.0f,0.0f,1.0f };
}
int Side_Other(int s){
	switch (s){
	case CUBE_SIDE_SOUTH: 	return CUBE_SIDE_NORTH;
	case CUBE_SIDE_EAST: 	return CUBE_SIDE_WEST;
	case CUBE_SIDE_NORTH: 	return CUBE_SIDE_SOUTH;
	case CUBE_SIDE_WEST: 	return CUBE_SIDE_EAST;
	case CUBE_SIDE_TOP: 	return CUBE_SIDE_BOTTOM;
	case CUBE_SIDE_BOTTOM: 	return CUBE_SIDE_TOP;
	}
	return CUBE_SIDE_SOUTH;
}
vec3d Neighbour_AllSide(int s){
	switch (s){
		case CUBE_SIDE_000: return (vec3d){-1.0f,-1.0f,-1.0f,1.0f };
		case CUBE_SIDE_001: return (vec3d){-1.0f,-1.0f, 0.0f,1.0f };
		case CUBE_SIDE_100: return (vec3d){ 0.0f,-1.0f,-1.0f,1.0f };
		case CUBE_SIDE_101: return (vec3d){ 0.0f,-1.0f, 0.0f,1.0f };
		case CUBE_SIDE_002: return (vec3d){-1.0f,-1.0f, 1.0f,1.0f };
		case CUBE_SIDE_102: return (vec3d){ 0.0f,-1.0f, 1.0f,1.0f };
		case CUBE_SIDE_200: return (vec3d){ 1.0f,-1.0f,-1.0f,1.0f };
		case CUBE_SIDE_201: return (vec3d){ 1.0f,-1.0f, 0.0f,1.0f };
		case CUBE_SIDE_202: return (vec3d){ 1.0f,-1.0f, 1.0f,1.0f };
		case CUBE_SIDE_010: return (vec3d){-1.0f, 0.0f,-1.0f,1.0f };
		case CUBE_SIDE_011: return (vec3d){-1.0f, 0.0f, 0.0f,1.0f };
		case CUBE_SIDE_110: return (vec3d){ 0.0f, 0.0f,-1.0f,1.0f };
		case CUBE_SIDE_111: return (vec3d){ 0.0f, 0.0f, 0.0f,1.0f };
		case CUBE_SIDE_012: return (vec3d){-1.0f, 0.0f, 1.0f,1.0f };
		case CUBE_SIDE_112: return (vec3d){ 0.0f, 0.0f, 1.0f,1.0f };
		case CUBE_SIDE_210: return (vec3d){ 1.0f, 0.0f,-1.0f,1.0f };
		case CUBE_SIDE_211: return (vec3d){ 1.0f, 0.0f, 0.0f,1.0f };
		case CUBE_SIDE_212: return (vec3d){ 1.0f, 0.0f, 1.0f,1.0f };
		case CUBE_SIDE_020: return (vec3d){-1.0f, 1.0f,-1.0f,1.0f };
		case CUBE_SIDE_021: return (vec3d){-1.0f, 1.0f, 0.0f,1.0f };
		case CUBE_SIDE_120: return (vec3d){ 0.0f, 1.0f,-1.0f,1.0f };
		case CUBE_SIDE_121: return (vec3d){ 0.0f, 1.0f, 0.0f,1.0f };
		case CUBE_SIDE_022: return (vec3d){-1.0f, 1.0f, 1.0f,1.0f };
		case CUBE_SIDE_122: return (vec3d){ 0.0f, 1.0f, 1.0f,1.0f };
		case CUBE_SIDE_220: return (vec3d){ 1.0f, 1.0f,-1.0f,1.0f };
		case CUBE_SIDE_221: return (vec3d){ 1.0f, 1.0f, 0.0f,1.0f };
		case CUBE_SIDE_222: return (vec3d){ 1.0f, 1.0f, 1.0f,1.0f };
	}
	return (vec3d){ 0.0f,0.0f,0.0f,1.0f };
}

typedef struct mesh{
	Vector tris;
} mesh;

typedef struct mat4x4{
	float m[4][4];
} mat4x4;

mat4x4 Matrix_Null()
{
	return (mat4x4){{
		{ 0.0f,0.0f,0.0f,0.0f },
		{ 0.0f,0.0f,0.0f,0.0f },
		{ 0.0f,0.0f,0.0f,0.0f },
		{ 0.0f,0.0f,0.0f,0.0f }
	}};
}
vec3d Matrix_MultiplyVector(mat4x4 m, vec3d i)
{
	vec3d v;
	v.x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + i.w * m.m[3][0];
	v.y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + i.w * m.m[3][1];
	v.z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + i.w * m.m[3][2];
	v.w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + i.w * m.m[3][3];
	return v;
}
mat4x4 Matrix_MakeIdentity()
{
	mat4x4 matrix = Matrix_Null();
	matrix.m[0][0] = 1.0f;
	matrix.m[1][1] = 1.0f;
	matrix.m[2][2] = 1.0f;
	matrix.m[3][3] = 1.0f;
	return matrix;
}
mat4x4 Matrix_MakeRotationX(float fAngleRad)
{
	mat4x4 matrix = Matrix_Null();
	matrix.m[0][0] = 1.0f;
	matrix.m[1][1] = cosf(fAngleRad);
	matrix.m[1][2] = sinf(fAngleRad);
	matrix.m[2][1] = -sinf(fAngleRad);
	matrix.m[2][2] = cosf(fAngleRad);
	matrix.m[3][3] = 1.0f;
	return matrix;
}
mat4x4 Matrix_MakeRotationY(float fAngleRad)
{
	mat4x4 matrix = Matrix_Null();
	matrix.m[0][0] = cosf(fAngleRad);
	matrix.m[0][2] = sinf(fAngleRad);
	matrix.m[2][0] = -sinf(fAngleRad);
	matrix.m[1][1] = 1.0f;
	matrix.m[2][2] = cosf(fAngleRad);
	matrix.m[3][3] = 1.0f;
	return matrix;
}
mat4x4 Matrix_MakeRotationZ(float fAngleRad)
{
	mat4x4 matrix = Matrix_Null();
	matrix.m[0][0] = cosf(fAngleRad);
	matrix.m[0][1] = sinf(fAngleRad);
	matrix.m[1][0] = -sinf(fAngleRad);
	matrix.m[1][1] = cosf(fAngleRad);
	matrix.m[2][2] = 1.0f;
	matrix.m[3][3] = 1.0f;
	return matrix;
}
mat4x4 Matrix_MakeTranslation(float x, float y, float z)
{
	mat4x4 matrix = Matrix_Null();
	matrix.m[0][0] = 1.0f;
	matrix.m[1][1] = 1.0f;
	matrix.m[2][2] = 1.0f;
	matrix.m[3][3] = 1.0f;
	matrix.m[3][0] = x;
	matrix.m[3][1] = y;
	matrix.m[3][2] = z;
	return matrix;
}
mat4x4 Matrix_MakeProjection(float fFovDegrees, float fAspectRatio, float fNear, float fFar)
{
	float fFovRad = 1.0f / tanf(fFovDegrees * 0.5f / 180.0f * 3.14159f);
	mat4x4 matrix = Matrix_Null();
	matrix.m[0][0] = fAspectRatio * fFovRad;
	matrix.m[1][1] = fFovRad;
	matrix.m[2][2] = fFar / (fFar - fNear);
	matrix.m[3][2] = (-fFar * fNear) / (fFar - fNear);
	matrix.m[2][3] = 1.0f;
	matrix.m[3][3] = 0.0f;
	return matrix;
}
mat4x4 Matrix_MultiplyMatrix(mat4x4 m1, mat4x4 m2)
{
	mat4x4 matrix = Matrix_Null();
	for (int c = 0; c < 4; c++)
		for (int r = 0; r < 4; r++)
			matrix.m[r][c] = m1.m[r][0] * m2.m[0][c] + m1.m[r][1] * m2.m[1][c] + m1.m[r][2] * m2.m[2][c] + m1.m[r][3] * m2.m[3][c];
	return matrix;
}

float dist(vec3d plane_p, vec3d plane_n,vec3d p){
	vec3d n = vec3d_Normalise(p);
	return (plane_n.x * p.x + plane_n.y * p.y + plane_n.z * p.z - vec3d_DotProduct(plane_n, plane_p));
}

int Triangle_ClipAgainstPlane(vec3d plane_p, vec3d plane_n, triangle in_tri, triangle* out_tri1, triangle* out_tri2){
	// Make sure plane normal is indeed normal
	plane_n = vec3d_Normalise(plane_n);
	
	// Create two temporary storage arrays to classify points either side of plane
	// If distance sign is positive, point lies on "inside" of plane
	vec3d* inCUBE_SIDE_points[3];  int nInsidePointCount = 0;
	vec3d* outCUBE_SIDE_points[3]; int nOutsidePointCount = 0;
	// Get signed distance of each point in triangle to plane
	float d0 = dist(plane_p,plane_n,in_tri.p[0]);
	float d1 = dist(plane_p,plane_n,in_tri.p[1]);
	float d2 = dist(plane_p,plane_n,in_tri.p[2]);
	if (d0 >= 0) { inCUBE_SIDE_points[nInsidePointCount++] = &in_tri.p[0]; }
	else { outCUBE_SIDE_points[nOutsidePointCount++] = &in_tri.p[0]; }
	if (d1 >= 0) { inCUBE_SIDE_points[nInsidePointCount++] = &in_tri.p[1]; }
	else { outCUBE_SIDE_points[nOutsidePointCount++] = &in_tri.p[1]; }
	if (d2 >= 0) { inCUBE_SIDE_points[nInsidePointCount++] = &in_tri.p[2]; }
	else { outCUBE_SIDE_points[nOutsidePointCount++] = &in_tri.p[2]; }
	// Now classify triangle points, and break the input triangle into 
	// smaller output triangles if required. There are four possible
	// outcomes...
	if (nInsidePointCount == 0)
	{
		// All points lie on the outside of plane, so clip whole triangle
		// It ceases to exist
		return 0; // No returned triangles are valid
	}
	if (nInsidePointCount == 3)
	{
		// All points lie on the inside of plane, so do nothing
		// and allow the triangle to simply pass through
		*out_tri1 = in_tri;
		return 1; // Just the one returned original triangle is valid
	}
	if (nInsidePointCount == 1 && nOutsidePointCount == 2)
	{
		// Triangle should be clipped. As two points lie outside
		// the plane, the triangle simply becomes a smaller triangle
		// Copy appearance info to new triangle
		out_tri1->c =  in_tri.c;
		// The inside point is valid, so keep that...
		out_tri1->p[0] = *inCUBE_SIDE_points[0];
		// but the two new points are at the locations where the 
		// original sides of the triangle (lines) intersect with the plane
		out_tri1->p[1] = vec3d_IntersectPlane(plane_p, &plane_n, *inCUBE_SIDE_points[0], *outCUBE_SIDE_points[0]);
		out_tri1->p[2] = vec3d_IntersectPlane(plane_p, &plane_n, *inCUBE_SIDE_points[0], *outCUBE_SIDE_points[1]);
		return 1; // Return the newly formed single triangle
	}
	if (nInsidePointCount == 2 && nOutsidePointCount == 1)
	{
		// Triangle should be clipped. As two points lie inside the plane,
		// the clipped triangle becomes a "quad". Fortunately, we can
		// represent a quad with two new triangles
		// Copy appearance info to new triangles
		out_tri1->c =  in_tri.c;
		out_tri2->c =  in_tri.c;
		// The first triangle consists of the two inside points and a new
		// point determined by the location where one side of the triangle
		// intersects with the plane
		out_tri1->p[0] = *inCUBE_SIDE_points[0];
		out_tri1->p[1] = *inCUBE_SIDE_points[1];
		out_tri1->p[2] = vec3d_IntersectPlane(plane_p, &plane_n, *inCUBE_SIDE_points[0], *outCUBE_SIDE_points[0]);
		// The second triangle is composed of one of he inside points, a
		// new point determined by the intersection of the other side of the 
		// triangle and the plane, and the newly created point above
		out_tri2->p[0] = *inCUBE_SIDE_points[1];
		out_tri2->p[1] = out_tri1->p[2];
		out_tri2->p[2] = vec3d_IntersectPlane(plane_p, &plane_n, *inCUBE_SIDE_points[1], *outCUBE_SIDE_points[0]);
		return 2; // Return two newly formed triangles which form a quad
	}
}

mat4x4 Matrix_PointAt(vec3d pos, vec3d target, vec3d up)
{
	// Calculate new forward direction
	vec3d newForward = vec3d_Sub(target, pos);
	newForward = vec3d_Normalise(newForward);
	// Calculate new Up direction
	vec3d a = vec3d_Mul(newForward, vec3d_DotProduct(up,newForward));
	vec3d newUp = vec3d_Sub(up, a);
	newUp = vec3d_Normalise(newUp);
	// New Right direction is easy, its just cross product
	vec3d newRight = vec3d_CrossProduct(newUp, newForward);
	// Construct Dimensioning and Translation Matrix	
	mat4x4 matrix;
	matrix.m[0][0] = newRight.x;	matrix.m[0][1] = newRight.y;	matrix.m[0][2] = newRight.z;	matrix.m[0][3] = 0.0f;
	matrix.m[1][0] = newUp.x;		matrix.m[1][1] = newUp.y;		matrix.m[1][2] = newUp.z;		matrix.m[1][3] = 0.0f;
	matrix.m[2][0] = newForward.x;	matrix.m[2][1] = newForward.y;	matrix.m[2][2] = newForward.z;	matrix.m[2][3] = 0.0f;
	matrix.m[3][0] = pos.x;			matrix.m[3][1] = pos.y;			matrix.m[3][2] = pos.z;			matrix.m[3][3] = 1.0f;
	return matrix;
}
mat4x4 Matrix_QuickInverse(mat4x4 m) // Only for Rotation/Translation Matrices
{
	mat4x4 matrix;
	matrix.m[0][0] = m.m[0][0]; matrix.m[0][1] = m.m[1][0]; matrix.m[0][2] = m.m[2][0]; matrix.m[0][3] = 0.0f;
	matrix.m[1][0] = m.m[0][1]; matrix.m[1][1] = m.m[1][1]; matrix.m[1][2] = m.m[2][1]; matrix.m[1][3] = 0.0f;
	matrix.m[2][0] = m.m[0][2]; matrix.m[2][1] = m.m[1][2]; matrix.m[2][2] = m.m[2][2]; matrix.m[2][3] = 0.0f;
	matrix.m[3][0] = -(m.m[3][0] * matrix.m[0][0] + m.m[3][1] * matrix.m[1][0] + m.m[3][2] * matrix.m[2][0]);
	matrix.m[3][1] = -(m.m[3][0] * matrix.m[0][1] + m.m[3][1] * matrix.m[1][1] + m.m[3][2] * matrix.m[2][1]);
	matrix.m[3][2] = -(m.m[3][0] * matrix.m[0][2] + m.m[3][1] * matrix.m[1][2] + m.m[3][2] * matrix.m[2][2]);
	matrix.m[3][3] = 1.0f;
	return matrix;
}


typedef struct Camera {
	mat4x4 matProj;
	vec3d vCamera;
	vec3d vVelocity;
	vec3d vLookDir;
	float fYaw;	
	float fPitch;
	vec3d vLength;
} Camera;

Camera Camera_New(vec3d p){
	Camera c;
	c.matProj = Matrix_MakeProjection(90.0f,(float)GetHeight() / (float)GetWidth(),0.1f,1000.0f);
	c.vCamera = p;
	c.vVelocity = (vec3d){ 0.0f,0.0f,0.0f,1.0f };
	c.vLookDir = (vec3d){ 0.0f,0.0f,0.0f,1.0f };
	c.fYaw = 0.0f;
	c.fPitch = 0.0f;
	c.vLength = (vec3d){ 0.0f,0.0f,0.0f,1.0f };
	return c;
}

#endif