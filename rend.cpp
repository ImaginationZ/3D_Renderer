#include	"stdafx.h"
#include	"stdio.h"
#define _USE_MATH_DEFINES
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"
#include	"disp.h"
#include <algorithm>

typedef struct {
	float startX;
	float endX;
	float startY;
	float endY;
} BoundingBox;

typedef struct {
	GzCoord screen;
	GzCoord modelPos;
	GzCoord modelNormal;
	GzCoord worldPos;
	GzCoord worldNormal;
	GzCoord image;
	GzCoord normal;
	GzTextureIndex uv;
	GzColor color;
} TriVertex;

typedef struct {
	TriVertex start;
	TriVertex end;
	float A;
	float B;
	float C;
	bool isColored;
} TriEdge;

typedef struct {
	float A;
	float B;
	float C;
	float D;
} Plane;

int RenderTriangle(GzRender *render, GzCoord* modelSpaceVerticies, GzCoord* modelSpaceNormals, GzTextureIndex* uvList);
int MultiplyVectorByCoefficient(GzCoord vec, GzCoord coef, GzCoord product);
int GetIdentityMatrix(GzMatrix mat);
int ScaleAllMatrixTerms(GzMatrix mat, float scaleFactor);
int PushImageMatrix(GzRender *render, GzMatrix matrix);
int PushNormMatrix(GzRender *render, GzMatrix matrix);
int PushWorldMatrix(GzRender *render, GzMatrix matrix);
int GetNormalPlane(TriEdge* edges, int vectorComponent, Plane* plane);
float GetPerspectiveFactor(float screenSpaceZ);
int GetReflectionAcrossNormal(GzRender* render, GzCoord worldPos, GzCoord worldNorm, GzCoord* reflectionOutput);
int GetWorldPosPlane(TriEdge* edges, int vectorComponent, Plane* plane);
int GetWorldNormalPlane(TriEdge* edges, int vectorComponent, Plane* plane);
void getSortedEdges(GzRender* render, GzCoord* screenVerticies, GzCoord* modelVerticies, GzCoord* modelNormals, GzTextureIndex* uvList, TriEdge* edges);
int GetBasePlane(TriEdge* edges, float* valuesToInterpolate, Plane* plane);

float DegreesToRadians(float degrees) {
	return degrees * M_PI / 180.0;
}

int GzRotXMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along x axis
// Pass back the matrix using mat value

	float radians = DegreesToRadians(degree);

	mat[0][0] = 1;
	mat[0][1] = 0;
	mat[0][2] = 0;
	mat[0][3] = 0;

	mat[1][0] = 0;
	mat[1][1] = cos(radians);
	mat[1][2] = -sin(radians);
	mat[1][3] = 0;

	mat[2][0] = 0;
	mat[2][1] = sin(radians);
	mat[2][2] = cos(radians);
	mat[2][3] = 0;

	mat[3][0] = 0;
	mat[3][1] = 0;
	mat[3][2] = 0;
	mat[3][3] = 1;

	return GZ_SUCCESS;
}


int GzRotYMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along y axis
// Pass back the matrix using mat value

	float radians = DegreesToRadians(degree);

	mat[0][0] = cos(radians);
	mat[0][1] = 0;
	mat[0][2] = sin(radians);
	mat[0][3] = 0;

	mat[1][0] = 0;
	mat[1][1] = 1;
	mat[1][2] = 0;
	mat[1][3] = 0;

	mat[2][0] = -sin(radians);
	mat[2][1] = 0;
	mat[2][2] = cos(radians);
	mat[2][3] = 0;

	mat[3][0] = 0;
	mat[3][1] = 0;
	mat[3][2] = 0;
	mat[3][3] = 1;

	return GZ_SUCCESS;
}


int GzRotZMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along z axis
// Pass back the matrix using mat value

	float radians = DegreesToRadians(degree);

	mat[0][0] = cos(radians);
	mat[0][1] = -sin(radians);
	mat[0][2] = 0;
	mat[0][3] = 0;

	mat[1][0] = sin(radians);
	mat[1][1] = cos(radians);
	mat[1][2] = 0;
	mat[1][3] = 0;

	mat[2][0] = 0;
	mat[2][1] = 0;
	mat[2][2] = 1;
	mat[2][3] = 0;

	mat[3][0] = 0;
	mat[3][1] = 0;
	mat[3][2] = 0;
	mat[3][3] = 1;

	return GZ_SUCCESS;
}


int GzTrxMat(GzCoord translate, GzMatrix mat)
{
// Create translation matrix
// Pass back the matrix using mat value

	mat[0][0] = 1;
	mat[0][1] = 0;
	mat[0][2] = 0;
	mat[0][3] = translate[X];

	mat[1][0] = 0;
	mat[1][1] = 1;
	mat[1][2] = 0;
	mat[1][3] = translate[Y];

	mat[2][0] = 0;
	mat[2][1] = 0;
	mat[2][2] = 1;
	mat[2][3] = translate[Z];

	mat[3][0] = 0;
	mat[3][1] = 0;
	mat[3][2] = 0;
	mat[3][3] = 1;

	return GZ_SUCCESS;
}


int GzScaleMat(GzCoord scale, GzMatrix mat)
{
// Create scaling matrix
// Pass back the matrix using mat value

	mat[0][0] = scale[X];
	mat[0][1] = 0;
	mat[0][2] = 0;
	mat[0][3] = 0;

	mat[1][0] = 0;
	mat[1][1] = scale[Y];
	mat[1][2] = 0;
	mat[1][3] = 0;

	mat[2][0] = 0;
	mat[2][1] = 0;
	mat[2][2] = scale[Z];
	mat[2][3] = 0;

	mat[3][0] = 0;
	mat[3][1] = 0;
	mat[3][2] = 0;
	mat[3][3] = 1;

	return GZ_SUCCESS;
}

int GetIdentityMatrix(GzMatrix mat) {
	mat[0][0] = 1;
	mat[0][1] = 0;
	mat[0][2] = 0;
	mat[0][3] = 0;

	mat[1][0] = 0;
	mat[1][1] = 1;
	mat[1][2] = 0;
	mat[1][3] = 0;

	mat[2][0] = 0;
	mat[2][1] = 0;
	mat[2][2] = 1;
	mat[2][3] = 0;

	mat[3][0] = 0;
	mat[3][1] = 0;
	mat[3][2] = 0;
	mat[3][3] = 1;

	return GZ_SUCCESS;
}

int GzMatrixMultiply(GzMatrix mat0, GzMatrix mat1, GzMatrix output) {

	for (int row = 0; row < 4; row++) {
		for (int col = 0; col < 4; col++) {
			output[row][col] = 0;
			for (int i = 0; i < 4; i++) {
				output[row][col] += mat0[row][i] * mat1[i][col];
			}
		}
	}

	return GZ_SUCCESS;
}

int ScaleAllMatrixTerms(GzMatrix mat, float scaleFactor) {

	for (int row = 0; row < 4; row++) {
		for (int col = 0; col < 4; col++) {
			mat[row][col] *= scaleFactor;
		}
	}

	return GZ_SUCCESS;
}

int VectorAdd(GzCoord a, GzCoord b, GzCoord diff) {

	diff[X] = a[X] + b[X];
	diff[Y] = a[Y] + b[Y];
	diff[Z] = a[Z] + b[Z];

	return GZ_SUCCESS;
}

int VectorSubtract(GzCoord a, GzCoord b, GzCoord diff) {

	diff[X] = a[X] - b[X];
	diff[Y] = a[Y] - b[Y];
	diff[Z] = a[Z] - b[Z];

	return GZ_SUCCESS;
}

int VectorNormalize(GzCoord vector) {

	float length = sqrt(pow(vector[X], 2) + pow(vector[Y], 2) + pow(vector[Z], 2));
	vector[X] /= length;
	vector[Y] /= length;
	vector[Z] /= length;

	return GZ_SUCCESS;
}

int VectorScale(GzCoord a, float term, GzCoord product) {

	product[X] = a[X] * term;
	product[Y] = a[Y] * term;
	product[Z] = a[Z] * term;

	return GZ_SUCCESS;
}

float VectorDotProduct(GzCoord a, GzCoord b) {
	return a[X] * b[X] + a[Y] * b[Y] + a[Z] * b[Z];
}

int VectorCrossProduct(GzCoord a, GzCoord b, GzCoord product) {

	product[X] = a[Y] * b[Z] - a[Z] * b[Y];
	product[Y] = a[Z] * b[X] - a[X] * b[Z];
	product[Z] = a[X] * b[Y] - a[Y] * b[X];

	return GZ_SUCCESS;
}

int VectorNegate(GzCoord a, GzCoord result) {

	result[X] = -a[X];
	result[Y] = -a[Y];
	result[Z] = -a[Z];

	return GZ_SUCCESS;
}

int CalculateXsp(GzRender *render, GzDisplay *display) {

	render->Xsp[0][0] = display->xres / 2;
	render->Xsp[0][1] = 0;
	render->Xsp[0][2] = 0;
	render->Xsp[0][3] = display->xres / 2;

	render->Xsp[1][0] = 0;
	render->Xsp[1][1] = -(display->yres / 2);
	render->Xsp[1][2] = 0;
	render->Xsp[1][3] = display->yres / 2;

	render->Xsp[2][0] = 0;
	render->Xsp[2][1] = 0;
	render->Xsp[2][2] = INT_MAX;
	render->Xsp[2][3] = 0;

	render->Xsp[3][0] = 0;
	render->Xsp[3][1] = 0;
	render->Xsp[3][2] = 0;
	render->Xsp[3][3] = 1;

	return GZ_SUCCESS;
}

int CalculateXpi(GzRender *render) {
	float fovInRadians = DegreesToRadians(render->camera.FOV);
	float projectionTerm = tan(fovInRadians / 2);
	
	render->camera.Xpi[0][0] = 1;
	render->camera.Xpi[0][1] = 0;
	render->camera.Xpi[0][2] = 0;
	render->camera.Xpi[0][3] = 0;

	render->camera.Xpi[1][0] = 0;
	render->camera.Xpi[1][1] = 1;
	render->camera.Xpi[1][2] = 0;
	render->camera.Xpi[1][3] = 0;

	render->camera.Xpi[2][0] = 0;
	render->camera.Xpi[2][1] = 0;
	render->camera.Xpi[2][2] = projectionTerm;
	render->camera.Xpi[2][3] = 0;

	render->camera.Xpi[3][0] = 0;
	render->camera.Xpi[3][1] = 0;
	render->camera.Xpi[3][2] = projectionTerm;
	render->camera.Xpi[3][3] = 1;

	return GZ_SUCCESS;
}

int CalculateXiw(GzRender *render) {

	GzCoord camZ;
	VectorSubtract(render->camera.lookat, render->camera.position, camZ);
	VectorNormalize(camZ);

	GzCoord camY, zAndUp;
	float zScale = VectorDotProduct(render->camera.worldup, camZ);
	VectorScale(camZ, zScale, zAndUp);
	VectorSubtract(render->camera.worldup, zAndUp, camY);
	VectorNormalize(camY);

	GzCoord camX;
	VectorCrossProduct(camY, camZ, camX);
	VectorNormalize(camX);

	render->camera.Xiw[0][0] = camX[X];
	render->camera.Xiw[0][1] = camX[Y];
	render->camera.Xiw[0][2] = camX[Z];
	render->camera.Xiw[0][3] = -VectorDotProduct(camX, render->camera.position);

	render->camera.Xiw[1][0] = camY[X];
	render->camera.Xiw[1][1] = camY[Y];
	render->camera.Xiw[1][2] = camY[Z];
	render->camera.Xiw[1][3] = -VectorDotProduct(camY, render->camera.position);

	render->camera.Xiw[2][0] = camZ[X];
	render->camera.Xiw[2][1] = camZ[Y];
	render->camera.Xiw[2][2] = camZ[Z];
	render->camera.Xiw[2][3] = -VectorDotProduct(camZ, render->camera.position);

	render->camera.Xiw[3][0] = 0;
	render->camera.Xiw[3][1] = 0;
	render->camera.Xiw[3][2] = 0;
	render->camera.Xiw[3][3] = 1;

	return GZ_SUCCESS;
}

int CalculateXforms(GzRender *render) {

	int status = GZ_SUCCESS;

	status |= CalculateXpi(render);
	status |= CalculateXiw(render);

	return status;
}

int BuildXformsStack(GzRender *render) {

	int status = GZ_SUCCESS;

	GzMatrix identity;
	GetIdentityMatrix(identity);

	status |= PushImageMatrix(render, render->Xsp);
	status |= PushNormMatrix(render, identity);
	status |= PushWorldMatrix(render, identity);
	render->matlevel++;

	status |= PushImageMatrix(render, render->camera.Xpi);
	status |= PushNormMatrix(render, identity);
	status |= PushWorldMatrix(render, identity);
	render->matlevel++;

	status |= PushImageMatrix(render, render->camera.Xiw);
	status |= PushNormMatrix(render, render->camera.Xiw);
	status |= PushWorldMatrix(render, identity);
	render->matlevel++;

	return status;
}

int MultiplyVector3ByMatrix4(GzCoord vertex, GzMatrix matrix, GzCoord result) {

	float vector4d[4], result4d[4];
	vector4d[X] = vertex[X];
	vector4d[Y] = vertex[Y];
	vector4d[Z] = vertex[Z];
	vector4d[3] = 1;

	for (int row = 0; row < 4; row++) {
		result4d[row] = 0;
		for (int i = 0; i < 4; i++) {
			result4d[row] += matrix[row][i] * vector4d[i];
		}
	}
	
	result[X] = result4d[X] / result4d[3];
	result[Y] = result4d[Y] / result4d[3];
	result[Z] = result4d[Z] / result4d[3];

	return GZ_SUCCESS;
}

//----------------------------------------------------------
// Begin main functions

int GzNewRender(GzRender **render, GzDisplay *display)
{
/* 
- malloc a renderer struct
- span interpolator needs pointer to display for pixel writes
*/
	*render = new GzRender();
	(*render)->display = display;

	(*render)->camera.position[X] = DEFAULT_IM_X;      
  	(*render)->camera.position[Y] = DEFAULT_IM_Y;
  	(*render)->camera.position[Z] = DEFAULT_IM_Z;

  	(*render)->camera.lookat[X] = 0.0;
  	(*render)->camera.lookat[Y] = 0.0;
  	(*render)->camera.lookat[Z] = 0.0;

  	(*render)->camera.worldup[X] = 0.0;
  	(*render)->camera.worldup[Y] = 1.0;
  	(*render)->camera.worldup[Z] = 0.0;

	(*render)->camera.FOV = DEFAULT_FOV;

	CalculateXsp(*render, display);

	(*render)->interp_mode = GZ_FLAT;

	(*render)->aaShiftX = 0;
	(*render)->aaShiftY = 0;

	return GZ_SUCCESS;
}


int GzFreeRender(GzRender *render)
{
/* 
-free all renderer resources
*/
	free(render);
	render->matlevel = 0;

	return GZ_SUCCESS;
}


int GzBeginRender(GzRender	*render)
{
/* 
- set up for start of each frame - init frame buffer
*/
	int status = GZ_SUCCESS;

	status |= GzInitDisplay(render->display);

	status |= CalculateXforms(render);

	status |= BuildXformsStack(render);

	return status;
}


int GzPutCamera(GzRender *render, GzCamera *camera)
{
/*
- overwrite renderer camera structure with new camera definition
*/
	memcpy(&(render->camera), camera, sizeof(GzCamera));

	return GZ_SUCCESS;	
}

int GzPushMatrix(GzRender *render, GzMatrix	matrix)
{
/*
- push a matrix onto the Ximage stack
- check for stack overflow
*/
	if (render->matlevel >= MATLEVELS) {
		return GZ_FAILURE;
	}

	PushImageMatrix(render, matrix);
	PushNormMatrix(render, matrix);
	PushWorldMatrix(render, matrix);
	render->matlevel++;

	return GZ_SUCCESS;
}

int GzPopMatrix(GzRender *render)
{
/*
- pop a matrix off the Ximage stack
- check for stack underflow
*/
	if (render->matlevel == 0) {
		return GZ_FAILURE;
	}
	render->matlevel--;

	return GZ_SUCCESS;
}

int PushImageMatrix(GzRender *render, GzMatrix matrix)
{
	if (render->matlevel == 0) {
		memcpy(render->Ximage[0], matrix, sizeof(GzMatrix));
	} else {
		GzMatrix top;
		GzMatrixMultiply(render->Ximage[render->matlevel - 1], matrix, top);
		memcpy(render->Ximage[render->matlevel], top, sizeof(GzMatrix));
	}

	return GZ_SUCCESS;
}

int PushNormMatrix(GzRender *render, GzMatrix matrix)
{
	//remove translation
	matrix[0][3] = 0;
	matrix[1][3] = 0;
	matrix[2][3] = 0;

	//ensure unitary rotation
	float scaleFactor = 1 / sqrt(matrix[0][0] * matrix[0][0] + matrix[0][1] * matrix[0][1] + matrix[0][2] * matrix[0][2]);
	ScaleAllMatrixTerms(matrix, scaleFactor);

	if (render->matlevel == 0) {
		memcpy(render->Xnorm[0], matrix, sizeof(GzMatrix));
	} else {
		GzMatrix top;
		GzMatrixMultiply(render->Xnorm[render->matlevel - 1], matrix, top);
		memcpy(render->Xnorm[render->matlevel], top, sizeof(GzMatrix));
	}

	return GZ_SUCCESS;
}

int PushWorldMatrix(GzRender *render, GzMatrix matrix)
{
	if (render->matlevel == 0) {
		memcpy(render->Xworld[0], matrix, sizeof(GzMatrix));
	} else {
		GzMatrix top;
		GzMatrixMultiply(render->Xworld[render->matlevel - 1], matrix, top);
		memcpy(render->Xworld[render->matlevel], top, sizeof(GzMatrix));
	}

	return GZ_SUCCESS;
}

int GzPutAttribute(GzRender	*render, int numAttributes, GzToken	*nameList, 
	GzPointer *valueList) /* void** valuelist */
{
/*
- set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
- later set shaders, interpolaters, texture maps, and lights
*/
	for (int i = 0; i < numAttributes; i++) {
		switch(nameList[i]) {
		case GZ_RGB_COLOR:
			GzColor* color;
			color = static_cast<GzColor*>(valueList[i]);
			render->flatcolor[RED] = (*color)[RED];
			render->flatcolor[GREEN] = (*color)[GREEN];
			render->flatcolor[BLUE] = (*color)[BLUE];
			break;
		case GZ_AMBIENT_LIGHT:
			render->ambientlight = *(static_cast<GzLight*>(valueList[i]));
			break;
		case GZ_AMBIENT_COEFFICIENT:
			GzColor* ka;
			ka = static_cast<GzColor*>(valueList[i]);
			render->Ka[RED] = (*ka)[RED];
			render->Ka[GREEN] = (*ka)[GREEN];
			render->Ka[BLUE] = (*ka)[BLUE];
			break;
		case GZ_DIFFUSE_COEFFICIENT:
			GzColor* kd;
			kd = static_cast<GzColor*>(valueList[i]);
			render->Kd[RED] = (*kd)[RED];
			render->Kd[GREEN] = (*kd)[GREEN];
			render->Kd[BLUE] = (*kd)[BLUE];
			break;
		case GZ_SPECULAR_COEFFICIENT:
			GzColor* ks;
			ks = static_cast<GzColor*>(valueList[i]);
			render->Ks[RED] = (*ks)[RED];
			render->Ks[GREEN] = (*ks)[GREEN];
			render->Ks[BLUE] = (*ks)[BLUE];
			break;
		case GZ_DISTRIBUTION_COEFFICIENT:
			render->spec = *((float*)valueList[i]);
			break;
		case GZ_DIRECTIONAL_LIGHT:
			render->lights[render->numlights] = *(static_cast<GzLight*>(valueList[i]));
			render->numlights++;
			break;
		case GZ_INTERPOLATE:
			render->interp_mode = *(static_cast<int*>(valueList[i]));
			break;
		case GZ_TEXTURE_MAP:
			render->tex_fun = static_cast<GzTexture>(valueList[i]);
			break;
		case GZ_AASHIFTX:
			render->aaShiftX = *((float*)valueList[i]);
			break;
		case GZ_AASHIFTY:
			render->aaShiftY = *((float*)valueList[i]);
			break;
		case GZ_REFLECTIVE:
			render->isReflective = valueList[i];
			break;
		case GZ_CUBE_MAP:
			render->cubetex_fun = static_cast<GzCubeMap>(valueList[i]);
			break;
		}
	}

	return GZ_SUCCESS;
}

/* convert float color to GzIntensity short */
short	ctoi(float color)
{
  return(short)((int)(color * ((1 << 12) - 1)));
}

int sortVerticiesByX(const void* a, const void* b) {
	return ((float*)a)[X] - ((float*)b)[X];
}

int sortVerticiesByY(const void* a, const void* b) {
	TriVertex* vertA = (TriVertex*)(a);
	TriVertex* vertB = (TriVertex*)(b);
	return vertA->screen[Y] - vertB->screen[Y];
}

void getSortedEdges(GzRender* render, GzCoord* screenVerticies, GzCoord* modelVerticies, GzCoord* modelNormals, GzTextureIndex* uvList, TriEdge* edges) {

	GzCoord imageVerticies[3], imageNormals[3], worldVerticies[3], worldNormals[3];
	for (int i = 0; i < 3; i++) {
		MultiplyVector3ByMatrix4(modelVerticies[i], render->Xnorm[render->matlevel - 1], imageVerticies[i]);
		MultiplyVector3ByMatrix4(modelNormals[i], render->Xnorm[render->matlevel - 1], imageNormals[i]);
		MultiplyVector3ByMatrix4(modelVerticies[i], render->Xworld[render->matlevel - 1], worldVerticies[i]);
		MultiplyVector3ByMatrix4(modelNormals[i], render->Xworld[render->matlevel - 1], worldNormals[i]);
	}

	//set various coords (screen, model, world, image, uv) for each vertex
	TriVertex triVerticies[3];
	for (int i = 0; i < 3; i++) {
		memcpy(triVerticies[i].screen, screenVerticies[i], sizeof(GzCoord));
		memcpy(triVerticies[i].modelPos, modelVerticies[i], sizeof(GzCoord));
		memcpy(triVerticies[i].modelNormal, modelNormals[i], sizeof(GzCoord));
		memcpy(triVerticies[i].worldPos, worldVerticies[i], sizeof(GzCoord));
		memcpy(triVerticies[i].worldNormal, worldNormals[i], sizeof(GzCoord));
		memcpy(triVerticies[i].image, imageVerticies[i], sizeof(GzCoord));
		memcpy(triVerticies[i].normal, imageNormals[i], sizeof(GzCoord));
		memcpy(triVerticies[i].uv, uvList[i], sizeof(GzTextureIndex));
	}

	std::qsort(triVerticies, 3, sizeof(TriVertex), sortVerticiesByY);

	//set which edges should be colored
	float slopeOfOppositeEdge = (triVerticies[2].screen[Y] - triVerticies[0].screen[Y]) / (triVerticies[2].screen[X] - triVerticies[0].screen[X]);
	float midVertexX = ((triVerticies[1].screen[Y] - triVerticies[0].screen[Y]) / slopeOfOppositeEdge) + triVerticies[0].screen[X];
	
	if (midVertexX > triVerticies[1].screen[X]) {
		edges[0].start = triVerticies[0];
		edges[0].end = triVerticies[1];
		edges[0].isColored = true;

		edges[1].start = triVerticies[1];
		edges[1].end = triVerticies[2];
		edges[1].isColored = true;

		edges[2].start = triVerticies[2];
		edges[2].end = triVerticies[0];
		edges[2].isColored = false;
	} else {
		edges[0].start = triVerticies[0];
		edges[0].end = triVerticies[2];
		edges[0].isColored = true;

		edges[1].start = triVerticies[2];
		edges[1].end = triVerticies[1];
		edges[1].isColored = false;

		edges[2].start = triVerticies[1];
		edges[2].end = triVerticies[0];
		edges[2].isColored = false;
	}
	
	//compute A, B, C terms for LEE
	for (int i = 0; i < 3; i++) {
		float dX = edges[i].end.screen[X] - edges[i].start.screen[X];
		float dY = edges[i].end.screen[Y] - edges[i].start.screen[Y];
		float startX = edges[i].start.screen[X];
		float startY = edges[i].start.screen[Y];
		edges[i].A = dY;
		edges[i].B = -dX;
		edges[i].C = (dX * startY) - (dY * startX);
	}
}

BoundingBox* getBoundingBox(GzCoord* verticies) {
	BoundingBox* bb = new BoundingBox();
	
	GzCoord verticiesCopy[3];
	for(int i = 0; i < 3; i++) {
		memcpy(verticiesCopy[i], verticies[i], sizeof(GzCoord));
	}

	qsort(verticiesCopy, 3, sizeof(GzCoord), sortVerticiesByY);
	bb->startY = floor(verticiesCopy[0][Y]);
	bb->endY = ceil(verticiesCopy[2][Y]);

	qsort(verticiesCopy, 3, sizeof(GzCoord), sortVerticiesByX);
	bb->startX = floor(verticiesCopy[0][X]);
	bb->endX = ceil(verticiesCopy[2][X]);

	return bb;
}

float getInterpolatedZ(TriEdge* edges, int x, int y) {
	
	float x0 = edges[0].end.screen[X] - edges[0].start.screen[X];
	float y0 = edges[0].end.screen[Y] - edges[0].start.screen[Y];
	float z0 = edges[0].end.screen[Z] - edges[0].start.screen[Z];
	float x1 = edges[1].end.screen[X] - edges[1].start.screen[X];
	float y1 = edges[1].end.screen[Y] - edges[1].start.screen[Y];
	float z1 = edges[1].end.screen[Z] - edges[1].start.screen[Z];

	//cross product
	float A = (y0 * z1) - (z0 * y1);
	float B = (z0 * x1) - (x0 * z1);
	float C = (x0 * y1) - (y0 * x1);
	float D = -(A * edges[0].start.screen[X] + B * edges[0].start.screen[Y] + C * edges[0].start.screen[Z]);

	return (A * x + B * y + D) / (-C);
}

int GetColorPlane(TriEdge* edges, int colorChannel, Plane* plane) {

	float valuesToInter[4] = {
		edges[0].start.color[colorChannel],
		edges[0].end.color[colorChannel],
		edges[1].start.color[colorChannel],
		edges[1].end.color[colorChannel]
	};
	return GetBasePlane(edges, valuesToInter, plane);
}

int GetNormalPlane(TriEdge* edges, int vectorComponent, Plane* plane) {
	
	float valuesToInter[4] = {
		edges[0].start.normal[vectorComponent],
		edges[0].end.normal[vectorComponent],
		edges[1].start.normal[vectorComponent],
		edges[1].end.normal[vectorComponent]
	};
	return GetBasePlane(edges, valuesToInter, plane);
}

int GetWorldNormalPlane(TriEdge* edges, int vectorComponent, Plane* plane) {

	float valuesToInter[4] = {
		edges[0].start.worldNormal[vectorComponent],
		edges[0].end.worldNormal[vectorComponent],
		edges[1].start.worldNormal[vectorComponent],
		edges[1].end.worldNormal[vectorComponent]
	};
	return GetBasePlane(edges, valuesToInter, plane);
}

int GetWorldPosPlane(TriEdge* edges, int vectorComponent, Plane* plane) {
	float valuesToInter[4] = {
		edges[0].start.worldPos[vectorComponent],
		edges[0].end.worldPos[vectorComponent],
		edges[1].start.worldPos[vectorComponent],
		edges[1].end.worldPos[vectorComponent]
	};
	return GetBasePlane(edges, valuesToInter, plane);
}

int GetUVPlane(TriEdge* edges, int vectorComponent, Plane* plane) {

	float zValues[4];
	zValues[0] = getInterpolatedZ(edges, edges[0].start.screen[X], edges[0].start.screen[Y]);
	zValues[1] = getInterpolatedZ(edges, edges[0].end.screen[X], edges[0].end.screen[Y]);
	zValues[2] = getInterpolatedZ(edges, edges[1].start.screen[X], edges[1].start.screen[Y]);
	zValues[3] = getInterpolatedZ(edges, edges[1].end.screen[X], edges[1].end.screen[Y]);

	float valuesToInter[4] = {
		(edges[0].start.uv[vectorComponent] / GetPerspectiveFactor(zValues[0])),
		(edges[0].end.uv[vectorComponent] / GetPerspectiveFactor(zValues[1])),
		(edges[1].start.uv[vectorComponent] / GetPerspectiveFactor(zValues[2])),
		(edges[1].end.uv[vectorComponent] / GetPerspectiveFactor(zValues[3]))
	};
	return GetBasePlane(edges, valuesToInter, plane);
}

//valuesToInterpolate = [valueAtEdge0Start, valueAtEdge0End, valueAtEdge1Start, valueAtEdge1End]  
int GetBasePlane(TriEdge* edges, float* valuesToInterpolate, Plane* plane) {
	GzCoord vectA, vectB, result;

	vectA[X] = edges[0].end.screen[X] - edges[0].start.screen[X];
	vectA[Y] = edges[0].end.screen[Y] - edges[0].start.screen[Y];
	vectA[Z] = valuesToInterpolate[1] - valuesToInterpolate[0];
	vectB[X] = edges[1].end.screen[X] - edges[1].start.screen[X];
	vectB[Y] = edges[1].end.screen[Y] - edges[1].start.screen[Y];
	vectB[Z] = valuesToInterpolate[3] - valuesToInterpolate[2];

	VectorCrossProduct(vectA, vectB, result);

	plane->A = result[X];
	plane->B = result[Y];
	plane->C = result[Z];
	plane->D = -(plane->A * edges[0].end.screen[X] + plane->B * edges[0].end.screen[Y] + plane->C * valuesToInterpolate[1]);
	
	return GZ_SUCCESS;
}

float SolveForPlaneZ(Plane plane, int x, int y) {
	return (plane.A * x + plane.B * y + plane.D) / (-plane.C);
}

bool IsTriangleHidden(GzRender* render, GzCoord* screenSpaceVerticies) {

	bool hidden = false;

	hidden |= (screenSpaceVerticies[0][X] < 0 && screenSpaceVerticies[1][X] < 0 && screenSpaceVerticies[2][X] < 0);
	hidden |= (screenSpaceVerticies[0][X] >= render->display->xres && screenSpaceVerticies[1][X] >= render->display->xres && screenSpaceVerticies[2][X] >= render->display->xres);

	hidden |= (screenSpaceVerticies[0][Y] < 0 && screenSpaceVerticies[1][Y] < 0 && screenSpaceVerticies[2][Y] < 0);
	hidden |= (screenSpaceVerticies[0][Y] >= render->display->yres && screenSpaceVerticies[1][Y] >= render->display->yres && screenSpaceVerticies[2][Y] >= render->display->yres);

	return hidden;
}

int GzPutTriangle(GzRender *render, int	numParts, GzToken *nameList,
	GzPointer *valueList)
{

	GzCoord* modelSpaceVerticies = (GzCoord*)malloc(sizeof (GzCoord));
	GzCoord* modelSpaceNormals = (GzCoord*)malloc(sizeof (GzCoord));
	GzTextureIndex* uvList = (GzTextureIndex*)malloc(sizeof (GzTextureIndex));

	for (int i = 0; i < numParts; i++) {
		switch(nameList[i]) {
		case GZ_NULL_TOKEN:
			break;
		case GZ_POSITION:
			modelSpaceVerticies = static_cast<GzCoord*>(valueList[i]);
			break;
		case GZ_NORMAL:
			modelSpaceNormals = static_cast<GzCoord*>(valueList[i]);
			break;
		case GZ_TEXTURE_INDEX:
			uvList = static_cast<GzTextureIndex*>(valueList[i]);
			break;
		}
	}

	RenderTriangle(render, modelSpaceVerticies, modelSpaceNormals, uvList);

	return GZ_SUCCESS;
}

int GetReflectionAcrossNormal(GzRender* render, GzCoord worldPos, GzCoord worldNorm, GzCoord* reflectionOutput) {

	VectorNormalize(worldNorm);

	GzCoord incident;
	VectorSubtract(worldPos, render->camera.position, incident);
	VectorNormalize(incident);

	float nDotE = VectorDotProduct(worldNorm, incident);

	GzCoord reflection;
	VectorScale(worldNorm, 2 * nDotE, reflection);
	VectorSubtract(incident, reflection, reflection);
	VectorNormalize(reflection);

	memcpy(reflectionOutput, reflection, sizeof(GzCoord));

	return GZ_SUCCESS;
}

int GetColorAtNormal(GzRender* render, GzCoord normal, GzColor textureColor, GzColor* normalColor) {
	
	GzColor color = {0, 0, 0};
	
	GzCoord eyeVector;
	eyeVector[X] = 0;
	eyeVector[Y] = 0;
	eyeVector[Z] = -1;

	GzCoord diffuseSum;
	diffuseSum[RED] = diffuseSum[GREEN] = diffuseSum[BLUE] = 0;
	GzCoord specularSum;
	specularSum[RED] = specularSum[GREEN] = specularSum[BLUE] = 0;

	VectorNormalize(normal);

	for (int lightIndex = 0; lightIndex < render->numlights; lightIndex++) {
		float nDotL = VectorDotProduct(normal, render->lights[lightIndex].direction);
		float nDotE = VectorDotProduct(normal, eyeVector);

		GzCoord vertexNormal;

		if (nDotL > 0 && nDotE > 0) {
			memcpy(vertexNormal, normal, sizeof(GzCoord));
		} else if (nDotL < 0 && nDotE < 0) {
			VectorScale(normal, -1, vertexNormal);
			nDotL *= -1;
			nDotE *= -1;
		} else {
			//light and eye are on different sides, skip
			continue;
		}

		GzCoord reflection;
		VectorScale(normal, 2 * nDotL, reflection);
		VectorSubtract(reflection, render->lights[lightIndex].direction, reflection);
		VectorNormalize(reflection);

		float rDotE = VectorDotProduct(reflection, eyeVector);
		if (rDotE < 0) {
			rDotE = 0;
		}

		GzCoord specularVector;
		VectorScale(render->lights[lightIndex].color, pow(rDotE, render->spec), specularVector);
		VectorAdd(specularSum, specularVector, specularSum);

		GzCoord diffuseVector;
		VectorScale(render->lights[lightIndex].color, nDotL, diffuseVector);
		VectorAdd(diffuseSum, diffuseVector, diffuseSum);
	}

	GzColor ambient, diffuse, specular;

	if (textureColor) {
		//textured Phong
		MultiplyVectorByCoefficient(render->ambientlight.color, textureColor, ambient);
		MultiplyVectorByCoefficient(diffuseSum, textureColor, diffuse);
		MultiplyVectorByCoefficient(specularSum, render->Ks, specular);
	} else if (render->tex_fun == 0) {
		//no texture
		MultiplyVectorByCoefficient(render->ambientlight.color, render->Ka, ambient);
		MultiplyVectorByCoefficient(diffuseSum, render->Kd, diffuse);
		MultiplyVectorByCoefficient(specularSum, render->Ks, specular);
	} else {
		//textured Gouraud
		memcpy(ambient, render->ambientlight.color, sizeof(GzColor));
		memcpy(diffuse, diffuseSum, sizeof(GzColor));
		memcpy(specular, specularSum, sizeof(GzColor));
	}

	VectorAdd(color, ambient, color);
	VectorAdd(color, diffuse, color);
	VectorAdd(color, specular, color);

	for (int i = 0; i < 3; i++) {
		if (color[i] > 1) {
			color[i] = 1;
		}
	}

	memcpy(normalColor, color, sizeof(GzColor));

	return GZ_SUCCESS;
}

int RenderTriangle(GzRender *render, GzCoord* modelSpaceVerticies, GzCoord* modelSpaceNormals, GzTextureIndex* uvList) {

	GzCoord* screenSpaceVerticies = (GzCoord*)malloc(3 * sizeof(GzCoord));

	bool skipTriangle = false;
	for (int i = 0; i < 3; i++) {
		MultiplyVector3ByMatrix4(modelSpaceVerticies[i], render->Ximage[render->matlevel - 1], screenSpaceVerticies[i]);
		
		//handle AA offset
		screenSpaceVerticies[i][X] -= render->aaShiftX;
		screenSpaceVerticies[i][Y] -= render->aaShiftY;
		
		if (screenSpaceVerticies[i][Z] < 0) {
			skipTriangle = true;
			break;
		}
	}

	skipTriangle |= IsTriangleHidden(render, screenSpaceVerticies);

	if (skipTriangle) {
		return GZ_SUCCESS;
	}
		
	TriEdge edges[3];
	getSortedEdges(render, screenSpaceVerticies, modelSpaceVerticies, modelSpaceNormals, uvList, edges);
	BoundingBox* box = getBoundingBox(screenSpaceVerticies);

	Plane colorPlanes[3];
	Plane normalPlanes[3];
	Plane uvPlanes[2];
	switch(render->interp_mode) {
	case GZ_COLOR:

		for (int i = 0; i < 3; i++) {
			GzColor vertexColor;
			GetColorAtNormal(render, edges[i].end.normal, 0, &vertexColor);
			memcpy(edges[i].end.color, vertexColor, sizeof(GzColor));
		}
		memcpy(edges[1].start.color, edges[0].end.color, sizeof(GzColor));
		memcpy(edges[2].start.color, edges[1].end.color, sizeof(GzColor));
		memcpy(edges[0].start.color, edges[2].end.color, sizeof(GzColor));

		GetColorPlane(edges, RED, &colorPlanes[RED]);
		GetColorPlane(edges, GREEN, &colorPlanes[GREEN]);
		GetColorPlane(edges, BLUE, &colorPlanes[BLUE]);
		break;
	case GZ_NORMALS:
		GetNormalPlane(edges, X, &normalPlanes[X]);
		GetNormalPlane(edges, Y, &normalPlanes[Y]);
		GetNormalPlane(edges, Z, &normalPlanes[Z]);
		break;
	}

	//world planes used for environment mapping
	Plane worldPosPlanes[3];
	Plane worldNormalPlanes[3];
	for (int coordIndex = 0; coordIndex < 3; coordIndex++) {
		GetWorldPosPlane(edges, coordIndex, &worldPosPlanes[coordIndex]);
		GetWorldNormalPlane(edges, coordIndex, &worldNormalPlanes[coordIndex]);
	}

	if (uvList) {
		GetUVPlane(edges, U, &uvPlanes[U]);
		GetUVPlane(edges, V, &uvPlanes[V]);
	}

	for (int y = box->startY; y <= box->endY; y++) {
		for (int x = box->startX; x <= box->endX; x++) {

			GzIntensity r, g, b, a;
			GzDepth bufferZ;
					
			if (GzGetDisplay(render->display, x, y, &r, &g, &b, &a, &bufferZ)) {
				//pixel is outside viewing area
				continue;
			}

			bool colorPixel = true;
			for (int edgeIndex = 0; edgeIndex < 3; edgeIndex++) {
				TriEdge edge = edges[edgeIndex];
				float sign = edge.A * x + edge.B * y + edge.C;

				//pixel is on edge
				if (abs(sign) < 0.000001) {
					colorPixel = edge.isColored;
					break;
				}

				//pixel is outside triangle
				if (sign < 0) {
					colorPixel = false;
					break;
				}
			}

			if (colorPixel) {
				float newZ = getInterpolatedZ(edges, x, y);
				
				if (newZ < bufferZ) {
					GzColor color, textureColor;

					GzCoord interUV;
					switch (render->interp_mode) {
					case GZ_COLOR:

						//base color
						for (int colorChannel = 0; colorChannel < 3; colorChannel++) {
							color[colorChannel] = SolveForPlaneZ(colorPlanes[colorChannel], x, y);
						}

						if (render->tex_fun) {
							//texture color
							interUV[U] = SolveForPlaneZ(uvPlanes[U], x, y) * GetPerspectiveFactor(newZ);
							interUV[V] = SolveForPlaneZ(uvPlanes[V], x, y) * GetPerspectiveFactor(newZ);
							render->tex_fun(interUV[U], interUV[V], textureColor);
							MultiplyVectorByCoefficient(color, textureColor, color);
						}

						break;
					case GZ_NORMALS:

						//base color
						GzCoord interNormal, interWorldPos, interWorldNorm;
						for (int vertexComp = 0; vertexComp < 3; vertexComp++) {
							interNormal[vertexComp] = SolveForPlaneZ(normalPlanes[vertexComp], x, y);
							interWorldPos[vertexComp] = SolveForPlaneZ(worldPosPlanes[vertexComp], x, y);
							interWorldNorm[vertexComp] = SolveForPlaneZ(worldNormalPlanes[vertexComp], x, y);
						}

						if (render->isReflective) {
							//compute reflection into environment map
							GzCoord reflection;
							GetReflectionAcrossNormal(render, interWorldPos, interWorldNorm, &reflection);

							render->cubetex_fun(reflection, textureColor);
							GetColorAtNormal(render, interNormal, textureColor, &color);
						} else if (render->tex_fun) {
							//texture color
							interUV[U] = SolveForPlaneZ(uvPlanes[U], x, y) * GetPerspectiveFactor(newZ);
							interUV[V] = SolveForPlaneZ(uvPlanes[V], x, y) * GetPerspectiveFactor(newZ);
							render->tex_fun(interUV[U], interUV[V], textureColor);

							GetColorAtNormal(render, interNormal, textureColor, &color);
						} else {
							GetColorAtNormal(render, interNormal, 0, &color);
						}

						break;
					default:
						color[RED] = render->flatcolor[RED];
						color[GREEN] = render->flatcolor[GREEN];
						color[BLUE] = render->flatcolor[BLUE];
						break;
					}

					GzPutDisplay(render->display, x, y, ctoi(color[RED]), ctoi(color[GREEN]), ctoi(color[BLUE]), 1, newZ);
				}
			}
		}
	}

	return GZ_SUCCESS;
}

int MultiplyVectorByCoefficient(GzCoord vec, GzCoord coef, GzCoord product) {

	for (int i = 0; i < 3; i++) {
		product[i] = vec[i] * coef[i];
	}

	return GZ_SUCCESS;
}

float GetPerspectiveFactor(float screenSpaceZ) {
	return (screenSpaceZ / (INT_MAX - screenSpaceZ)) + 1;
}
