// ApplicationFinal.cpp: implementation of the Application class.
//
//////////////////////////////////////////////////////////////////////

/*
* application code for final CS580 Project
*/

#include "stdafx.h"
#include "CS580HW.h"
#include "ApplicationFinal.h"
#include "Gz.h"
#include "disp.h"
#include "rend.h"
#include <functional>

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

#define OUTFILE "output.ppm"
#define EM_OUTFILE	"em_output.ppm"

#define MODEL_COUNT	2
const char* ModelFiles[MODEL_COUNT] = {
	"ppot.asc",
	"plane.asc"
};
bool IsModelReflective[MODEL_COUNT] = {
	false,
	false
};
bool IsModelRefractive[MODEL_COUNT] = { //NOTE in order to make a model transparant, both refractive and reflective 
										//indices must be true for it. Therefore, if you set a material to be transparant, its reflective index, will be changed to true, automatically later.
	true,
	false
};

#define TEX_NONE 0
#define TEX_FILE 1
#define TEX_PROC 2
const int ModelTextures[MODEL_COUNT] = {
	TEX_NONE,
	TEX_PROC
};

const float AAFilter[AAKERNEL_SIZE][3] = // X-shift, Y-shift, weight
{
	-0.52, 0.38, 0.128, 		0.41, 0.56, 0.119,		0.27, 0.08, 0.294,
	-0.17, -0.29, 0.249,		0.58, -0.55, 0.104,		-0.31, -0.71, 0.106
};

extern int tex_fun(float u, float v, GzColor color); /* image texture function */
extern int ptex_fun(float u, float v, GzColor color); /* procedural texture function */
extern int cubetex_fun(GzCoord reflection, GzColor color); /* environment map texture function */

void shade(GzCoord norm, GzCoord color);

void FinishAA(GzDisplay** aaDisplays, GzDisplay* outputDisplay);

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ApplicationFinal::ApplicationFinal()
{

}

ApplicationFinal::~ApplicationFinal()
{
	Clean();
}

int ApplicationFinal::Initialize()
{
	GzCamera	camera;  
	int		    xRes, yRes;	/* display parameters */ 

	GzToken		nameListShader[9]; 	    /* shader attribute names */
	GzPointer   valueListShader[9];		/* shader attribute pointers */
	GzToken     nameListLights[10];		/* light info */
	GzPointer   valueListLights[10];
	GzToken     nameListAA[2];
	GzPointer   valueListAA[2];
	int			shaderType, interpStyle;
	float		specpower;
	int		status; 

	status = 0;

	/* 
	* Allocate memory for user input
	*/
	m_pUserInput = new GzInput;

	/* 
	* initialize the display and the renderer 
	*/ 
	m_nWidth = 512;		// frame buffer and display width
	m_nHeight = 512;    // frame buffer and display height

	status |= GzNewFrameBuffer(&m_pFrameBuffer, m_nWidth, m_nHeight);

	//setup final display
	status |= GzNewDisplay(&m_pDisplay, m_nWidth, m_nHeight);
	status |= GzGetDisplayParams(m_pDisplay, &xRes, &yRes);
	status |= GzInitDisplay(m_pDisplay); 

	//malloc AA structs
	m_pAADisplays = (GzDisplay**) malloc(AAKERNEL_SIZE * sizeof (GzDisplay*));
	m_pAARenders = (GzRender**) malloc(AAKERNEL_SIZE * sizeof (GzRender*));

	for (int aaPass = 0; aaPass < AAKERNEL_SIZE; aaPass++) {

		status |= GzNewDisplay(&(m_pAADisplays[aaPass]), m_nWidth, m_nHeight);
		status |= GzGetDisplayParams(m_pAADisplays[aaPass], &xRes, &yRes);

		status |= GzNewRender(&(m_pAARenders[aaPass]), m_pAADisplays[aaPass]);

#if 0 	/* set up app-defined camera if desired, else use camera defaults */
		camera.position[X] = 0;
		camera.position[Y] = 0;
		camera.position[Z] = -6;

		camera.lookat[X] = 0;
		camera.lookat[Y] = 0;
		camera.lookat[Z] = 0;

		camera.worldup[X] = 0.0;
		camera.worldup[Y] = 1.0;
		camera.worldup[Z] = 0.0;

		camera.FOV = 63.7;              /* degrees */

		status |= GzPutCamera(m_pAARenders[aaPass], &camera);
#endif 
		camera.position[X] = 0;
		camera.position[Y] = 5;
		camera.position[Z] = -7;

		camera.lookat[X] = 0;
		camera.lookat[Y] = 0;
		camera.lookat[Z] = 0;

		camera.worldup[X] = 0.0;
		camera.worldup[Y] = 1.0;
		camera.worldup[Z] = 0.0;

		camera.FOV = 63.7;              /* degrees */

		status |= GzPutCamera(m_pAARenders[aaPass], &camera); 

		/* Start Renderer */
		status |= GzBeginRender(m_pAARenders[aaPass]);

		/* Light */
		GzLight	light1 = { {-0.7071, 0.7071, 0}, {0.5, 0.5, 0.9} };
		GzLight	light2 = { {0, -0.7071, -0.7071}, {0.9, 0.2, 0.3} };
		GzLight	light3 = { {0.7071, 0.0, -0.7071}, {0.2, 0.7, 0.3} };
		GzLight	ambientlight = { {0, 0, 0}, {0.3, 0.3, 0.3} };

		/* Material property */
		GzColor specularCoefficient = { 0.3, 0.3, 0.3 };
		GzColor ambientCoefficient = { 0.1, 0.1, 0.1 };
		GzColor diffuseCoefficient = {0.7, 0.7, 0.7};

		/* 
		renderer is ready for frame --- define lights and shader at start of frame 
		*/

		/*
		* Tokens associated with light parameters
		*/
		nameListLights[0] = GZ_DIRECTIONAL_LIGHT;
		valueListLights[0] = (GzPointer)&light1;
		nameListLights[1] = GZ_DIRECTIONAL_LIGHT;
		valueListLights[1] = (GzPointer)&light2;
		nameListLights[2] = GZ_DIRECTIONAL_LIGHT;
		valueListLights[2] = (GzPointer)&light3;
		status |= GzPutAttribute(m_pAARenders[aaPass], 3, nameListLights, valueListLights);

		nameListLights[0] = GZ_AMBIENT_LIGHT;
		valueListLights[0] = (GzPointer)&ambientlight;
		status |= GzPutAttribute(m_pAARenders[aaPass], 1, nameListLights, valueListLights);

		/*
		* Tokens associated with shading 
		*/
		nameListShader[0]  = GZ_DIFFUSE_COEFFICIENT;
		valueListShader[0] = (GzPointer)diffuseCoefficient;

		/* 
		* Select either GZ_COLOR or GZ_NORMALS as interpolation mode  
		*/
		nameListShader[1]  = GZ_INTERPOLATE;
		interpStyle = GZ_NORMALS;         /* Phong shading */
		valueListShader[1] = (GzPointer)&interpStyle;

		nameListShader[2]  = GZ_AMBIENT_COEFFICIENT;
		valueListShader[2] = (GzPointer)ambientCoefficient;
		nameListShader[3]  = GZ_SPECULAR_COEFFICIENT;
		valueListShader[3] = (GzPointer)specularCoefficient;
		nameListShader[4]  = GZ_DISTRIBUTION_COEFFICIENT;
		specpower = 32;
		valueListShader[4] = (GzPointer)&specpower;

		status |= GzPutAttribute(m_pAARenders[aaPass], 5, nameListShader, valueListShader);

		nameListAA[0]  = GZ_AASHIFTX;
		valueListAA[0] = (GzPointer)&(AAFilter[aaPass][X]);
		nameListAA[1]  = GZ_AASHIFTY;
		valueListAA[1] = (GzPointer)&(AAFilter[aaPass][Y]);

		status |= GzPutAttribute(m_pAARenders[aaPass], 2, nameListAA, valueListAA);

		GzMatrix	scale = 
		{ 
			1.0,	0.0,	0.0,	0.0, 
			0.0,	1.0,	0.0,	0.0, 
			0.0,	0.0,	1.0,	0.0, 
			0.0,	0.0,	0.0,	1.0 
		};

		GzMatrix rotateY, rotateX;
		GzRotYMat(0, rotateY);
		GzRotXMat(0, rotateX);

		status |= GzPushMatrix(m_pAARenders[aaPass], scale);  
		status |= GzPushMatrix(m_pAARenders[aaPass], rotateY); 
		status |= GzPushMatrix(m_pAARenders[aaPass], rotateX); 
	}

	if (status) exit(GZ_FAILURE); 

	if (status) 
		return(GZ_FAILURE); 
	else 
		return(GZ_SUCCESS); 
}

int ApplicationFinal::Render() 
{
	GzToken		nameListTriangle[3]; 	/* vertex attribute names */
	GzPointer	valueListTriangle[3]; 	/* vertex attribute pointers */
	GzCoord		vertexList[3];	/* vertex position coordinates */ 
	GzCoord		normalList[3];	/* vertex normals */ 
	GzTextureIndex  	uvList[3];		/* vertex texture map indices */ 
	GzToken		nameListShader[9]; 	    /* shader attribute names */
	GzPointer   valueListShader[9];
	char		dummy[256]; 
	int			status; 

	//clean before re-rendering
	for (int aaPass = 0; aaPass < AAKERNEL_SIZE; aaPass++) {
		GzInitDisplay(m_pAADisplays[aaPass]);
	}
	GzInitDisplay(m_pDisplay);

	GzDisplay* emDisplay;
	int		    xRes, yRes;	/* display parameters */ 
	int emWidth = m_nWidth * 4;
	int emHeight = m_nHeight * 3;
	status |= GzNewDisplay(&emDisplay, emWidth, emHeight);
	status |= GzGetDisplayParams(emDisplay, &xRes, &yRes);
	status |= GzInitDisplay(emDisplay);

	//TODO: dummy array will be replaced with real environment map

#if 0
	float dummyEmColors[6][7] = {
		m_nWidth, m_nWidth * 2, 0, m_nWidth, 4000, 2000, 0,
		0, m_nWidth, m_nWidth, m_nWidth * 2, 2000, 2000, 0,
		m_nWidth, m_nWidth * 2, m_nWidth, m_nWidth * 2, 2000, 2000, 2000,
		m_nWidth, m_nWidth * 2, m_nWidth * 2, m_nWidth * 3, 0, 2000, 0,
		m_nWidth * 2, m_nWidth * 3, m_nWidth, m_nWidth * 2, 2000, 4000, 2000,
		m_nWidth * 3, m_nWidth * 4, m_nWidth, m_nWidth * 2, 4000, 0, 4000
	};
	for (int faceIndex = 0; faceIndex < 6; faceIndex++) {
		for (int x = dummyEmColors[faceIndex][0]; x < dummyEmColors[faceIndex][1]; x++) {
			for (int y = dummyEmColors[faceIndex][2]; y < dummyEmColors[faceIndex][3]; y++) {
				GzPutDisplay(emDisplay, x, y, 
					dummyEmColors[faceIndex][4], 
					dummyEmColors[faceIndex][5], 
					dummyEmColors[faceIndex][6], 1, INT_MAX);
			}
		}
	}

	FILE *em_outfile;
	if( (em_outfile  = fopen( EM_OUTFILE , "wb" )) == NULL )
	{
		AfxMessageBox( "The environment map output file was not opened\n" );
		return GZ_FAILURE;
	}
	GzFlushDisplay2File(em_outfile, emDisplay);
#endif
	/* 
	* Tokens associated with triangle vertex values 
	*/ 
	nameListTriangle[0] = GZ_POSITION; 
	nameListTriangle[1] = GZ_NORMAL; 
	nameListTriangle[2] = GZ_TEXTURE_INDEX;  

	for (int fileIndex = 0; fileIndex < MODEL_COUNT; fileIndex++) {
		// I/O File open
		FILE *infile;
		if( (infile  = fopen( ModelFiles[fileIndex], "r" )) == NULL )
		{
			AfxMessageBox( "The input file was not opened\n" );
			return GZ_FAILURE;
		}

		//set texture per model
		float refractionIndex = 1.2;

		nameListShader[0]  = GZ_REFLECTIVE;
		nameListShader[2]  = GZ_REFRACTION_INDEX;
		valueListShader[2] = (GzPointer)&refractionIndex;
		nameListShader[3]  = GZ_REFRACTIVE;

		if (IsModelReflective[fileIndex]) {
			valueListShader[0] = (GzPointer)true;
			nameListShader[1]  = GZ_CUBE_MAP;
			valueListShader[1] = (GzPointer)(cubetex_fun);
			valueListShader[3] = (GzPointer)true;
		} else if (ModelTextures[fileIndex] == TEX_FILE) {
			valueListShader[0] = (GzPointer)false;
			nameListShader[1]  = GZ_TEXTURE_MAP;
			valueListShader[1] = (GzPointer)(tex_fun);
			valueListShader[3] = (GzPointer)false;
		} else if (ModelTextures[fileIndex] == TEX_PROC) {
			valueListShader[0] = (GzPointer)false;
			nameListShader[1]  = GZ_TEXTURE_MAP;
			valueListShader[1] = (GzPointer)(ptex_fun);
			valueListShader[3] = (GzPointer)false;
		} else {
			valueListShader[0] = (GzPointer)false;
			nameListShader[1]  = GZ_TEXTURE_MAP;
			valueListShader[1] = 0;
			valueListShader[3] = (GzPointer)false;
		}

		for (int aaPass = 0; aaPass < AAKERNEL_SIZE; aaPass++) {
			status |= GzPutAttribute(m_pAARenders[aaPass], 4, nameListShader, valueListShader);
		}

		/* 
		* Walk through the list of triangles, set color 
		* and render each triangle 
		*/ 
		while( fscanf(infile, "%s", dummy) == 1) { 	/* read in tri word */
			fscanf(infile, "%f %f %f %f %f %f %f %f", 
				&(vertexList[0][0]), &(vertexList[0][1]),  
				&(vertexList[0][2]), 
				&(normalList[0][0]), &(normalList[0][1]), 	
				&(normalList[0][2]), 
				&(uvList[0][0]), &(uvList[0][1]) ); 
			fscanf(infile, "%f %f %f %f %f %f %f %f", 
				&(vertexList[1][0]), &(vertexList[1][1]), 	
				&(vertexList[1][2]), 
				&(normalList[1][0]), &(normalList[1][1]), 	
				&(normalList[1][2]), 
				&(uvList[1][0]), &(uvList[1][1]) ); 
			fscanf(infile, "%f %f %f %f %f %f %f %f", 
				&(vertexList[2][0]), &(vertexList[2][1]), 	
				&(vertexList[2][2]), 
				&(normalList[2][0]), &(normalList[2][1]), 	
				&(normalList[2][2]), 
				&(uvList[2][0]), &(uvList[2][1]) ); 

			/* 
			* Set the value pointers to the first vertex of the 	
			* triangle, then feed it to the renderer 
			* NOTE: this sequence matches the nameList token sequence
			*/ 
			valueListTriangle[0] = (GzPointer)vertexList; 
			valueListTriangle[1] = (GzPointer)normalList; 
			valueListTriangle[2] = (GzPointer)uvList;

			for (int aaPass = 0; aaPass < AAKERNEL_SIZE; aaPass++) {
				GzPutTriangle(m_pAARenders[aaPass], 3, nameListTriangle, valueListTriangle); 
			}
		}

		if( fclose( infile ) )
			AfxMessageBox( "The input file was not closed\n" );
	}

	FinishAA(m_pAADisplays, m_pDisplay);

	FILE *outfile;
	if( (outfile  = fopen( OUTFILE , "wb" )) == NULL )
	{
		AfxMessageBox( "The output file was not opened\n" );
		return GZ_FAILURE;
	}

	GzFlushDisplay2File(outfile, m_pDisplay); 	/* write out or update display to file*/
	GzFlushDisplay2FrameBuffer(m_pFrameBuffer, m_pDisplay);	// write out or update display to frame buffer

	/* 
	* Close file
	*/ 

	if( fclose( outfile ) )
		AfxMessageBox( "The output file was not closed\n" );

	return GZ_SUCCESS;
}

void FinishAA(GzDisplay** aaDisplays, GzDisplay* outputDisplay) 
{
	GzIntensity tempR, tempG, tempB, tempA, finalR, finalG, finalB, finalA;
	//can ignore Z
	GzDepth bufferZ;

	for (int x = 0; x < outputDisplay->xres; x++) {
		for (int y = 0; y < outputDisplay->yres; y++) {

			finalR = finalG = finalB = finalA = 0;
			for (int aaPass = 0; aaPass < AAKERNEL_SIZE; aaPass++) {
				GzGetDisplay(aaDisplays[aaPass], x, y, &tempR, &tempG, &tempB, &tempA, &bufferZ);
				finalR += tempR * AAFilter[aaPass][AAFILTER_WEIGHT];
				finalG += tempG * AAFilter[aaPass][AAFILTER_WEIGHT];
				finalB += tempB * AAFilter[aaPass][AAFILTER_WEIGHT];
				finalA += tempA * AAFilter[aaPass][AAFILTER_WEIGHT];
			}
			GzPutDisplay(outputDisplay, x, y, finalR, finalG, finalB, finalA, 0);
		}
	}
}

int ApplicationFinal::SetReflective(bool isReflective) 
{
	IsModelReflective[0] = isReflective;
	return GZ_SUCCESS;
}

int ApplicationFinal::Clean()
{
	/* 
	* Clean up and exit 
	*/ 
	int	status = 0; 
	
	status |= GzFreeDisplay(m_pDisplay);
	status |= GzFreeTexture();

	if (status) 
		return(GZ_FAILURE); 
	else 
		return(GZ_SUCCESS);
}



