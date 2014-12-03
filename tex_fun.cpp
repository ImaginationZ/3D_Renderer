/* Texture functions for cs580 GzLib	*/
#include    "stdafx.h" 
#include	"stdio.h"
#include	"Gz.h"
#include	"disp.h"
#include	"math.h"
#include	<functional>

GzColor	*image=NULL;
GzColor	*cubeMap=NULL;
int xs, ys;
int cubeWidth;
int reset = 1;
int environmentReset = 1;

#define	TEX_ARRAY(x,y)	(x+(y*xs))
#define	CUBE_ARRAY(x,y)	(x+(y*cubeWidth*4))

float scaleToPositiveRange(float input);

/* Image texture function */
int tex_fun(float u, float v, GzColor color)
{
  unsigned char		pixel[3];
  unsigned char     dummy;
  char  		foo[8];
  int   		i, j;
  FILE			*fd;

  if (reset) {          /* open and load texture file */
    fd = fopen ("texture", "rb");
    if (fd == NULL) {
      fprintf (stderr, "texture file not found\n");
      exit(-1);
    }
    fscanf (fd, "%s %d %d %c", foo, &xs, &ys, &dummy);
    image = (GzColor*)malloc(sizeof(GzColor)*(xs+1)*(ys+1));
    if (image == NULL) {
      fprintf (stderr, "malloc for texture image failed\n");
      exit(-1);
    }

    for (i = 0; i < xs*ys; i++) {	/* create array of GzColor values */
      fread(pixel, sizeof(pixel), 1, fd);
      image[i][RED] = (float)((int)pixel[RED]) * (1.0 / 255.0);
      image[i][GREEN] = (float)((int)pixel[GREEN]) * (1.0 / 255.0);
      image[i][BLUE] = (float)((int)pixel[BLUE]) * (1.0 / 255.0);
      }

    reset = 0;          /* init is done */
	fclose(fd);
  }

/* bounds-test u,v to make sure nothing will overflow image array bounds */
/* determine texture cell corner values and perform bilinear interpolation */
/* set color to interpolated GzColor value and return */

  if (u < 0) {
	  u = 0;
  }
  if (u > 1) {
	  u = 1;
  }
  if (v < 0) {
	  v = 0;
  }
  if (v > 1) {
	  v = 1;
  }
  
  float textureX = u * (xs - 1);
  float textureY = v * (ys - 1);

  GzColor colorA, colorB, colorC, colorD;
  memcpy(colorA, image[(int)TEX_ARRAY(floor(textureX), floor(textureY))], sizeof(GzColor));
  memcpy(colorB, image[(int)TEX_ARRAY(ceil(textureX), floor(textureY))], sizeof(GzColor));
  memcpy(colorC, image[(int)TEX_ARRAY(ceil(textureX), ceil(textureY))], sizeof(GzColor));
  memcpy(colorD, image[(int)TEX_ARRAY(floor(textureX), ceil(textureY))], sizeof(GzColor));

  float s = textureX - floor(textureX);
  float t = textureY - floor(textureY);
  for (int colorComp = 0; colorComp < 3; colorComp++) {
	  color[colorComp] = s * t * colorC[colorComp] + 
		  (1 - s) * t * colorD[colorComp] + 
		  s * (1 - t) * colorB[colorComp] + 
		  (1 - s) * (1 - t) * colorA[colorComp];
  }

  return GZ_SUCCESS;
}

/* Procedural texture function */
int ptex_fun(float u, float v, GzColor color)
{
	GzColor foreColor;
	foreColor[RED] = 0.9;
	foreColor[GREEN] = 0.9;
	foreColor[BLUE] = 0.9;

	GzColor backgroundColor;
	backgroundColor[RED] = 0.1;
	backgroundColor[GREEN] = 0.1;
	backgroundColor[BLUE] = 0.1;

	bool isColored = false;


#if 1
	//used for stripe pattern
	int bandCount = 4;
	float bandWidth = .01;
	for (int i = 0; i < bandCount; i++) {
		float bandLocation = (1.0 / bandCount) * i;
		if ((u > bandLocation - bandWidth && u < bandLocation + bandWidth) ||
			(v > bandLocation - bandWidth && v < bandLocation + bandWidth)) {
			isColored = true;
			break;
		}
	}
#else
	int bandCount = 5;
	float bandWidth = .05;
	for (int i = 0; i < bandCount; i++) {
		for (int j = 0; j < bandCount; j++) {
			float bandA = ((1.0 / bandCount) * i) + bandWidth;
			float bandB = ((1.0 / bandCount) * j) + bandWidth;
			if (sqrt(pow(u - bandA, 2) + pow(v - bandB, 2)) < bandWidth) {
				isColored = true;
				break;
			}
		}
		if (isColored) {
			break;
		}
	}
#endif

	if (isColored) {
		memcpy(color, foreColor, sizeof(GzColor));
	} else {
		memcpy(color, backgroundColor, sizeof(GzColor));
	}

	return GZ_SUCCESS;
}

int cubetex_fun(GzCoord reflection, GzColor color) {

	unsigned char		pixel[3];
	unsigned char     dummy;
	char  		foo[8];
	int   		i, j;
	int			texWidth, texHeight;
	FILE			*fd;

	/* open and load environment map file */
	if (environmentReset) {
		fd = fopen ("em_skybox.ppm", "rb");
		if (fd == NULL) {
			fprintf (stderr, "environment map texture file (em_output.ppm) not found\n");
			exit(-1);
		}
		fscanf (fd, "%s %d %d %c", foo, &texWidth, &texHeight, &dummy);
		cubeWidth = texWidth / 4; //four cube faces across
		cubeMap = (GzColor*)malloc(sizeof(GzColor)*(texWidth+1)*(texHeight+1));
		if (cubeMap == NULL) {
			fprintf (stderr, "malloc for texture image failed\n");
			exit(-1);
		}

		for (i = 0; i < texWidth*texHeight; i++) {	/* create array of GzColor values */
			fread(pixel, sizeof(pixel), 1, fd);
			cubeMap[i][RED] = (float)((int)pixel[RED]) * (1.0 / 255.0);
			cubeMap[i][GREEN] = (float)((int)pixel[GREEN]) * (1.0 / 255.0);
			cubeMap[i][BLUE] = (float)((int)pixel[BLUE]) * (1.0 / 255.0);
		}

		environmentReset = 0;
		fclose(fd);
	}

	//check which cube face this reflection maps to
	float absX = abs(reflection[X]);
	float absY = abs(reflection[Y]);
	float absZ = abs(reflection[Z]);

	float faceOffsetX, faceOffsetY;
	float textureX, textureY;
	
	if (absX > absY && absX > absZ) {
		textureX = (reflection[Z] / reflection[X]);
		textureY = (reflection[Y] / reflection[X]);
		if (reflection[X] >= 0) {
			faceOffsetX = cubeWidth * 2;
			faceOffsetY = cubeWidth;
		} else {
			faceOffsetX = 0;
			faceOffsetY = cubeWidth;
		}
	} else if (absY > absX && absY > absZ) {
		textureX = (reflection[X] / reflection[Y]);
		textureY = (reflection[Z] / reflection[Y]);
		if (reflection[Y] >= 0) {
			faceOffsetX = cubeWidth;
			faceOffsetY = 0;
		} else {
			faceOffsetX = cubeWidth;
			faceOffsetY = cubeWidth * 2;
		}
	} else {
		textureX = (reflection[X] / reflection[Z]);
		textureY = (reflection[Y] / reflection[Z]);
		if (reflection[Z] >= 0) {
			faceOffsetX = cubeWidth;
			faceOffsetY = cubeWidth;
		} else {
			faceOffsetX = cubeWidth * 3;
			faceOffsetY = cubeWidth;
		}
	}
	textureX = (scaleToPositiveRange(textureX) * (cubeWidth - 1)) + faceOffsetX;
	textureY = (scaleToPositiveRange(textureY) * (cubeWidth - 1)) + faceOffsetY;
	
	GzColor colorA, colorB, colorC, colorD;
	memcpy(colorA, cubeMap[(int)CUBE_ARRAY(floor(textureX), floor(textureY))], sizeof(GzColor));
	memcpy(colorB, cubeMap[(int)CUBE_ARRAY(ceil(textureX), floor(textureY))], sizeof(GzColor));
	memcpy(colorC, cubeMap[(int)CUBE_ARRAY(ceil(textureX), ceil(textureY))], sizeof(GzColor));
	memcpy(colorD, cubeMap[(int)CUBE_ARRAY(floor(textureX), ceil(textureY))], sizeof(GzColor));

	float s = textureX - floor(textureX);
	float t = textureY - floor(textureY);
	for (int colorComp = 0; colorComp < 3; colorComp++) {
		color[colorComp] = s * t * colorC[colorComp] + 
			(1 - s) * t * colorD[colorComp] + 
			s * (1 - t) * colorB[colorComp] + 
			(1 - s) * (1 - t) * colorA[colorComp];
	}

	return GZ_SUCCESS;
}

//scale input from range [-1, 1] to [0, 1]
float scaleToPositiveRange(float input) {
	float inMin = -1;
	float inMax = 1;
	float outMin = 0;
	float outMax = 1;
	return ((outMax - outMin) * (input - inMin) / (inMax - inMin)) + outMin;
}

/* Free texture memory */
int GzFreeTexture()
{
	if(image!=NULL)
		free(image);
	return GZ_SUCCESS;
}

