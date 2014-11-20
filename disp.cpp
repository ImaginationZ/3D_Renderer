/*   CS580 HW   */
#include    "stdafx.h"  
#include	"Gz.h"
#include	"disp.h"
#include <iostream>
#include <string>
#include <limits.h>

//C doesn't have byte type, so we use char array
int GzNewFrameBuffer(char** framebuffer, int width, int height)
{
/* create a framebuffer:
 -- allocate memory for framebuffer : (sizeof)GzPixel x width x height
 -- pass back pointer 
 -- NOTE: this function is optional and not part of the API, but you may want to use it within the display function.
*/
	*framebuffer = (char*)malloc(sizeof(GzPixel) * width * height);
	
	return GZ_SUCCESS;
}

int GzNewDisplay(GzDisplay	**display, int xRes, int yRes)
{
/* create a display:
  -- allocate memory for indicated resolution
  -- pass back pointer to GzDisplay object in display
*/
	*display = new GzDisplay();
	(*display)->xres = xRes;
	(*display)->yres = yRes;

	char* fbuf;
	GzNewFrameBuffer(&fbuf, xRes, yRes);
	(*display)->fbuf = (GzPixel*)fbuf;

	return GZ_SUCCESS;
}


int GzFreeDisplay(GzDisplay	*display)
{
/* clean up, free memory */

	free(display->fbuf);
	free(display);

	return GZ_SUCCESS;
}


int GzGetDisplayParams(GzDisplay *display, int *xRes, int *yRes)
{
/* pass back values for a display */

	*xRes = display->xres;
	*yRes = display->yres;

	return GZ_SUCCESS;
}


int GzInitDisplay(GzDisplay	*display)
{
/* set everything to some default values - start a new frame */

	int xRes = display->xres;
	int yRes = display->yres;

	for (int y = 0; y < yRes; y++) {
		for (int x = 0; x < xRes; x++) {
			GzPutDisplay(display, x, y, 2000, 2000, 4000, 1, INT_MAX);
		}
	}

	return GZ_SUCCESS;
}


GzIntensity snapToRange(GzIntensity val) 
{
	if (val < 0) {
		return 0;
	} else if (val > 4095) {
		return 4095;
	} else {
		return val;
	}
}


int GzPutDisplay(GzDisplay *display, int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
/* write pixel values into the display */

	if (i >= display->xres || i < 0 || j >= display->yres || j < 0) {
		return GZ_SUCCESS;
	}

	int index = ARRAY(i, j);
	GzPixel px = {snapToRange(r), snapToRange(g), snapToRange(b), a, z};
	*(display->fbuf + index) = px;

	return GZ_SUCCESS;
}


int GzGetDisplay(GzDisplay *display, int i, int j, GzIntensity *r, GzIntensity *g, GzIntensity *b, GzIntensity *a, GzDepth *z)
{
	/* pass back pixel value in the display */
	if (i >= display->xres || i < 0 || j >= display->yres || j < 0) {
		return GZ_FAILURE;
	}

	GzPixel pixel = display->fbuf[ARRAY(i, j)];
	*r = pixel.red;
	*g = pixel.green;
	*b = pixel.blue;
	*a = pixel.alpha;
	*z = pixel.z;

	return GZ_SUCCESS;
}


int GzFlushDisplay2File(FILE* outfile, GzDisplay *display)
{

	/* write pixels to ppm file -- "P6 %d %d 255\r" */

	fprintf(outfile, "P6 %d %d 255\r", display->xres, display->yres);
	for (int y = 0; y < display->yres; y++) {
		for (int x = 0; x < display->xres; x++) {
			int displayIndex = ARRAY(x, y);
			
			GzPixel* thisPixel = display->fbuf + displayIndex;
			
			char r = thisPixel->red >> 4;
			char g = thisPixel->green >> 4;
			char b = thisPixel->blue >> 4;
			fwrite(&r, 1, 1, outfile);
			fwrite(&g, 1, 1, outfile);
			fwrite(&b, 1, 1, outfile);
		}
	}
	
	return GZ_SUCCESS;
}

int GzFlushDisplay2FrameBuffer(char* framebuffer, GzDisplay *display)
{

	/* write pixels to framebuffer: 
		- Put the pixels into the frame buffer
		- Caution: store the pixel to the frame buffer as the order of blue, green, and red 
		- Not red, green, and blue !!!
	*/

	for (int y = 0; y < display->yres; y++) {
		for (int x = 0; x < display->xres; x++) {
			int displayIndex = ARRAY(x, y);
			int frameIndex = displayIndex * sizeof(char) * 3;
			
			GzPixel* thisPixel = display->fbuf + displayIndex;
		
			*(framebuffer + frameIndex) = (char)(thisPixel->blue >> 4);
			*(framebuffer + frameIndex + 1) = (char)(thisPixel->green >> 4);
			*(framebuffer + frameIndex + 2) = (char)(thisPixel->red >> 4);
		}
	}

	return GZ_SUCCESS;
}