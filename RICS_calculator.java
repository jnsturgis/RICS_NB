/**
*	@file RICS_calculator.java
*
*	@brief This file implments the RICS_Calculator class.
*
*	Calculation of an auto correlation function from a series of images. This file
* contains a single class that implements the calculation of the autocorrelation
* function in java, interfaced for use in the ImageJ image manipulation suite.
* This autocorrelation function can then be used for a RICS (Raster Image
* Correlation Spectrum) analysis of particle number brightness and mobility. The
* RICS method was described by Michelle Digman et al. 2005 in Measuring Fast
* Dynamics in Solutions and Cells with a Laser Scanning Microscope Biophysical
* Journal 89: 1317-1327.
*
*	@author James N Sturgis
*	@date 2009-2022
*/

import ij.*;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.*;
import ij.plugin.filter.PlugInFilter;
import ij.gui.GenericDialog;
import ij.io.OpenDialog;
import ij.io.Opener;
import java.util.regex.*;

/**
*	@brief A plugin for ImageJ to calculate an image auto-correlation function.
*
*	This class is part of my implementation of RICS analysis for ImageJ.
*	This class takes a stack of images as input and calculates the spatial
* auto-correlation function.
*	\f{
*		ACF(x,y) = \frac{\langle I(\zeta ,\psi)\timesI(x+\zeta ,y+\psi )\rangle _{\zeta ,\psi}}
*					{\langle I(\zeta ,\psi)\rangle _{\zeta ,\psi}}
*	\f}
*	The resulting correlation can then be used to analyse the microscope PSF and
*	the diffusion of fluorescent objects currently implemented by the RICS_fit
*	class.
*
*	@todo The initialization section should use a dialog to allow several
*		  different but related calculations intensity-correlation or fluctuation-
*		  correlation, cross-correlation or auto-correlation, with or without
*		  baseline subtraction and allow user specification of the size of the
*     correlation map.
*
*	@todo Sanity checks and error management should be implemented.
*
*	@todo Improve display of equations in Ref_Manual.
*
*	@todo Introduction in the User_Manual explaining equations, logic and
*	      acquisition parameters with examples.
*/

public class RICS_calculator implements PlugInFilter {
  private ImagePlus image; /**< Place to strore the image initialized in setup(). */

  public int setup(String arg, ImagePlus image) {
  /**
  *
  *	@brief	Initialization function.
  *
  *	@param arg String passed from ImageJ
  *	@param image The image that is concerned, this is saved.
  *	@return Flags determining the types of objects that can be treated
  * STACK_REQUIRED + DOES_ALL
  *
  *	This function is used by ImageJ to determine what types of image that can be
  * treated, for this plug in a stack of images is requred (STACK_REQUIRED), but
  * the type of image (8bit, 16bit, grey scale or color) is not important
  * (DOES_ALL). The function is also used to initialise the image with the
  * appropriate structure.
  */
    this.image = image;
    return STACK_REQUIRED + DOES_ALL;
  }

  public void run(ImageProcessor ip) {
  /**
  *
  *	@brief The main function of the RICS_calculator plugin.
  *
  *	This function has three parts an initialization section, where various
  *	sizes are calculated, a section where the mean pixel intensity (the
  *	reference for the correlation calculation) is computed and then the last
  *	section calculates in real space the correlation map, finally the program
  * creates an image of the auto-correlation map and displays it.
  *
  *	The calculated values for the mean intensity are:
  *		\f{equation}
  *			\langle I\rangle = \frac{\sum_{pixels} I_{image} }{n_{pixel}}
  *		\f}
  *	with the sum made over all pixels in all images of the stack.
  *
  *	For the correlation function at an offset of \f$\delta \f$ the calculation
  *	is:
  *		\f{equation}
  *			C(\delta )= \frac{\sum_{pixels}(I(0)_{ref}-<I>)(I(\delta)_{image}-<I>)}
  *					{n_{pixel}}
  *		\f}
  *	With the sum made over all pairs possible ensuring a constant number of
  *	pixels for all different offsets. This is thus the image height and width
  *	reduced by the size of the correlation map.
  */

	/*************************************
	**									                **
	**		Initialization section.		    **
	**									                **
	*************************************/

    int	width  = image.getWidth(); 			// Data on current image.
    int height = image.getHeight();
    int Nimage = image.getStack().getSize();

    int	max_offset = 16;					// Space for calculation.
    float [][] correlation = new float[2*max_offset+1][2*max_offset+1];
    float	count = (width - 2 * max_offset)*(height - 2*max_offset) * Nimage;


	/*************************************
	**									                **
	**	Calculation of Mean Intensity.  **
	**									                **
	*************************************/

    double	mean  = 0.0;

    for (int i=1; i <= Nimage; i=i+1) {
    	ImageProcessor imgtest = image.getStack().getProcessor(i);
    	for( int x_pos = 0; x_pos < width; x_pos++){
			  for( int y_pos = 0; y_pos < height; y_pos++){
	  			mean += imgtest.getPixel( x_pos, y_pos );
			  }
      }
    }
    mean /= width * height * Nimage;
    IJ.log( "Mean intensity is "+mean );
    double ref = mean;

	/*************************************
	**									                **
	**  Calculation of correlation map. **
	**									                **
	*************************************/

    for (int i=1; i <= Nimage; i=i+1) {
    	ImageProcessor imgtest = image.getStack().getProcessor(i);
    	ImageProcessor imgref  = imgtest;
    	for (int x_pos = max_offset; x_pos < width - max_offset; x_pos++ ) {
			  for (int y_pos = max_offset; y_pos < height - max_offset; y_pos++ ) {
	  			for (int x_offset = -max_offset; x_offset <= max_offset; x_offset++ ) {
	    			for (int y_offset = -max_offset; y_offset <= max_offset; y_offset++ ) {
	    				double pixeltest =
	    					imgtest.getPixel( x_pos + x_offset, y_pos + y_offset) - ref;
						  double pixelref  =
							  imgref.getPixel( x_pos, y_pos ) - ref;
						  correlation[x_offset + max_offset][y_offset + max_offset] +=
							  (pixeltest * pixelref)/(ref*ref*count);
					  }
				  }
			  }
    	}
    }

	/*************************************
	**									                **
	**	  Display of result and exit.   **
	**									                **
	*************************************/

    ImageProcessor myprocessor = new FloatProcessor( correlation );
    image = new ImagePlus("Stack Correlation", myprocessor );
    image.show();
    image.updateAndDraw();

  }
}
