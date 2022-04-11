/**
*	@file NB_calculator.java
*
*	@brief This file implments the NB_Calculator class
*
*	Number and brightness analysis for ImageJ. The file contains a single public
* class that implements number and brightness analysis in java. This method was
* described by Michelle Digman et al. 2008 in "Mapping the Number of Molecules
* and Brightness in the Laser Scanning Microscope" Biophysical Journal 94:
* 2320-2332. The method relies on analysing the noise/variation in a series of
* images of the same field.
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
*	@brief A plugin class for ImageJ to implement number and brightness analysis.
*
*	This class is an implementation of Number and Brightness analysis for ImageJ.
*	The initial input is a stack of images. For details on how these images should
*	be acquired see the Digman article or the User manual. The primary analyses
* calculates for each pixel an estimation of the number and brightness of the
* fluorescent objects contributing to the image.
*
*	The primary analysis calculates the following:
*	\f{eqnarray*}
*		Mean&=&\langle k\rangle _{pixel}
*				= \frac{ \sum_{images} k_{pixel} }{N_{images}}\\
*		Variance&=&\sigma^{2}
*				= \frac{ \sum_{images} (k_{pixel} -\langle k\rangle _{pixel} )^{2} }{N_{images}}\\
*		Number&=&\langle N\rangle
*				= \frac{\langle k\rangle ^{2}}{\sigma^{2}}\\
*		Brightness&=&B = \frac{\langle k\rangle }{\langle N\rangle }
*				= \frac{\sigma^{2}}{\langle k\rangle }\\
*		Corrected Number&=&\langle N\rangle _{corr}
*				= \frac{\langle k>^{2}}{(\sigma^{2}-\langle k\rangle }\\
*		Corrected Brightness&=&B_{corr}
*				= \frac{\sigma^{2}-\langle k\rangle }{\langle k\rangle }\\
*	\f}
*
*	The stack resulting from this primary analysis can then be used for secondary
*	analyses to examine correlations between the parameters, the distribution of
* different values and then used to select different regions.
*
*	@todo Add possibilities for secondary analysis (see Gratton lecture)
*	This should include at least table like output and graphs for passage to a
*	graphing package and selection of regions in the average image based on
*	number and brightness criteria. Implementation of this should also allow
*	secondary analysis to be performed on precomputed stacks that can be read in.
* Perhaps this needs an NB_analysis class rather than an extension to the
* calculator class.
*
*	@todo Improve display of equations in Ref_Manual.
*
*	@todo Introduction in the User_Manual explaining equations, logic and
* acquisition parameters with examples.
*/
public class NB_calculator implements PlugInFilter {
  private ImagePlus image; /**< Storage for the image to work on */

  public int setup(String arg, ImagePlus image) {
  /**
  *
  *	@brief	Initialization function.
  *
  *	@param arg   String passed from ImageJ.
  *	@param image The image that is concerned, this is saved in the private
  *              member data image.
  *	@return      Flags determining the types of objects that can be treated
  *              STACK_REQUIRED + DOES_ALL
  *
  *	This function is used by ImageJ to determine what types of image that can be
  * treated, for this plugin a stack of images is requred (STACK_REQUIRED), but
  * the type of image (8bit, 16bit,	grey scale or color) is not important
  * (DOES_ALL). The function is also used to initialise	the image with the
  * appropriate structure.
  **/
	this.image = image;					       // Initialize the image member data.
	return STACK_REQUIRED + DOES_ALL;  // We need a stack of images to work on.
  }

  public void run(ImageProcessor ip) {
  /**
  *
  *	@brief The main function of the NB_calculator plugin.
  *
  *	This function after the initialization section, first calculates the sum and
  * the sum of squared intensities, then calculates the Mean, Variance, Number,
  * Brightness, Corrected number and Corrected brightness images, finally the
  * different images are packed into a stack and displayed.
  */

	/*************************************
	**									                **
	**		Initialization section.		    **
	**									                **
	*************************************/

	int	width  = image.getWidth();					      // Size of the images to process
	int height = image.getHeight();					      // width and height and number
	int Nimage = image.getStack().getSize();		  // of images in the stack.
                                                // Space for output images
	float [][] frame1 = new float[width][height];	// average,
	float [][] frame2 = new float[width][height];	// variance,
	float [][] frame3 = new float[width][height];	// number,
	float [][] frame4 = new float[width][height];	// brightness
	float [][] frame5 = new float[width][height];	// corrected number
	float [][] frame6 = new float[width][height];	// and corrected brightness

	/*************************************
	**									                **
	**		Calculation of Sums.  		    **
	**									                **
	*************************************/

    for (int i=1; i <= Nimage; i=i+1) {				  // Loop over images in stack
      ImageProcessor imgtest = image.getStack().getProcessor(i);
      for( int x_pos = 0; x_pos < width; x_pos++){ // Loop over pixels
        for( int y_pos = 0; y_pos < height; y_pos++){
          frame1[x_pos][y_pos] += imgtest.getPixelValue(x_pos,y_pos);
          frame2[x_pos][y_pos] += imgtest.getPixelValue(x_pos,y_pos)
                                * imgtest.getPixelValue(x_pos,y_pos);
        }
      }
    }

	/*************************************
	**									                **
	**	Calculation of Derived images.  **
	**									                **
	*************************************/

    for( int x_pos = 0; x_pos < width; x_pos++){	// Loop over pixels doing calculations
    	for( int y_pos = 0; y_pos < height; y_pos++){
        frame1[x_pos][y_pos] /= Nimage;
        frame2[x_pos][y_pos] /= Nimage;
        frame2[x_pos][y_pos] -= frame1[x_pos][y_pos] * frame1[x_pos][y_pos];
        frame3[x_pos][y_pos] = (frame1[x_pos][y_pos] * frame1[x_pos][y_pos])
                               / frame2[x_pos][y_pos];
        frame4[x_pos][y_pos] =  frame2[x_pos][y_pos] / frame1[x_pos][y_pos];
        frame5[x_pos][y_pos] = (frame1[x_pos][y_pos] * frame1[x_pos][y_pos])
                              / (frame2[x_pos][y_pos] - frame1[x_pos][y_pos]);
        if( frame5[x_pos][y_pos] < 0.0 ) frame5[x_pos][y_pos] = (float)0.0;
        frame6[x_pos][y_pos] = (frame2[x_pos][y_pos] - frame1[x_pos][y_pos])
                             /  frame1[x_pos][y_pos];
        if( frame6[x_pos][y_pos] < 0.0 ) frame6[x_pos][y_pos] = (float)0.0;
      }
    }

	/*************************************
	**									                **
	**	  Display of result and exit.   **
	**									                **
	*************************************/
												// Put arrays into Image processors
    ImageProcessor proc1 = new FloatProcessor( frame1 );
    ImageProcessor proc2 = new FloatProcessor( frame2 );
    ImageProcessor proc3 = new FloatProcessor( frame3 );
    ImageProcessor proc4 = new FloatProcessor( frame4 );
    ImageProcessor proc5 = new FloatProcessor( frame5 );
    ImageProcessor proc6 = new FloatProcessor( frame6 );
    											// Create an Image stack
    ImageStack stack = new ImageStack( width, height );

    stack.addSlice( "Average", proc1 );			// Add Image processors to the stack
    stack.addSlice( "Variance", proc2 );
    stack.addSlice( "Raw Number", proc3 );
    stack.addSlice( "Raw Brightness", proc4 );
    stack.addSlice( "Number", proc5 );
    stack.addSlice( "Brightness", proc6 );

    											// Create an Image plus for the result
    											// reusing initial object.
    image = new ImagePlus("N & B analysis", stack );
    image.show();								// put it on the screen,
    image.updateAndDraw();						// and make sure it is visible.
  }
}
