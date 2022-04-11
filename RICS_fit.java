/**
*	@file RICS_fit.java
*
*	This file implements tha RICS_fit class.
*
*	@author James N Sturgis
*	@date 2009-2015
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
*	@brief A plugin class for ImageJ to estimate diffusion and PSF from RICS data.
*
*	This class is part of my implementation of RICS analysis for ImageJ.
*	This class takes as input an image created by the RICS_calculator class that
*	contains an image auto-correlation map. This map, when calculated on an image
*	obtained for RICS analysis using fluctuation correlations represents the point
*	spread function (PSF) of the microscope as distorted by movement of objects during
*	image acquisition. The aim of this program is to deconvolute this image into
*	two parts the PSF and the fluorescence correlation spectrum (FCS) associated with
*	the movement induced distortion.
*
*	This is described by the equation:
*	\f{equation}
*		Autocorrelation(x,y) = PSF(x,y) x FCS( delay(x,y) )
*	\f}
*	Where the \f$ delay() \f$ function describes the time lag between the different
*	positions in the autocorrelation function and the central point. This in turn
*	depends on the scanning parameters.
*
*	This deconvolution can be used in two ways, first to estimate the PSF by
*	providing information on the movement of the fluorophores during the
*	acquisition, or second to estimate the FCS by providing information on
*	the PSF. The deconvolved functions can be exported for further treatment the
*	PSF as an image and the FCS as a table.
*
*	The program also includes simple analyses to fit basic models to the estimated
*	PSF or FCS. The fitting is performed to the Autocorrelation data adjusting either
*	the PSF model or the FCS model as appropriate. The basic PSF model is a symmetric
*	two dimensionnal gaussian:
*	\f{equation}
		PSF(x,y) = exp^{-\frac{(x-x_0)^2 
					+ (y-y_0)^2 *Dx*Dx}{(1+\frac{\tau}{\tau_d})(Wxy*Wxy)}}
	\f}
*	Where \f$ x_0 \f$ is the x origin, \f$ y_0 \f$ is the y origin, Dx is the
*	size of a pixel in microns, and Wxy is the beam waist in microns. The
*	delay function is given by:
*	\f{equation}
		delay(x,y) = |x-x_0| * D_{pixel} + |y-y_0| * D_{line}
	\f}
*	Where \f$ D_{pixel} \f$ is the pixel step time and \f$ D_{line} \f$ is the line
*	step time. Finally the basic FCS model is a simple 3 dimensionnal diffusion
*	given by:
*	\f{equation}
		FCS(\tau) = \frac{1}{(1+\frac{\tau}{\tau_d})*\sqrt(1+Z^{2}*\frac{\tau}{\tau_d})}
	\f}
*	Where Z is the ratio of the beam waist in the x,y plane to that in the z direction, and
*	\f$ \tau_d \f$ is the characteristic diffusion time.
*/

public class RICS_fit implements PlugInFilter {
  private ImagePlus image;	/**< Place to store the image initialized in setup(). */

  public int setup(String arg, ImagePlus image) {
/**
*
*	@brief	Initialization function.
*
*	This function is used by ImageJ to determine the types of image that can be treated.
*	@param arg String passed from ImageJ
*	@param image The image that is concerned, this is saved.
*	@return DOES_32 Flags determining the types of object treated, in this case float images.
*
*/
    this.image = image;
    return DOES_32;	// We need an image of the auto-correlation function to fit.
  }

  public void run(ImageProcessor ip) {
/**
*
*	@brief The main function of the RICS_calculator plugin.
*
*	This function first in the initialization section calls the fitting
*	dialog so the user can enter parameters and decide what calculations to
*	to perform. Then second the program as necessary loads an additional image or
*	calculates one. Third the program calculates by division an estimation of
*	the missing  part of the deconvolution. If requested by the user this
*	estimation is then saved. Finally if requested the model function is
*	adjusted to fit the data as well a possible.
*
*	The parameters array is used to pass information between the different routines
*	and is initialized with a set of default values. The order of the parameters
*	in the array and their default values are:
*		- Concentration of objects 1.0 nM.
*		- Molecular brightness in 1.0 Mcps at psf optimum.
*		- Diffusion coefficient \f$ 0.0 \mu m^2 sec^{-1} \f$.
*		- Beam waist in x and y directions in \f$ 0.4 \mu m \f$.
*		- Beam waist/height in z direction in \f$ 1.2 \mu m \f$.
*		- Pixel size in \f$ 0.44 \mu m \f$.
*		- Pixel dwell time in \f$ 8.0 \mu sec \f$.
*		- Line scane time in \f$ 4.0 msec \f$.
*		- A float cast from an integer containing the different
*		  flag values (Load PSF Calc FCS * 64 + Calc PSF calc FCS * 32 +
*		  Calc FCS calc PSF * 16 + Save PSF * 8 + Save FCS * 4 +
*		  Fit PSF * 2 + Fit FCS ). Default 0, do nothing!
*		- Image size (assumed to be square )
*
*	@todo	documentation of fcs output.
*	@todo	handle precision beter with error propagation and real estimates.
*	@todo	only output half of fcs using weighted average values
*	@todo	documentation of fitting output.
*/

	/*************************************
	**									**
	**		Initialization section.		**
	**									**
	*************************************/

	int		size = (int) ip.getWidth();			// ACF size
    int		di = (size - 1)/2;		            // Position of midpoint.
    int     i, j;                               // Counters

	float [] parameters = new float[10];		// Array to hold parameters
	parameters[0] = (float) 1.0;				// Concentration
	parameters[1] = (float) 1.0;				// Molecular brightness
	parameters[2] = (float) 0.0;				// Diffusion coefficient
	parameters[3] = (float) 0.4;				// Beam waist
	parameters[4] = (float) 1.2;				// Beam height
	parameters[5] = (float) 0.044;				// Pixel size
	parameters[6] = (float) 8.0;				// Pixel time
	parameters[7] = (float) 4.0;				// Line time
	parameters[8] = (float) 0.0;				// Flags
	parameters[9] = (float) size;	            // Size of ACF

    GenericDialog	gd = new GenericDialog("RICS fitting");
    FIT_dialog( gd, parameters );				// Handle dialog and decide what to do.
    if( gd.wasCanceled() ) return;

    float [][] psf   = new float[size][size];   // Array for calculated psf data
    float [][] delay = new float[size][size];   // Array for pixel delays
    float [][] fcs   = new float[size][size];   // Array for calculated fcs data
	float [][] acf   = new float[size][size];   // Array for input acf
	float [][] resid = new float[size][size];   // Array for calculated residuals

    for( i = 0; i < size; i++ )					// Extract the ACF from the image
        for( j = 0; j < size; j++ ){
            acf[i][j]  = ip.getPixelValue(i,j);
        }

	int	flags = (int) parameters[8];
	double adjust[] = new double[4];
	adjust[0] = 4.0*parameters[2]/(parameters[5]*parameters[5]);
												// rate constant in pixel^2 sec-1
	adjust[1] = parameters[3]/parameters[5];	// Waist in pixels
	adjust[2] = parameters[4]/parameters[3];	// Height/Width ratio
	adjust[3] = 1.0;							// g(0)


	delay_fill( delay, parameters );
	fcs_fill( fcs, size, delay, adjust );
	psf_fill( psf, size, delay, adjust );
	
	fcs[di][di] = (float)0.0;					// Ignore central point contains all
	psf[di][di] = (float)0.0;					// Uncorrelated stuff.
	acf[di][di] = (float)0.0;

	float sum_A  = (float)0.0;					// Estimate scale_factor so sum of
	float sum_FP = (float)0.0;					// residuals is 0.0
	
    for( i = 0; i < size; i++ )
        for( j = 0; j < size; j++ ){
        	sum_A  += acf[i][j];
        	sum_FP += psf[i][j] * fcs[i][j];
        }

	float scale_factor = sum_A/sum_FP;

	adjust[3] *= scale_factor;					// g(0)
	
	IJ.log("Flags ="+flags );

	/*************************************
	**									**
	**	Decide what to do based on the  **
	**  flags set by the dialog. 		**
	**	    0x0001	fit FCS				**
	**	    0x0002	fit PSF				**
	**  Otherwise use model based on	**
	**	input parameters only without	**
	**	any adjustment/fitting.			**
	**									**
	*************************************/

	if ((flags & 0x0002 ) != 0 ){				// Adjust to fit diffusion and g(0)
		my_curve_fit( size, delay, acf, 0, adjust );
	} else if ((flags & 0x0001 ) != 0 ){		// Adjust to fit beam waist W and g(0)
		my_curve_fit( size, delay, acf, 1, adjust );
	} else {									// No adjustment required.
	}
	parameters[2] = (float)(adjust[0] * parameters[5] * parameters[5] / 4.0);
	parameters[3] = (float)(adjust[1] * parameters[5]);
	parameters[4] = (float)(adjust[2] * parameters[3]);
	scale_factor = (float)adjust[3];
	
	fcs_fill( fcs, size, delay, adjust );
	psf_fill( psf, size, delay, adjust );
	fcs[di][di] = (float)0.0;					// Ignore central point contains all
	psf[di][di] = (float)0.0;					// Uncorrelated stuff.
	
	float	sum = (float)0.0;	
	for( i = 0; i < size; i++ )					// Put g(0) into fcs part,
	  for( j = 0; j < size; j++ ){				// and calculate residuals.
        resid[i][j] = acf[i][j] - psf[i][j] * fcs[i][j];
        sum += resid[i][j] * resid[i][j];
      }

	/*************************************
	**									**
	**	  Report results as requested.	**
	**									**
	*************************************/

    IJ.log("Residual is          \t="+sum );
    IJ.log("Diffusion coefficient\t="+parameters[2]+"um²/sec" );
	IJ.log("Beam height          \t="+parameters[4]+"um" );
	IJ.log("Beam waist           \t="+parameters[3]+"um" );
	IJ.log("G(0)                 \t="+scale_factor );
	
	/*************************************
	**									**
	**	Display of residuals and the	**
	**	two parts of the acf and exit.  **
	**									**
	*************************************/

    ImageProcessor myprocessor = new FloatProcessor( resid );
    image = new ImagePlus("Residuals", myprocessor );
    image.show();
    image.updateAndDraw();
    
    ImageProcessor myprocessor2 = new FloatProcessor( fcs );
    image = new ImagePlus("FCS part", myprocessor2 );
    image.show();
    image.updateAndDraw();
    
    ImageProcessor myprocessor3 = new FloatProcessor( psf );
    image = new ImagePlus("PSF part", myprocessor3 );
    image.show();
    image.updateAndDraw();
    
  }

  private void	FIT_dialog(GenericDialog gd, float [] parameters )
  {
/**
*	@brief	Member to handle dialog box
*
*	@param gd an ImageJ generic dialog object.
*	@param parameters an array holding the various parameter values.
*	this is used to pass the default values and return the input
*	values. The order of the parameters in the array is:
*		- Concentration of molecules in nM.
*		- Molecular brightness in Mcps at psf optimum.
*		- Diffusion coefficient \f$ \mu m^2 sec^{-1} \f$.
*		- Beam waist in x and y directions in \f$ \mu m \f$.
*		- Beam waist/height in z direction in \f$ \mu m \f$.
*		- Pixel size in \f$ \mu m \f$.
*		- Pixel dwell time in \f$ \mu sec \f$.
*		- Line scane time in \f$ msec \f$.
*		- A float cast from an integer containing the different
*		  flag values (Load PSF Calc FCS * 64 + Calc PSF calc FCS * 32 +
*		  Calc FCS calc PSF * 16 + Save PSF * 8 + Save FCS * 4 +
*		  Fit PSF * 2 + Fit FCS ).
*	Undobtedly there is a more elegant way of doing this.
*
*	This routine creates the dialog box filling the fields with the
*	initial values from the parameters array, it then displays the
*	dialog box and on exit fills the parameters array as appropriate.
*/
	/*************************************
	**									**
	**	  Build the dialog box with		**
	**	  the default parameters.		**
	**									**
	*************************************/

    gd.addMessage("Object properties");
    gd.addNumericField("Concentration :", parameters[0], 3, 7, "nM" );
    gd.addNumericField("Molecular brightness :", parameters[1], 0, 7, "Mcps" );
    gd.addNumericField("Diffusion coef. :", parameters[2], 1, 7, "um²/sec");

    gd.addMessage("Microscope Optical Properties");
    gd.addNumericField("Beam waist  :", parameters[3], 3, 7, "um" );
    gd.addNumericField("Beam height :", parameters[4], 3, 7, "um" );

    gd.addMessage("Microscope Acquisition Properties");
    gd.addNumericField("Pixel size :", parameters[5], 3, 7, "um" );
    gd.addNumericField("Pixel time :", parameters[6], 0, 7, "usec" );
    gd.addNumericField("Line time  :", parameters[7], 0, 7, "msec"  );

	gd.addMessage("What do you want to do?" );

	gd.addCheckbox("Calculate beam properties and fit diffusion", 0==1 );
	gd.addCheckbox("Calculate diffusion and fit beam properties", 0==1 );

	/*************************************
	**									**
	**	  Allow user interaction.		**
	**									**
	*************************************/

    gd.showDialog();

    if( gd.wasCanceled() ) return;

	/*************************************
	**									**
	**	  Extract the users values.		**
	**									**
	*************************************/

    parameters[0]  = (float) gd.getNextNumber();
    parameters[1]  = (float) gd.getNextNumber();
    parameters[2]  = (float) gd.getNextNumber();
    parameters[3]  = (float) gd.getNextNumber();
    parameters[4]  = (float) gd.getNextNumber();
    parameters[5]  = (float) gd.getNextNumber();
    parameters[6]  = (float) gd.getNextNumber();
    parameters[7]  = (float) gd.getNextNumber();

    int flags = 0;

    flags +=  2 * (gd.getNextBoolean()?1:0);
    flags +=  1 * (gd.getNextBoolean()?1:0);
	parameters[8] = (float)flags;

    return;
  }

  private void  psf_fill( float[][] psf, int size, float[][] delay, double[] adjust )
  {
  /**
  *
  *		@brief Fill an array with a theroetical psf.
  *
  *		@param psf An array to fill with the psf.
  *		@param parameters The set of parameters for the calculation.
  *
  *		The equation for the theoretical psf is taken from Digman et al. 2005, note
  *		the intrusion of the diffusion into the PSF part of the acf as we might
  *		follow a moving particle.
  */
    int   i, j;
    int   di = (size - 1)/2;

	double tau;

    for( i= -di; i <= di; i++ )
      for( j= -di; j <= di; j++ ){
        tau   = delay[i+di][j+di] * adjust[0];
        psf[i+di][j+di] = (float)Math.exp( -(j*j+i*i)/(tau+adjust[1]*adjust[1]) );
      }
    psf[di][di] = (float)0.0;
  }

  private void	fcs_fill( float[][] fcs, int size, float[][] delay, double[] adjust )
  {
  /**
  *
  *		@brief Fill an array with a theroetical raster time delay correlation function.
  *
  *		@param fcs An array to fill with the fcs.
  *		@param parameters The set of parameters for the calculation.
  *
  *		The routine extracts from the parameter array the size of the
  *		destination array, the beam waist and height in \f$ \mu m \f$,
  *		the diffusion constant in \f$ \mu m^2 sec^{-1} \f$, and the
  *		pixel dwell time and line delay time in \f$ \mu sec \f$ and msec
  *		respectively. These values are then used to calculate a theoretical
  *		correlation function based on the time delay between points the
  *		the diffusion constant and the psf size.
  *
  */
    int 	i, j;
    int	di = (size - 1)/2;

	double tau;

	for( i= -di; i <= di; i++ ){
	  for( j= -di; j <= di; j++ ){
    	tau   = delay[i+di][j+di] * adjust[0] / (adjust[1]*adjust[1]);
    	fcs[i+di][j+di]  = (float)( adjust[3] / ((1.0 + tau )
    						* Math.sqrt(1.0 + tau /( adjust[2]*adjust[2] ))));
      }
    }
    fcs[di][di] = (float)0.0;
  }

  private void	delay_fill( float[][] delay, float[] parameters )
  {
  /**
  *
  *		@brief Fill an array with raster time delays based on scan parameters.
  *
  *		@param delay An array to fill with the timing info.
  *		@param parameters The set of parameters for the calculation.
  *
  *		The routine extracts from the parameter array the size of the
  *		destination array, the pixel dwell time and line delay time in
  *		\f$ \mu sec \f$ and msec respectively. These values are then used
  *		to calculate the delay associated with each pixel in seconds.
  *
  */
      int i, j;
	  int	size = (int) parameters[9];
      int	di = (size - 1)/2;
      float	Dtx = parameters[6];
      float	Dty = parameters[7];

	  for( i= -di; i <= di; i++ ){
	  	for( j= -di; j <= di; j++ ){
            delay[i+di][j+di] = (float)Math.abs( Dtx * i * 1E-6 + Dty * j * 1E-3 );
        }
      }
  }

  private boolean my_curve_fit( int size, float[][] delays,
  		float[][] ydata, int which, double[] adjust )
  {
/**
*
*	@brief Curve fitting routine adjusts parameters to fit function to ydata.
*
*/
  	double	fit_param[][] = new double[4][5];
  	double	the_param[]   = new double[4];
  	float   fit_psf[][]   = new float[size][size];
  	float	fit_fcs[][]	  = new float[size][size];
  	
  	int		i,j,k,l;					// General counters
	int		di = (size - 1)/2;
	String	str = "";

  	/***********************************
  	*                                  *
  	*   Initialize fill the fit_param  *
  	*                                  *
  	***********************************/

	if( adjust[which] == 0.0 ) adjust[which] = 1.0;

	for(j = 0; j < 3; j++ )
	  for( i = 0; i < 4; i++ )
		fit_param[j][i+1] = adjust[i];
	
	fit_param[2][which + 1] *= 2.0;
	fit_param[0][which + 1] *= 0.5;

	double sum_y  = 0.0;
	double sum_pf = 0.0;
	double sum_ee = 0.0;
	double scale_factor;
	double error = 0.0;
	double eps_fit;
	double temp;
	
	for( i = 0; i < size; i++ )
	  for(j = 0; j < size; j++ )
	  	sum_y += ydata[i][j];
										// Calculate scale_factors
	for( j = 0; j <= 3; j++ ){			// Calculate sum of square errors
	  for( i = 0; i < 4; i++ ) the_param[i] = fit_param[j][i+1];
	  
	  psf_fill( fit_psf, size, delays, the_param );
	  fcs_fill( fit_fcs, size, delays, the_param );
	  sum_pf = 0.0;
	  for( i = 0; i < size; i++ )
		for(k = 0; k < size; k++ ) sum_pf += fit_psf[i][k]*fit_fcs[i][k];
		
	  scale_factor = sum_y / sum_pf;
	  fit_param[j][4] *= scale_factor;
	  sum_ee = 0.0;
	  for( i = 0; i < size; i++ )
		for(k = 0; k < size; k++ ){
		  error = ydata[i][k] - scale_factor * fit_psf[i][k]*fit_fcs[i][k];
	  	  sum_ee += error * error;
	  	}
	  fit_param[j][0] = sum_ee;
	}
										// Sort the initial guesses not elegant.
	if( fit_param[1][0] < fit_param[0][0] ){
	  for( i = 0; i < 5; i++ ){
	  	fit_param[3][i] = fit_param[0][i];
	  	fit_param[0][i] = fit_param[1][i];
	  	fit_param[1][i] = fit_param[3][i];
	  }
	}
	if( fit_param[2][0] < fit_param[1][0] ){
	  for( i = 0; i < 5; i++ ){
	  	fit_param[3][i] = fit_param[1][i];
	  	fit_param[1][i] = fit_param[2][i];
	  	fit_param[2][i] = fit_param[3][i];
	  }
	}
	if( fit_param[1][0] < fit_param[0][0] ){
	  for( i = 0; i < 5; i++ ){
	  	fit_param[3][i] = fit_param[0][i];
	  	fit_param[0][i] = fit_param[1][i];
	  	fit_param[1][i] = fit_param[3][i];
	  }
	}
	
	for( i = 0; i < 3; i++ ){
	  str = "SSE ="+fit_param[i][0]+", with";
	  for( j = 0; j < 4; j++ ) 
	  	str = str+" "+fit_param[i][j+1];
	  IJ.log( str );
	}
	
	eps_fit = Math.abs((fit_param[2][which + 1] - fit_param[0][which + 1])/fit_param[1][which + 1]);
	IJ.log("Size is "+eps_fit );

  	/***********************************
  	*                                  *
  	*   Aim for minimum and improve    *
  	*                                  *
  	***********************************/
  	int	step = 0;
  	int	max_step = 100;
	double eps = 1E-6;
	
	double	a, b;
	
  	while((step < max_step) & (eps_fit > eps)){
  	  step++;
										// Quadratic fit to three points and find minimum.
	  a  = fit_param[0][0]*(fit_param[1][which+1]-fit_param[2][which+1]);
	  a += fit_param[1][0]*(fit_param[2][which+1]-fit_param[0][which+1]);
	  a += fit_param[2][0]*(fit_param[0][which+1]-fit_param[1][which+1]);

	  b  = fit_param[0][0]*(fit_param[2][which+1]*fit_param[2][which+1]
	  		-fit_param[1][which+1]*fit_param[1][which+1]);
	  b += fit_param[1][0]*(fit_param[0][which+1]*fit_param[0][which+1]
	  		-fit_param[2][which+1]*fit_param[2][which+1]);
	  b += fit_param[2][0]*(fit_param[1][which+1]*fit_param[1][which+1]
	  		-fit_param[0][which+1]*fit_param[0][which+1]);

	  fit_param[3][which+1] = -b / (2.0*a);
										// Calculate for the new set of values
	  for( i = 0; i < 4; i++ ) the_param[i] = fit_param[3][i+1];
	  psf_fill( fit_psf, size, delays, the_param );
	  fcs_fill( fit_fcs, size, delays, the_param );
	  sum_pf = 0.0;
	  for( i = 0; i < size; i++ )
		for(k = 0; k < size; k++ )
		  sum_pf += fit_psf[i][k]*fit_fcs[i][k];
	  scale_factor = sum_y / sum_pf;
	  fit_param[3][4] *= scale_factor;
	  sum_ee = 0.0;
	  for( i = 0; i < size; i++ )
		for(k = 0; k < size; k++ ){
		  error = ydata[i][k] - scale_factor * fit_psf[i][k]*fit_fcs[i][k];
	  	  sum_ee += error * error;
	  	}
	  fit_param[3][0] = sum_ee;
										// Sort the results
	  for( j = 2; (j >= 0 && (fit_param[j+1][0] < fit_param[j][0])); j-- )
		for( i = 0; i < 5; i++ ){
		  temp = fit_param[j][i];
		  fit_param[j][i] = fit_param[j+1][i];
		  fit_param[j+1][i] = temp;	 
		}

	  eps_fit = Math.abs((fit_param[2][which + 1] - fit_param[0][which + 1])/fit_param[1][which + 1]);
	  IJ.log("Size is "+eps_fit );
  	}

  	/***********************************
  	*                                  *
  	*   Return best result or error    *
  	*                                  *
  	***********************************/

  	for( i = 0; i<4; i++ )
  		adjust[i] = fit_param[0][i+1];

	IJ.log("Optimisation stopped after "+step+" moves." );
	for( i = 0; i < 4; i++ ){
	  str = "SSE ="+fit_param[i][0]+", with";
	  for( j = 0; j < 4; j++ ) 
	  	str = str+" "+fit_param[i][j+1];
	  IJ.log( str );
	}
  	return (step >= max_step); 
  }
}

