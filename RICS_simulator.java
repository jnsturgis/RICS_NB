/**
*	@file RICS_simulator.java
*
*	This file implments the RICS_Simulator class.
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
*	@brief A plugin class for ImageJ to simulate the results of a RICS experiment.
*
*	This class is part of my implementation of RICS analysis for ImageJ. Using a
*	description of a microscope and sample it will produce a stack of images
*	that could result from such an experiment.
*
*	@todo	Improve code commenting and documentation.
*/
public class RICS_simulator implements PlugInFilter {

  private ImagePlus image; 			/**< Place to store the image initialized in
  										setup(). 									*/
  private	double	Conc = 1.0 ;	/**< Object concentration in nM (default 1.0)	*/
  private	double	Bril = 1.0;		/**< Object brilliance at focus in Mcps
  										(Default 1.0)								*/
  private	double	Diff = 0.0;		/**< Object diffusion coefficient in \f$\mu m^{2}
                                        sec^{-1}\f$ (default 0.0).					*/

  private	double	W0   = 0.400;	/**< Beam half waist \f$\mu m\f$ (default 0.4).	*/
  private	double	Z0   = 1.200;	/**< Beam half height \f$\mu m\f$ (default 1.2).*/

  private	double	Dtx  = 8.0;		/**< Delay between pixels in \f$\mu sec\f$
  										(default 8.0).								*/
  private	double	Dty  = 4.0;		/**< Delay between lines in msec (default 4.0).	*/
  private	double	Dx   = 0.044;	/**< Pixel size in \f$\mu m\f$ (default 0.044).	*/

  private	int	Nx   = 128;			/**< Image width (default 128.					*/
  private	int	Ny   = 128;			/**< Image height (default 128).				*/
  private	int	Nt   = 16;			/**< Number of frames (default 16).				*/

  private	double	Noise= 0.0;		/**< Detector noise cps (default 0.0).			*/

/**
*
*	@brief	Initialization function.
*
*	This function is used by ImageJ to determine the types of image that can be
* treated.
*	@param arg String passed from ImageJ
*	@param image The image that is concerned, this is saved.
*	@return NO_IMAGE_REQUIRED Flags determining the types of object treated.
*
*/
  public int setup(String arg, ImagePlus image) {
    this.image = image;
    return NO_IMAGE_REQUIRED; //We dont care about the current image
  }


/**
*
*	@brief The main function of the RICS_simulator plugin.
*
*	This function first calls the dialog box to allow the user to describe the
*	simuation to perform and set the different parameters, then various derived
*	parameters are calculated and space allocated, then in a loop over the images
*	the different images are generated, finally the resulting stack is shown.
*
*	@todo improve the code in the loop and verify that the time step is right.
*/

  public void run(ImageProcessor ip) {

	/*************************************
	**									                **
	**		User input into dialog.		    **
	**									                **
	*************************************/

    GenericDialog	gd = new GenericDialog("RICS simulator");
    RICS_dialog( gd );
    if( gd.wasCanceled() ) return;

	/*************************************
	**									                **
	**		Initialization of derived	    **
	**			parameter values.		        **
	**									                **
	*************************************/

    int 	Nz   = (int) Math.floor(1.5 * Z0 / Dx);
    										                // Simulation box height in pixels
    double	Nav  = Conc * 0.6022 * Nx * Ny * 2.0 * Nz * Dx * Dx * Dx;
    										                // Number of particles in box
    double	w2 = (W0 / Dx)*(W0/Dx);			// Variance of PSF in grid units
    double	z2 = (Z0 / Dx)*(Z0/Dx);			// and in the z axis
    double	prob = Bril * Dtx;				  // Average intensity per pixel with an
    										                // object in the focus in photons.
    double	step_per_pix = Dtx * 6E-6 * Diff / (Dx*Dx);
    										                // MSD in grid units during acquisition
    										                // of a pixel
    double line_delay = (Dty * 1000.0)/Dtx - Nx;
      										              // How many pixel dwell times between lines?
    ImageStack mystack = new ImageStack(Nx, Ny);
    										                // Make a new stack for the data.

	/*************************************
	**									                **
	**		Calculation loop.			        **
	**									                **
	*************************************/

    for (int f=0;f<Nt;f++) {				    // Loop over images;
    	IJ.showProgress(f,Nt);
    	ImageProcessor theslice = new ShortProcessor(Nx,Ny);
    	int Npart = Math.round((float)(Nav + gauss_dev() * Math.sqrt(Nav)));
											                  // How many particles in this frame ?
		IJ.log("Window "+f+"has"+Npart+" particles.");

		for( int p = 0; p < Npart; p++ ){	  // Loop over particles...
											                  // put particle 'p' at a random location
											                  // in the simulation box.
			int	x = (int) Math.floor(lin_dev( Nx ));
			int	y = (int) Math.floor(lin_dev( Ny ));
			int	z = (int) Math.floor(lin_dev( 2.0 * Nz ) - Nz);

			for( int i = 0; i < Nx; i++ ){	  // Add images of this particle to frame
	  	  for( int j = 0; j < Ny; j++ ){
										                   	// Calculate probability of excitation
										                  	// for the object
	    		double dx2 = ((x-i)*(x-i))/w2;
	    		double dy2 = ((y-j)*(y-j))/w2;
	    		double dz2 = ( z * z )/z2;
	    		double irf = Math.exp(-(dx2+dy2+dz2)*2.0);
	    	                								// Add photons in pixel
	    		int	Npix = poisson_dev( irf * prob );
	    		theslice.putPixel(j,i,theslice.getPixel(j,i)+Npix);

	    		int Nstep = poisson_dev( step_per_pix );
	    		for( int l = 0; l < Nstep; l++){		// Now move the object
	      		double R_num = lin_dev(1.0);
	      		if      ( R_num < 1.0/6.0 ) { x += 1; }	// Choose a direction
	      		else if ( R_num < 2.0/6.0 ) { x -= 1; }
	      		else if ( R_num < 3.0/6.0 ) { y += 1; }
	      		else if ( R_num < 4.0/6.0 ) { y -= 1; }
	      		else if ( R_num < 5.0/6.0 ) { z += 1; }
	      		else                        { z -= 1; }
	      		if( x < 0   ) x += Nx;	// Periodic conditions
	      		if( y < 0   ) y += Ny;
	      		if( x > Nx ) x -= Nx;
	      		if( y > Ny ) y -= Ny;
	      		if( z < Nz  ) z += 2 * Nz;
	      		if( z > Nz  ) z -= 2 * Nz;
	    		}
	  		}
	  		if( line_delay > 0.0 ){
	    		int Nstep = poisson_dev( step_per_pix * line_delay );
	    		for( int l = 0; l < Nstep; l++){		// Now move the object
	      		double R_num = lin_dev(1.0);
	      		if      ( R_num < 1.0/6.0 ) { x += 1; }	// Choose a direction
	      		else if ( R_num < 2.0/6.0 ) { x -= 1; }
	      		else if ( R_num < 3.0/6.0 ) { y += 1; }
	      		else if ( R_num < 4.0/6.0 ) { y -= 1; }
	      		else if ( R_num < 5.0/6.0 ) { z += 1; }
	      		else                        { z -= 1; }
	      		if( x < 0   ) x += Nx;	// Periodic conditions
	      		if( y < 0   ) y += Ny;
	      		if( x > Nx ) x -= Nx;
	      		if( y > Ny ) y -= Ny;
	      		if( z < Nz  ) z += 2 * Nz;
	      		if( z > Nz  ) z -= 2 * Nz;
	    		}
	  		}
			}
    }
    for( int i = 0; i < Nx; i++ ){			// Add noise
		  for( int j = 0; j < Ny; j++ ){
	  		double av = Noise * Dtx / 1000000.0;
	  		int	Npix = poisson_dev( av );
	  		theslice.putPixel(j,i,theslice.getPixel(j,i)+Npix);
			}
    }
    mystack.addSlice("", theslice );		//Add processor into the stack and show it
  }

	/*************************************
	**									                **
	**	  Display of result and exit.   **
	**									                **
	*************************************/

  image = new ImagePlus("RICS simulation", mystack );
  image.show();
  image.updateAndDraw();
}

/**
*	@brief Display and handle data input dialog.
*/

private void RICS_dialog( GenericDialog gd )
{

  gd.addMessage("Object properties");
  gd.addNumericField("Concentration :", Conc, 3, 7, "nM" );
  gd.addNumericField("Molecular brightness :", Bril, 0, 7, "Mcps" );
  gd.addNumericField("Diffusion coef. :", Diff, 1, 7, "um²/sec");

  gd.addMessage("Microscope properties");
  gd.addNumericField("Beam waist :", W0, 3, 7, "um" );
  gd.addNumericField("Beam height :", Z0, 3, 7, "um" );

  gd.addNumericField("Pixel size :", Dx, 3, 7, "um" );

  gd.addNumericField("Pixel time :", Dtx, 0, 7, "usec" );
  gd.addNumericField("Line time :", Dty, 0, 7, "msec"  );

  gd.addNumericField("Image width :", Nx, 0, 4, "pixels"  );
  gd.addNumericField("Image Height :", Ny, 0, 4, "pixels" );
  gd.addNumericField("Stack size :", Nt, 0, 4, "frames" );

  gd.addNumericField("Detector noise :", Noise, 1, 5, "cps" );

  gd.showDialog();

  if( gd.wasCanceled() ) return;

  Conc  = (double) gd.getNextNumber();
  Bril  = (double) gd.getNextNumber();
  Diff  = (double) gd.getNextNumber();
  W0    = (double) gd.getNextNumber();
  Z0    = (double) gd.getNextNumber();
  Dx    = (double) gd.getNextNumber();
  Dtx   = (double) gd.getNextNumber();
  Dty   = (double) gd.getNextNumber();
  Nx    =    (int) gd.getNextNumber();
  Ny    =    (int) gd.getNextNumber();
  Nt    =    (int) gd.getNextNumber();
  Noise = (double) gd.getNextNumber();

  return;
}

private	double cof[] = {76.18009172947146, -86.50532032941677,
	24.01409824083091, -1.231739572450155,
	0.1208650973866179e-2, -0.5395239384953e-5};

/**
*	@brief Define the ln gamma function
*
*	Based on numerical Recipies.
*/

	private double	gammln( double xx )
	{
		double	x, y, tmp, ser;
		int	j;

		y = x = xx;
		tmp = x + 5.5;
		tmp -= (x+0.5)*Math.log(tmp);
		ser = 1.000000000190015;
		for(j=0;j<6;j++) ser += cof[j]/++y;
		return -tmp+Math.log(2.5066282746310005*ser/x);
	}

	/**
	*	@brief return a random number between 0 and max.
	*
	*	@param double max the largest number that will be returned.
	*
	*/
	private	double	lin_dev( double max )
	{
		return max * Math.random();
	}

	private	int	Spare = 0;				/*< Is there a spare value in for gauss_dev()? */
	private	double	SValue = 0.0;		/*< The spare value to return. */

	/**
	*	@brief Calculate a gaussian deviate (mean 0.0, standard deviation 1.0).
	*
	*	This code is adapted from Numerical Recipies and uses the Box-Muller method.
	*/
	private	double	gauss_dev()
	{
		double	result;
		double	V1, V2, r;

		if( Spare == 0 ){
			do {
				V1 = lin_dev(2.0) - 1.0;
				V2 = lin_dev(2.0) - 1.0;
				r = V1 * V1 + V2 * V2;
			} while( r > 1.0 );
			double fac = Math.sqrt( - 2.0 * Math.log(r)/r );
			SValue = V1 * fac;
			result = V2 * fac;
			Spare = 1;
		} else {
			result = SValue;
			Spare = 0;
		}
		return result;
	}

	private double	sq, alxm, g, oldm = (-1.0);

	/**
	*	@brief Calculate a possion deviate.
	*
	*	@param mean the desired average number returned.
	*
	*	This code is adapted from Numerical Recipies.
	*/
	private	int	poisson_dev( double mean )
	{
		double		em, t, y;

		if( mean < 12.0 ){
			if( mean != oldm ){
				oldm = mean;
				g = Math.exp( -mean );
			}
			em = -1.0;
			t = 1.0;
			do {
				++em;
				t *= lin_dev(1.0);
			} while ( t > g );
		} else {
			if( mean != oldm ){
				oldm = mean;
				sq = Math.sqrt( 2.0 * mean );
				alxm = Math.log( mean );
				g = mean * alxm - gammln( mean + 1.0 );
			}
			do {
				do {
					y = Math.tan( Math.PI * lin_dev( 1.0 ));
					em = sq * y + mean;
				} while ( em < 0.0 );
				em = Math.floor (em );
				t = 0.9 * (1.0+y*y)*Math.exp(em*alxm-gammln(em+1.0)-g);
			} while( lin_dev(1.0) > t );
		}
		return (int) em;
	}
}
