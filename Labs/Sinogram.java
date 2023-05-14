import java.util.Arrays ;

import java.awt.* ;
import javax.swing.* ;


public class Sinogram {

    static final int N = 512 ;

    static final int CELL_SIZE = 1 ;

    static final double SCALE = 0.0045 ;  // think of a better way to
                                          // parametrize this later...

    static final int CUTOFF = N/4 ;  // in ramp filter

    static final float GREY_SCALE_LO = 0.95f, GREY_SCALE_HI = 1.05f ;
        // Clipping, for display only.  See for example Figure 1 in:
        //    http://bigwww.epfl.ch/thevenaz/shepplogan/

    public static void main(String [] args) {

        long startTime = System.currentTimeMillis();


        double [] [] density = new double [N] [N] ;

        for(int i = 0 ; i < N ; i++) {
            double x = SCALE * (i - N/2) ;
            for(int j = 0 ; j < N ; j++) {
                double y = SCALE * (j - N/2) ;

                density [i] [j] = sheppLoganPhantom(x, y) ;
            }
        } 

        DisplayDensity display1 =
                new DisplayDensity(density, N, "Source Model",
                                GREY_SCALE_LO, GREY_SCALE_HI) ;

        // Radon tranform of density (as measured by detectors):

        double [] [] sinogram = new double [N] [N] ;

        for(int iTheta = 0 ; iTheta < N ; iTheta++) {
            double theta = (Math.PI * iTheta) / N ;
            double cos = Math.cos(theta) ;
            double sin = Math.sin(theta) ;
            for(int iR = 0 ; iR < N ; iR++) {
                double r = SCALE * (iR - N/2) ;
                double sum = 0 ;
                for(int iS = 0 ; iS < N ; iS++) {
                    double s = SCALE * (iS - N/2) ;
                    double x = r * cos + s * sin ;
                    double y = r * sin - s * cos ;
                    sum += sheppLoganPhantom(x, y) ;
                }
                sinogram [iTheta] [iR] = sum ;
            }
        }

        // DisplayDensity display2 = new DisplayDensity(sinogram, N, "Sinogram") ;

        // inferred integral of density points (actually sum of density
        // points, here) for laternormalization of reconstruction

        double normDensity = norm1(sinogram [0]) ;


        // ... Insert sinogram filtering code here! ...
        // long startTime = System.currentTimeMillis();


        // double [] [] sinogramFTRe = new double [N] [N],
        //             sinogramFTIm = new double [N] [N] ;
        
        // double [] [] sinogramFTRe2 = new double [N] [N],
        //             sinogramFTIm2 = new double [N] [N] ;

        double [] [] sinogramFTRe3 = new double [N] [N],
                    sinogramFTIm3 = new double [N] [N] ; 

        for(int iTheta = 0 ; iTheta < N ; iTheta++) {
            for(int iR = 0 ; iR < N ; iR++) {
                // sinogramFTRe [iTheta] [iR] = sinogram [iTheta] [iR] ;
                // sinogramFTRe2 [iTheta] [iR] = sinogram [iTheta] [iR] ;
                sinogramFTRe3 [iTheta] [iR] = sinogram [iTheta] [iR] ;


            }
        }

        for(int iTheta = 0; iTheta < N; iTheta++) {
            // FFT.fft1d(sinogramFTRe[iTheta], sinogramFTIm[iTheta], 1);
            // FFT.fft1d(sinogramFTRe2[iTheta], sinogramFTIm2[iTheta], 1);
            FFT.fft1d(sinogramFTRe3[iTheta], sinogramFTIm3[iTheta], 1);


        }

        // DisplaySinogramFT display3 =
        //         new DisplaySinogramFT(sinogramFTRe, sinogramFTIm, N,
        //                             "Sinogram radial Fourier Transform") ;

        // // RAMP FILTER
        // for(int iTheta = 0 ; iTheta < N ; iTheta++) {
        //     for(int iK = 0 ; iK < N ; iK++) {
        //         int kSigned = iK <= N/2 ? iK : iK - N ;
        //         sinogramFTRe2[iTheta][iK] *= Math.abs(kSigned);
        //         sinogramFTIm2[iTheta][iK] *= Math.abs(kSigned);
        //     }
        // }

        // RAM LAK FILTER
        // for(int iTheta = 0 ; iTheta < N ; iTheta++) {
        //     for(int iK = 0 ; iK < N ; iK++) {
        //         int kSigned = iK <= N/2 ? iK : iK - N ;
        //         double ramlak = Math.abs(kSigned) > CUTOFF ? 0 : Math.abs(kSigned);
        //         sinogramFTRe[iTheta][iK] *= ramlak;
        //         sinogramFTIm[iTheta][iK] *= ramlak;
        //     }
        // }

        // LOW COSINE FILTER
        for(int iTheta = 0 ; iTheta < N ; iTheta++) {
            for(int iK = 0 ; iK < N ; iK++) {
                int kSigned = iK <= N/2 ? iK : iK - N ;
                sinogramFTRe3[iTheta][iK] *= Math.abs(kSigned);
                sinogramFTIm3[iTheta][iK] *= Math.abs(kSigned);
                sinogramFTIm3[iTheta][iK] *= Math.cos((Math.PI * kSigned) / 2 * CUTOFF);
                sinogramFTRe3[iTheta][iK] *= Math.cos((Math.PI * kSigned) / 2 * CUTOFF);
                if (Math.abs(kSigned) > CUTOFF) {
                    sinogramFTIm3[iTheta][iK] = 0;
                    sinogramFTRe3[iTheta][iK] = 0;
                }
            }
        }



        for(int iTheta = 0; iTheta < N; iTheta++) {
            // FFT.fft1d(sinogramFTRe[iTheta], sinogramFTIm[iTheta], -1);
            // FFT.fft1d(sinogramFTRe2[iTheta], sinogramFTIm2[iTheta], -1);
            FFT.fft1d(sinogramFTRe3[iTheta], sinogramFTIm3[iTheta], -1);


        }


        // DisplayDensity display4 =
        //         new DisplayDensity(sinogramFTRe, N, "Ram Lak Filtered sinogram") ;

        // DisplayDensity display9 =
        //         new DisplayDensity(sinogramFTRe2, N, "Ramp Filter Filtered sinogram") ;

        DisplayDensity display10 =
                new DisplayDensity(sinogramFTRe3, N, "Low-Pass Cosine Filtered sinogram") ;

        // double [] [] backProjection = new double [N] [N] ;
        // double [] [] backProjection2 = new double [N] [N] ;
        double [] [] backProjection3 = new double [N] [N] ;


        // backProject(backProjection, sinogramFTRe) ;
        // backProject(backProjection2, sinogramFTRe2) ;
        backProject(backProjection3, sinogramFTRe3) ;



        // Normalize reconstruction, to have same sum as inferred for
        // original density

        // double factor = normDensity / norm2(backProjection) ;
        // double factor2 = normDensity / norm2(backProjection2) ;
        double factor3 = normDensity / norm2(backProjection3) ;


        for(int i = 0 ; i < N ; i++) {
            for(int j = 0 ; j < N ; j++) {
                // backProjection [i] [j] *= factor ;
                // backProjection2 [i] [j] *= factor2 ;
                backProjection3 [i] [j] *= factor3 ;

            }
        }

        // DisplayDensity display5 =
        //         new DisplayDensity(backProjection, N,
        //                         "RAM-LAK FILTER", GREY_SCALE_LO, GREY_SCALE_HI) ;
                                
        // DisplayDensity display6 =
        //         new DisplayDensity(backProjection2, N,
        //                         "RAMP FILTER", GREY_SCALE_LO, GREY_SCALE_HI) ;

        DisplayDensity display7 =
                new DisplayDensity(backProjection3, N,
                                "LOW COSINE FILTER", GREY_SCALE_LO, GREY_SCALE_HI) ;
        
        long endTime = System.currentTimeMillis();

        System.out.println("Calculation completed in " +
                            (endTime - startTime) + " milliseconds");
    }

    static void backProject(double [] [] projection, double [] [] sinogram) {

        // Back Projection operation on sinogram

        for(int i = 0 ; i < N ; i++) {
            double x = SCALE * (i - N/2) ;
            for(int j = 0 ; j < N ; j++) {
                double y = SCALE * (j - N/2) ;

                double sum = 0 ;
                for(int iTheta = 0 ; iTheta < N ; iTheta++) {
                    double theta = (Math.PI * iTheta) / N ;
                    double cos = Math.cos(theta) ;
                    double sin = Math.sin(theta) ;

                    double r = x * cos + y * sin ;

                    double rBox = N/2 + r/SCALE ;

                    if(rBox < 0) continue ;  // assume centred object, with
                                             // support radius < N/2

                    int iR = (int) rBox ; 

                    double offR = rBox - iR ;
                    int iPlusR = iR + 1 ; 

                    if(iPlusR >= N) continue ;  // ditto.

                    // linear interpolation
                    double sinogramVal =
                            (1 - offR) * sinogram [iTheta] [iR] +
                            offR * sinogram [iTheta] [iPlusR] ;
                    sum += sinogramVal ;
                }
                projection [i] [j] = sum ;
            }
        }
    }

    // Shepp-Logan Phantom:
    //
    //   https://en.wikipedia.org/wiki/Shepp%E2%80%93Logan_phantom

    static final Ellipse [] sheppLoganEllipses = {
        new Ellipse(0.0, 0.0, 0.69, 0.92, 0, 2.0),
        new Ellipse(0.0, -0.0184, 0.6624, 0.874, 0, -0.98),
        new Ellipse(0.22, 0, 0.11, 0.31, -18.0, -0.02),
        new Ellipse(-0.22, 0, 0.16, 0.41, 18.0, -0.02),
        new Ellipse(0, 0.35, 0.21, 0.25, 0, 0.01),
        new Ellipse(0, 0.1, 0.046, 0.046, 0, 0.01),
        new Ellipse(0, -0.1, 0.046, 0.046, 0, 0.01),
        new Ellipse(-0.08, -0.605, 0.046, 0.023, 0, 0.01),
        new Ellipse(0, -0.605, 0.023, 0.023, 0, 0.01),
        new Ellipse(0.06, -0.605, 0.023, 0.046, 0, 0.01),
    } ;

    static double sheppLoganPhantom (double x, double y) {

        double total = 0 ;
        for(Ellipse ellipse : sheppLoganEllipses) {
            total += ellipse.localDensity(x, y) ;
        }
        return total ;
    }

    static class Ellipse {

        double centreX ;
        double centreY ;
        double major ;
        double minor ;
        double theta ;
        double density ;
        double cos, sin ;

        Ellipse(double centreX, double centreY,
                double major, double minor, double theta, double density) {

            this.centreX = centreX ;
            this.centreY = centreY ;
            this.major = major ;
            this.minor = minor ;
            this.theta = theta ;
            if(theta == 0) {
                cos = 1 ;
                sin = 0 ;
            }
            else {
                double rad = Math.PI * theta / 180 ;
                cos = Math.cos(rad) ; 
                sin = Math.sin(rad) ; 
            }
            this.density = density;
        }

        double localDensity(double x, double y) {

            double xOff, yOff ;
            xOff = x - centreX ;
            yOff = y - centreY ;

            double xRot, yRot ;
            if(theta == 0) {
                xRot = xOff ;
                yRot = yOff ;
            }
            else {
                // Rotate so x/y aligned with major/minor axes.
                xRot = cos * xOff - sin * yOff ;
                yRot = sin * xOff + cos * yOff ;
            }
            double xNorm = xRot / major ;
            double yNorm = yRot / minor ;
            if(xNorm * xNorm + yNorm * yNorm < 1) {
                return density ;
            }
            else {
                return 0 ;
            }
        }
    }

    static double norm1(double [] density) {

        double norm = 0 ;
        for(int i = 0 ; i < N ; i++) {
            norm += density [i] ;
        }
        return norm ;
    }

    static double norm2(double [] [] density) {

        double norm = 0 ;
        for(int i = 0 ; i < N ; i++) {
            for(int j = 0 ; j < N ; j++) {
                if(density [i] [j] > 0) {
                    norm += density [i] [j] ;
                }
            }
        }
        return norm ;
    }

}
