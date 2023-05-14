public class FFTImageFiltering {

    public static int N = 256 ;

    public static void main(String [] args) throws Exception {

        // Calculate startTime
        long startTime = System.currentTimeMillis();


        double [] [] X = new double [N] [N] ;
        ReadPGM.read(X, "C:\\Users\\Kyle\\Documents\\University\\Masters\\Scientific Computing\\Labs\\wolf.pgm", N) ;

        DisplayDensity display =
                new DisplayDensity(X, N, "Original Image") ;

        // create array for in-place FFT, and copy original data to it
        double [] [] CRe = new double [N] [N], CIm = new double [N] [N] ;
        for(int k = 0 ; k < N ; k++) {
            for(int l = 0 ; l < N ; l++) {
                CRe [k] [l] = X [k] [l] ;
            }
        }

        fft2d(CRe, CIm, 1) ;  // Fourier transform

        Display2dFT display2 =
                new Display2dFT(CRe, CIm, N, "Discrete FT") ;

        

        // create array for in-place inverse FFT, and copy FT to it
        double [] [] reconRe = new double [N] [N],
                     reconIm = new double [N] [N] ;
        for(int k = 0 ; k < N ; k++) {
            for(int l = 0 ; l < N ; l++) {
                reconRe [k] [l] = CRe [k] [l] ;
                reconIm [k] [l] = CIm [k] [l] ;
            }
        }

        fft2d(reconRe, reconIm, -1) ;  // Inverse Fourier transform

        DisplayDensity display3 =
                new DisplayDensity(reconRe, N, "Reconstructed Image") ;

        // End benchmark
        long endTime = System.currentTimeMillis();

        System.out.println("Calculation completed in " +
                            (endTime - startTime) + " milliseconds");
    }

    static void transpose(double [] [] a) {

        int cutoff = N/8;

        for(int k = 0 ; k < N ; k++) {
            int kSigned = k <= N/2 ? k : k - N ;
            for(int l = 0 ; l < N ; l++) {
                int lSigned = l <= N/2 ? l : l - N ;
                if(Math.abs(kSigned) < cutoff || Math.abs(lSigned) < cutoff) {
                    a [k] [l] = 0 ;
                }
            }
        }

        for(int i = 0 ; i < cutoff ; i++) {
            for(int j = 0 ; j < i ; j++) {
                 a [i] [j] = a [j] [i];
            }
        }
    }

    static void fft2d(double [] [] re, double [] [] im, int isgn) {

        // For simplicity, assume square arrays

        for (int i = 0; i < N; i++){
            FFT.fft1d(re[i], im[i], isgn);
        }

        transpose(re) ;
        transpose(im) ;

        for (int j = 0; j < N; j++){
            FFT.fft1d(re[j], im[j], isgn);
        }

        transpose(re) ;
        transpose(im) ;
    }

}