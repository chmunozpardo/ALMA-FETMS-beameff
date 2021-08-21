/*******************************************************************************
* ALMA - Atacama Large Millimeter Array
* (c) Associated Universities Inc., 2015
*
*This library is free software; you can redistribute it and/or
*modify it under the terms of the GNU Lesser General Public
*License as published by the Free Software Foundation; either
*version 2.1 of the License, or (at your option) any later version.
*
*This library is distributed in the hope that it will be useful,
*but WITHOUT ANY WARRANTY; without even the implied warranty of
*MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*Lesser General Public License for more details.
*
*You should have received a copy of the GNU Lesser General Public
*License along with this library; if not, write to the Free Software
*Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
*
*
*/

#include "FitAmplitude.h"
#include "ScanData.h"
#include "ScanDataRaster.h"
#include <math.h>

extern "C" {
    #include "nr/nr.h"
}

namespace BeamFitting {

    // forward declare fitness functions:
    double function_amp(double p[]);
    void dfunction_amp(double p[], double df[]);

    const int nTerms_m(6);         ///< terms of amplitude fit model
    ScanData *fitAmpScan(NULL);    ///< Scan we are currently fitting against
    static float azActual;         ///< actual (center of mass) pointing in Az
    static float elActual;         ///< actual (center of mass) pointing in El

    void FitAmplitude(ScanData *currentScan, float _azActual, float _elActual) {
        fitAmpScan = currentScan;
        azActual = _azActual;
        elActual = _elActual;

        float ftol = pow(10, -1.0 * fitAmpScan -> getBand());
        int iter_amp;
        double fret_amp;
        double p[nTerms_m + 1];
        void (*dfunctionamp)(double p[], double df[]);
        double (*functionamp)(double p[]);

        functionamp = &function_amp;
        dfunctionamp = &dfunction_amp;

        //Initial guess values may be substituted later
        //p1=amp
        //p2=width (deg)
        //p3=u_offset (deg)
        //p4=v_offset (deg)
        //p5=D_0-90
        //p6=D_45-135
        p[0] = 0.0;
        p[1] = 1.0;
        p[2] = 3.0;
        p[3] = 0.01 * azActual;
        p[4] = 0.01 * elActual;
        p[5] = 0.0;
        p[6] = 0.0;

        frprmn(p, nTerms_m, ftol, &iter_amp, &fret_amp, functionamp, dfunctionamp);

        fitAmpScan -> setAmpFitResults(p[1], p[2], p[3], p[4], p[5], p[6]);
    }

    void dfunction_amp(double p[], double df[]) {
    /* This function should compute the gradients of the chi-squared function,
     * which are stored in the array "df". Since it is not a analytic function,
     * we must compute the partial derivatives numerically, which is done using:
     *    d(chi-square) / dp[j] = *chiSquare(p[j]+del) - chiSquare(p[j])) / del
     *
     */
        int i,j;
        double par[nTerms_m + 1];
        double delta = 0.01, del;

        const ScanDataRaster *pscan = fitAmpScan -> getFFScan();
        if (!pscan)
            return;

        for (j = 1; j <= nTerms_m; j++) {
            /* set all the parameters back to the current values */
            for (i = 1; i <= nTerms_m; i++) {
                par[i] = p[i];
            }
            /* apply a small offset to the parameter being adjusted this time through the loop */
            if (fabs(par[j]) > (delta / 100000)) {
                del = delta * par[j];
            } else {
                /* this takes care of the unique case where the initial guess is zero */
                del = delta;
            }
            par[j] += del;
            df[j] = ((pscan -> calcAmplitudeEfficiency(par, azActual, elActual))
                  -  (pscan -> calcAmplitudeEfficiency(p, azActual, elActual))) / del;
        }
    }


    double function_amp(double p[]) {
    /* This function should return the chi-squared value for the current model
     * parameters that are passed in as an argument.  In this case, we want
     * to minimize the phase "inefficiency", so we simply calculate that.
     */
        const ScanDataRaster *pscan = fitAmpScan -> getFFScan();
        if (!pscan)
            return 0.0;

        return (pscan -> calcAmplitudeEfficiency(p, azActual, elActual));
    }

}; //namespace
