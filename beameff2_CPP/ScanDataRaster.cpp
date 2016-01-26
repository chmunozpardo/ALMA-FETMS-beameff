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

#include "ScanDataRaster.h"
#include <iostream>
#include <stdio.h>
#include <math.h>
using namespace std;

void ScanDataRaster::clear() {
    size_m = 0;
    startRow_m = 0;
    NSIDateTime_m = "";

    xArray_m.clear();
    yArray_m.clear();
    ampArray_m.clear();
    phiArray_m.clear();
    EArray_m.clear();
    RadiusArray_m.clear();
    MaskArray_m.clear();

    memset(&results_m, 0, sizeof(AnalyzeResults));
    results_m.maxAmp = -999;
}

bool ScanDataRaster::loadFile(const std::string filename, const std::string delim, bool rotate) {
    // Start with empty arrays:
    clear();

    // Open the file:
    FILE *f = fopen(filename.c_str(), "r");
    if (!f) {
        cout << "ERROR: ScanDataRaster::loadFile(): Could not open '" << filename << "'" << endl;
        return false;
    }

    // Make scan string for valid data rows:
    char scanStr[20];
    char *ptr;
    sprintf(scanStr, "%s%s%s%s%s%s%s", "%f", delim.c_str(), "%f", delim.c_str(), "%f", delim.c_str(), "%f");

    // Targets for sscanf:
    float x, y, amp, phi;

    // File is open.  Start reading rows:
    const int LINESIZE = 1024;
    char buf[LINESIZE + 1];
    unsigned rowIndex = 0;
    int scanRet = 0;

    // Loop to read all lines:
    while ((ptr = fgets(buf, LINESIZE, f)) != NULL) {
        rowIndex++;
        // check for NSI format file marker indicating start of data:
        if (!startRow_m && strstr(ptr, "line:")) {
            startRow_m = rowIndex;      // data starts on 0-indexed row of the file.

        // check for NSI format date/time string:
        } else if (NSIDateTime_m.empty() && strstr(ptr, "date/time:")) {
            // Parse the datetime from the listing file:
            ptr = strtok(ptr, ":");
            ptr = strtok(NULL, ",");
            NSIDateTime_m = ptr;        // save the NSI date/time string.

        // normal data line or header line:
        } else {
            scanRet = sscanf(ptr, scanStr, &x, &y, &amp, &phi);
            if (rotate) {
                x = -x;
                y = -y;
                phi = -phi;
            }

            // Only store rows with exactly four values succefully parsed:
            if (scanRet == 4) {
                xArray_m.push_back(x);
                yArray_m.push_back(y);
                ampArray_m.push_back(amp);
                phiArray_m.push_back(phi);
                size_m++;

            // Report error only if we are in the body of the file:
            } else if (startRow_m || size_m) {
                cout << "ERROR: ScanDataRaster::loadFile(): sscanf returned " << scanRet << " on line " << rowIndex << "." << endl;
                return false;
            }
        }
    }
    return true;
}

bool ScanDataRaster::saveFile(const std::string filename, float copolPeakAmp,
                              const std::string &delim, const std::string &lineterm) const
{
    // delete the file if it already exists:
    remove(filename.c_str());

    // open the output file:
    FILE *f = fopen(filename.c_str(), "w");
    if (!f) {
        cout << "ERROR: ScanDataRaster::saveFile(): Could not open '" << filename << "'" << endl;
        return false;
    }

    float x, y, amp, phase; // individual pixels in file
    float lastX(0), lastY(0);     // use these to determine if the file is X-major or Y-major.
    char textline[500];

    // iterators for main loop:
    vector<float>::const_iterator itX = xArray_m.begin();
    vector<float>::const_iterator itY = yArray_m.begin();
    vector<float>::const_iterator itAmp = ampArray_m.begin();
    vector<float>::const_iterator itPhi = phiArray_m.begin();

    // to detect first and second time in loop:
    bool first = true;
    bool second = false;
    bool Xmajor = true;

    while (itX != xArray_m.end()) {
        x = *itX;
        y = *itY;
        amp = *itAmp - copolPeakAmp;    // normalize to the copol peak
        phase = *itPhi;

        // if first iteration, save the last value of x and y
        if (first) {
            lastX = x;
            lastY = y;
            first = false;
            second = true;

        // if second iteration: figure out if X-major or Y-major:
        } else if (second) {
            if (x != lastX)
                Xmajor = false;
            else if (y != lastY)
                Xmajor = true;
            second = false;

        // if the major axis has incremented, put an extra newline in the file:
        } else if ((Xmajor && (x != lastX)) || (!Xmajor && (y != lastY))) {
            lastX = x;
            lastY = y;
            fputs(lineterm.c_str(), f);
        }

        // write a line:
        sprintf(textline, "%f%s%f%s%f%s%f%s",
                          x, delim.c_str(),
                          y, delim.c_str(),
                          amp, delim.c_str(),
                          phase, lineterm.c_str());
        fputs(textline, f);

        // increment iterators:
        itX++;
        itY++;
        itAmp++;
        itPhi++;
    }
    fclose(f);
    return true;
}

void ScanDataRaster::calcStepSize() {
    // try just looking at the first two elements:
    results_m.xStepSize = fabs(xArray_m[1] - xArray_m[0]);

    // if that didn't work, loop to find the first jump, if any:
    if (results_m.xStepSize == 0.0) {
        vector<float>::const_iterator it = xArray_m.begin();
        float x0 = xArray_m[0];
        while (results_m.xStepSize == 0.0 && it != xArray_m.end()) {
            results_m.xStepSize = fabs(*it - x0);
            it++;
        }
    }

    // try just looking at the first two elements:
    results_m.yStepSize = fabs(yArray_m[1] - yArray_m[0]);

    // if that didn't work, loop to find the first jump, if any:
    if (results_m.yStepSize == 0.0) {
        vector<float>::const_iterator it = yArray_m.begin();
        float y0 = yArray_m[0];
        while (results_m.yStepSize == 0.0 && it != yArray_m.end()) {
            results_m.yStepSize = fabs(*it - y0);
            it++;
        }
    }
}

void ScanDataRaster::calcPeakAndPhase() {
    // reset peak variables:
    results_m.maxAmp = -999;
    results_m.phaseAtPeak = 0;

    // loop to find peak power and phase at peak:
    vector<float>::const_iterator itAmp = ampArray_m.begin();
    vector<float>::const_iterator itPhi = phiArray_m.begin();
    float amp, phi;

    while (itAmp != ampArray_m.end()) {
        amp = *itAmp;
        phi = *itPhi;
        if (amp > results_m.maxAmp) {
            results_m.maxAmp = amp;
            results_m.phaseAtPeak = phi;
        }
        itAmp++;
        itPhi++;
    }
}

void ScanDataRaster::calcCenterOfMass() {
    results_m.xCenterOfMass = 0;
    results_m.yCenterOfMass = 0;

    // loop to accumulate center of mass:
    vector<float>::const_iterator itX = xArray_m.begin();
    vector<float>::const_iterator itY = yArray_m.begin();
    vector<float>::const_iterator itAmp = ampArray_m.begin();

    float threshold = results_m.maxAmp - 10;  // points above peak - 10 dB.
    float x, y, amp;
    double sumXAmp = 0.0;
    double sumYAmp = 0.0;
    double sumAmp = 0.0;

    while (itX != xArray_m.end()) {
        x = *itX;
        y = *itY;
        amp = *itAmp;

        // Accumulate center of mass:
        if (amp > threshold) {
            // convert to linear units:
            amp = pow(10.0, amp / 10.0);
            // accumulate X and Y weighted sums:
            sumXAmp += x * amp;
            sumYAmp += y * amp;
            // accumulate overall sum:
            sumAmp += amp;
        }
        itX++;
        itY++;
        itAmp++;
    }
    results_m.xCenterOfMass = sumXAmp / sumAmp;
    results_m.yCenterOfMass = sumYAmp / sumAmp;
}

void ScanDataRaster::subtractForAttenuator(float ifAttenDiff) {
   vector<float>::iterator itAmp = ampArray_m.begin();
   while (itAmp != ampArray_m.end()) {
       (*itAmp) -= ifAttenDiff;
       itAmp++;
   }
}

void ScanDataRaster::analyzeBeam(float azNominal, float elNominal, float subreflectorRadius, float copolPeakAmp) {

    // values to use for computing the mask:
    float inner = subreflectorRadius - (getStepSize() / 2.0);
    float outer = subreflectorRadius + (getStepSize() / 2.0);

    // variables used in the loop:
    float x, y, azOnSub, elOnSub, amp, E, radius, mask, EOnSec, powerOnSec;

    // empty the output arrays:
    EArray_m.clear();
    MaskArray_m.clear();
    RadiusArray_m.clear();

    // reset the beam statistics:
    results_m.sumMask = 0.0;
    results_m.sumE = 0.0;
    results_m.sumSqE = 0.0;
    results_m.sumEOnSec = 0.0;
    results_m.sumSqEOnSec = 0.0;
    results_m.sumPowerOnSec = 0.0;
    results_m.sumSqPowerOnSec = 0.0;
    results_m.sumEOnEdge = 0.0;
    results_m.edgeNumPts = 0;

    // loop to calculate E, mask, etc.:
    vector<float>::const_iterator itX = xArray_m.begin();
    vector<float>::const_iterator itY = yArray_m.begin();
    vector<float>::const_iterator itAmp = ampArray_m.begin();

    while (itX != xArray_m.end()) {
        x = *itX;
        y = *itY;
        amp = *itAmp;

        // Compute E: array of electric field voltages from 2D pattern, normalized to peak = 1.0 of copol scan:
        E = pow(10.0, (amp - copolPeakAmp) / 20.0);
        EArray_m.push_back(E);

        // angles relative to nominal pointing:
        azOnSub = x - azNominal;
        elOnSub = y - elNominal;

        // radius angle relative to secondary reflector center:
        radius = sqrt(pow(azOnSub, 2.0) + pow(elOnSub, 2.0));
        RadiusArray_m.push_back(radius);

        // Compute mask of secondary reflector:
        if (radius > outer) {
            // zero outside the subreflectorRadius angle
            mask = 0.0;

        } else if (radius < inner) {
            // one inside the subreflector
            mask = 1.0;

        } else {
            // at points on the edge, a linear taper between 0 and 1:
            mask = (outer - radius) / getStepSize();
        }
        MaskArray_m.push_back(mask);

        // accumulate sum of the mask:
        results_m.sumMask += mask;

        // accumulate sum and sum of squares of the electric field voltage:
        results_m.sumE += E;
        results_m.sumSqE += pow(E, 2.0);

        // accumulate sum and sum of squares of voltage on the secondary:
        EOnSec = mask * E;
        results_m.sumEOnSec += EOnSec;
        results_m.sumSqEOnSec += pow(EOnSec, 2.0);

        // accumulate sum and sum of squares of power on the secondary:
        powerOnSec = mask * pow(E, 2.0);
        results_m.sumPowerOnSec += powerOnSec;
        results_m.sumSqPowerOnSec += pow(powerOnSec, 2.0);

        // accumulate the sum and count of voltages in the subreflector edge region:
        if (pow(radius - subreflectorRadius, 2.0) < pow(getStepSize(), 2.0)) {
            results_m.sumEOnEdge += E;
            results_m.edgeNumPts++;
        }
        itX++;
        itY++;
        itAmp++;
    }
}

float ScanDataRaster::calcPhaseEfficiency(float p[], float azNominal, float elNominal) const {
    float x, y, radiusSq, E, mask, phaseRad;                        // values from data arrays
    float phi_fit, phi_err, eta_phase;                              // calculated values
    double costerm(0.0), sinterm(0.0), normalizationFactor(0.0);    // accumulate fit errors in the loop

    unsigned i;
    for (i = 0; i < size_m; i++) {
        x = xArray_m[i] - azNominal;
        y = yArray_m[i] - elNominal;
        radiusSq = pow(RadiusArray_m[i], 2.0);
        E = EArray_m[i];
        mask = MaskArray_m[i];
        phaseRad = phiArray_m[i] * M_PI / 180.0;

        phi_fit = p[1] * x + p[2] * y + p[3] * radiusSq / 2.0;
        phi_err = phaseRad - phi_fit;
        costerm += mask * E * cos(phi_err);
        sinterm += mask * E * sin(phi_err);
        normalizationFactor += mask * E;
    }
    eta_phase = (pow(costerm, 2.0) + pow(sinterm, 2.0)) / pow(normalizationFactor, 2.0);
    return eta_phase;
}

float ScanDataRaster::calcAmplitudeEfficiency(float p[], float azActual, float elActual) const {
    float amp_mod_term, amp_diff_term;
    float x, y, E, mask;
    double amp_fit = 0.0;

    unsigned i;
    for (i = 0; i < size_m; i++) {
        x = xArray_m[i] - azActual;
        y = yArray_m[i] - elActual;
        E = EArray_m[i];
        mask = MaskArray_m[i];

        if (mask > 0.0) {
            amp_mod_term = p[1] * (
                1 - pow(
                    (
                        1 - exp(
                            -(1 / pow(p[2], 2.0)) * (
                                pow((x - p[3]), 2.0) +
                                pow((y - p[4]), 2.0) -
                                p[5] * (pow((x - p[3]), 2.0) -
                                pow((y - p[4]), 2.0)) -
                                2.0 * p[6] * (x - p[3]) * (y - p[4])
                            )
                        )
                    ), 2.0
                )
            );

            amp_diff_term = mask * (E - amp_mod_term);
            amp_fit += pow(amp_diff_term, 2.0);
        }
    }
    return static_cast<float>(amp_fit);
}

void ScanDataRaster::print(int _indent, unsigned headRows, unsigned tailRows) const {
    string indent(_indent, ' ');

    cout << indent << "size_m = " << size_m << endl;
    cout << indent << "startRow_m = " << startRow_m << endl;
    cout << indent << "NSIDateTime_m = '" << NSIDateTime_m << "'" << endl;

    results_m.print(_indent);

    if (size_m) {
        cout << indent << "  x, y, amp, phi:" << endl;
        unsigned index;
        for (index = 0; index < headRows; index++) {
            cout << indent << "  ";
            printRow(index);
        }
        cout << indent << "  ..." << endl;
        for (index = size_m - tailRows; index < size_m; index++) {
            cout << indent << "  ";
            printRow(index);
        }
    }
}

void ScanDataRaster::printRow(unsigned index) const {
    cout << xArray_m[index] << ", " << yArray_m[index] << ", "
         << ampArray_m[index] << ", " << phiArray_m[index] << endl;
}



