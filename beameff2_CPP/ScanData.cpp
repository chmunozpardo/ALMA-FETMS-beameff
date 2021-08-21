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

#include "ScanData.h"
#include "ScanDataRaster.h"
#include "ALMAConstants.h"
#include "iniparser.h"
#include <math.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

ScanData::ScanData()
  : NF_m(NULL),
    FF_m(NULL),
    NF2_m(NULL),
    FF2_m(NULL)
{
    clear();
}

ScanData::~ScanData() {
    clear();
}

void ScanData::clear() {
    freeArrays();
    inputSection_m.clear();
    measDateTime_m.clear();
    NSIDateTime_m.clear();
    filenameNF_m.clear();
    filenameFF_m.clear();
    filenameNF2_m.clear();
    filenameFF2_m.clear();

    scanId_m = 0;
    scanType_m = UNDEFINED;
    rfGHz_m = 0.0;
    pol_m = sb_m = -1;
    tilt_m = 0;
    nominalFocusZ_m = 0;
    ifAtten_m = 0;

    startrowNF_m = startrowFF_m = startrowNF2_m = startrowFF2_m = 0;
}

void ScanData::freeArrays() {
    delete NF_m;
    delete FF_m;
    delete NF2_m;
    delete FF2_m;
    NF_m = FF_m = NF2_m = FF2_m = NULL;
}

bool ScanData::loadFromIni(const dictionary *dict, const std::string inputSection) {

    inputSection_m = inputSection;

    string sectionKey;
    string value;

    // band_m
    sectionKey = inputSection_m;
    sectionKey += ":band";
    band_m = iniparser_getint(dict, sectionKey.c_str(), 0);

    // scanId_m
    sectionKey = inputSection_m;
    sectionKey += ":scan_id";
    scanId_m = iniparser_getint(dict, sectionKey.c_str(), 0);

    // measDateTime_m
    sectionKey = inputSection_m;
    sectionKey += ":ts";
    measDateTime_m = iniparser_getstring(dict, sectionKey.c_str(), "");

    // scanType_m
    sectionKey = inputSection_m;
    sectionKey += ":type";
    value = iniparser_getstring(dict, sectionKey.c_str(), "");

    if (value == "copol" || value == "COPOL")
        scanType_m = ScanData::COPOL;
    else if (value == "xpol" || value == "XPOL")
        scanType_m = ScanData::XPOL;
    else if (value == "dual" || value == "DUAL")
        scanType_m = ScanData::DUAL;
    else if (value == "copol180" || value == "COPOL180")
        scanType_m = ScanData::COPOL180;
    else
        scanType_m = ScanData::UNDEFINED;

    //is it a 4545 scan?
    sectionKey = inputSection_m;
    sectionKey += ":4545_scan";
    value = iniparser_getstring(dict, sectionKey.c_str(), "0");

    if (value == "1" || value[0] == 'T' || value[0] == 't' || value[0] == 'Y' || value[0] == 'y')
        scanType_m = ScanData::DUAL;

    //pol_m
    sectionKey = inputSection_m;
    sectionKey += ":pol";
    pol_m = iniparser_getint(dict, sectionKey.c_str(), 0);
    if (pol_m < 0 || pol_m > 1)
        pol_m = -1; // undefined

    //rfGHz_m
    sectionKey = inputSection_m;
    sectionKey += ":f";
    rfGHz_m = iniparser_getdouble(dict, sectionKey.c_str(), 0.0);

    //sb_m
    sectionKey = inputSection_m;
    sectionKey += ":sb";
    sb_m = iniparser_getint(dict, sectionKey.c_str(), 0);
    if (sb_m < 1 || sb_m > 2) {
        sb_m = 2;
        cout << "ScanData::loadFromIni('" << inputSection << "'): sb must be 1 or 2.  Assuming 2=LSB." << endl;
    }

    //tilt_m
    sectionKey = inputSection_m;
    sectionKey += ":tilt";
    tilt_m = iniparser_getint(dict, sectionKey.c_str(), 0);

    //probe Z distance
    sectionKey = inputSection_m;
    sectionKey += ":zdistance";
    nominalFocusZ_m = static_cast<float>(iniparser_getdouble(dict, sectionKey.c_str(), 0.0));
    //If zero or not provided, assume the FETMS default of 200 mm:
    if (nominalFocusZ_m == 0.0) {
        nominalFocusZ_m = -200.0;
        cout << "zdistance not specified in input. Using " << nominalFocusZ_m << " mm<br>" << endl;
    }
    //Start the nominal focus is negative relative to the probe (but may get inverted in loadListings):
    nominalFocusZ_m = -nominalFocusZ_m;

    //ifAtten_m
    sectionKey = inputSection_m;
    sectionKey += ":ifatten";
    ifAtten_m = static_cast<float>(iniparser_getdouble(dict, sectionKey.c_str(), 0.0));

    //filenameNF_m
    sectionKey = inputSection_m;
    sectionKey += ":nf";
    filenameNF_m = iniparser_getstring(dict, sectionKey.c_str(), "");

    //filenameFF_m
    sectionKey = inputSection_m;
    sectionKey += ":ff";
    filenameFF_m = iniparser_getstring(dict, sectionKey.c_str(), "");

    //filenameNF2_m
    sectionKey = inputSection_m;
    sectionKey += ":nf2";
    filenameNF2_m = iniparser_getstring(dict, sectionKey.c_str(), "");

    //filenameFF2_m
    sectionKey = inputSection_m;
    sectionKey += ":ff2";
    filenameFF2_m = iniparser_getstring(dict, sectionKey.c_str(), "");

    //startrowNF_m
    sectionKey = inputSection_m;
    sectionKey += ":nf_startrow";
    startrowNF_m = static_cast<unsigned>(iniparser_getint(dict, sectionKey.c_str(), 0));

    //startrowFF_m
    sectionKey = inputSection_m;
    sectionKey += ":ff_startrow";
    startrowNF_m = static_cast<unsigned>(iniparser_getint(dict, sectionKey.c_str(), 0));

    //startrowNF2_m
    sectionKey = inputSection_m;
    sectionKey += ":nf2_startrow";
    startrowNF_m = static_cast<unsigned>(iniparser_getint(dict, sectionKey.c_str(), 0));

    //startrowFF2_m
    sectionKey = inputSection_m;
    sectionKey += ":ff2_startrow";
    startrowNF_m = static_cast<unsigned>(iniparser_getint(dict, sectionKey.c_str(), 0));

    if (band_m < 1 || band_m > 10) {
        cout << "ScanData::loadFromIni('" << inputSection << "'): band must be in 1-10." << endl;
        return false;
    }
    if (scanType_m == ScanData::UNDEFINED) {
        cout << "ScanData::loadFromIni('" << inputSection << "'): scanType is UNDEFINED." << endl;
        return false;
    }
    if (pol_m == -1) {
        cout << "ScanData::loadFromIni('" << inputSection << "'): pol must be 0 or 1." << endl;
        return false;
    }
    return true;
}

bool ScanData::loadListings(const std::string &delim, ALMAConstants::InvertPhaseOptions invertPhaseOption) {
    // make sure there is not already raster data allocated:
    freeArrays();
    const bool rotateNF = false;    // never rotate NF scans
    bool rotateFF = false;          // switch to rotate FF scans
    bool invertPhase = false;       // switch to invert phase

    // determine whether to invert phase and rotate FF scans:
    switch (invertPhaseOption) {
    case ALMAConstants::INVERT_LSB:
        rotateFF = invertPhase = (sb_m == 2);
        break;

    case ALMAConstants::INVERT_USB:
        rotateFF = invertPhase = (sb_m == 1);
        break;

    case ALMAConstants::INVERT_ALL:
        rotateFF = invertPhase = true;
        break;

    case ALMAConstants::INVERT_NONE:
    default:
        break;
    }

    if (rotateFF)
        cout << "Pol " << pol_m << " " << getScanTypeString() << ": rotating scans.<br>" << endl;

    if (invertPhase)
        // we now expect the nominal focus to be at positive Z:
        nominalFocusZ_m = -nominalFocusZ_m;

    // create raster objects and load each file which is specified:
    if (!filenameNF_m.empty()) {
        NF_m = new ScanDataRaster;
        if (!NF_m -> loadFile(filenameNF_m, delim, rotateNF, invertPhase))
            return false;
        startrowNF_m = NF_m -> getStartRow();
        NSIDateTime_m = NF_m -> getNSIDateTime();
    }
    if (!filenameFF_m.empty()) {
        FF_m = new ScanDataRaster;
        if (!FF_m -> loadFile(filenameFF_m, delim, rotateFF, invertPhase))
            return false;
        startrowFF_m = FF_m -> getStartRow();
    }
    if (!filenameNF2_m.empty()) {
        NF2_m = new ScanDataRaster;
        if (!NF2_m -> loadFile(filenameNF2_m, delim, rotateNF, invertPhase))
            return false;
        startrowNF2_m = NF2_m -> getStartRow();
    }
    if (!filenameFF2_m.empty()) {
        FF2_m = new ScanDataRaster;
        if (!FF2_m -> loadFile(filenameFF2_m, delim, rotateFF, invertPhase))
            return false;
        startrowFF2_m = FF2_m -> getStartRow();
    }
    return true;
}

void ScanData::calcCenterOfMass() {
    if (FF_m)
        FF_m -> calcCenterOfMass();
}

void ScanData::getFFCenterOfMass(float &xCenter, float &yCenter) const {
    xCenter = yCenter = 0.0;
    if (FF_m)
        FF_m -> getCenterOfMass(xCenter, yCenter);
}

float ScanData::getFFTotalPower() const {
    if (FF_m)
        return FF_m -> getTotalPower();
    else
        return 0.0;
}
void ScanData::combineDualZScans() {
    if (NF_m && NF2_m) {
        combineDualZScans_impl(*NF_m, *NF2_m);
        delete NF2_m;
        NF2_m = NULL;
    }
    if (FF_m && FF2_m) {
        combineDualZScans_impl(*FF_m, *FF2_m);
        delete FF2_m;
        FF2_m = NULL;
    }
}

void ScanData::combineDualZScans_impl(ScanDataRaster &z1, ScanDataRaster &z2) {
    float amp1, phase1, amp2, phase2;

    z1.calcPeakAndPhase();
    z2.calcPeakAndPhase();

    amp1 = z1.getPeak(&phase1);
    amp2 = z2.getPeak(&phase2);

    // Calculate whether to add or subtract z2 based on peak amps and phases found:
    float intensity1 = pow(10.0, (amp1 / 20.0));
    float z1real = intensity1 * cos(phase1);
    float z1imag = intensity1 * sin(phase1);
    float intensity2 = pow(10.0, (amp2 / 20.0));
    float z2real = intensity2 * cos(phase2);
    float z2imag = intensity2 * sin(phase2);
    float magZPLUS = sqrt(pow(z1real - z2imag, 2.0) + pow(z1imag + z1real, 2.0));
    float magZMINUS = sqrt(pow(z1real + z2imag, 2.0) + pow(z1imag - z2real, 2.0));

    bool addingZ;
    if (magZPLUS >= magZMINUS)
        addingZ = true;  //z1+iz2
    else
        addingZ = false; //z1-iz2

    // get the raw implitude and phase raster data:
    const vector<float> &z1Amps = z1.getAmpArray();
    const vector<float> &z1Phis = z1.getPhaseArray();
    const vector<float> &z2Amps = z2.getAmpArray();
    const vector<float> &z2Phis = z2.getPhaseArray();

    // will compute new combined rasters here:
    vector<float> newAmps;
    vector<float> newPhis;

    // loop on all four arrays:
    vector<float>::const_iterator itA1 = z1Amps.begin();
    vector<float>::const_iterator itP1 = z1Phis.begin();
    vector<float>::const_iterator itA2 = z2Amps.begin();
    vector<float>::const_iterator itP2 = z2Phis.begin();

    while (itA1 != z1Amps.end()) {
        amp1 = *itA1;
        phase1 = *itP1;
        amp2 = *itA2;
        phase2 = *itP2;

        // Compute new output values from values read:
        intensity1 = pow(10.0 , (amp1 / 20.0));
        z1real = intensity1 * cos(phase1);
        z1imag = intensity1 * sin(phase1);
        intensity2 = pow(10.0 , (amp2 / 20.0));
        z2real = intensity2 * cos(phase2);
        z2imag = intensity2 * sin(phase2);

        if (addingZ) {
           z1real = 0.5 * (z1real - z2imag);
           z1imag = 0.5 * (z1imag + z2real);
        } else {
           z1real = 0.5 * (z1real + z2imag);
           z1imag = 0.5 * (z1imag - z2real);
        }
        amp1 = 20.0 * log10(sqrt(pow(z1real, 2.0) + pow(z1imag, 2.0)));
        phase1 = atan2(z1imag, z1real);  // radians
        newAmps.push_back(amp1);
        newPhis.push_back(phase1);

        itA1++;
        itP1++;
        itA2++;
        itP2++;
    }
    // replace the z1 amplitude and phase arrays with the new combined results:
    z1.replaceAmpAndPhase(newAmps, newPhis);
}

void ScanData::analyzeBeams(ALMAConstants::PointingOptions pointingOption,
                            float azPointing, float elPointing,
                            float copolPeakAmp, bool doUnwrapPhase)
{
    if (NF_m) {
        NF_m -> calcStepSize();
        NF_m -> calcPeakAndPhase();
    }
    if (FF_m) {
        FF_m -> calcStepSize();
        FF_m -> calcPeakAndPhase();
        // for copol scans, compute the unwrapped phase, needed for finding the phase center:

        // if the pointing to use for calculation is not already specified from another beam,
        if (azPointing == 0.0 && elPointing == 0.0) {
            // use the ACTUAL pointing of this beam:
            if (pointingOption == ALMAConstants::ACTUAL)
                FF_m -> getCenterOfMass(azPointing, elPointing);
            // or use the nominal pointing for the requested pointingOption:
            else
                ALMAConstants::getNominalAngles(band_m, pointingOption, azPointing, elPointing);
        }

        // get the secondary reflector radius for the selected pointing option:
        float subreflectorRadius = ALMAConstants::getSubreflectorRadius(pointingOption);

        // Use the copolPeakAmp for normalizing if specified. Otherwise normalize using this beam's peak:
        float peak = (copolPeakAmp != 0.0) ? copolPeakAmp : FF_m -> getPeak();

        // Calculate sums to compute efficiencies using the pointing and radius determined above:
        FF_m -> analyzeBeam(azPointing, elPointing, subreflectorRadius, peak, doUnwrapPhase);
    }
}

std::string ScanData::getScanTypeString() const {
    string ret("UNDEFINED");
    switch (scanType_m) {
    case ScanData::COPOL:
        ret = "copol";
        break;

    case ScanData::XPOL:
        ret = "xpol";
        break;

    case ScanData::COPOL180:
        ret = "copol180";
        break;

    case ScanData::DUAL:
        ret = "4545_scan";
        break;

    case ScanData::UNDEFINED:
    default:
        break;
    }
    return ret;
}

float ScanData::getFFPeak(float *phase) const {
    return FF_m ? FF_m -> getPeak(phase) : 0.0;
}

float ScanData::getNFPeak(float *phase) const {
    return NF_m ? NF_m -> getPeak(phase) : 0.0;
}

void ScanData::setPhaseFitResults(float deltaX, float deltaY, float deltaZ, float etaPhase) {
    if (FF_m) {
        AnalyzeResults &FFRes = FF_m -> useAnalyzeResults();
        FFRes.deltaX = deltaX;
        FFRes.deltaY = deltaY;
        FFRes.deltaZ = deltaZ;
        FFRes.etaPhase = etaPhase;
    }
}

void ScanData::setAmpFitResults(float ampFitAmp, float ampFitWidthDeg, float ampFitUOffDeg,
                                float ampFitVOffDeg, float ampFitD_0_90, float ampFitD_45_135)
{
    if (FF_m) {
        AnalyzeResults &FFRes = FF_m -> useAnalyzeResults();
        FFRes.ampFitAmp = ampFitAmp;
        FFRes.ampFitWidthDeg = ampFitWidthDeg;
        FFRes.ampFitUOffDeg = ampFitUOffDeg;
        FFRes.ampFitVOffDeg = ampFitVOffDeg;
        FFRes.ampFitD_0_90 = ampFitD_0_90;
        FFRes.ampFitD_45_135 = ampFitD_45_135;
    }
}

void ScanData::subtractForAttenuator(float ifAttenDiff) {
    if (FF_m)
        FF_m -> subtractForAttenuator(ifAttenDiff);
    if (NF_m)
        NF_m -> subtractForAttenuator(ifAttenDiff);
}

void ScanData::print(int _indent) const {
    string indent(_indent, ' ');

    cout << indent << "inputSection_m = '" << inputSection_m << "'" << endl;
    cout << indent << "scanId_m = " << scanId_m << endl;
    cout << indent << "measDateTime_m = '" << measDateTime_m << "'" << endl;
    cout << indent << "NSIDateTime_m = '" << NSIDateTime_m << "'" << endl;

    cout << indent << "filenameNF_m = '" << filenameNF_m << "'" << endl;
    cout << indent << "filenameFF_m = '" << filenameFF_m << "'" << endl;
    cout << indent << "filenameNF2_m = '" << filenameNF2_m << "'" << endl;
    cout << indent << "filenameFF2_m = '" << filenameFF2_m << "'" << endl;

    cout << indent << "scanType_m = " << getScanTypeString() << endl;
    cout << indent << "rfGHz_m = " << rfGHz_m << endl;
    cout << indent << "pol_m = " << pol_m << endl;
    cout << indent << "sb_m = " << sb_m << endl;
    cout << indent << "nominalFocusZ_m = " << nominalFocusZ_m << endl;
    cout << indent << "ifAtten_m = " << ifAtten_m << endl;

    cout << indent << "startrowNF_m = " << startrowNF_m << endl;
    cout << indent << "startrowFF_m = " << startrowFF_m << endl;
    cout << indent << "startrowNF2_m = " << startrowNF2_m << endl;
    cout << indent << "startrowFF2_m = " << startrowFF2_m << endl;

    if (FF_m) {
        cout << indent << "ScanDataRaster(FF):" << endl;
        FF_m -> print(_indent + 2, 3, 2);
    }
    if (NF_m) {
        cout << indent << "ScanDataRaster(NF):" << endl;
        NF_m -> print(_indent + 2, 3, 2);
    }
    if (FF2_m) {
        cout << indent << "ScanDataRaster(FF2):" << endl;
        FF2_m -> print(_indent + 2, 3, 2);
    }
    if (NF2_m) {
        cout << indent << "ScanDataRaster(NF2):" << endl;
        NF2_m -> print(_indent + 2, 3, 2);
    }
}

