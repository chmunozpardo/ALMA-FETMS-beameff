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

#include "ScanSet.h"
#include "ScanData.h"
#include "ScanDataRaster.h"
#include "ALMAConstants.h"
#include "SWVersion.h"
#include "FitPhase.h"
#include "FitAmplitude.h"
#include "iniparser.h"
#include "stringConvert.h"
#include <math.h>
#include <iostream>
#include <cstdlib>
using namespace std;

ScanSet::ScanSet(unsigned scanSet)
  : pointingOption_m(ALMAConstants::NOMINAL),
    CopolPol0_m(NULL),
    XpolPol0_m(NULL),
    CopolPol1_m(NULL),
    XpolPol1_m(NULL),
    Copol180_m(NULL)
{
    clear();
    scanSet_m = scanSet;
}

ScanSet::~ScanSet() {
    clear();
}

#define WHACK(P) { if (P) delete P; P = NULL; }

void ScanSet::clear() {
    scanSet_m = 0;
    scanSetId_m = 0;
    FEConfigId_m = 0;
    TestDataHeaderId_m = 0;
    band_m = 0;
    pointingOption_m = ALMAConstants::NOMINAL;
    WHACK(CopolPol0_m)
    WHACK(XpolPol0_m)
    WHACK(CopolPol1_m)
    WHACK(XpolPol1_m)
    WHACK(Copol180_m)
}

void ScanSet::setDatabaseKeys(unsigned scanSetId, unsigned FEConfigId, unsigned TestDataHeaderId) {
    scanSetId_m = scanSetId;
    FEConfigId_m = FEConfigId;
    TestDataHeaderId_m = TestDataHeaderId;
}

bool ScanSet::loadScan(const dictionary *dict, const std::string inputSection,
                       const std::string delim, ALMAConstants::InvertPhaseOptions invertPhaseOption)
{
    bool error = false;

    ScanData *scan = new ScanData;
    // load the contents of the section:
    if (!scan -> loadFromIni(dict, inputSection)) {
        cout << "ERROR: ScanSet::loadScan(): loadFromIni failed." << endl;
        error = true;

    // load the listing files:
    } else if (!scan -> loadListings(delim, invertPhaseOption)) {
        cout << "ERROR: ScanSet::loadScan(): loadListings failed." << endl;
        error = true;

    } else {
        // get the band number from the first scan:
        if (band_m == 0)
            band_m = scan -> getBand();

        // combine dual Z scans:
        scan -> combineDualZScans();

        // calculate center of mass for copol scans:
        if (scan -> getScanType() == ScanData::COPOL || scan -> getScanType() == ScanData::COPOL180)
            scan -> calcCenterOfMass();

        // determine which slot to store the scan object in:
        // check whether this is a 180-degree copol scan for phase center correction (can be either pol):
        if (scan -> getScanType() == ScanData::COPOL180) {
            if (Copol180_m) {
                cout << "ERROR: ScanSet::loadScan(): already a scan loaded for 180 degree copol." << endl;
                error = true;
            } else
                Copol180_m = scan;

        // check for pol 0 scans:
        } else if (scan -> getPol() == 0) {
            if (scan -> getScanType() == ScanData::COPOL) {
                if (CopolPol0_m) {
                    cout << "ERROR: ScanSet::loadScan(): already a scan loaded for copol pol0." << endl;
                    error = true;
                } else
                    CopolPol0_m = scan;
            }
            if (scan -> getScanType() == ScanData::XPOL) {
                if (XpolPol0_m) {
                    cout << "ERROR: ScanSet::loadScan(): already a scan loaded for xpol pol0." << endl;
                    error = true;
                } else
                    XpolPol0_m = scan;
            }

        // check for pol 1 scans:
        } else if (scan -> getPol() == 1) {
            if (scan -> getScanType() == ScanData::COPOL) {
                if (CopolPol1_m) {
                    cout << "ERROR: ScanSet::loadScan(): already a scan loaded for copol pol1." << endl;
                    error = true;
                } else
                    CopolPol1_m = scan;
            }
            if (scan -> getScanType() == ScanData::XPOL) {
                if (XpolPol1_m) {
                    cout << "ERROR: ScanSet::loadScan(): already a scan loaded for xpol pol1." << endl;
                    error = true;
                } else
                    XpolPol1_m = scan;
            }
        }
    }
    if (error) {
        delete scan;
        return false;
    } else
        return true;
}

bool ScanSet::calcEfficiencies(ALMAConstants::PointingOptions pointingOption, bool unwrapPhaseOption) {
    bool ret = true;

    // save the pointing and unwrap phase options specified:
    pointingOption_m = pointingOption;
    unwrapPhaseOption_m = unwrapPhaseOption;

    // calculate for pol0:
    ret = ret && calcEfficiencies_impl(CopolPol0_m, XpolPol0_m, Pol0Eff_m);

    // calculate for pol1:
    ret = ret && calcEfficiencies_impl(CopolPol1_m, XpolPol1_m, Pol1Eff_m);

    // calculate average Z offset of both scans:
    if (CopolPol0_m && CopolPol1_m) {
        const ScanDataRaster *pol0FF = CopolPol0_m -> getFFScan();
        const ScanDataRaster *pol1FF = CopolPol1_m -> getFFScan();
        if (pol0FF && pol1FF) {
            const AnalyzeResults &pol0 = pol0FF -> getAnalyzeResults();
            const AnalyzeResults &pol1 = pol1FF -> getAnalyzeResults();
            CombinedEff_m.nominal_z_offset = (pol0.deltaZ + pol1.deltaZ) / 2.0;
        }
    }

    // calculate phase 2 for pol0:
    ret = ret && calcEfficiencies_impl2(CopolPol0_m, XpolPol0_m, Pol0Eff_m);

    // calculate phase 2 for pol1:
    ret = ret && calcEfficiencies_impl2(CopolPol1_m, XpolPol1_m, Pol1Eff_m);

    // calculate beam squint:
    ret = ret && calcSquint_impl();
    return ret;
}

bool ScanSet::calcEfficiencies_impl(ScanData *copolScan, ScanData *xpolScan, EfficiencyData::OnePol &eff) {
    if (!copolScan)
        return false;

    // find peak and atten differences:
    if (xpolScan) {
        // Find the difference in copol and xpol IF attenuation:
        eff.ifatten_diff = copolScan -> getIfAtten() - xpolScan -> getIfAtten();
        // Adjust the xpol data to account for the difference:
        xpolScan -> subtractForAttenuator(eff.ifatten_diff);
    }

    float azPointing(0), elPointing(0), copolPeakAmp(0);

    // analyze copol, returns azPointing, elPointing for the selected pointingOption:
    if (analyzeCopol_impl(copolScan, azPointing, elPointing)) {
        // save the copol peak power to use for normalizing xpol beams below.
        copolPeakAmp = copolScan -> getFFPeak();
    }

    // analyze xpol, normalized by copolPeakAmp:
    // if using the ACTUAL pointing option, azActual and elActual have already been substituted above.
    if (xpolScan)
        xpolScan -> analyzeBeams(pointingOption_m, azPointing, elPointing, copolPeakAmp);

    return true;
}

bool ScanSet::calcEfficiencies_impl2(const ScanData *copolScan, const ScanData *xpolScan, EfficiencyData::OnePol &eff) {
    float M, psi_o, psi_m, plateFactor, dishDiameter;
    //, lambda, delta, beta;
    ALMAConstants::getAntennaParameters(pointingOption_m, M, psi_o, psi_m, plateFactor, dishDiameter);

    if (copolScan && xpolScan) {
        // calculate amplitude peak differences between copol and xpol:
        eff.peak_diff_FF = copolScan -> getFFPeak() - xpolScan -> getFFPeak();
        eff.peak_diff_NF = copolScan -> getNFPeak() - xpolScan -> getNFPeak();

        // get the analysis results to compute the rest of the efficiencies:
        const ScanDataRaster *copolFF = copolScan -> getFFScan();
        const ScanDataRaster *xpolFF = xpolScan -> getFFScan();

        if (copolFF && xpolFF) {
            const AnalyzeResults &copol = copolFF -> getAnalyzeResults();
            const AnalyzeResults &xpol = xpolFF -> getAnalyzeResults();

            // TICRA method for computing spillover:
            eff.eta_spill_co_cross = (copol.sumPowerOnSec + xpol.sumPowerOnSec) / (copol.sumSqE + xpol.sumSqE);

            // TICRA method ratio of copol beam to copol+xpol on secondary:
            eff.eta_pol_on_secondary = (copol.sumPowerOnSec) / (copol.sumPowerOnSec + xpol.sumPowerOnSec);

            // Total polarization*spillover efficiency:
            eff.eta_pol_spill = eff.eta_spill_co_cross * eff.eta_pol_on_secondary;

            // Spillover efficiency, using the 'alternative' definition from R.Hills paper,
            // is the ratio of copol power on the secondary to total copol power:
            eff.eta_spillover = copol.sumPowerOnSec / copol.sumSqE;

            // Polarization efficiency using the 'alternative' definition from R.Hills paper,
            // is the ratio of total copol power to total copol+xpol power, NOT masked for the secondary:
            eff.eta_pol = copol.sumSqE / (copol.sumSqE + xpol.sumSqE);

            // Taper efficiency (Amplitude efficiency in R.Hills paper) is the ratio of the illumination
            // of the secondary to a uniform illumination:
            eff.eta_taper = pow(copol.sumEOnSec, 2.0) / (copol.sumMask * copol.sumPowerOnSec);

            // Total efficiency excluding polarization and defocus:
            eff.deltaX = copol.deltaX;
            eff.deltaY = copol.deltaY;
            eff.deltaZ = copol.deltaZ;
            eff.eta_phase = copol.etaPhase;
            eff.eta_tot_np = eff.eta_phase * eff.eta_spillover * eff.eta_taper;

            // Illumination efficiency is the product of taper and spillover efficiency:
            eff.eta_illumination = eff.eta_taper * eff.eta_spillover;

            // Power in dB at edge of the secondary is computed from the average voltage in the edge region:
            eff.edge_taper_db = 20 * log10(copol.sumEOnEdge / copol.edgeNumPts);

            // Total efficiency excluding defocus:
            eff.eta_tot_nd = eff.eta_tot_np * eff.eta_pol;

            // Defocus efficiency:
            eff.eta_defocus = 1 - 0.3 * pow(
                    (copolScan -> getKWaveNumber() * (eff.deltaZ - CombinedEff_m.nominal_z_offset) / 1000) *
                    pow(32.0, -2.0), 2.0);

            // Total efficiency including defocus:
            eff.total_aperture_eff = eff.eta_tot_nd * eff.eta_defocus;

            return true;
        }
    }
    return false;
}

bool ScanSet::analyzeCopol_impl(ScanData *copolScan, float &azPointing, float &elPointing) {
    if (!copolScan)
        return false;

    // accumulate values needed for efficiency calculation, using the selected pointing option:
    copolScan -> analyzeBeams(pointingOption_m, 0.0, 0.0, 0.0, unwrapPhaseOption_m);

    // get the actual pointing angles for THIS BEAM:
    float azActual, elActual;
    copolScan -> getFFCenterOfMass(azActual, elActual);

    // pointing we will use for FitPhase which may be different:
    float azPhaseFit, elPhaseFit;

    // If the pointing option is ACTUAL then find the weighted average pointing of pol0 and po1 beams:
    if (pointingOption_m == ALMAConstants::ACTUAL) {
        float az0, el0, az1, el1;
        CopolPol0_m -> getFFCenterOfMass(az0, el0);
        CopolPol1_m -> getFFCenterOfMass(az1, el1);
        float totalPower0 = CopolPol0_m -> getFFTotalPower();
        float totalPower1 = CopolPol1_m -> getFFTotalPower();
        float weight0 = (totalPower0 / (totalPower0 + totalPower1));
        azPhaseFit = (az0 * weight0) + (az1 * (1.0 - weight0));
        elPhaseFit = (el0 * weight0) + (el1 * (1.0 - weight0));
        // but return the pointing of THIS BEAM:
        azPointing = azActual;
        elPointing = elActual;
        cout << "FitPhase using average of pol0, pol1 pointings: Az=" << azPhaseFit << " El=" << elPhaseFit << "<br>" << endl;
    // otherwise use the nominal pointing for the specified band and pointing option:
    } else {
        ALMAConstants::getNominalAngles(band_m, pointingOption_m, azPointing, elPointing);
        azPhaseFit = azPointing;
        elPhaseFit = elPointing;
        cout << "FitPhase using nominal pointing: Az=" << azPhaseFit << " El=" << elPhaseFit << "<br>" << endl;
    }

    // perform phase fit using the selected pointing angle:
    BeamFitting::FitPhase(copolScan, azPhaseFit, elPhaseFit);

    // perform amplitude fit always using the ACTUAL pointing of THIS BEAM.
    BeamFitting::FitAmplitude(copolScan, azActual, elActual);
    return true;
}

bool ScanSet::calcSquint_impl() {

    // check that both pol data sets are available:
    if (!CopolPol0_m && CopolPol1_m)
        return false;   // error.

    // get the analysis results for both polarizations:
    const ScanDataRaster *pol0FF = CopolPol0_m -> getFFScan();
    const ScanDataRaster *pol1FF = CopolPol1_m -> getFFScan();
    const AnalyzeResults &pol0 = pol0FF -> getAnalyzeResults();
    const AnalyzeResults &pol1 = pol1FF -> getAnalyzeResults();

    // if the 180-degree scan is present, perform analysis:
    int pol180 = -1;    // -1=UNDEFINED pol
    const ScanDataRaster *copol180FF = NULL;
    if (Copol180_m) {
        // perform the beam analysis on the 180-degree scan:
        float azPointing, elPointing;
        analyzeCopol_impl(Copol180_m, azPointing, elPointing);
        // get the polarization of the 180-degree scan:
        pol180 = Copol180_m -> getPol();
        // get the FF data set:
        copol180FF = Copol180_m -> getFFScan();
    }

    // Get the phase centers of the 0 and 180 deg co-pol patterns for one the first polarization,
    // and the 90 deg co-pol pattern for the other polarization.
    // Call the (x,y) values: (x0,y0), (x180,y180), and (x90,y90), respectively.

    float x0(0.0), y0(0.0), x90(0.0), y90(0.0), x180(0.0), y180(0.0);

    switch (pol180) {
    case 0:
        x0 = pol0.deltaX;       // the 180 deg scan is pol0
        y0 = pol0.deltaY;       // so (X0,Y0) are pol0
        x90 = pol1.deltaX;      // and (X90,Y90) are pol1.
        y90 = pol1.deltaY;
        break;

    case 1:
        x0 = pol1.deltaX;       // the 180 deg scan is pol1
        y0 = pol1.deltaY;       // so (X0,y0) are pol1
        x90 = pol0.deltaX;      // and (X90,Y90) are pol0.
        y90 = pol0.deltaY;
        break;

    case -1:
    default:
        // Default case we will calculate with no 180-degree scan, so no corrections.
        break;
    }

    float x_diff(0.0), y_diff(0.0);

    // no 180-deg scan so skip correction:
    if (!copol180FF) {
        // uncorrected difference between phase centers:
        x_diff = pol0.deltaX - pol1.deltaX;
        y_diff = pol0.deltaY - pol1.deltaY;
        CombinedEff_m.dist_between_centers_mm = fabs((sqrt(pow(x_diff, 2.0) + pow(y_diff, 2.0))));
        // No corrections:
        CombinedEff_m.x_corr = 0.0;
        CombinedEff_m.y_corr = 0.0;
        // Assign with no corrections applied:
        Pol0Eff_m.correctedX = pol0.deltaX;
        Pol0Eff_m.correctedY = pol0.deltaY;
        Pol1Eff_m.correctedX = pol1.deltaX;
        Pol1Eff_m.correctedY = pol1.deltaY;
        CombinedEff_m.corrected_pol = -1;

    // only perform correction if there is a 180 deg scan:
    } else {
        // Algorithm to correct phase centers described here:
        // https://safe.nrao.edu/wiki/bin/view/ALMA/BeamSquintFromSingleScan#correctionProcedure

        const AnalyzeResults &copol180 = copol180FF -> getAnalyzeResults();
        x180 = copol180.deltaX;    // Get (x180, y180)
        y180 = copol180.deltaY;

        // uncorrected difference between phase centers:
        x_diff = x0 - x90;
        y_diff = y0 - y90;

        float dX = x180 - x0;
        float dY = y180 - y0;
        float abs_diff = fabs(dX) - fabs(dY);

        // "There are only two possible corrections:
        //  depending on whether the relative scan angle of the second polarization was +90 or -90.
        //  To determine whether it was +90 or -90, compare the signs of x_diff, y_diff, and abs_diff.
        //  If all three of them are negative or exactly two of them are positive, then it was +90.
        //  Otherwise it was -90."

        bool positive90 = false;
        if ((x_diff * y_diff * abs_diff) < 0)
            positive90 = true;

        //compute x_corr and y_corr:
        if (positive90) {
            CombinedEff_m.x_corr = ((-1.0 * dX) - dY) / 2.0;
            CombinedEff_m.y_corr = (dX - dY) / 2.0;

        } else {
            CombinedEff_m.x_corr = ((-1.0 * dX) + dY) / 2.0;
            CombinedEff_m.y_corr = ((-1.0 * dX) - dY) / 2.0;
        }

        // apply corrections to x90 and y90 for distance and squint calculations:
        x90 += CombinedEff_m.x_corr;
        y90 += CombinedEff_m.y_corr;

        // calc corrected distance between beam centers:
        CombinedEff_m.dist_between_centers_mm = fabs((sqrt(pow(x0 - x90, 2.0) + pow(y0 - y90, 2.0))));

        // save the corrected phase centers:
        switch (pol180) {
        case 0:
            Pol0Eff_m.correctedX = x0;
            Pol0Eff_m.correctedY = y0;
            Pol1Eff_m.correctedX = x90;
            Pol1Eff_m.correctedY = y90;
            CombinedEff_m.corrected_pol = 1;
            break;
        case 1:
            Pol0Eff_m.correctedX = x90;
            Pol0Eff_m.correctedY = y90;
            Pol1Eff_m.correctedX = x0;
            Pol1Eff_m.correctedY = y0;
            CombinedEff_m.corrected_pol = 0;
            break;
        default:
            // impossible case
            break;
        }
    }

    // Calculate squint using the phase centers:
    float M, psi_o, psi_m, plateFactor, dishDiameter;
    ALMAConstants::getAntennaParameters(pointingOption_m, M, psi_o, psi_m, plateFactor, dishDiameter);

    // Calculate squint:
    CombinedEff_m.squint_arcseconds = CombinedEff_m.dist_between_centers_mm * plateFactor;

    // Calculate squint in percentage of units of FWHM of the beam:

    // lambda in mm, hence 1.0E6, not 1.0E9:
    float lambda = ALMAConstants::c / (CopolPol0_m -> getRFGhz() * 1.0E6);

    // 1.15 is the coefficient to multiply by lambda/D to get FWHM where D=diameter of primary mirror in mm.
    // 60.0 * 60.0 converts to arcseconds.
    CombinedEff_m.squint_percent = (100.0 * CombinedEff_m.squint_arcseconds) / (1.15 * lambda * ALMAConstants::RAD_TO_DEG * 60.0 * 60.0 / dishDiameter);

    cout << "squint_percent = " << CombinedEff_m.squint_percent << endl;

    return true;
}

bool ScanSet::writeOutputFile(dictionary *dict, const std::string outputFilename) {
    string section;

    bool ret = true;

    // Save the software version to the settings section:
    updateDictionary(dict, "settings", "software_version", BEAMEFF_SW_VERSION_STRING.c_str());

    // make a scanset_id or RF specific results section for results which apply to both polarizations:
    section = "results_";
    if (scanSetId_m) {
        section += "ssid";
        section += to_string(scanSetId_m);
    } else {
        section += "rf";
        section += to_string(CopolPol0_m -> getRFGhz());
    }
    iniparser_set(dict, section.c_str(), NULL);

    // Pointing angles plot:
    updateDictionary(dict, section, "pointingangles", CombinedEff_m.FFPointingPlot);
    // Also put the pointing angles plot in "settings" for legacy LV app:
    updateDictionary(dict, "settings", "pointingangles", CombinedEff_m.FFPointingPlot);

    // Z offset and squint:
    updateDictionary(dict, section, "nominal_z_offset",     to_string(CombinedEff_m.nominal_z_offset, std::fixed, 6));
    updateDictionary(dict, section, "squint_percent",       to_string(CombinedEff_m.squint_percent, std::fixed, 2));
    updateDictionary(dict, section, "squint_arcseconds",    to_string(CombinedEff_m.squint_arcseconds, std::fixed, 6));

    // Phase center correction:
    updateDictionary(dict, section, "x_corr",    to_string(CombinedEff_m.x_corr, std::fixed, 3));
    updateDictionary(dict, section, "y_corr",    to_string(CombinedEff_m.y_corr, std::fixed, 3));
    updateDictionary(dict, section, "corrected_pol",    to_string(CombinedEff_m.corrected_pol));
    updateDictionary(dict, section, "dist_between_centers_mm",    to_string(CombinedEff_m.dist_between_centers_mm, std::fixed, 3));

    // output the results for pol0 copol:
    if (CopolPol0_m) {
        section = CopolPol0_m -> getInputSection();
        ret = ret && updateCopolSection(dict, section, CopolPol0_m, Pol0Eff_m);
    }
    // output the results for pol0 xpol:
    if (XpolPol0_m) {
        section = XpolPol0_m -> getInputSection();
        ret = ret && updateXpolSection(dict, section, XpolPol0_m, Pol0Eff_m);
    }
    // output the results for pol1 copol:
    if (CopolPol1_m) {
        section = CopolPol1_m -> getInputSection();
        ret = ret && updateCopolSection(dict, section, CopolPol1_m, Pol1Eff_m);
    }
    // output the results for pol1 xpol:
    if (XpolPol1_m) {
        section = XpolPol1_m -> getInputSection();
        ret = ret && updateXpolSection(dict, section, XpolPol1_m, Pol1Eff_m);
    }

    // TODO: also write out what we calculated for the copol180 section.

    // write out the whole dictionary to the output file:
    FILE *f = fopen(outputFilename.c_str(), "w");
    if (!f)
        ret = false;
    else {
        iniparser_dump_ini(dict, f);
        fclose(f);
    }
    return ret;
}

bool ScanSet::updateCopolSection(dictionary *dict, const std::string &section, const ScanData *copol, const EfficiencyData::OnePol &eff) {
    if (!copol)
        return false;

    // get the nominal angles, just for the output file:
    // Override ACTUAL with NOMINAL but otherwise use the requested option.
    float azNominal, elNominal;
    ALMAConstants::PointingOptions option = (pointingOption_m == ALMAConstants::ACTUAL)
                                          ? ALMAConstants::NOMINAL
                                          : pointingOption_m;
    ALMAConstants::getNominalAngles(band_m, option, azNominal, elNominal);

    // get the analysis results
    const ScanDataRaster *pFF = copol -> getFFScan();
    const ScanDataRaster *pNF = copol -> getNFScan();
    if (!pFF || !pNF)
        return false;

    const AnalyzeResults &FF = pFF -> getAnalyzeResults();
    const AnalyzeResults &NF = pNF -> getAnalyzeResults();

    updateDictionary(dict, section, "zdistance",            to_string(copol -> getNominalFocusZ(), std::fixed, 2));

    updateDictionary(dict, section, "ifatten_diff",         to_string(eff.ifatten_diff, std::fixed, 3));
    updateDictionary(dict, section, "eta_spillover",        to_string(eff.eta_spillover, std::fixed, 6));
    updateDictionary(dict, section, "eta_taper",            to_string(eff.eta_taper, std::fixed, 6));
    updateDictionary(dict, section, "eta_illumination",     to_string(eff.eta_illumination, std::fixed, 6));
    updateDictionary(dict, section, "delta_x",              to_string(eff.deltaX, std::fixed, 6));
    updateDictionary(dict, section, "delta_y",              to_string(eff.deltaY, std::fixed, 6));
    updateDictionary(dict, section, "delta_z",              to_string(eff.deltaZ, std::fixed, 6));
    updateDictionary(dict, section, "corrected_x",          to_string(eff.correctedX, std::fixed, 6));
    updateDictionary(dict, section, "corrected_y",          to_string(eff.correctedY, std::fixed, 6));
    updateDictionary(dict, section, "eta_phase",            to_string(eff.eta_phase, std::fixed, 6));
    updateDictionary(dict, section, "edge_dB",              to_string(eff.edge_taper_db, std::fixed, 3));
    updateDictionary(dict, section, "eta_tot_np",           to_string(eff.eta_tot_np, std::fixed, 6));
    updateDictionary(dict, section, "eta_pol",              to_string(eff.eta_pol, std::fixed, 6));
    updateDictionary(dict, section, "eta_tot_nd",           to_string(eff.eta_tot_nd, std::fixed, 6));
    updateDictionary(dict, section, "eta_defocus",          to_string(eff.eta_defocus, std::fixed, 6));
    updateDictionary(dict, section, "defocus_efficiency",   to_string(eff.eta_defocus, std::fixed, 6));  // 2.0.9 keep old key name for now.

//    updateDictionary(dict, section, "shift_from_focus_mm",  to_string(eff.shift_from_focus_mm, std::fixed, 6));
//    updateDictionary(dict, section, "subreflector_shift_mm",to_string(eff.subreflector_shift_mm, std::fixed, 6));
//    updateDictionary(dict, section, "defocus_efficiency_due_to_moving_the_subreflector",
//                                                            to_string(eff.defocus_efficiency, std::fixed, 6));
    updateDictionary(dict, section, "total_aperture_eff",   to_string(eff.total_aperture_eff, std::fixed, 6));

    updateDictionary(dict, section, "ff_xcenter",           to_string(FF.xCenterOfMass, std::fixed, 6));
    updateDictionary(dict, section, "ff_ycenter",           to_string(FF.yCenterOfMass, std::fixed, 6));
    updateDictionary(dict, section, "nf_xcenter",           to_string(NF.xCenterOfMass, std::fixed, 6));
    updateDictionary(dict, section, "nf_ycenter",           to_string(NF.yCenterOfMass, std::fixed, 6));
    updateDictionary(dict, section, "az_nominal",           to_string(azNominal, std::fixed, 6));
    updateDictionary(dict, section, "el_nominal",           to_string(elNominal, std::fixed, 6));

    updateDictionary(dict, section, "max_ff_amp_db",        to_string(FF.maxAmp, std::fixed, 3));
    updateDictionary(dict, section, "max_nf_amp_db",        to_string(NF.maxAmp, std::fixed, 3));
    updateDictionary(dict, section, "ampfit_amp",           to_string(FF.ampFitAmp, std::fixed, 6));
    updateDictionary(dict, section, "ampfit_width_deg",     to_string(FF.ampFitWidthDeg, std::fixed, 6));
    updateDictionary(dict, section, "ampfit_u_off_deg",     to_string(FF.ampFitUOffDeg, std::fixed, 6));
    updateDictionary(dict, section, "ampfit_v_off_deg",     to_string(FF.ampFitVOffDeg, std::fixed, 6));
    updateDictionary(dict, section, "ampfit_d_0_90",        to_string(FF.ampFitD_0_90, std::fixed, 6));
    updateDictionary(dict, section, "ampfit_d_45_135",      to_string(FF.ampFitD_45_135, std::fixed, 6));

    updateDictionary(dict, section, "ampfit_d_45_135",      to_string(FF.ampFitD_45_135, std::fixed, 6));

    updateDictionary(dict, section, "plot_copol_ffamp",     eff.FFCopolAmpPlot);
    updateDictionary(dict, section, "plot_copol_ffphase",   eff.FFCopolPhasePlot);
    if (unwrapPhaseOption_m)
        updateDictionary(dict, section, "plot_copol_ffphase_unwrap",    eff.FFCopolPhaseUnwrappedPlot);
    updateDictionary(dict, section, "plot_copol_ffphase_model",     eff.FFCopolPhaseModelPlot);
    updateDictionary(dict, section, "plot_copol_nfamp",     eff.NFCopolAmpPlot);
    updateDictionary(dict, section, "plot_copol_nfphase",   eff.NFCopolPhasePlot);

    // for backwards compatibility with LabVIEW wrapper, put squint results in the pol1, copol section:
    if (copol -> getPol() == 1) {
        updateDictionary(dict, section, "squint_percent",       to_string(CombinedEff_m.squint_percent, std::fixed, 2));
        updateDictionary(dict, section, "squint_arcseconds",    to_string(CombinedEff_m.squint_arcseconds, std::fixed, 6));
    }
    return true;
}

bool ScanSet::updateXpolSection(dictionary *dict, const std::string &section, const ScanData *xpol, const EfficiencyData::OnePol &eff) {
    if (!xpol)
        return false;

    const ScanDataRaster *pFF = xpol -> getFFScan();
    if (!pFF)
        return false;

    const AnalyzeResults &FF = pFF -> getAnalyzeResults();

    updateDictionary(dict, section, "zdistance",            to_string(xpol -> getNominalFocusZ(), std::fixed, 2));
    updateDictionary(dict, section, "eta_spill_co_cross",   to_string(eff.eta_spill_co_cross, std::fixed, 6));
    updateDictionary(dict, section, "eta_pol_on_secondary", to_string(eff.eta_pol_on_secondary, std::fixed, 6));
    updateDictionary(dict, section, "eta_pol_spill",        to_string(eff.eta_pol_spill, std::fixed, 6));
    updateDictionary(dict, section, "max_dbdifference",     to_string(eff.peak_diff_FF, std::fixed, 3));
    updateDictionary(dict, section, "max_ff_amp_db",        to_string(FF.maxAmp, std::fixed, 3));

    updateDictionary(dict, section, "plot_xpol_ffamp",      eff.FFXpolAmpPlot);
    updateDictionary(dict, section, "plot_xpol_ffphase",    eff.FFXpolPhasePlot);
    updateDictionary(dict, section, "plot_xpol_nfamp",      eff.NFXpolAmpPlot);
    updateDictionary(dict, section, "plot_xpol_nfphase",    eff.NFXpolPhasePlot);
    return true;
}

void ScanSet::updateDictionary(dictionary *dict, const std::string &section, const std::string &key, const std::string &value) {
    string sectionKey(section);
    sectionKey += ":";
    sectionKey += key;
    iniparser_set(dict, sectionKey.c_str(), value.c_str());
}


bool ScanSet::makePlots(const std::string &outputDirectory, const std::string &gnuplotPath) {
    string fileNameFF, fileNameNF;
    bool ret = true;

    float azPointing(0), elPointing(0);
    ALMAConstants::getNominalAngles(band_m, pointingOption_m, azPointing, elPointing);

    float copolPeakAmpFF(0), copolPeakAmpNF(0);

    // Make the pol0 copol plots:
    if (CopolPol0_m) {

        // get the actual pointing angles for copol pol0:
        if (pointingOption_m == ALMAConstants::ACTUAL)
            CopolPol0_m -> getFFCenterOfMass(azPointing, elPointing);

        // get the peak amplitude for normalization:
        copolPeakAmpFF = CopolPol0_m -> getFFPeak();
        copolPeakAmpNF = CopolPol0_m -> getNFPeak();

        // write out the NF and FF data files for gnuplot to use:
        if (writePlotDataFiles(outputDirectory, CopolPol0_m, fileNameFF, fileNameNF, copolPeakAmpFF, copolPeakAmpNF)) {
            // make the Farfield Amplitude plot:
            ret = ret && makeOnePlot(outputDirectory, gnuplotPath, fileNameFF, CopolPol0_m,
                                     false, false, false, azPointing, elPointing, Pol0Eff_m.FFCopolAmpPlot);
            // make the Farfield Phase plot:
            ret = ret && makeOnePlot(outputDirectory, gnuplotPath, fileNameFF, CopolPol0_m,
                                     false, true, false, azPointing, elPointing, Pol0Eff_m.FFCopolPhasePlot);
            // make the Unwrapped Farfield Phase plot:
            if (unwrapPhaseOption_m)
                ret = ret && makeOnePlot(outputDirectory, gnuplotPath, fileNameFF, CopolPol0_m,
                                         false, true, true, azPointing, elPointing, Pol0Eff_m.FFCopolPhaseUnwrappedPlot);

            // make the Farfield PhaseFit model plot:
            ret = ret && makePhaseFitPlot(outputDirectory, gnuplotPath, CopolPol0_m, azPointing, elPointing, Pol0Eff_m.FFCopolPhaseModelPlot);

            // make the Nearfield Amplitude plot:
            ret = ret && makeOnePlot(outputDirectory, gnuplotPath, fileNameNF, CopolPol0_m,
                                     true, false, false, azPointing, elPointing, Pol0Eff_m.NFCopolAmpPlot);
            // make the Nearfield Phase plot:
            ret = ret && makeOnePlot(outputDirectory, gnuplotPath, fileNameNF, CopolPol0_m,
                                     true, true, false, azPointing, elPointing, Pol0Eff_m.NFCopolPhasePlot);
            // delete the tempory data files:
            remove(fileNameFF.c_str());
            remove(fileNameNF.c_str());
        }
    }

    // Make the pol0 xpol plots:
    if (XpolPol0_m) {

        // Assume we are using the copol pointing angles and peak amplitude from above.

        // write out the NF and FF data files for gnuplot to use:
        if (writePlotDataFiles(outputDirectory, XpolPol0_m, fileNameFF, fileNameNF, copolPeakAmpFF, copolPeakAmpNF)) {
            // make the Farfield Amplitude plot:
            ret = ret && makeOnePlot(outputDirectory, gnuplotPath, fileNameFF, XpolPol0_m,
                                     false, false, false, azPointing, elPointing, Pol0Eff_m.FFXpolAmpPlot);
            // make the Farfield Phase plot:
            ret = ret && makeOnePlot(outputDirectory, gnuplotPath, fileNameFF, XpolPol0_m,
                                     false, true, false, azPointing, elPointing, Pol0Eff_m.FFXpolPhasePlot);
            // make the Nearfield Amplitude plot:
            ret = ret && makeOnePlot(outputDirectory, gnuplotPath, fileNameNF, XpolPol0_m,
                                     true, false, false, azPointing, elPointing, Pol0Eff_m.NFXpolAmpPlot);
            // make the Nearfield Phase plot:
            ret = ret && makeOnePlot(outputDirectory, gnuplotPath, fileNameNF, XpolPol0_m,
                                     true, true, false, azPointing, elPointing, Pol0Eff_m.NFXpolPhasePlot);
            // delete the tempory data files:
            remove(fileNameFF.c_str());
            remove(fileNameNF.c_str());
        }
    }
    // reset the pointing assumption to nominal:
    ALMAConstants::getNominalAngles(band_m, pointingOption_m, azPointing, elPointing);

    // reset the amplitude peak registers to 0:
    copolPeakAmpFF = copolPeakAmpNF = 0.0;

    // Make the pol1 copol plots:
    if (CopolPol1_m) {

        // get the actual pointing angles for copol pol1:
        if (pointingOption_m == ALMAConstants::ACTUAL)
            CopolPol1_m -> getFFCenterOfMass(azPointing, elPointing);

        // get the peak amplitude for normalization:
        copolPeakAmpFF = CopolPol1_m -> getFFPeak();
        copolPeakAmpNF = CopolPol1_m -> getNFPeak();

        // write out the NF and FF data files for gnuplot to use:
        if (writePlotDataFiles(outputDirectory, CopolPol1_m, fileNameFF, fileNameNF, copolPeakAmpFF, copolPeakAmpNF)) {
            // make the Farfield Amplitude plot:
            ret = ret && makeOnePlot(outputDirectory, gnuplotPath, fileNameFF, CopolPol1_m,
                                     false, false, false, azPointing, elPointing, Pol1Eff_m.FFCopolAmpPlot);
            // make the Farfield Phase plot:
            ret = ret && makeOnePlot(outputDirectory, gnuplotPath, fileNameFF, CopolPol1_m,
                                     false, true, false, azPointing, elPointing, Pol1Eff_m.FFCopolPhasePlot);

            // make the Unwrapped Farfield Phase plot:
            if (unwrapPhaseOption_m)
                ret = ret && makeOnePlot(outputDirectory, gnuplotPath, fileNameFF, CopolPol1_m,
                                         false, true, true, azPointing, elPointing, Pol1Eff_m.FFCopolPhaseUnwrappedPlot);

            // make the Farfield PhaseFit model plot:
            ret = ret && makePhaseFitPlot(outputDirectory, gnuplotPath, CopolPol1_m, azPointing, elPointing, Pol1Eff_m.FFCopolPhaseModelPlot);

            // make the Nearfield Amplitude plot:
            ret = ret && makeOnePlot(outputDirectory, gnuplotPath, fileNameNF, CopolPol1_m,
                                     true, false, false, azPointing, elPointing, Pol1Eff_m.NFCopolAmpPlot);
            // make the Nearfield Phase plot:
            ret = ret && makeOnePlot(outputDirectory, gnuplotPath, fileNameNF, CopolPol1_m,
                                     true, true, false, azPointing, elPointing, Pol1Eff_m.NFCopolPhasePlot);
            // delete the tempory data files:
            remove(fileNameFF.c_str());
            remove(fileNameNF.c_str());
        }
    }
    // Make the pol1 xpol plots:
    if (XpolPol1_m) {

        // Assume we are using the copol pointing angles and peak amplitude from above.

        // write out the NF and FF data files for gnuplot to use:
        if (writePlotDataFiles(outputDirectory, XpolPol1_m, fileNameFF, fileNameNF, copolPeakAmpFF, copolPeakAmpNF)) {
            // make the Farfield Amplitude plot:
            ret = ret && makeOnePlot(outputDirectory, gnuplotPath, fileNameFF, XpolPol1_m,
                                     false, false, false, azPointing, elPointing, Pol1Eff_m.FFXpolAmpPlot);
            // make the Farfield Phase plot:
            ret = ret && makeOnePlot(outputDirectory, gnuplotPath, fileNameFF, XpolPol1_m,
                                     false, true, false, azPointing, elPointing, Pol1Eff_m.FFXpolPhasePlot);
            // make the Nearfield Amplitude plot:
            ret = ret && makeOnePlot(outputDirectory, gnuplotPath, fileNameNF, XpolPol1_m,
                                     true, false, false, azPointing, elPointing, Pol1Eff_m.NFXpolAmpPlot);
            // make the Nearfield Phase plot:
            ret = ret && makeOnePlot(outputDirectory, gnuplotPath, fileNameNF, XpolPol1_m,
                                     true, true, false, azPointing, elPointing, Pol1Eff_m.NFXpolPhasePlot);
            // delete the tempory data files:
            remove(fileNameFF.c_str());
            remove(fileNameNF.c_str());
        }
    }
    // Make the pointing angles plot:
    ret = ret && makePointingAnglesPlot(outputDirectory, gnuplotPath, CombinedEff_m.FFPointingPlot);

    return ret;
}

bool ScanSet::writePlotDataFiles(const std::string &outputDirectory, const ScanData *scan,
                                 std::string &fileNameFF, std::string &fileNameNF,
                                 float copolPeakAmpFF, float copolPeakAmpNF)
{
    fileNameFF.clear();
    fileNameNF.clear();

    if (!scan)
        return false;

    bool ret = true;

    // get the scan raster objects:
    const ScanDataRaster *pFF = scan -> getFFScan();
    const ScanDataRaster *pNF = scan -> getNFScan();

    if (pFF) {
        fileNameFF = outputDirectory;
        fileNameFF += "ffdata_temp.txt";
        ret = ret && pFF -> saveFile(fileNameFF, copolPeakAmpFF);
    }
    if (pNF) {
        fileNameNF = outputDirectory;
        fileNameNF += "nfdata_temp.txt";
        ret = ret && pNF -> saveFile(fileNameNF, copolPeakAmpNF);
    }
    return ret;
}

bool ScanSet::makeOnePlot(const std::string &outputDirectory, const std::string &gnuplotPath,
                          const std::string &dataFilename, const ScanData *scan,
                          bool nf, bool phase, bool unwrapped,
                          float azPointing, float elPointing, std::string &fileNamePlot)
{
    fileNamePlot.clear();
    if (!scan)
        return false;

    fileNamePlot = outputDirectory;
    fileNamePlot += "band";
    fileNamePlot += to_string(scan -> getBand());
    fileNamePlot += "_pol";
    fileNamePlot += to_string(scan -> getPol());
    fileNamePlot += "_";
    fileNamePlot += to_string(scan -> getScanTypeString());
    fileNamePlot += "_";
    fileNamePlot += to_string(scan -> getRFGhz());
    fileNamePlot += "GHz";
    fileNamePlot += (nf) ? "_nf" : "_ff";
    fileNamePlot += (phase) ? "phase" : "amp";
    fileNamePlot += (unwrapped) ? "_unw" : "";
    fileNamePlot += "_tilt";
    fileNamePlot += to_string(scan -> getTilt());
    fileNamePlot += "_scanset";
    fileNamePlot += to_string(scanSetId_m);
    fileNamePlot += ".png";
    // delete the plot file if it already exists:
    remove(fileNamePlot.c_str());

    string fileNameCmd = outputDirectory;
    fileNameCmd += "gnuplot_cmd";
    fileNameCmd += (nf) ? "_nf" : "_ff";
    fileNameCmd += (phase) ? "phase" : "amp";
    fileNameCmd += (unwrapped) ? "_unw" : "";
    fileNameCmd += ".txt";
    // delete the command file if it already exists:
    remove(fileNameCmd.c_str());

    FILE *f = fopen(fileNameCmd.c_str(), "w");
    if (!f)
        return false;

    string title = "Band ";
    title += to_string(scan -> getBand());
    title += ", pol ";
    title += to_string(scan -> getPol());
    title += " ";
    title += scan -> getScanTypeString();
    title += ", RF ";
    title += to_string(scan -> getRFGhz());
    title += " GHz, tilt ";
    title += to_string(scan -> getTilt());
    title += " deg";

    // newline constant.   TODO: See if just \n works on Windows and Linux:
    const char *NL = "\n";

    // 500x500 pixel .png files.  crop=no blank space around plot:
    fprintf(f, "set terminal png size 500, 500 crop%s", NL);
    fprintf(f, "set output '%s'%s", fileNamePlot.c_str(), NL);
    fprintf(f, "set title '%s'%s", title.c_str(), NL);
    // X and Y axis labels:
    fprintf(f, "set xlabel '%s'%s", (nf) ? "X(m)" : "Az(deg)", NL);
    fprintf(f, "set ylabel '%s'%s", (nf) ? "Y(m)" : "El(deg)", NL);
    // Palette limited to -50 to 0:
    fprintf(f, "set palette model RGB defined (-50 'purple', -40 'blue', -30 'green', -20 'yellow', -10 'orange', 0 'red')%s", NL);
    // Label for legend:
    fprintf(f, "set cblabel '%s %s'%s", (nf) ? "Nearfield" : "Farfield", (phase) ? "Phase (deg)" : "Amplitude (dB)", NL);
    // Set the range of values which are colored using the current palette:
    if (!phase)
        fprintf(f, "set cbrange [-50:0]%s", NL);
    else if (!unwrapped)
        fprintf(f, "set cbrange [-180:180]%s", NL);
    // Top-down view:
    fprintf(f, "set view 0,0%s", NL);
    // Palette mapped 3d style for splot:
    fprintf(f, "set pm3d map%s", NL);
    // Force a square plot within the canvas:
    fprintf(f, "set size square%s", NL);

    // For FF plots use parametric mode in degrees to draw the subreflector circle:
    if (!nf) {
        // for FF plots, set the resolution of the subreflector circle:
        fprintf(f, "set isosamples 13%s", NL);
        // draw parametric:
        fprintf(f, "set parametric%s", NL);
        fprintf(f, "set angles degrees%s", NL);
        fprintf(f, "set urange [0:360]%s", NL);
    }

    // Set X and Y major tics:
    if (nf)
        fprintf(f, "set xtics rotate%s", NL);
    else {
        fprintf(f, "set xtics rotate 2%s", NL);
        fprintf(f, "set ytics 2%s", NL);
    }

    // Add the database keys label:
    string label;
    getDatabaseKeysLabel(label, scan -> getScanId());
    if (!label.empty())
        fprintf(f, "set label '%s' at screen 0.01, 0.05%s", label.c_str(), NL);

    // Add the measurmement info label:
    getMeasInfoLabel(label, *scan);
    fprintf(f, "set label '%s' at screen 0.01, 0.02%s", label.c_str(), NL);

    // plot from the data file.  Column 3 is amplitude, 4 is phase, 5 is unwrapped phase:
    fprintf(f, "splot '%s' using 1:2:%s title ''", dataFilename.c_str(), (phase) ? ((unwrapped) ? "5" : "4") : "3");

    // for FF plots, draw a circle for the subreflector position:
    if (!nf) {
        float subreflectorRadius = ALMAConstants::getSubreflectorRadius(pointingOption_m);
        int color = (phase) ? -111 : -40;  // blue
        // draw three circles close together to make lines appear thick:
        // splot doesn't seem to support linewidth - the resolution is controlled by the input data grid.
        fprintf(f, ",%f + %.2f*cos(u),%f + %.2f*sin(u),%d notitle",
                   azPointing, subreflectorRadius, elPointing, subreflectorRadius, color);
        fprintf(f, ",%f + %.2f*cos(u * 0.98),%f + %.2f*sin(u * 0.98),%d notitle",
                   azPointing, subreflectorRadius, elPointing, subreflectorRadius, color);
        fprintf(f, ",%f + %.2f*cos(u * 0.96),%f + %.2f*sin(u * 0.96),%d notitle",
                   azPointing, subreflectorRadius, elPointing, subreflectorRadius, color);
    }
    fprintf(f, "%s", NL);
    fclose(f);

    // call gnuplot:
    string gnuplotCommand = gnuplotPath;
    gnuplotCommand += " ";
    gnuplotCommand += fileNameCmd;
    std::system(gnuplotCommand.c_str());
    // delete the temporary command file:
    remove(fileNameCmd.c_str());
    return true;
}

bool ScanSet::makePhaseFitPlot(const std::string &outputDirectory, const std::string &gnuplotPath,
                               const ScanData *scan,
                               float azPointing, float elPointing, std::string &fileNamePlot)
{
    fileNamePlot.clear();
    if (!scan)
        return false;

    const ScanDataRaster *pFF = scan -> getFFScan();
    if (!pFF)
        return false;

    const AnalyzeResults &res = pFF -> getAnalyzeResults();

    fileNamePlot = outputDirectory;
    fileNamePlot += "band";
    fileNamePlot += to_string(scan -> getBand());
    fileNamePlot += "_pol";
    fileNamePlot += to_string(scan -> getPol());
    fileNamePlot += "_";
    fileNamePlot += to_string(scan -> getScanTypeString());
    fileNamePlot += "_";
    fileNamePlot += to_string(scan -> getRFGhz());
    fileNamePlot += "GHz";
    fileNamePlot += "_phasefit";
    fileNamePlot += "_tilt";
    fileNamePlot += to_string(scan -> getTilt());
    fileNamePlot += "_scanset";
    fileNamePlot += to_string(scanSetId_m);
    fileNamePlot += ".png";
    // delete the plot file if it already exists:
    remove(fileNamePlot.c_str());

    string fileNameCmd = outputDirectory;
    fileNameCmd += "gnuplot_cmd_phasefit.txt";
    // delete the command file if it already exists:
    remove(fileNameCmd.c_str());

    FILE *f = fopen(fileNameCmd.c_str(), "w");
    if (!f)
        return false;

    string title = "Band ";
    title += to_string(scan -> getBand());
    title += ", pol ";
    title += to_string(scan -> getPol());
    title += " ";
    title += scan -> getScanTypeString();
    title += ", RF ";
    title += to_string(scan -> getRFGhz());
    title += " GHz, tilt ";
    title += to_string(scan -> getTilt());
    title += " deg";

    // 500x500 pixel .png files.  crop=no blank space around plot:
    fprintf(f, "set terminal png size 500, 500 crop\r\n");
    fprintf(f, "set output '%s'\r\n", fileNamePlot.c_str());
    fprintf(f, "set title '%s'\r\n", title.c_str());
    // X and Y axis labels:
    fprintf(f, "set xlabel 'Az(deg)'\r\n");
    fprintf(f, "set ylabel 'El(deg)'\r\n");
    // Using the reverse range of colors from phase fit plots so we can visually compare the peaks:
    fprintf(f, "set palette model RGB defined (-50 'black', -49.9 'red', -40 'orange', -30 'yellow', -20 'green', -10 'blue', 0 'purple')\r\n");
    // Scale it to this many unwrapped degrees:
    fprintf(f, "set cbrange [-4000:]\r\n");
    // Label for legend:
    fprintf(f, "set cblabel 'Fitted FF Phase (deg)'\r\n");
    // Top-down view:
    fprintf(f, "set view 0,0\r\n");
    // Palette mapped 3d style for splot:
    fprintf(f, "set pm3d map\r\n");
    // Force a square plot within the canvas:
    fprintf(f, "set size square\r\n");

    // calculate in rad:
    fprintf(f, "set angles radians\r\n");

    // Set the X and Y resolution of the plot:
    fprintf(f, "set isosamples 201,201\r\n");

    fprintf(f, "set xtics rotate 2\r\n");
    fprintf(f, "set ytics 2\r\n");
    fprintf(f, "set xrange [-16:16]\r\n");
    fprintf(f, "set yrange [-16:16]\r\n");

    // Add the database keys label:
    string label;
    getDatabaseKeysLabel(label, scan -> getScanId());
    if (!label.empty())
        fprintf(f, "set label '%s' at screen 0.01, 0.05\r\n", label.c_str());

    // Add the measurmement info label:
    getMeasInfoLabel(label, *scan);
    fprintf(f, "set label '%s' at screen 0.01, 0.02\r\n", label.c_str());

    // convert mm to radians:
    float k = scan -> getKWaveNumber();    // rad/m
    float deltaX = res.deltaX * 0.001 * k; // rad
    float deltaY = res.deltaY * 0.001 * k;
    float deltaZ = res.deltaZ * 0.001 * k;

    // our plot axes are degrees so define function to convert to radians:
    fprintf(f, "torad(a) = a * pi / 180\r\n");

    // define functions to offset the phase center by the nominal pointing:
    fprintf(f, "azoffset(az) = az - %f\r\n", azPointing);
    fprintf(f, "eloffset(el) = el - %f\r\n", elPointing);

//  phaseFit = p[1]*sin(Az)*cos(El) + p[2]*sin(El) + p[3]*cos(Az)*cos(El);
    string func = "phase(x, y) = ";
    func += to_string(deltaX);
    func += " * sin(torad(azoffset(x)))*cos(torad(eloffset(y))) + ";
    func += to_string(deltaY);
    func += " * sin(torad(eloffset(y))) + ";
    func += to_string(deltaZ);
    func += " * cos(torad(azoffset(x)))*cos(torad(eloffset(y)))";
    fprintf(f, "%s\r\n", func.c_str());

    // convert to degrees for plotting:
    fprintf(f, "todeg(a) = a * 180 / pi\r\n");
    // a mod b:
    fprintf(f, "mod(a,b) = a - (b * floor(a / b))\r\n");
    // plot in (-inf:0]:
    fprintf(f, "splot todeg(phase(x, y)) - todeg(phase(%f, %f)) title ''\r\n", azPointing, elPointing);
    fclose(f);

    // call gnuplot:
    string gnuplotCommand = gnuplotPath;
    gnuplotCommand += " ";
    gnuplotCommand += fileNameCmd;
    std::system(gnuplotCommand.c_str());
    // delete the temporary command file:
    remove(fileNameCmd.c_str());
    return true;
}

bool ScanSet::makePointingAnglesPlot(const std::string &outputDirectory, const std::string &gnuplotPath,
                                     std::string &fileNamePlot)
{
    fileNamePlot.clear();
    if (!CopolPol0_m || !CopolPol1_m)
        return false;

    fileNamePlot = outputDirectory;
    fileNamePlot += "band";
    fileNamePlot += to_string(CopolPol0_m -> getBand());
    if (scanSetId_m) {
        fileNamePlot += "_scanset";
        fileNamePlot += to_string(scanSetId_m);
    } else {
        fileNamePlot += "_rf";
        fileNamePlot += to_string(CopolPol0_m -> getRFGhz());
    }

    fileNamePlot += "_pointingangles";
    fileNamePlot += ".png";
    // delete the plot file if it already exists:
    remove(fileNamePlot.c_str());

    string fileNameCmd = outputDirectory;
    fileNameCmd += "gnuplot_cmd.txt";
    // delete the command file if it already exists:
    remove(fileNameCmd.c_str());

    FILE *f = fopen(fileNameCmd.c_str(), "w");
    if (!f)
        return false;

    string title = "Band ";
    title += to_string(CopolPol0_m -> getBand());
    title += " Pointing Angles, RF ";
    title += to_string(CopolPol0_m -> getRFGhz());
    title += " GHz, tilt ";
    title += to_string(CopolPol0_m -> getTilt());
    title += " deg";

    // always show nominal pointing in this plot, even if efficiencies are using ACTUAL:
    float azNominal, elNominal, az0, el0, az1, el1;
    ALMAConstants::PointingOptions option = (pointingOption_m == ALMAConstants::ACTUAL)
                                          ? ALMAConstants::NOMINAL
                                          : pointingOption_m;
    ALMAConstants::getNominalAngles(band_m, option, azNominal, elNominal);

    float subreflectorRadius = ALMAConstants::getSubreflectorRadius(pointingOption_m);

    // 800x800 pixel .png files.  crop=no blank space around plot:
    fprintf(f, "set terminal png size 800, 800 crop\r\n");
    fprintf(f, "set output '%s'\r\n", fileNamePlot.c_str());
    fprintf(f, "set title '%s'\r\n", title.c_str());

    // X and Y axis labels:
    fprintf(f, "set xlabel 'Az(deg)'\r\n");
    fprintf(f, "set ylabel 'El(deg)'\r\n");

    // use parametric mode in degrees to draw the subreflector circle:
    fprintf(f, "set parametric\r\n");

    // border is moved inward to make room for the key:
    fprintf(f, "set key outside\r\n");

    // x and y range large enough for the subreflector circle:
    fprintf(f, "set xrange [%.2f:%.2f]\r\n", azNominal - (0.5 * subreflectorRadius) - 2, azNominal + (0.5 * subreflectorRadius) + 2);
    fprintf(f, "set yrange [%.2f:%.2f]\r\n", elNominal - (0.5 * subreflectorRadius) - 2, elNominal + (0.5 * subreflectorRadius) + 2);

    // Force a square plot within the canvas:
    fprintf(f, "set size square\r\n");

    // Add the database keys label:
    string label;
    getDatabaseKeysLabel(label);
    if (!label.empty())
        fprintf(f, "set label '%s' at screen 0.01, 0.11\r\n", label.c_str());

    // Add the measurement info label:
    getMeasInfoLabel(label, *CopolPol0_m);
    fprintf(f, "set label '%s' at screen 0.01, 0.09\r\n", label.c_str());

    // Add the pointing angles label
    CopolPol0_m -> getFFCenterOfMass(az0, el0);
    CopolPol1_m -> getFFCenterOfMass(az1, el1);
    label = "pol 0  az: ";
    label += to_string(az0, std::fixed, 3);
    label += "   el: ";
    label += to_string(el0, std::fixed, 3);
    fprintf(f, "set label '%s' at screen 0.70, 0.60\r\n", label.c_str());
    label = "pol 1  az: ";
    label += to_string(az1, std::fixed, 3);
    label += "   el: ";
    label += to_string(el1, std::fixed, 3);
    fprintf(f, "set label '%s' at screen 0.70, 0.57\r\n", label.c_str());

    // plot the subreflector circle:
    fprintf(f, "plot [0:2*pi] %f+%.2f*sin(t),%f+%.2f*cos(t) title 'subreflector'",
                azNominal, subreflectorRadius, elNominal, subreflectorRadius);

    // plot the nominal pointing angle:
    fprintf(f, ", %f,%f with points lw 1 pt 1 ", azNominal, elNominal);

    if (pointingOption_m == ALMAConstants::ACA7METER)
        fprintf(f, "title 'ACA 7m nominal pointing'");
    else
        fprintf(f, "title 'nominal pointing angle'");

    fprintf(f, ", %f,%f with points lw 1 pt 4 title 'pol 0'", az0, el0);
    fprintf(f, ", %f,%f with points lw 1 pt 3 title 'pol 1'", az1, el1);
    fprintf(f, "\r\n");

    fclose(f);

    // call gnuplot:
    string gnuplotCommand = gnuplotPath;
    gnuplotCommand += " ";
    gnuplotCommand += fileNameCmd;
    std::system(gnuplotCommand.c_str());
    // delete the temporary command file:
    remove(fileNameCmd.c_str());
    return true;
}

const std::string &ScanSet::getMeasInfoLabel(std::string &toFill, const ScanData &scan) const {
    // measurement date and software version:
    toFill = "";
    toFill += "MeasDate: ";
    toFill += scan.getMeasDateTime();
    toFill += ", BeamEff v";
    toFill += BEAMEFF_SW_VERSION_STRING;
    return toFill;
}

const std::string &ScanSet::getDatabaseKeysLabel(std::string &toFill, unsigned scanId) const {
    // database keys string...
    toFill = "";
    // Add either the test data header or the ScanSet id:
    if (TestDataHeaderId_m) {
        toFill += "TDH=" + to_string(TestDataHeaderId_m);
    } else if (scanSetId_m) {
        toFill += "ScanSet=" + to_string(scanSetId_m);
    }
    // Add the ScanDetails id if defined:
    if (scanId) {
        if (!toFill.empty())
            toFill += ", ";
        toFill += "ScanDetails=" + to_string(scanId);
    }
    // Add the FE config id if defined:
    if (FEConfigId_m) {
        if (!toFill.empty())
            toFill += ", ";
        toFill += "FEConfig=" + to_string(FEConfigId_m);
    }
    return toFill;
}

void ScanSet::print(int _indent) {
    string indent(_indent, ' ');
    cout << indent << "ScanSet(" << scanSet_m << "):" << endl;
    cout << indent << "band_m = " << band_m << endl;

    if (CopolPol0_m || XpolPol0_m) {
        cout << indent << "EfficiencyData(pol0):" << endl;
        Pol0Eff_m.print(_indent + 2);
    }
    cout << indent << "ScanData(pol0 copol): ";
    if (!CopolPol0_m)
        cout << indent << "NULL" << endl;
    else {
        cout << indent << endl;
        CopolPol0_m -> print(_indent + 2);
    }
    cout << indent << "ScanData(pol0 xpol): ";
    if (!XpolPol0_m)
        cout << indent << "NULL" << endl;
    else {
        cout << indent << endl;
        XpolPol0_m -> print(_indent + 2);
    }

    if (CopolPol1_m || XpolPol1_m) {
        cout << indent << "EfficiencyData(pol1):" << endl;
        Pol1Eff_m.print(_indent + 2);
    }
    cout << indent << "ScanData(pol1 copol): ";
    if (!CopolPol1_m)
        cout << indent << "NULL" << endl;
    else {
        cout << indent << endl;
        CopolPol1_m -> print(_indent + 2);
    }
    cout << indent << "ScanData(pol1 xpol): ";
    if (!XpolPol1_m)
        cout << indent << "NULL" << endl;
    else {
        cout << indent << endl;
        XpolPol1_m -> print(_indent + 2);
    }
    cout << indent << "ScanData(180): ";
    if (!Copol180_m)
        cout << indent << "NULL" << endl;
    else {
        cout << indent << endl;
        Copol180_m -> print(_indent + 2);
    }

    cout << indent << "EfficiencyData(combined pols): " << endl;
    CombinedEff_m.print(_indent + 2);
}



