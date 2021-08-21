#ifndef SCANDATARASTER_H_
#define SCANDATARASTER_H_
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

#include "AnalyzeResults.h"
#include <string>
#include <vector>

/// Contains the raw and intermediate-processed x, y, amplitude, and phase from a single scan.
class ScanDataRaster {
public:
    ScanDataRaster()
      { clear(); }
    ///< default construct

    ~ScanDataRaster()
      {}
    ///< destructor

    void clear();
    ///< reset the object to its just-constructed state.

    bool loadFile(const std::string filename, const std::string delim, bool rotate, bool invertPhase);
    ///< load a file from disk, including parsing for special NSI text file markers.
    ///< @param filename to load
    ///< @param delim: delimiter to use when parsing the file
    ///< @param rotate: if true, beam XY coordinates are rotated 180 degrees.
    ///< @param invertPhase: if true, beam sign is inverted on phase data.
    ///< @return true if no errors loading and parsing the file

    bool saveFile(const std::string filename, float copolPeakAmp = 0.0,
                  const std::string &delim = "\t", const std::string &lineterm = "\r\n") const;
    ///< save a file to disk
    ///< @param filename to create or over-write
    ///< @param copolPeakAmp: If non-zero, use to normalize power levels before saving
    ///< @param delim: delimiter to use between columns in the file.
    ///< @param lineterm: newline string to put at the end of each row.
    ///< @return true if no errors saving the file.

    unsigned getStartRow() const
      { return startRow_m; }
    ///< @return the first row containing raster data in the file.

    const std::string &getNSIDateTime() const
      { return NSIDateTime_m; }
    ///< @return the NSI date and time string parsed from the file, if any.

    const std::vector<float> &getAmpArray() const
      { return ampArray_m; }
    ///< @return const access to the amplitude array.

    const std::vector<float> &getPhaseArray() const
      { return phiArray_m; }
    ///< @return const access to the phase array as radians.

    void replaceAmpAndPhase(const std::vector<float> &newAmpArray, const std::vector<float> &newPhaseArray)
      { ampArray_m = newAmpArray;
        phiArray_m = newPhaseArray; }
    ///< replace the contents of the amplitude and phase arrays.
    ///< @param newAmpArray new array of amplitudes
    ///< @param newPhaseArray new array of phases in radians.

    void calcStepSize();
    ///< calculate the X and Y step sizes
    ///< calculate the X and Y dimensions

    float getStepSize() const
      { return (results_m.xStepSize != 0.0) ? results_m.xStepSize : results_m.yStepSize; }
    ///< @return the scan step size.  It is assumed to be the same for X and Y.
    ///< Uses the Y step size if there are no X steps, e.g. a vertical cut.

    void calcPeakAndPhase();
    ///< calculate beam peak amplitude and the phase at the peak.
    ///< results are stored internally and can be retrieved with getPeak()

    float getPeak(float *phase = NULL) const
      { if (phase)
          *phase = results_m.phaseAtPeak;
        return results_m.maxAmp; }
    ///< @return the peak amplitude seen and optionally the phase at the peak.
    ///< @param phase: if not NULL return the phase in radians at the peak.

    void calcCenterOfMass();
    ///< calculate the center of mass of the beam.
    ///< results are stored and can be retrieved with getCenterOfMass

    void getCenterOfMass(float &xCenter, float &yCenter) const
    ///< retrieve the caculated center of mass of the beam.
    ///< @param xCenter will contain the x center of mass
    ///< @param yCenter will contain the y center of mass
      { xCenter = results_m.xCenterOfMass;
        yCenter = results_m.yCenterOfMass; }

    float getTotalPower() const
      { return results_m.totalPower; }
    ///< @return the total linear power found in the Center Of Mass calculation

    void subtractForAttenuator(float ifAttenDiff);
    ///< adjust all points for the difference in attenuation.
    ///< @param ifAttenDiff: copol - xpol attenuator setting in dB.

    void analyzeBeam(float azNominal, float elNominal, float subreflectorRadius, float copolPeakAmp = 0.0, bool doUnwrapPhase = false);
    ///< calculate statistics of the beam used for efficiency calculation.
    ///< This function only makes sense to use for far-field scans where x=Az and y=El.
    ///< @param azNominal: the nominal azimuth pointing to use for the subreflector mask.
    ///< @param elNominal: the nominal elevation pointing to use for the subreflector mask.
    ///< @param subreflectorRadius: in degrees.
    ///< @param copolPeakAmp: If non-zero, use to normalize the FF beam (to normalize xpol to copol peak.)
    ///< @param doUnwrapPhase: If true, unwrap the phase for copol beams

    bool unwrapPhase();
    ///< unwraps the phase array

    double calcPhaseEfficiency(double p[], float azNominal, float elNominal,
                              bool approx = false, bool reduceSubreflector = false) const;
    ///< calculate phase efficiency for the given phase center model p[].
    ///< called from within FitPhase()
    ///< @param p: array giving phase phase center guess in radians {delta_x, delta_y, delta_z}.
    ///< @param azNominal, elNominal: angle to center of the subreflector in degrees
    ///< @param approx:  if true, use approximate equation for fitPhase else use exact equation.
    ///< @param reduceSubreflector:  if true, reduce the subreflector radius to calculate over.
    ///< @return phase efficiency in 0-1.

    double calcAmplitudeEfficiency(double p[], float azActual, float elActual) const;
    ///< calculate amplitude efficiency at the given amplitude fit model p[].
    ///< called from within FitAmplitude()
    ///< @param p: array giving amplitude fit model
    ///< @param azActual: actual (center of mass) pointing in Az
    ///< @param elActual: actual (center of mass) pointing in El
    ///< @return amplitude efficiency in 0-1.

    const AnalyzeResults &getAnalyzeResults() const
      { return results_m; }
    ///< @return const analysis results

    AnalyzeResults &useAnalyzeResults()
      { return results_m; }
    ///< @return non-const analysis results

    void print(int indent = 0, unsigned headRows = 0, unsigned tailRows = 0) const;
    ///< output object state
    ///< @param indent number of spaces
    ///< @param headRows number of rows from head to show.
    ///< @param tailRows number of rows from tail to show.

private:
    // metadata about the file and loaded data:
    unsigned size_m;               ///< size of our arrays
    unsigned startRow_m;           ///< row on which actual raster data starts
    std::string NSIDateTime_m;          ///< timestamp from NSI format text file.

    // raw raster data arrays:
    std::vector<float> xArray_m;            ///< x or az coordinates
    std::vector<float> yArray_m;            ///< y or el coordinates
    std::vector<float> ampArray_m;          ///< amplitudes in dB
    std::vector<float> phiArray_m;          ///< phases in radians
    std::vector<float> phiArrayUnwrapped_m; ///< unwrapped phases
    std::vector<float> *phiArrayPhaseFit_mp; ///< pointer to which phiArray to use for phasefit
    std::vector<float> EArray_m;            ///< electric field voltages
    std::vector<float> RadiusArray_m;       ///< distance from nominal beam center
    std::vector<float> MaskArray_m;         ///< subreflector mask
    std::vector<float> MaskArrayReduced_m;  ///< reduced subreflector mask option phase

    AnalyzeResults results_m;               ///< results from calcPeakAndPhase(), calcCenterOfMass(), and analyzeBeam()

    void printRow(unsigned index) const;
    ///< print a single row.

    void calcPeakAndPhase_impl(float &maxAmp, float &phaseAtPeak,
                               const std::vector<float> &ampArray,
                               const std::vector<float> &phiArray) const;
    ///< private implementation
};

#endif /* SCANDATARASTER_H_ */
