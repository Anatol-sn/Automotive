#ifndef GEOMAGHEADER_H
#define GEOMAGHEADER_H
#pragma once
#include <emmintrin.h>

struct MagneticModelCoeff
{
    double Main_Field_Coeff_G;
    double Main_Field_Coeff_H;
    double Secular_Var_Coeff_G;
    double Secular_Var_Coeff_H;
};

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
// MagneticModelCoeff 2025/01/01 - 2029/12/31
const MagneticModelCoeff g_WMM_COF[] =
{
    //    2025.0            WMM-2025     11/13/2024
        {     0.0,       0.0,        0.0,        0.0},
        {-29351.8,       0.0,       12.0,        0.0},
        { -1410.8,    4545.4,        9.7,      -21.5},
        { -2556.6,       0.0,      -11.6,        0.0},
        {  2951.1,   -3133.6,       -5.2,      -27.7},
        {  1649.3,    -815.1,       -8.0,      -12.1},
        {  1361.0,       0.0,       -1.3,        0.0},
        { -2404.1,     -56.6,       -4.2,        4.0},
        {  1243.8,     237.5,        0.4,       -0.3},
        {   453.6,    -549.5,      -15.6,       -4.1},
        {   895.0,       0.0,       -1.6,        0.0},
        {   799.5,     278.6,       -2.4,       -1.1},
        {    55.7,    -133.9,       -6.0,        4.1},
        {  -281.1,     212.0,        5.6,        1.6},
        {    12.1,    -375.6,       -7.0,       -4.4},
        {  -233.2,       0.0,        0.6,        0.0},
        {   368.9,      45.4,        1.4,       -0.5},
        {   187.2,     220.2,        0.0,        2.2},
        {  -138.7,    -122.9,        0.6,        0.4},
        {  -142.0,      43.0,        2.2,        1.7},
        {    20.9,     106.1,        0.9,        1.9},
        {    64.4,       0.0,       -0.2,        0.0},
        {    63.8,     -18.4,       -0.4,        0.3},
        {    76.9,      16.8,        0.9,       -1.6},
        {  -115.7,      48.8,        1.2,       -0.4},
        {   -40.9,     -59.8,       -0.9,        0.9},
        {    14.9,      10.9,        0.3,        0.7},
        {   -60.7,      72.7,        0.9,        0.9},
        {    79.5,       0.0,       -0.0,        0.0},
        {   -77.0,     -48.9,       -0.1,        0.6},
        {    -8.8,     -14.4,       -0.1,        0.5},
        {    59.3,      -1.0,        0.5,       -0.8},
        {    15.8,      23.4,       -0.1,        0.0},
        {     2.5,      -7.4,       -0.8,       -1.0},
        {   -11.1,     -25.1,       -0.8,        0.6},
        {    14.2,      -2.3,        0.8,       -0.2},
        {    23.2,       0.0,       -0.1,        0.0},
        {    10.8,       7.1,        0.2,       -0.2},
        {   -17.5,     -12.6,        0.0,        0.5},
        {     2.0,      11.4,        0.5,       -0.4},
        {   -21.7,      -9.7,       -0.1,        0.4},
        {    16.9,      12.7,        0.3,       -0.5},
        {    15.0,       0.7,        0.2,       -0.6},
        {   -16.8,      -5.2,       -0.0,        0.3},
        {     0.9,       3.9,        0.2,        0.2},
        {     4.6,       0.0,       -0.0,        0.0},
        {     7.8,     -24.8,       -0.1,       -0.3},
        {     3.0,      12.2,        0.1,        0.3},
        {    -0.2,       8.3,        0.3,       -0.3},
        {    -2.5,      -3.3,       -0.3,        0.3},
        {   -13.1,      -5.2,        0.0,        0.2},
        {     2.4,       7.2,        0.3,       -0.1},
        {     8.6,      -0.6,       -0.1,       -0.2},
        {    -8.7,       0.8,        0.1,        0.4},
        {   -12.9,      10.0,       -0.1,        0.1},
        {    -1.3,       0.0,        0.1,        0.0},
        {    -6.4,       3.3,        0.0,        0.0},
        {     0.2,       0.0,        0.1,       -0.0},
        {     2.0,       2.4,        0.1,       -0.2},
        {    -1.0,       5.3,       -0.0,        0.1},
        {    -0.6,      -9.1,       -0.3,       -0.1},
        {    -0.9,       0.4,        0.0,        0.1},
        {     1.5,      -4.2,       -0.1,        0.0},
        {     0.9,      -3.8,       -0.1,       -0.1},
        {    -2.7,       0.9,       -0.0,        0.2},
        {    -3.9,      -9.1,       -0.0,       -0.0},
        {     2.9,       0.0,        0.0,        0.0},
        {    -1.5,       0.0,       -0.0,       -0.0},
        {    -2.5,       2.9,        0.0,        0.1},
        {     2.4,      -0.6,        0.0,       -0.0},
        {    -0.6,       0.2,        0.0,        0.1},
        {    -0.1,       0.5,       -0.1,       -0.0},
        {    -0.6,      -0.3,        0.0,       -0.0},
        {    -0.1,      -1.2,       -0.0,        0.1},
        {     1.1,      -1.7,       -0.1,       -0.0},
        {    -1.0,      -2.9,       -0.1,        0.0},
        {    -0.2,      -1.8,       -0.1,        0.0},
        {     2.6,      -2.3,       -0.1,        0.0},
        {    -2.0,       0.0,        0.0,        0.0},
        {    -0.2,      -1.3,        0.0,       -0.0},
        {     0.3,       0.7,       -0.0,        0.0},
        {     1.2,       1.0,       -0.0,       -0.1},
        {    -1.3,      -1.4,       -0.0,        0.1},
        {     0.6,      -0.0,       -0.0,       -0.0},
        {     0.6,       0.6,        0.1,       -0.0},
        {     0.5,      -0.1,       -0.0,       -0.0},
        {    -0.1,       0.8,        0.0,        0.0},
        {    -0.4,       0.1,        0.0,       -0.0},
        {    -0.2,      -1.0,       -0.1,       -0.0},
        {    -1.3,       0.1,       -0.0,        0.0},
        {    -0.7,       0.2,       -0.1,       -0.1}
};

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
/*These error values are the NCEI error model
 *
 */
#define WMM_UNCERTAINTY_F 145
#define WMM_UNCERTAINTY_H 128
#define WMM_UNCERTAINTY_X 131
#define WMM_UNCERTAINTY_Y 94
#define WMM_UNCERTAINTY_Z 157
#define WMM_UNCERTAINTY_I 0.21
#define WMM_UNCERTAINTY_D_OFFSET 0.26
#define WMM_UNCERTAINTY_D_COEF 5625


#ifndef M_PI
#define M_PI    ((2)*(acos(0.0)))
#endif

#define RAD2DEG(rad)    ((rad)*(180.0L/M_PI))
#define DEG2RAD(deg)    ((deg)*(M_PI/180.0L))

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
struct MAGtype_MagneticModel 
{
    double  m_dEditionDate;
    double  m_dEpoch; /*Base time of Geomagnetic model epoch (yrs)*/
    char    m_szModelName[32];

    std::vector <double> m_vdMain_Field_Coeff_G;      /* C - Gauss coefficients of main geomagnetic model (nT) Index is (n * (n + 1) / 2 + m) */
    std::vector <double> m_vdMain_Field_Coeff_H;      /* C - Gauss coefficients of main geomagnetic model (nT) */
    std::vector <double> m_vdSecular_Var_Coeff_G;     /* CD - Gauss coefficients of secular geomagnetic model (nT/yr) */
    std::vector <double> m_vdSecular_Var_Coeff_H;     /* CD - Gauss coefficients of secular geomagnetic model (nT/yr) */

    int m_nMax; /* Maximum degree of spherical harmonic model */
    int m_nMaxSecVar; /* Maximum degree of spherical harmonic secular model */
    int m_iSecularVariationUsed; /* Whether or not the magnetic secular variation vector will be needed by program*/
    double m_dCoefficientFileEndDate;

    MAGtype_MagneticModel(double epoch = 2025) : m_dEpoch (epoch)
    { 
        m_dCoefficientFileEndDate = m_dEpoch + 5;
        strcpy_s(m_szModelName, "WMM-2020");  
        m_nMax = 12;
        m_nMaxSecVar = 12;

        m_dEditionDate = 0;
        m_iSecularVariationUsed = 0;
    }
 
    ~MAGtype_MagneticModel() 
    {
        m_vdMain_Field_Coeff_G.clear();
        m_vdMain_Field_Coeff_H.clear();
        m_vdSecular_Var_Coeff_G.clear();
        m_vdSecular_Var_Coeff_H.clear();
    }
};

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
struct MAGtype_Ellipsoid 
{
    double a; /*semi-major axis of the ellipsoid*/
    double b; /*semi-minor axis of the ellipsoid*/
    double fla; /* flattening */
    double epssq; /*first eccentricity squared */
    double eps; /* first eccentricity */
    double re; /* mean radius of  ellipsoid*/
};

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
struct MAGtype_CoordGeodetic
{
    double lambda; /* longitude */
    double phi; /* geodetic latitude */
    double HeightAboveEllipsoid; /* height above the ellipsoid (HaE) */
    double HeightAboveGeoid; /* (height above the EGM96 geoid model ) */
    int UseGeoid;
};

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
struct MAGtype_CoordSpherical 
{
    double lambda; /* longitude*/
    double phig; /* geocentric latitude*/
    double r; /* distance from the center of the ellipsoid*/
};

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
struct  MAGtype_Date{
    uint32_t    m_uYear;
    uint32_t    m_uMonth;
    uint32_t    m_uDay;
    double      m_dDecimalYear; /* decimal years */
};

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
struct MAGtype_GeoMagneticElements{
    double m_dDeclination; /* 1. Angle between the magnetic field vector and true north, positive east*/
    double m_dIncl; /*2. Angle between the magnetic field vector and the horizontal plane, positive down*/
    double m_dF; /*3. Magnetic Field Strength*/
    double m_dH; /*4. Horizontal Magnetic Field Strength*/
    double m_dX; /*5. Northern component of the magnetic field vector*/
    double m_dY; /*6. Eastern component of the magnetic field vector*/
    double m_dZ; /*7. Downward component of the magnetic field vector*/
    double m_dGV; /*8. The Grid Variation*/
    double m_dDecldot; /*9. Yearly Rate of change in declination*/
    double m_dIncldot; /*10. Yearly Rate of change in inclination*/
    double m_dFdot; /*11. Yearly rate of change in Magnetic field strength*/
    double m_dHdot; /*12. Yearly rate of change in horizontal field strength*/
    double m_dXdot; /*13. Yearly rate of change in the northern component*/
    double m_dYdot; /*14. Yearly rate of change in the eastern component*/
    double m_dZdot; /*15. Yearly rate of change in the downward component*/
    double m_dGVdot; /*16. Yearly rate of change in grid variation*/
};

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
struct MAGtype_LegendreFunction
{
    std::vector <double>    m_vdPcup;       /* Legendre Function */
    std::vector <double>    m_vdDerPcup;      /* Derivative of Legendre fcn */

    ~MAGtype_LegendreFunction()
    {
        m_vdPcup.clear();
        m_vdDerPcup.clear();
    };
};

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
struct  MAGtype_SphericalHarmonicVariables 
{
    std::vector <double>    m_vdRelativeRadiusPower; /* [earth_reference_radius_km / sph. radius ]^n  */
    std::vector <double>    m_vdCos_mlambda; /*cp(m)  - cosine of (m*spherical coord. longitude)*/
    std::vector <double>    m_vdSin_mlambda; /* sp(m)  - sine of (m*spherical coord. longitude) */

    ~MAGtype_SphericalHarmonicVariables()
    {
        m_vdRelativeRadiusPower.clear();
        m_vdCos_mlambda.clear();
        m_vdSin_mlambda.clear();
    }
};

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
struct MAGtype_MagneticResults 
{
    double Bx; /* North */
    double By; /* East */
    double Bz; /* Down */
};

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
class CDeclinationProc
{
    MAGtype_Ellipsoid     m_Ellip;
    MAGtype_MagneticModel m_MagneticModel;
    
   
    HRESULT zMAG_robustReadMagModels();
    void    zMAG_SetDefaults();
    HRESULT zMAG_DateToYear(MAGtype_Date* CalendarDate);
    void    zMAG_GeodeticToSpherical(MAGtype_Ellipsoid* pEllip, MAGtype_CoordGeodetic CoordGeodetic, MAGtype_CoordSpherical* CoordSpherical);
    void    zMAG_TimelyModifyMagneticModel(MAGtype_Date UserDate, MAGtype_MagneticModel& rMagneticModel, MAGtype_MagneticModel& rTimedMagneticModel);
    HRESULT zMAG_Geomag(MAGtype_Ellipsoid* pEllip, MAGtype_CoordSpherical CoordSpherical, MAGtype_CoordGeodetic CoordGeodetic, MAGtype_MagneticModel* TimedMagneticModel, MAGtype_GeoMagneticElements* GeoMagneticElements);
    void    zMAG_ComputeSphericalHarmonicVariables(MAGtype_Ellipsoid* pEllip, MAGtype_CoordSpherical CoordSpherical, int nMax, MAGtype_SphericalHarmonicVariables* pSphVariables);
    HRESULT zMAG_AssociatedLegendreFunction(MAGtype_CoordSpherical CoordSpherical, int nMax, MAGtype_LegendreFunction* LegendreFunction);
    HRESULT zMAG_PcupLow(std::vector<double>& rvPcup, std::vector<double>& rvdPcup, double x, int nMax);
    HRESULT zMAG_PcupHigh(std::vector<double>& rvPcup, std::vector<double>& rvdPcup, double x, int nMax);
    HRESULT zMAG_Summation(MAGtype_LegendreFunction* pLegendreFunction, MAGtype_MagneticModel* pMagneticModel, MAGtype_SphericalHarmonicVariables SphVariables, MAGtype_CoordSpherical CoordSpherical, MAGtype_MagneticResults* pMagneticResults);
    HRESULT zMAG_SummationSpecial(MAGtype_MagneticModel* pMagneticModel, MAGtype_SphericalHarmonicVariables SphVariables, MAGtype_CoordSpherical CoordSpherical, MAGtype_MagneticResults* pMagneticResults);
    HRESULT zMAG_SecVarSummation(MAGtype_LegendreFunction* pLegendreFunction, MAGtype_MagneticModel* pMagneticModel, MAGtype_SphericalHarmonicVariables SphVariables, MAGtype_CoordSpherical CoordSpherical, MAGtype_MagneticResults* pMagneticResults);
    HRESULT zMAG_SecVarSummationSpecial(MAGtype_MagneticModel* pMagneticModel, MAGtype_SphericalHarmonicVariables SphVariables, MAGtype_CoordSpherical CoordSpherical, MAGtype_MagneticResults* pMagneticResults);
    void    zMAG_RotateMagneticVector(MAGtype_CoordSpherical CoordSpherical, MAGtype_CoordGeodetic CoordGeodetic, MAGtype_MagneticResults MagneticResultsSph, MAGtype_MagneticResults* MagneticResultsGeo);
    void    zMAG_CalculateGeoMagneticElements(MAGtype_MagneticResults* pMagneticResultsGeo, MAGtype_GeoMagneticElements* pGeoMagneticElements);
    void    zMAG_CalculateSecularVariationElements(MAGtype_MagneticResults MagneticVariation, MAGtype_GeoMagneticElements* pMagneticElements);
    void    zMAG_WMMErrorCalc(double H, MAGtype_GeoMagneticElements* Uncertainty);
public:
    CDeclinationProc() { Init(); }
    HRESULT Init() 
    {
        HRESULT hr = zMAG_robustReadMagModels();
        if (FAILED(hr))
            return hr;

        zMAG_SetDefaults();

        return S_OK;
    }

    //----------------------------------------------------------------------------------------------
    // https://www.ncei.noaa.gov/products/world-magnetic-model
    // The model, associated software, and documentation are distributed by NCEI on behalf of the National Geospatial - Intelligence Agency(NGA).
    // The model is produced at 5 - year intervals, with the current model expiring on December 31, 2024.
    // North latitude positive, South negative, East longitude positive, West negative. 
    // All degrees are decimal. Height above WGS-84 ellipsoid in km
    // For Check https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml#igrfwmm 
    HRESULT GetDeclination(double dLatitude, double dLongitude, double dHeight, double& rdDeclination, double& rdDeclinationError, SYSTEMTIME* pSt = NULL);
};

#endif /*GEOMAGHEADER_H*/