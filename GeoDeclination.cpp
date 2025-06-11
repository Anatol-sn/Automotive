#include "stdafx.h"
#include "GeoDeclination.h"


//!!!! Поки немає нових таблиц магнитного склонення ми беремо за останню валідну дату - 2024/12/31 (update - вже є)
//#define _USE_OLD_DECLINATION_TABLE_ 1

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//typedef std::unique_ptr< double > TDoubleArrayPtr;
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
HRESULT CDeclinationProc::zMAG_robustReadMagModels()
{
    int num_terms = m_MagneticModel.m_nMax * (m_MagneticModel.m_nMax + 1) / 2 + m_MagneticModel.m_nMax;

    assert(num_terms + 1 == _countof(g_WMM_COF));

    m_MagneticModel.m_vdMain_Field_Coeff_G.resize(num_terms + 1);
    m_MagneticModel.m_vdMain_Field_Coeff_H.resize(num_terms + 1);
    m_MagneticModel.m_vdSecular_Var_Coeff_G.resize(num_terms + 1);
    m_MagneticModel.m_vdSecular_Var_Coeff_H.resize(num_terms + 1);

    for each (auto & rCoef in g_WMM_COF)
    {
        m_MagneticModel.m_vdMain_Field_Coeff_G.push_back(rCoef.Main_Field_Coeff_G);
        m_MagneticModel.m_vdMain_Field_Coeff_H.push_back(rCoef.Main_Field_Coeff_H);
        m_MagneticModel.m_vdSecular_Var_Coeff_G.push_back(rCoef.Secular_Var_Coeff_G);
        m_MagneticModel.m_vdSecular_Var_Coeff_H.push_back(rCoef.Secular_Var_Coeff_H);
    }

    return S_OK;
}

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
void CDeclinationProc::zMAG_SetDefaults()
{
     /* Sets default values for WMM subroutines. */
     /* Sets WGS-84 parameters */

     m_Ellip.a = 6378.137; /*semi-major axis of the ellipsoid in */
     m_Ellip.b = 6356.7523142; /*semi-minor axis of the ellipsoid in */
     m_Ellip.fla = 1 / 298.257223563; /* flattening */
     m_Ellip.eps = sqrt(1 - (m_Ellip.b * m_Ellip.b) / (m_Ellip.a * m_Ellip.a)); /*first eccentricity */
     m_Ellip.epssq = (m_Ellip.eps * m_Ellip.eps); /*first eccentricity squared */
     m_Ellip.re = 6371.2; /* Earth's radius */
}

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
HRESULT CDeclinationProc::zMAG_DateToYear(MAGtype_Date* CalendarDate )

/* Converts a given calendar date into a decimal year,
it also outputs an error string if there is a problem
INPUT  CalendarDate  Pointer to the  data  structure with the following elements
                        int	Year;
                        int	Month;
                        int	Day;
                        double DecimalYear;      decimal years
OUTPUT  CalendarDate  Pointer to the  data  structure with the following elements updated
                        double DecimalYear;      decimal years
                Error	pointer to an error string
CALLS : none

 */
{
    //M_LOGGER_EXTERNALSTR(0, NULL);

    int temp = 0; /*Total number of days */
    int MonthDays[13];
    int ExtraDay = 0;
    int i;
    if (CalendarDate->m_uMonth == 0)
    {
        CalendarDate->m_dDecimalYear = CalendarDate->m_uYear;
        return TRUE;
    }
    if ((CalendarDate->m_uYear % 4 == 0 && CalendarDate->m_uYear % 100 != 0) || CalendarDate->m_uYear % 400 == 0)
        ExtraDay = 1;
    MonthDays[0] = 0;
    MonthDays[1] = 31;
    MonthDays[2] = 28 + ExtraDay;
    MonthDays[3] = 31;
    MonthDays[4] = 30;
    MonthDays[5] = 31;
    MonthDays[6] = 30;
    MonthDays[7] = 31;
    MonthDays[8] = 31;
    MonthDays[9] = 30;
    MonthDays[10] = 31;
    MonthDays[11] = 30;
    MonthDays[12] = 31;

    /******************Validation********************************/
    if (CalendarDate->m_uMonth <= 0 || CalendarDate->m_uMonth > 12)
    {
        wprintf_s(L"Error: The Month entered is invalid, valid months are '1 to 12'");
        return E_FAIL;
    }
    if (CalendarDate->m_uDay <= 0 || CalendarDate->m_uDay > MonthDays[CalendarDate->m_uMonth])
    {
        wprintf_s(L"The number of days in month % d is % d. Error: The day entered is invalid", CalendarDate->m_uMonth, MonthDays[CalendarDate->m_uMonth]);
        return E_FAIL;
    }
    /****************Calculation of t***************************/
    for (i = 1; i <= CalendarDate->m_uMonth; i++)
        temp += MonthDays[i - 1];
    temp += CalendarDate->m_uDay;
    CalendarDate->m_dDecimalYear = CalendarDate->m_uYear + (temp - 1) / (365.0 + ExtraDay);

    return S_OK;
} /*MAG_DateToYear*/


//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
void CDeclinationProc::zMAG_GeodeticToSpherical(MAGtype_Ellipsoid* pEllip, MAGtype_CoordGeodetic CoordGeodetic, MAGtype_CoordSpherical* CoordSpherical)

/* Converts Geodetic coordinates to Spherical coordinates

  INPUT   Ellip  data  structure with the following elements
                        double a; semi-major axis of the ellipsoid
                        double b; semi-minor axis of the ellipsoid
                        double fla;  flattening
                        double epssq; first eccentricity squared
                        double eps;  first eccentricity
                        double re; mean radius of  ellipsoid

                CoordGeodetic  Pointer to the  data  structure with the following elements updates
                        double lambda; ( longitude )
                        double phi; ( geodetic latitude )
                        double HeightAboveEllipsoid; ( height above the WGS84 ellipsoid (HaE) )
                        double HeightAboveGeoid; (height above the EGM96 Geoid model )

 OUTPUT		CoordSpherical 	Pointer to the data structure with the following elements
                        double lambda; ( longitude)
                        double phig; ( geocentric latitude )
                        double r;  	  ( distance from the center of the ellipsoid)

CALLS : none

 */
{
    double CosLat, SinLat, rc, xp, zp; /*all local variables */

    /*
     ** Convert geodetic coordinates, (defined by the WGS-84
     ** reference ellipsoid), to Earth Centered Earth Fixed Cartesian
     ** coordinates, and then to spherical coordinates.
     */

    CosLat = cos(DEG2RAD(CoordGeodetic.phi));
    SinLat = sin(DEG2RAD(CoordGeodetic.phi));

    /* compute the local radius of curvature on the WGS-84 reference ellipsoid */

    rc = pEllip->a / sqrt(1.0 - pEllip->epssq * SinLat * SinLat);

    /* compute ECEF Cartesian coordinates of specified point (for longitude=0) */

    xp = (rc + CoordGeodetic.HeightAboveEllipsoid) * CosLat;
    zp = (rc * (1.0 - pEllip->epssq) + CoordGeodetic.HeightAboveEllipsoid) * SinLat;

    /* compute spherical radius and angle lambda and phi of specified point */

    CoordSpherical->r = sqrt(xp * xp + zp * zp);
    CoordSpherical->phig = RAD2DEG(asin(zp / CoordSpherical->r)); /* geocentric latitude */
    CoordSpherical->lambda = CoordGeodetic.lambda; /* longitude */
}/*MAG_GeodeticToSpherical*/

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
void CDeclinationProc::zMAG_TimelyModifyMagneticModel(MAGtype_Date UserDate, MAGtype_MagneticModel& rMagneticModel, MAGtype_MagneticModel& rTimedMagneticModel)

/* Time change the Model coefficients from the base year of the model using secular variation coefficients.
Store the coefficients of the static model with their values advanced from epoch t0 to epoch t.
Copy the SV coefficients.  If input "t�" is the same as "t0", then this is merely a copy operation.
If the address of "TimedMagneticModel" is the same as the address of "MagneticModel", then this procedure overwrites
the given item "MagneticModel".

INPUT: UserDate
           MagneticModel
OUTPUT:TimedMagneticModel
CALLS : none
 */
{
    int n, m, index, a, b;
    rTimedMagneticModel.m_dEditionDate = rMagneticModel.m_dEditionDate;
    rTimedMagneticModel.m_dEpoch = rMagneticModel.m_dEpoch;
    rTimedMagneticModel.m_nMax = rMagneticModel.m_nMax;
    rTimedMagneticModel.m_nMaxSecVar = rMagneticModel.m_nMaxSecVar;
    a = rTimedMagneticModel.m_nMaxSecVar;
    b = (a * (a + 1) / 2 + a);
    strcpy_s(rTimedMagneticModel.m_szModelName, rMagneticModel.m_szModelName);
    for (n = 1; n <= rMagneticModel.m_nMax; n++)
    {
        for (m = 0; m <= n; m++)
        {
            index = (n * (n + 1) / 2 + m);
            if (index <= b)
            {
                rTimedMagneticModel.m_vdMain_Field_Coeff_H[index] = rMagneticModel.m_vdMain_Field_Coeff_H[index] + (UserDate.m_dDecimalYear - rMagneticModel.m_dEpoch) * rMagneticModel.m_vdSecular_Var_Coeff_H[index];
                rTimedMagneticModel.m_vdMain_Field_Coeff_G[index] = rMagneticModel.m_vdMain_Field_Coeff_G[index] + (UserDate.m_dDecimalYear - rMagneticModel.m_dEpoch) * rMagneticModel.m_vdSecular_Var_Coeff_G[index];
                rTimedMagneticModel.m_vdSecular_Var_Coeff_H[index] = rMagneticModel.m_vdSecular_Var_Coeff_H[index]; /* We need a copy of the secular var coef to calculate secular change */
                rTimedMagneticModel.m_vdSecular_Var_Coeff_G[index] = rMagneticModel.m_vdSecular_Var_Coeff_G[index];
            }
            else
            {
                rTimedMagneticModel.m_vdMain_Field_Coeff_H[index] = rMagneticModel.m_vdMain_Field_Coeff_H[index];
                rTimedMagneticModel.m_vdMain_Field_Coeff_G[index] = rMagneticModel.m_vdMain_Field_Coeff_G[index];
            }
        }
    }
} /* MAG_TimelyModifyMagneticModel */

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
HRESULT CDeclinationProc::zMAG_Geomag(MAGtype_Ellipsoid* pEllip, MAGtype_CoordSpherical CoordSpherical, MAGtype_CoordGeodetic CoordGeodetic,
    MAGtype_MagneticModel* pTimedMagneticModel, MAGtype_GeoMagneticElements* GeoMagneticElements)
    /*
    The main subroutine that calls a sequence of WMM sub-functions to calculate the magnetic field elements for a single point.
    The function expects the model coefficients and point coordinates as input and returns the magnetic field elements and
    their rate of change. Though, this subroutine can be called successively to calculate a time series, profile or grid
    of magnetic field, these are better achieved by the subroutine MAG_Grid.

    INPUT: Ellip
                  CoordSpherical
                  CoordGeodetic
                  TimedMagneticModel

    OUTPUT : GeoMagneticElements

    CALLS:  	MAG_AllocateLegendreFunctionMemory(NumTerms);  ( For storing the ALF functions )
                         MAG_ComputeSphericalHarmonicVariables( Ellip, CoordSpherical, TimedMagneticModel->nMax, &SphVariables); (Compute Spherical Harmonic variables  )
                         MAG_AssociatedLegendreFunction(CoordSpherical, TimedMagneticModel->nMax, LegendreFunction);  	Compute ALF
                         MAG_Summation(LegendreFunction, TimedMagneticModel, SphVariables, CoordSpherical, &MagneticResultsSph);  Accumulate the spherical harmonic coefficients
                         MAG_SecVarSummation(LegendreFunction, TimedMagneticModel, SphVariables, CoordSpherical, &MagneticResultsSphVar); Sum the Secular Variation Coefficients
                         MAG_RotateMagneticVector(CoordSpherical, CoordGeodetic, MagneticResultsSph, &MagneticResultsGeo); Map the computed Magnetic fields to Geodetic coordinates
                         MAG_CalculateGeoMagneticElements(&MagneticResultsGeo, GeoMagneticElements);   Calculate the Geomagnetic elements
                         MAG_CalculateSecularVariationElements(MagneticResultsGeoVar, GeoMagneticElements); Calculate the secular variation of each of the Geomagnetic elements

     */
{
    MAGtype_LegendreFunction               LegendreFunction;
    MAGtype_SphericalHarmonicVariables     SphVariables;
    int NumTerms;
    MAGtype_MagneticResults MagneticResultsSph, MagneticResultsGeo, MagneticResultsSphVar, MagneticResultsGeoVar;

    NumTerms = ((pTimedMagneticModel->m_nMax + 1) * (pTimedMagneticModel->m_nMax + 2) / 2);
    
    //pLegendreFunction = MAG_AllocateLegendreFunctionMemory(NumTerms); /* For storing the ALF functions */
    LegendreFunction.m_vdPcup.resize(NumTerms + 1);
    LegendreFunction.m_vdDerPcup.resize(NumTerms+1);

   
    //SphVariables = MAG_AllocateSphVarMemory(TimedMagneticModel->nMax);
    SphVariables.m_vdRelativeRadiusPower.resize(pTimedMagneticModel->m_nMax+1);
    SphVariables.m_vdCos_mlambda.resize(pTimedMagneticModel->m_nMax+1);
    SphVariables.m_vdSin_mlambda.resize(pTimedMagneticModel->m_nMax+1);

    zMAG_ComputeSphericalHarmonicVariables(pEllip, CoordSpherical, pTimedMagneticModel->m_nMax, &SphVariables); /* Compute Spherical Harmonic variables  */
    HRESULT hr; 
    hr = zMAG_AssociatedLegendreFunction(CoordSpherical, pTimedMagneticModel->m_nMax, &LegendreFunction); /* Compute ALF  */
    if (FAILED(hr))
        return hr;

    hr = zMAG_Summation(&LegendreFunction, pTimedMagneticModel, SphVariables, CoordSpherical, &MagneticResultsSph); /* Accumulate the spherical harmonic coefficients*/
    if (FAILED(hr))
        return hr;

    hr = zMAG_SecVarSummation(&LegendreFunction, pTimedMagneticModel, SphVariables, CoordSpherical, &MagneticResultsSphVar); /*Sum the Secular Variation Coefficients  */
    if (FAILED(hr))
        return hr;

    zMAG_RotateMagneticVector(CoordSpherical, CoordGeodetic, MagneticResultsSph, &MagneticResultsGeo); /* Map the computed Magnetic fields to Geodeitic coordinates  */
    zMAG_RotateMagneticVector(CoordSpherical, CoordGeodetic, MagneticResultsSphVar, &MagneticResultsGeoVar); /* Map the secular variation field components to Geodetic coordinates*/
    zMAG_CalculateGeoMagneticElements(&MagneticResultsGeo, GeoMagneticElements); /* Calculate the Geomagnetic elements, Equation 19 , WMM Technical report */
    zMAG_CalculateSecularVariationElements(MagneticResultsGeoVar, GeoMagneticElements); /*Calculate the secular variation of each of the Geomagnetic elements*/

    //MAG_FreeLegendreMemory(LegendreFunction);
    //MAG_FreeSphVarMemory(SphVariables);

    return S_OK;
} /*MAG_Geomag*/

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
void    CDeclinationProc::zMAG_ComputeSphericalHarmonicVariables(MAGtype_Ellipsoid* pEllip, MAGtype_CoordSpherical CoordSpherical, int nMax, MAGtype_SphericalHarmonicVariables* pSphVariables)

/* Computes Spherical variables
       Variables computed are (a/r)^(n+2), cos_m(lamda) and sin_m(lambda) for spherical harmonic
       summations. (Equations 10-12 in the WMM Technical Report)
       INPUT   Ellip  data  structure with the following elements
                             double a; semi-major axis of the ellipsoid
                             double b; semi-minor axis of the ellipsoid
                             double fla;  flattening
                             double epssq; first eccentricity squared
                             double eps;  first eccentricity
                             double re; mean radius of  ellipsoid
                     CoordSpherical 	A data structure with the following elements
                             double lambda; ( longitude)
                             double phig; ( geocentric latitude )
                             double r;  	  ( distance from the center of the ellipsoid)
                     nMax   integer 	 ( Maxumum degree of spherical harmonic secular model)\

     OUTPUT  SphVariables  Pointer to the   data structure with the following elements
             double RelativeRadiusPower[MAG_MAX_MODEL_DEGREES+1];   [earth_reference_radius_km  sph. radius ]^n
             double cos_mlambda[MAG_MAX_MODEL_DEGREES+1]; cp(m)  - cosine of (mspherical coord. longitude)
             double sin_mlambda[MAG_MAX_MODEL_DEGREES+1];  sp(m)  - sine of (mspherical coord. longitude)
     CALLS : none
 */
{
    double cos_lambda, sin_lambda;
    int m, n;
    cos_lambda = cos(DEG2RAD(CoordSpherical.lambda));
    sin_lambda = sin(DEG2RAD(CoordSpherical.lambda));
    /* for n = 0 ... model_order, compute (Radius of Earth / Spherical radius r)^(n+2)
    for n  1..nMax-1 (this is much faster than calling pow MAX_N+1 times).      */
    pSphVariables->m_vdRelativeRadiusPower[0] = (pEllip->re / CoordSpherical.r) * (pEllip->re / CoordSpherical.r);
    for (n = 1; n <= nMax; n++)
    {
        pSphVariables->m_vdRelativeRadiusPower[n] = pSphVariables->m_vdRelativeRadiusPower[n - 1] * (pEllip->re / CoordSpherical.r);
    }

    /*
     Compute cos(m*lambda), sin(m*lambda) for m = 0 ... nMax
           cos(a + b) = cos(a)*cos(b) - sin(a)*sin(b)
           sin(a + b) = cos(a)*sin(b) + sin(a)*cos(b)
     */
    pSphVariables->m_vdCos_mlambda[0] = 1.0;
    pSphVariables->m_vdSin_mlambda[0] = 0.0;

    pSphVariables->m_vdCos_mlambda[1] = cos_lambda;
    pSphVariables->m_vdSin_mlambda[1] = sin_lambda;
    for (m = 2; m <= nMax; m++)
    {
        pSphVariables->m_vdCos_mlambda[m] = pSphVariables->m_vdCos_mlambda[m - 1] * cos_lambda - pSphVariables->m_vdSin_mlambda[m - 1] * sin_lambda;
        pSphVariables->m_vdSin_mlambda[m] = pSphVariables->m_vdCos_mlambda[m - 1] * sin_lambda + pSphVariables->m_vdSin_mlambda[m - 1] * cos_lambda;
    }
} /*MAG_ComputeSphericalHarmonicVariables*/

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
HRESULT CDeclinationProc::zMAG_AssociatedLegendreFunction(MAGtype_CoordSpherical CoordSpherical, int nMax, MAGtype_LegendreFunction* LegendreFunction)

/* Computes  all of the Schmidt-semi normalized associated Legendre
functions up to degree nMax. If nMax <= 16, function MAG_PcupLow is used.
Otherwise MAG_PcupHigh is called.
INPUT  CoordSpherical 	A data structure with the following elements
                                                double lambda; ( longitude)
                                                double phig; ( geocentric latitude )
                                                double r;  	  ( distance from the center of the ellipsoid)
                nMax        	integer 	 ( Maxumum degree of spherical harmonic secular model)
                LegendreFunction Pointer to data structure with the following elements
                                                double *Pcup;  (  pointer to store Legendre Function  )
                                                double *dPcup; ( pointer to store  Derivative of Lagendre function )

OUTPUT  LegendreFunction  Calculated Legendre variables in the data structure

 */
{
    double sin_phi;
    HRESULT hr;

    sin_phi = sin(DEG2RAD(CoordSpherical.phig)); /* sin  (geocentric latitude) */

    if (nMax <= 16 || (1 - fabs(sin_phi)) < 1.0e-10) /* If nMax is less tha 16 or at the poles */
        hr = zMAG_PcupLow(LegendreFunction->m_vdPcup, LegendreFunction->m_vdDerPcup, sin_phi, nMax);
    else 
        hr = zMAG_PcupHigh(LegendreFunction->m_vdPcup, LegendreFunction->m_vdDerPcup, sin_phi, nMax);
    
    return hr;
} /*MAG_AssociatedLegendreFunction */

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
HRESULT CDeclinationProc::zMAG_PcupLow(std::vector<double>& rvPcup, std::vector<double>& rvdPcup, double x, int nMax)

/*   This function evaluates all of the Schmidt-semi normalized associated Legendre
        functions up to degree nMax.

        Calling Parameters:
                INPUT
                        nMax:	 Maximum spherical harmonic degree to compute.
                        x:		cos(colatitude) or sin(latitude).

                OUTPUT
                        Pcup:	A vector of all associated Legendgre polynomials evaluated at
                                        x up to nMax.
                   dPcup: Derivative of Pcup(x) with respect to latitude

        Notes: Overflow may occur if nMax > 20 , especially for high-latitudes.
        Use MAG_PcupHigh for large nMax.

   Written by Manoj Nair, June, 2009 . Manoj.C.Nair@Noaa.Gov.

  Note: In geomagnetism, the derivatives of ALF are usually found with
  respect to the colatitudes. Here the derivatives are found with respect
  to the latitude. The difference is a sign reversal for the derivative of
  the Associated Legendre Functions.
 */
{
    int n, m, index, index1, index2, NumTerms;
    double k, z;
    rvPcup[0] = 1.0;
    rvdPcup[0] = 0.0;
    /*sin (geocentric latitude) - sin_phi */
    z = sqrt((1.0 - x) * (1.0 + x));

    NumTerms = ((nMax + 1) * (nMax + 2) / 2);
    
    //std::vector SchmidtQuasiNorm;
    std::vector<double> vdSchmidtQuasiNorm;
    vdSchmidtQuasiNorm.resize(NumTerms + 1);

  
    /*	 First,	Compute the Gauss-normalized associated Legendre  functions*/
    for (n = 1; n <= nMax; n++)
    {
        for (m = 0; m <= n; m++)
        {
            index = (n * (n + 1) / 2 + m);
            if (n == m)
            {
                index1 = (n - 1) * n / 2 + m - 1;
                rvPcup[index] = z * rvPcup[index1];
                rvdPcup[index] = z * rvdPcup[index1] + x * rvPcup[index1];
            }
            else if (n == 1 && m == 0)
            {
                index1 = (n - 1) * n / 2 + m;
                rvPcup[index] = x * rvPcup[index1];
                rvdPcup[index] = x * rvdPcup[index1] - z * rvPcup[index1];
            }
            else if (n > 1 && n != m)
            {
                index1 = (n - 2) * (n - 1) / 2 + m;
                index2 = (n - 1) * n / 2 + m;
                if (m > n - 2)
                {
                    rvPcup[index] = x * rvPcup[index2];
                    rvdPcup[index] = x * rvdPcup[index2] - z * rvPcup[index2];
                }
                else
                {
                    k = (double)(((n - 1) * (n - 1)) - (m * m)) / (double)((2 * n - 1) * (2 * n - 3));
                    rvPcup[index] = x * rvPcup[index2] - k * rvPcup[index1];
                    rvdPcup[index] = x * rvdPcup[index2] - z * rvPcup[index2] - k * rvdPcup[index1];
                }
            }
        }
    }
    /* Compute the ration between the the Schmidt quasi-normalized associated Legendre
     * functions and the Gauss-normalized version. */

    vdSchmidtQuasiNorm[0] = 1.0;
    for (n = 1; n <= nMax; n++)
    {
        index = (n * (n + 1) / 2);
        index1 = (n - 1) * n / 2;
        /* for m = 0 */
        vdSchmidtQuasiNorm[index] = vdSchmidtQuasiNorm[index1] * (double)(2 * n - 1) / (double)n;

        for (m = 1; m <= n; m++)
        {
            index = (n * (n + 1) / 2 + m);
            index1 = (n * (n + 1) / 2 + m - 1);
            vdSchmidtQuasiNorm[index] = vdSchmidtQuasiNorm[index1] * sqrt((double)((n - m + 1) * (m == 1 ? 2 : 1)) / (double)(n + m));
        }
    }

    /* Converts the  Gauss-normalized associated Legendre
              functions to the Schmidt quasi-normalized version using pre-computed
              relation stored in the variable schmidtQuasiNorm */

    for (n = 1; n <= nMax; n++)
    {
        for (m = 0; m <= n; m++)
        {
            index = (n * (n + 1) / 2 + m);
            rvPcup[index] = rvPcup[index] * vdSchmidtQuasiNorm[index];
            rvdPcup[index] = -rvdPcup[index] * vdSchmidtQuasiNorm[index];
            /* The sign is changed since the new WMM routines use derivative with respect to latitude
            insted of co-latitude */
        }
    }

    return S_OK;
} /*MAG_PcupLow */

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
HRESULT CDeclinationProc::zMAG_PcupHigh(std::vector<double>& rvPcup, std::vector<double>& rvdPcup, double x, int nMax)

/*	This function evaluates all of the Schmidt-semi normalized associated Legendre
        functions up to degree nMax. The functions are initially scaled by
        10^280 sin^m in order to minimize the effects of underflow at large m
        near the poles (see Holmes and Featherstone 2002, J. Geodesy, 76, 279-299).
        Note that this function performs the same operation as MAG_PcupLow.
        However this function also can be used for high degree (large nMax) models.

        Calling Parameters:
                INPUT
                        nMax:	 Maximum spherical harmonic degree to compute.
                        x:		cos(colatitude) or sin(latitude).

                OUTPUT
                        Pcup:	A vector of all associated Legendgre polynomials evaluated at
                                        x up to nMax. The lenght must by greater or equal to (nMax+1)*(nMax+2)/2.
                  dPcup:   Derivative of Pcup(x) with respect to latitude

                CALLS : none
        Notes:

  Adopted from the FORTRAN code written by Mark Wieczorek September 25, 2005.

  Manoj Nair, Nov, 2009 Manoj.C.Nair@Noaa.Gov

  Change from the previous version
  The prevous version computes the derivatives as
  dP(n,m)(x)/dx, where x = sin(latitude) (or cos(colatitude) ).
  However, the WMM Geomagnetic routines requires dP(n,m)(x)/dlatitude.
  Hence the derivatives are multiplied by sin(latitude).
  Removed the options for CS phase and normalizations.

  Note: In geomagnetism, the derivatives of ALF are usually found with
  respect to the colatitudes. Here the derivatives are found with respect
  to the latitude. The difference is a sign reversal for the derivative of
  the Associated Legendre Functions.

  The derivatives can't be computed for latitude = |90| degrees.
 */
{
    double pm2, pm1, pmm, plm, rescalem, z, scalef;
    int k, kstart, m, n, NumTerms;
    NumTerms = ((nMax + 1) * (nMax + 2) / 2);

    if (fabs(x) == 1.0)
    {
        wprintf_s(L"Error in PcupHigh: derivative cannot be calculated at poles");
        return E_FAIL;
    }

    std::vector<double> pdaF1, pdaF2, pdaPreSqr;
    
    pdaF1.resize(NumTerms + 1);
    pdaF2.resize(NumTerms + 1);
    pdaPreSqr.resize(NumTerms + 1);

    scalef = 1.0e-280;

    for (n = 0; n <= 2 * nMax + 1; ++n)
    {
        pdaPreSqr[n] = sqrt((double)(n));
    }

    k = 2;

    for (n = 2; n <= nMax; n++)
    {
        k = k + 1;
       
        pdaF1[k] = (double)(2 * n - 1) / (double)(n);
        pdaF2[k] = (double)(n - 1) / (double)(n);
        for (m = 1; m <= n - 2; m++)
        {
            k = k + 1;
            pdaF1[k] = (double)(2 * n - 1) / pdaPreSqr[n + m] / pdaPreSqr[n - m];
            pdaF2[k] = pdaPreSqr[n - m - 1] * pdaPreSqr[n + m - 1] / pdaPreSqr[n + m] / pdaPreSqr[n - m];
        }
        k = k + 2;
    }

    /*z = sin (geocentric latitude) */
    z = sqrt((1.0 - x) * (1.0 + x));
    pm2 = 1.0;
    rvPcup[0] = 1.0;
    rvdPcup[0] = 0.0;
    if (nMax == 0)
        return E_FAIL;
    pm1 = x;
    rvPcup[1] = pm1;
    rvdPcup[1] = z;
    k = 1;

    for (n = 2; n <= nMax; n++)
    {
        k = k + n;
        plm = pdaF1[k] * x * pm1 - pdaF2[k] * pm2;
        rvPcup[k] = plm;
        rvdPcup[k] = (double)(n) * (pm1 - x * plm) / z;
        pm2 = pm1;
        pm1 = plm;
    }

    pmm = pdaPreSqr[2] * scalef;
    rescalem = 1.0 / scalef;
    kstart = 0;

    for (m = 1; m <= nMax - 1; ++m)
    {
        rescalem = rescalem * z;

        /* Calculate Pcup(m,m)*/
        kstart = kstart + m + 1;
        pmm = pmm * pdaPreSqr[2 * m + 1] / pdaPreSqr[2 * m];
        rvPcup[kstart] = pmm * rescalem / pdaPreSqr[2 * m + 1];
        rvdPcup[kstart] = -((double)(m)*x * rvPcup[kstart] / z);
        pm2 = pmm / pdaPreSqr[2 * m + 1];
        /* Calculate Pcup(m+1,m)*/
        k = kstart + m + 1;
        pm1 = x * pdaPreSqr[2 * m + 1] * pm2;
        rvPcup[k] = pm1 * rescalem;
        rvdPcup[k] = ((pm2 * rescalem) * pdaPreSqr[2 * m + 1] - x * (double)(m + 1) * rvPcup[k]) / z;
        /* Calculate Pcup(n,m)*/
        for (n = m + 2; n <= nMax; ++n)
        {
            k = k + n;
            plm = x * pdaF1[k] * pm1 - pdaF2[k] * pm2;
            rvPcup[k] = plm * rescalem;
            rvdPcup[k] = (pdaPreSqr[n + m] * pdaPreSqr[n - m] * (pm1 * rescalem) - (double)(n)*x * rvPcup[k]) / z;
            pm2 = pm1;
            pm1 = plm;
        }
    }

    /* Calculate Pcup(nMax,nMax)*/
    rescalem = rescalem * z;
    kstart = kstart + m + 1;
    pmm = pmm / pdaPreSqr[2 * nMax];
    rvPcup[kstart] = pmm * rescalem;
    rvdPcup[kstart] = -(double)(nMax)*x * rvPcup[kstart] / z;

    return S_OK;
} /* MAG_PcupHigh */

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
HRESULT  CDeclinationProc::zMAG_Summation(MAGtype_LegendreFunction* pLegendreFunction, MAGtype_MagneticModel* pMagneticModel, MAGtype_SphericalHarmonicVariables SphVariables,
                 MAGtype_CoordSpherical CoordSpherical, MAGtype_MagneticResults* pMagneticResults)
{
    /* Computes Geomagnetic Field Elements X, Y and Z in Spherical coordinate system using
    spherical harmonic summation.


    The vector Magnetic field is given by -grad V, where V is Geomagnetic scalar potential
    The gradient in spherical coordinates is given by:

                     dV ^     1 dV ^        1     dV ^
    grad V = -- r  +  - -- t  +  -------- -- p
                     dr       r dt       r sin(t) dp


    INPUT :  LegendreFunction
                    MagneticModel
                    SphVariables
                    CoordSpherical
    OUTPUT : MagneticResults

    CALLS : MAG_SummationSpecial

    Manoj Nair, June, 2009 Manoj.C.Nair@Noaa.Gov
     */
    int m, n, index;
    double cos_phi;
    pMagneticResults->Bz = 0.0;
    pMagneticResults->By = 0.0;
    pMagneticResults->Bx = 0.0;
    
    for (n = 1; n <= pMagneticModel->m_nMax; n++)
    {
        for (m = 0; m <= n; m++)
        {
            index = (n * (n + 1) / 2 + m);

            /*		    nMax  	(n+2) 	  n     m            m           m
                    Bz =   -SUM (a/r)   (n+1) SUM  [g cos(m p) + h sin(m p)] P (sin(phi))
                                    n=1      	      m=0   n            n           n  */
                                    /* Equation 12 in the WMM Technical report.  Derivative with respect to radius.*/
            pMagneticResults->Bz -= SphVariables.m_vdRelativeRadiusPower[n] *
                (pMagneticModel->m_vdMain_Field_Coeff_G[index] * SphVariables.m_vdCos_mlambda[m] +
                    pMagneticModel->m_vdMain_Field_Coeff_H[index] * SphVariables.m_vdSin_mlambda[m])
                * (double)(n + 1) * pLegendreFunction->m_vdPcup[index];

            /*		  1 nMax  (n+2)    n     m            m           m
                    By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
                               n=1             m=0   n            n           n  */
                               /* Equation 11 in the WMM Technical report. Derivative with respect to longitude, divided by radius. */
            pMagneticResults->By += SphVariables.m_vdRelativeRadiusPower[n] *
                (pMagneticModel->m_vdMain_Field_Coeff_G[index] * SphVariables.m_vdSin_mlambda[m] -
                    pMagneticModel->m_vdMain_Field_Coeff_H[index] * SphVariables.m_vdCos_mlambda[m])
                * (double)(m)*pLegendreFunction->m_vdPcup[index];
            /*		   nMax  (n+2) n     m            m           m
                    Bx = - SUM (a/r)   SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
                               n=1         m=0   n            n           n  */
                               /* Equation 10  in the WMM Technical report. Derivative with respect to latitude, divided by radius. */

            pMagneticResults->Bx -= SphVariables.m_vdRelativeRadiusPower[n] *
                (pMagneticModel->m_vdMain_Field_Coeff_G[index] * SphVariables.m_vdCos_mlambda[m] +
                    pMagneticModel->m_vdMain_Field_Coeff_H[index] * SphVariables.m_vdSin_mlambda[m])
                * pLegendreFunction->m_vdDerPcup[index];
        }
    }

    cos_phi = cos(DEG2RAD(CoordSpherical.phig));
    if (fabs(cos_phi) > 1.0e-10)
    {
        pMagneticResults->By = pMagneticResults->By / cos_phi;
    }
    else
        /* Special calculation for component - By - at Geographic poles.
         * If the user wants to avoid using this function,  please make sure that
         * the latitude is not exactly +/-90. An option is to make use the function
         * MAG_CheckGeographicPoles.
         */
    {
        HRESULT hr;
        hr = zMAG_SummationSpecial(pMagneticModel, SphVariables, CoordSpherical, pMagneticResults);
        if (FAILED(hr))
            return hr;
    }

    return S_OK;
}/*MAG_Summation */

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
HRESULT  CDeclinationProc::zMAG_SummationSpecial(MAGtype_MagneticModel* pMagneticModel, MAGtype_SphericalHarmonicVariables SphVariables, MAGtype_CoordSpherical CoordSpherical, MAGtype_MagneticResults* pMagneticResults)
    /* Special calculation for the component By at Geographic poles.
    Manoj Nair, June, 2009 manoj.c.nair@noaa.gov
    INPUT: MagneticModel
               SphVariables
               CoordSpherical
    OUTPUT: MagneticResults
    CALLS : none
    See Section 1.4, "SINGULARITIES AT THE GEOGRAPHIC POLES", WMM Technical report

     */
{
    int n, index;
    double k, sin_phi, * PcupS, schmidtQuasiNorm1, schmidtQuasiNorm2, schmidtQuasiNorm3;

    PcupS = (double*)malloc((pMagneticModel->m_nMax + 1) * sizeof(double));

    if (PcupS == 0)
        return E_OUTOFMEMORY;

    PcupS[0] = 1;
    schmidtQuasiNorm1 = 1.0;

    pMagneticResults->By = 0.0;
    sin_phi = sin(DEG2RAD(CoordSpherical.phig));

    for (n = 1; n <= pMagneticModel->m_nMax; n++)
    {

        /*Compute the ration between the Gauss-normalized associated Legendre
  functions and the Schmidt quasi-normalized version. This is equivalent to
  sqrt((m==0?1:2)*(n-m)!/(n+m!))*(2n-1)!!/(n-m)!  */
        index = (n * (n + 1) / 2 + 1);
        schmidtQuasiNorm2 = schmidtQuasiNorm1 * (double)(2 * n - 1) / (double)n;
        schmidtQuasiNorm3 = schmidtQuasiNorm2 * sqrt((double)(n * 2) / (double)(n + 1));
        schmidtQuasiNorm1 = schmidtQuasiNorm2;
        if (n == 1)
        {
            PcupS[n] = PcupS[n - 1];
        }
        else
        {
            k = (double)(((n - 1) * (n - 1)) - 1) / (double)((2 * n - 1) * (2 * n - 3));
            PcupS[n] = sin_phi * PcupS[n - 1] - k * PcupS[n - 2];
        }

        /*		  1 nMax  (n+2)    n     m            m           m
                By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
                           n=1             m=0   n            n           n  */
                           /* Equation 11 in the WMM Technical report. Derivative with respect to longitude, divided by radius. */

        pMagneticResults->By += SphVariables.m_vdRelativeRadiusPower[n] *
            (pMagneticModel->m_vdMain_Field_Coeff_G[index] * SphVariables.m_vdSin_mlambda[1] -
                pMagneticModel->m_vdMain_Field_Coeff_H[index] * SphVariables.m_vdCos_mlambda[1])
            * PcupS[n] * schmidtQuasiNorm3;
    }

    if (PcupS)
        free(PcupS);
    
    return S_OK;
}/*MAG_SummationSpecial */

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
HRESULT CDeclinationProc::zMAG_SecVarSummation(MAGtype_LegendreFunction* pLegendreFunction, MAGtype_MagneticModel* pMagneticModel, MAGtype_SphericalHarmonicVariables SphVariables, MAGtype_CoordSpherical CoordSpherical, MAGtype_MagneticResults* pMagneticResults)
{
    /*This Function sums the secular variation coefficients to get the secular variation of the Magnetic vector.
    INPUT :  LegendreFunction
                    MagneticModel
                    SphVariables
                    CoordSpherical
    OUTPUT : MagneticResults

    CALLS : MAG_SecVarSummationSpecial

     */
    
    int m, n, index;
    double cos_phi;
    pMagneticModel->m_iSecularVariationUsed = 1;
    pMagneticResults->Bz = 0.0;
    pMagneticResults->By = 0.0;
    pMagneticResults->Bx = 0.0;
    for (n = 1; n <= pMagneticModel->m_nMaxSecVar; n++)
    {
        for (m = 0; m <= n; m++)
        {
            index = (n * (n + 1) / 2 + m);

            /*		    nMax  	(n+2) 	  n     m            m           m
                    Bz =   -SUM (a/r)   (n+1) SUM  [g cos(m p) + h sin(m p)] P (sin(phi))
                                    n=1      	      m=0   n            n           n  */
                                    /*  Derivative with respect to radius.*/
            pMagneticResults->Bz -= SphVariables.m_vdRelativeRadiusPower[n] *
                (pMagneticModel->m_vdSecular_Var_Coeff_G[index] * SphVariables.m_vdCos_mlambda[m] +
                    pMagneticModel->m_vdSecular_Var_Coeff_H[index] * SphVariables.m_vdSin_mlambda[m])
                * (double)(n + 1) * pLegendreFunction->m_vdPcup[index];

            /*		  1 nMax  (n+2)    n     m            m           m
                    By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
                               n=1             m=0   n            n           n  */
                               /* Derivative with respect to longitude, divided by radius. */
            pMagneticResults->By += SphVariables.m_vdRelativeRadiusPower[n] *
                (pMagneticModel->m_vdSecular_Var_Coeff_G[index] * SphVariables.m_vdSin_mlambda[m] -
                    pMagneticModel->m_vdSecular_Var_Coeff_H[index] * SphVariables.m_vdCos_mlambda[m])
                * (double)(m)*pLegendreFunction->m_vdPcup[index];
            /*		   nMax  (n+2) n     m            m           m
                    Bx = - SUM (a/r)   SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
                               n=1         m=0   n            n           n  */
                               /* Derivative with respect to latitude, divided by radius. */

            pMagneticResults->Bx -= SphVariables.m_vdRelativeRadiusPower[n] *
                (pMagneticModel->m_vdSecular_Var_Coeff_G[index] * SphVariables.m_vdCos_mlambda[m] +
                    pMagneticModel->m_vdSecular_Var_Coeff_H[index] * SphVariables.m_vdSin_mlambda[m])
                * pLegendreFunction->m_vdDerPcup[index];
        }
    }
    cos_phi = cos(DEG2RAD(CoordSpherical.phig));
    if (fabs(cos_phi) > 1.0e-10)
    {
        pMagneticResults->By = pMagneticResults->By / cos_phi;
    }
    else
        /* Special calculation for component By at Geographic poles */
    {
        HRESULT hr;
        hr = zMAG_SecVarSummationSpecial(pMagneticModel, SphVariables, CoordSpherical, pMagneticResults);
        if (FAILED(hr))
            return hr;
    }
    
    return S_OK;
} /*MAG_SecVarSummation*/

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
HRESULT CDeclinationProc::zMAG_SecVarSummationSpecial(MAGtype_MagneticModel* pMagneticModel, MAGtype_SphericalHarmonicVariables SphVariables, MAGtype_CoordSpherical CoordSpherical, MAGtype_MagneticResults* pMagneticResults)
{
/*Special calculation for the secular variation summation at the poles.


INPUT: MagneticModel
           SphVariables
           CoordSpherical
OUTPUT: MagneticResults
CALLS : none


 */
    int n, index;
    double k, sin_phi, schmidtQuasiNorm1, schmidtQuasiNorm2, schmidtQuasiNorm3;
    
    std::vector<double> pdaPcupS;
        
    pdaPcupS.resize(pMagneticModel->m_nMaxSecVar + 1);

    pdaPcupS[0] = 1;
    schmidtQuasiNorm1 = 1.0;
    
    pMagneticResults->By = 0.0;
    sin_phi = sin(DEG2RAD(CoordSpherical.phig));
    
    for (n = 1; n <= pMagneticModel->m_nMaxSecVar; n++)
    {
        index = (n * (n + 1) / 2 + 1);
        schmidtQuasiNorm2 = schmidtQuasiNorm1 * (double)(2 * n - 1) / (double)n;
        schmidtQuasiNorm3 = schmidtQuasiNorm2 * sqrt((double)(n * 2) / (double)(n + 1));
        schmidtQuasiNorm1 = schmidtQuasiNorm2;
        if (n == 1)
        {
            pdaPcupS[n] = pdaPcupS[n - 1];
        }
        else
        {
            k = (double)(((n - 1) * (n - 1)) - 1) / (double)((2 * n - 1) * (2 * n - 3));
            pdaPcupS[n] = sin_phi * pdaPcupS[n - 1] - k * pdaPcupS[n - 2];
        }
    
        /*		  1 nMax  (n+2)    n     m            m           m
                By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
                           n=1             m=0   n            n           n  */
                           /* Derivative with respect to longitude, divided by radius. */
    
        pMagneticResults->By += SphVariables.m_vdRelativeRadiusPower[n] *
            (pMagneticModel->m_vdSecular_Var_Coeff_G[index] * SphVariables.m_vdSin_mlambda[1] -
                pMagneticModel->m_vdSecular_Var_Coeff_H[index] * SphVariables.m_vdCos_mlambda[1])
            * pdaPcupS[n] * schmidtQuasiNorm3;
    }
    
    return S_OK;
}/*SecVarSummationSpecial*/


//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
void    CDeclinationProc::zMAG_RotateMagneticVector(MAGtype_CoordSpherical CoordSpherical, MAGtype_CoordGeodetic CoordGeodetic, MAGtype_MagneticResults MagneticResultsSph, MAGtype_MagneticResults* pMagneticResultsGeo)
/* Rotate the Magnetic Vectors to Geodetic Coordinates
Manoj Nair, June, 2009 Manoj.C.Nair@Noaa.Gov
Equation 16, WMM Technical report

INPUT : CoordSpherical : Data structure MAGtype_CoordSpherical with the following elements
                        double lambda; ( longitude)
                        double phig; ( geocentric latitude )
                        double r;  	  ( distance from the center of the ellipsoid)

                CoordGeodetic : Data structure MAGtype_CoordGeodetic with the following elements
                        double lambda; (longitude)
                        double phi; ( geodetic latitude)
                        double HeightAboveEllipsoid; (height above the ellipsoid (HaE) )
                        double HeightAboveGeoid;(height above the Geoid )

                MagneticResultsSph : Data structure MAGtype_MagneticResults with the following elements
                        double Bx;      North
                        double By;      East
                        double Bz;      Down

OUTPUT: MagneticResultsGeo Pointer to the data structure MAGtype_MagneticResults, with the following elements
                        double Bx;      North
                        double By;      East
                        double Bz;      Down

CALLS : none

 */
{
    double Psi;
    /* Difference between the spherical and Geodetic latitudes */
    Psi = (M_PI / 180) * (CoordSpherical.phig - CoordGeodetic.phi);

    /* Rotate spherical field components to the Geodetic system */
    pMagneticResultsGeo->Bz = MagneticResultsSph.Bx * sin(Psi) + MagneticResultsSph.Bz * cos(Psi);
    pMagneticResultsGeo->Bx = MagneticResultsSph.Bx * cos(Psi) - MagneticResultsSph.Bz * sin(Psi);
    pMagneticResultsGeo->By = MagneticResultsSph.By;

} /*MAG_RotateMagneticVector*/

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
void    CDeclinationProc::zMAG_CalculateGeoMagneticElements(MAGtype_MagneticResults* pMagneticResultsGeo, MAGtype_GeoMagneticElements* pGeoMagneticElements)

/* Calculate all the Geomagnetic elements from X,Y and Z components
INPUT     MagneticResultsGeo   Pointer to data structure with the following elements
                        double Bx;    ( North )
                        double By;	  ( East )
                        double Bz;    ( Down )
OUTPUT    GeoMagneticElements    Pointer to data structure with the following elements
                        double Decl; (Angle between the magnetic field vector and true north, positive east)
                        double Incl; Angle between the magnetic field vector and the horizontal plane, positive down
                        double F; Magnetic Field Strength
                        double H; Horizontal Magnetic Field Strength
                        double X; Northern component of the magnetic field vector
                        double Y; Eastern component of the magnetic field vector
                        double Z; Downward component of the magnetic field vector
CALLS : none
 */
{
    pGeoMagneticElements->m_dX = pMagneticResultsGeo->Bx;
    pGeoMagneticElements->m_dY = pMagneticResultsGeo->By;
    pGeoMagneticElements->m_dZ = pMagneticResultsGeo->Bz;

    pGeoMagneticElements->m_dH = sqrt(pMagneticResultsGeo->Bx * pMagneticResultsGeo->Bx + pMagneticResultsGeo->By * pMagneticResultsGeo->By);
    pGeoMagneticElements->m_dF = sqrt(pGeoMagneticElements->m_dH * pGeoMagneticElements->m_dH + pMagneticResultsGeo->Bz * pMagneticResultsGeo->Bz);
    pGeoMagneticElements->m_dDeclination = RAD2DEG(atan2(pGeoMagneticElements->m_dY, pGeoMagneticElements->m_dX));
    pGeoMagneticElements->m_dIncl = RAD2DEG(atan2(pGeoMagneticElements->m_dZ, pGeoMagneticElements->m_dH));
} /*MAG_CalculateGeoMagneticElements */

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
void    CDeclinationProc::zMAG_CalculateSecularVariationElements(MAGtype_MagneticResults MagneticVariation, MAGtype_GeoMagneticElements* pMagneticElements)
/*This takes the Magnetic Variation in x, y, and z and uses it to calculate the secular variation of each of the Geomagnetic elements.
        INPUT     MagneticVariation   Data structure with the following elements
                                double Bx;    ( North )
                                double By;	  ( East )
                                double Bz;    ( Down )
        OUTPUT   MagneticElements   Pointer to the data  structure with the following elements updated
                        double Decldot; Yearly Rate of change in declination
                        double Incldot; Yearly Rate of change in inclination
                        double Fdot; Yearly rate of change in Magnetic field strength
                        double Hdot; Yearly rate of change in horizontal field strength
                        double Xdot; Yearly rate of change in the northern component
                        double Ydot; Yearly rate of change in the eastern component
                        double Zdot; Yearly rate of change in the downward component
                        double GVdot;Yearly rate of chnage in grid variation
        CALLS : none

 */
{
    pMagneticElements->m_dXdot = MagneticVariation.Bx;
    pMagneticElements->m_dYdot = MagneticVariation.By;
    pMagneticElements->m_dZdot = MagneticVariation.Bz;
    pMagneticElements->m_dHdot = (pMagneticElements->m_dX * pMagneticElements->m_dXdot + pMagneticElements->m_dY * pMagneticElements->m_dYdot) / pMagneticElements->m_dH; /* See equation 19 in the WMM technical report */
    pMagneticElements->m_dFdot = (pMagneticElements->m_dX * pMagneticElements->m_dXdot + pMagneticElements->m_dY * pMagneticElements->m_dYdot + pMagneticElements->m_dZ * pMagneticElements->m_dZdot) / pMagneticElements->m_dF;
    pMagneticElements->m_dFdot = (pMagneticElements->m_dX * pMagneticElements->m_dXdot + pMagneticElements->m_dY * pMagneticElements->m_dYdot + pMagneticElements->m_dZ * pMagneticElements->m_dZdot) / pMagneticElements->m_dF;
    pMagneticElements->m_dDecldot = 180.0 / M_PI * (pMagneticElements->m_dX * pMagneticElements->m_dYdot - pMagneticElements->m_dY * pMagneticElements->m_dXdot) / (pMagneticElements->m_dH * pMagneticElements->m_dH);
    pMagneticElements->m_dIncldot = 180.0 / M_PI * (pMagneticElements->m_dH * pMagneticElements->m_dZdot - pMagneticElements->m_dZ * pMagneticElements->m_dHdot) / (pMagneticElements->m_dF * pMagneticElements->m_dF);
    pMagneticElements->m_dGVdot = pMagneticElements->m_dDecldot;
} /*MAG_CalculateSecularVariationElements*/

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
void CDeclinationProc::zMAG_WMMErrorCalc(double H, MAGtype_GeoMagneticElements* Uncertainty)
{
    double decl_variable, decl_constant;
    Uncertainty->m_dF = WMM_UNCERTAINTY_F;
    Uncertainty->m_dH = WMM_UNCERTAINTY_H;
    Uncertainty->m_dX = WMM_UNCERTAINTY_X;
    Uncertainty->m_dZ = WMM_UNCERTAINTY_Z;
    Uncertainty->m_dIncl = WMM_UNCERTAINTY_I;
    Uncertainty->m_dY = WMM_UNCERTAINTY_Y;
    decl_variable = (WMM_UNCERTAINTY_D_COEF / H);
    decl_constant = (WMM_UNCERTAINTY_D_OFFSET);
    Uncertainty->m_dDeclination = sqrt(decl_constant * decl_constant + decl_variable * decl_variable);
    if (Uncertainty->m_dDeclination > 180) {
        Uncertainty->m_dDeclination = 180;
    }
}

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
HRESULT CDeclinationProc::GetDeclination(double dLatitude, double dLongitude, double dHeight, double& rdDeclination, double& rdDeclinationError, SYSTEMTIME* pSt)
{
    //M_LOGGER_EXTERNALSTR(0, NULL);

    MAGtype_MagneticModel TimedMagneticModel;
    ////mauto_ptr <memmgr_new<MAGtype_MagneticModel>> pTimedMagneticModel;
    //pTimedMagneticModel = std::make_unique<MAGtype_MagneticModel>(1);
    //if (!pTimedMagneticModel)
    //    return E_OUTOFMEMORY;
    
    int num_term = m_MagneticModel.m_nMax * (m_MagneticModel.m_nMax + 1) / 2 + m_MagneticModel.m_nMax;

    TimedMagneticModel.m_vdMain_Field_Coeff_G.resize(num_term + 1);
    TimedMagneticModel.m_vdMain_Field_Coeff_H.resize(num_term + 1);
    TimedMagneticModel.m_vdSecular_Var_Coeff_G.resize(num_term + 1);
    TimedMagneticModel.m_vdSecular_Var_Coeff_H.resize(num_term + 1);
    
    MAGtype_GeoMagneticElements     GeoMagneticElements, Errors;
    MAGtype_CoordGeodetic           CoordGeodetic;
    MAGtype_CoordSpherical          CoordSpherical;

    MAGtype_Date dtCurrDate;

    if (!pSt)
    {
        SYSTEMTIME st = { 0 };
        GetSystemTime(&st);

        dtCurrDate.m_uYear = st.wYear;
        dtCurrDate.m_uMonth = st.wMonth;
        dtCurrDate.m_uDay = st.wDay;
    }
    else
    {
        dtCurrDate.m_uYear = pSt->wYear;
        dtCurrDate.m_uMonth = pSt->wMonth;
        dtCurrDate.m_uDay = pSt->wDay;
    }

#ifdef _USE_OLD_DECLINATION_TABLE_ 
    if ( dtCurrDate.m_uYear < pTimedMagneticModel->m_dEpoch )
    {
        wprintf_s(L"Invalid date. Valid range: 2020/1/1 - 2024/12/31");
        return E_FAIL;
    }

    if (dtCurrDate.m_uYear >= pTimedMagneticModel->m_dCoefficientFileEndDate)
    {
        dtCurrDate.m_uYear = 2024;
        dtCurrDate.m_uMonth = 12;
        dtCurrDate.m_uDay = 31;
    }
#else
    if (dtCurrDate.m_uYear >= TimedMagneticModel.m_dCoefficientFileEndDate || dtCurrDate.m_uYear < TimedMagneticModel.m_dEpoch)
    {
        wprintf_s(L"Invalid date. Valid range: 2020/1/1 - 2024/12/31");
        return E_FAIL;
    }
#endif
    HRESULT hr;
    hr = zMAG_DateToYear(&dtCurrDate);
    if (FAILED(hr))
        return hr;

    CoordGeodetic.lambda    = dLongitude; /* longitude */
    CoordGeodetic.phi       = dLatitude; /* geodetic latitude */
    CoordGeodetic.HeightAboveEllipsoid = dHeight; /* height above the ellipsoid (HaE) */

    zMAG_GeodeticToSpherical(&m_Ellip, CoordGeodetic, &CoordSpherical);
    zMAG_TimelyModifyMagneticModel(dtCurrDate, m_MagneticModel, TimedMagneticModel); /* Time adjust the coefficients, Equation 19, WMM Technical report */

    hr = zMAG_Geomag(&m_Ellip, CoordSpherical, CoordGeodetic, &TimedMagneticModel, &GeoMagneticElements); /* Computes the geoMagnetic field elements and their time change*/
    if (FAILED(hr))
        return hr;

    zMAG_WMMErrorCalc(GeoMagneticElements.m_dH, &Errors);

    rdDeclination = GeoMagneticElements.m_dDeclination;
    rdDeclinationError = Errors.m_dDeclination;

    return S_OK;
}