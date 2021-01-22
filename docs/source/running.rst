.. include:: ./substitutions.rst

Running |name|
==============

The main function to compute the Mann-Kendall test with the desired prewhitening method, temporal
segmentation and confidence limit is:

.. code-block:: Matlab

     result=MK_tempAggr(data_tempAgg, resolution, varargin)

For example:

How to call the function:

.. code-block:: Matlab

   result=MK_tempAggr(data, 0.01);
   result=MK_tempAggr(data, 2, ‘PW_method’,’TFPW_WS’,’alpha_MK’,99,’alpha_ak’,95, ‘alpha_CL’,95, ‘alpha_Xhomo’,90);
    
An exemple without temporal aggregation:
    
.. code-block:: Matlab

   test_data=timetable(datetime(datevec([datenum([2000 01 01]):365:datenum([2010 01 01])])));
   test_data.param=[1.92 2.28 2.71 2.89 2.82 4.02 3.25 3.49 4.94 3.75 3.25]';
   test_result=MK_tempAggr(test_data,0.001)

   test_result = 

   struct with fields:

   P: 0.0077
   ss: 95
   slope: 0.1861
   UCL: 0.3052
   LCL: 0.1301

An exemple with a temporal segmentation into 4 quartals: 

.. code-block:: Matlab

    test_data_sea1=timetable(datetime(datevec([datenum([2000 01 01]):365:datenum([2018 01 01])])));
    test_data_sea2=timetable(datetime(datevec([datenum([2000 04 01]):365:datenum([2018 04 01])])));
    test_data_sea3=timetable(datetime(datevec([datenum([2000 07 01]):365:datenum([2018 07 01])])));
    test_data_sea4=timetable(datetime(datevec([datenum([2000 10 01]):365:datenum([2018 10 01])])));

    test_data_sea1.param=[ 3.20 2.92 3.95 1.80 2.45 2.70 2.22 2.10 2.08 2.21 1.93 2.15 2.03 1.82 1.94 2.24 1.67 1.34 1.61]';
    test_data_sea2.param=[3.92 2.99 4.37 3.04 3.12 4.07 3.91 3.42 2.94 3.14 2.53 2.80 2.98 2.86 3.22 2.31 2.03 1.59 2.14]';
    test_data_sea3.param=[4.56 4.13 4.31 1.83 3.22 5.06 4.39 4.13 4.06 3.20 4.01 3.62 3.78 3.61 3.42 3.65 2.39 3.01 3.03]';
    test_data_sea4.param=[4.22 4.78 2.96 3.23 2.82 2.96 3.12 3.49 2.73 2.61 3.00 2.66 3.49 2.58 2.32 2.10 2.38 2.29 2.07]';

    test_data_sea(1).obs=test_data_sea1;
    test_data_sea(2).obs=test_data_sea2;
    test_data_sea(3).obs=test_data_sea3;
    test_data_sea(4).obs=test_data_sea4;
    test_result_sea=MK_tempAggr(test_data_sea,0.001);

.. csv-table::     Results
   :header: "  ", "P", "ss", "Slope", "UCL", "LCL"
   :widths: auto
   
   "1st quartal", 0.0280, 95, "-0.0631", "-0.0272", "-0.0888"
   "2nd quartal", 0.0051, 95, "-0.0947", "-0.0452", "-0.1599"
   "3rd quartal", 0.0689, -1, "-0.0626", "-0.0168", "-0.1167"
   "4th quartal", 0.0064, 95, "-0.0626", "-0.0337", "-0.1170"
   "Year", 1.365e-6, 95, "-0.029", "-0.0304", "-0.1168"
   

Function description:
---------------------

The MK test and the Sen slope are applied on the given time granularity, temporal aggregation and
prewhitening method. Five prewhitening methods can be chosen, 3PW being the default option:

  - ``3PW`` (Collaud Coen et al., 2020): 3 prewhitening methods are applied (PW and TFPW_Y to determine the statistic significance (ss) of the MK test and the VCTFPW method to compute the Sen's slope
  -	``PW`` (prewhitened, Kulkarni and von Storch, 1995)
  -	``TFPW_Y`` (trend free PW,Yue et al., 2001)
  -	``TFPW_WS`` (trend free PW, Wang and Swail, 2001)
  -	``VCTFPW`` (variance corrected trend free PW, Wang et al., 2015)

For the PW, only ss autocorrelation are taken into account. The default ss for the MK test is taken
at 95% confidence limit. The default ss for upper and lower confidence limits is 90% of the all
intervals differences distribution. The default ss for the autocorrelation coefficient is 95%. The
default ss for the homogeneity test between temporal aggregation of the MK test is 90%.
If seasonal Mann-Kendall is applied, the yearly trend is assigned only if the results of the seasonal test are homogeneous. The default ss for the homogeneity test between temporal aggregation of the seasonal MK test is 90%.


INPUT:
******

  - ``data_tempAgg`` (n_arrays)= a unique timetable if MK test without temporal aggregation has to
    be used of structure of timetables if the seasonal MK test has to be applied. Each "season" is
    given by one timetable in the structure. The timetable should have only one field (apart the time).

  - ``resolution`` (float)= interval to determine the number of ties. It should be similar to the
    resolution of the instrument.

  -	Optional input ``varargin`` :

     * ``PW_method`` (string)=  used PW method (3PW; PW, TFPW_Y, TFPW_WS, VCTFPW). Default is 3PW.
     * ``alpha_MK`` (float)= confidence limit for Mk test in %. Default value is 95%
     * ``alpha_CL`` (float)= confidence limit for the confidence limits of the Sen's slope in %. Default
       value is 90%
     * ``alpha_Xhomo`` (float)= confidence limit for the homogeneity between seasons in %. Default value
       is 90%
     * ``alpha_ak`` (float)= confidence limit for the first lag autocorrelation in %. Default value is 95%


OUTPUT:
*******

result (table)= comprises the following fields:

    * ``result.P`` (float): probability for the statistical significance. If 3PW is applied,
      P= max(P_PW, P_TFPW_Y);
    * ``result.ss`` (float)= statistical significance:

        - ``alpha_MK`` if the test is ss at the alpha confidence level. Default=95%
        - ``0`` if the test is not ss at the alpha_MK confidence level
        - ``-1`` if the test is a TFPW_Y fals epositive at alpha_MK confidence level
        - ``-2`` if the test is a PW false positive at alpha_MK confidence level

    * ``result.slope`` (float): Sen's slope in units/y
    * ``result.UCL`` (float): upper confidence level in units/y
    * ``result.LCL`` (float): lower confidence level in units/y


**Sources:**

  - Collaud Coen et al., Effects of the prewhitening method, the time granularity and the time segmentation on the Mann-Kendall trend detection and the associated Sen's slope,  Atmos. Meas. Tech. Discuss,https://doi.org/10.5194/amt-2020-178, in review, 2020.

  - Sirois, A.: A brief and biased overview of time-series analysis of how to find that evasive trend, WMO/EMEP Workshop on Advanced Statistical Methods and Their Application to Air Quality Data Sets, Annex E., Global Atmosphere Watch No. 133, TD- No. 956, World Meteorological Organization, Geneva, Switzerland, 1998. annexe E, p. 26

  - Gilbert, R.: Statistical Methods for Environmental Pollution Monitoring, Van Nostrand Reinhold Company, New York, 1987.
    and the explanations about MULTMK/PARTMK de C. Libiseller
