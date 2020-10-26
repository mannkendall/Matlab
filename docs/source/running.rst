.. include:: ./substitutions.rst

Running |name|
==============

The main function to compute the Mann-Kendall test with the desired prewhitening method, temporal
segmentation and confidence limit is:

.. code-block:: Matlab

     result=MK_tempAggr(data_tempAgg, resolution, varargin)

For example:

.. code-block:: Matlab

    result=MK_tempAggr(data, 0.01);
    result=MK_tempAggr(data, 2, ‘PW_method’,’TFPW_WS’,’alpha_MK’,99,’alpha_ak’,95, ‘alpha_CL’,95, ‘alpha_Xhomo’,90);


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
