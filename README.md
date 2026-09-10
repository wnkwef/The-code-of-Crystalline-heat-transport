# The code of *Crystalline-heat-transport*

1. **PDF**: All data figures in this work were plotted using **OriginLab** software (Version 2024b; **OriginLab** Corporation, Northampton, MA, USA). Non-linear curve fitting for the thermal conductivity and specific heat analyses was performed using the *built-in fitting module* of **OriginLab**, with Equations (1) and (2) in the main text as the fitting functions.

   Crystal structure visualizations were generated using the open-source **VESTA** software (version 3.4.4).[^1]

   nPDF structural refinements were carried out with the open-source **PDFgui** package (version 3.1.0).[^2] Refinements under different symmetry constraints were implemented using the software’s native *symmetry constraints* feature. All refined parameters are presented in Tables S6–S8 of the *Supplementary Information*.

2. **Wigner limit**: Raman data were fitted with Lorentzian peaks, as presented in the *Supplementary Note 2*, before solving the lifetime of each Raman peak. Run **WignerLimitPlot.m** to plot temperature dependency of the peaks' lifetime. The file **LorentzianFitWL.m** would be automatically invoked.

3. **Raman_Polarimetry**: Run **MainPeaksEgTensorDecomposition.m** to plot the polarimetry in polar plots. Both **MainPeaksEgTensorInput.m** and **FitAgBgModel.m**  would be automatically invoked.

4. **Stokes_shift**: Run **EmExPlot.m** to draw T-dependent PL (photoluminescence) & PLE (photoluminescence excitation) spectra. Run **TempDepend.m** to draw T-dependent Stokes shift. **EmExData.m** to generate input variables would be automatically invoked.

5. **Transient_reflectance**: Run **TRPlot.m** to plot TR (transient reflectivity) and T-dependent fitting parameters. **TRFit.m** would be automatically invoked.

6. **Input files for .m code:**

   - **Raman_Parallel.mat & Raman_Cross.mat**: Raman data for carrying out **Wigner limit** and **Raman_Polarimetry**, measured at varied temperatures.

   - **Temperature_numeric.mat**: A numerical and a cell array recording temperatures.

   - **pks.mat**: Positions of Raman peaks in wavenumber for different temperatures.

   - **TRData.mat**: TR data to plot ultrafast dynamics at varied temperatures.

   - **SSData.mat**: PL and PLE data at varied temperatures.

8. **System Requirements:** The code has been fully tested on **MATLAB R2025a** under **Windows 11 (64-bit)**. No non-standard hardware required.

9. **Expected outputs** are MATLAB figures. **Runtimes** are a few minutes on a standard desktop computer.

10. **Instructions for Use**: Run scripts abovementioned in MATLAB. It will quantitatively reproduce our data.

    

    [^1]: Momma, K. & Izumi, F. VESTA 3 for three-dimensional visualization of crystal, volumetric and morphology data. *Journal of Applied Crystallography* **44**, 1272-1276 (2011).
    [^2]: Farrow, C. L. et al. PDFfit2 and PDFgui: computer programs for studying nanostructure in crystals. *J. Phys.:Condens. Matter* **19**, 335219 (2007).

    

    

    

    

    

