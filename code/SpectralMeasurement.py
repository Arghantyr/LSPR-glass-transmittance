class SpectralMeasurement:
    """
    Basic spectral measurements: absorbance, transmittance
    """
    
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import integrate
    from scipy import stats
    
    # Get the important physical constants
    pi = np.pi

    def __init__(self, name: str=None, matrix: OpticalMaterial=None, nanoparticle: Nanoparticle=None,
                 thickness: float=1):

        """
        Presets the simulation of a spectral measurement of the matrix-nanoparticle composite of known thickness.
        """

        self.name = name
        self.thickness = thickness
        self.matrix = matrix
        self.nanoparticle = nanoparticle
        self.matrix_refractive_index = matrix.refractive_index
        self.matrix_real_permittivity = matrix.permittivity.iloc[:,1]
        self.nanoparticle_permittivity = nanoparticle.material.permittivity
        self.nanoparticle_concentration = nanoparticle.concentration
        self.nanoparticle_size = nanoparticle.size
        self.nanoparticle_stdev = nanoparticle.stdev
        self.f_abs = None
        self.f_sca = None
        self.cross_sections = None

        # Reflectance
        L = 1 - self.matrix_refractive_index['n']
        M = 1 + self.matrix_refractive_index['n']
        reflectance = 100 * abs(L/M)**2
        self.mreflectance = pd.concat([self.matrix_refractive_index.iloc[:,0], reflectance], axis=1,
                                      keys=[self.matrix_refractive_index.columns[0], "Reflectance, %"])

        # Raise exception if datasets are not aligned, i.e.: wavelength subsets are not identical.

        set_difference = matrix.permittivity.iloc[:,0] - nanoparticle.material.permittivity.iloc[:,0]
        if set_difference.any() != 0:
            print("The Matrix and Nanoparticle wavelengths are different, which may result in errors. Align the data with 'align_datasets', 'even_xspacing' or other tool to avoid incorrect results.")

        # Absorbance

        # Define f_sca function

        def f_sca(eps_m, e_1, e_2, l_0, C_tot, t):

            """
            A function determining the position and shape of the scattering cross section.

            Parameters:
                eps_m:
                    Real permittivity of the matrix. [a.u]
                e_1:
                    Real permittivity of the material of the nanoparticles. [a.u.]
                e_2:
                    Imaginary permittivity of the material of the nanoparticles. [a.u.]
                l_0:
                    Wavelength. [µm]
                C_tot:
                    Total concentration of nanoparticles. [1/cm3]
                t:
                    Sample thickness. [µm]

            """

            c_sca = ((0.8686 * (pi ** 5)) * 1e-30) / 3
            L_sca = c_sca * C_tot * t * (((eps_m ** 2) / (l_0 ** 4)) * (((e_1 - eps_m) ** 2) + (e_2 ** 2)))
            M_sca = (((e_1 + 2*eps_m) ** 2) + (e_2 ** 2))
            f_sca = L_sca / M_sca
            return f_sca

        # Define f_abs function
        def f_abs(eps_m, e_1, e_2, l_0, C_tot, t):

            """
            A function determining the position and shape of the absorption cross section.

            Parameters:
                eps_m:
                    Real permittivity of the matrix. [a.u]
                e_1:
                    Real permittivity of the material of the nanoparticles. [a.u.]
                e_2:
                    Imaginary permittivity of the material of the nanoparticles. [a.u.]
                l_0:
                    Wavelength. [µm]
                C_tot:
                    Total concentration of nanoparticles. [1/cm3]
                t:
                    Sample thickness. [µm]
            """

            c_abs = 1.3029 * (pi ** 2) * 1e-21
            L_abs = c_abs * C_tot * t * (((eps_m ** 1.5) / l_0) * e_2)
            M_abs = (((e_1 + 2*eps_m) ** 2) + (e_2 ** 2))
            f_abs = L_abs / M_abs
            return f_abs

        # Define integrals 1 and 2 - I_1
        def I_1(n, mean, sd):

            """
            Integrated n-th order power of the normally distributed diameter as part
            of the absorption and scattering cross sections of the Mie dipol approximation.

            Parameters:
                n:
                    order of the exponentiation
                mean:
                    mean size of the nanoparticle size distribution in [nm]
                sd:
                    standard deviation of the nanoparticle distribution in [nm]

            Returns:
                Integral of product:
                    normal distribution N(mean, sd)   X   d**n
            """

            integrand_1 = lambda x: (x**n)*stats.norm.pdf(x, loc=mean, scale=sd)
            integrand_2 = lambda x: stats.norm.pdf(x, loc=mean, scale=sd)
            return integrate.quad(integrand_1, 0, np.inf)[0] / integrate.quad(integrand_2, 0, np.inf)[0]

        n_glass = self.matrix_refractive_index['n']
        k_glass = self.matrix_refractive_index['k']
        eps_m = self.matrix_real_permittivity
        l_0 = self.matrix_refractive_index.iloc[:,0]
        e_1 = self.nanoparticle_permittivity.iloc[:,1]
        e_2 = self.nanoparticle_permittivity.iloc[:,2]
        C_tot = self.nanoparticle_concentration
        mean = self.nanoparticle_size
        sd = self.nanoparticle_stdev
        t = self.thickness

        # Reflection losses - approximation using matrix only
        mrefl_correction = 1 - 0.01 * self.mreflectance.iloc[:,1]

        # Define A_glass function
        L1 = k_glass * t
        A_matrix = 5.4576 * (L1 / l_0)

        # A_nano - inner function
        A_nano = f_abs(eps_m, e_1, e_2, l_0, C_tot, t) * I_1(3, mean, sd) + f_sca(eps_m, e_1, e_2, l_0, C_tot, t) * I_1(6, mean, sd)

        A_comp = A_matrix + A_nano

        self.Abs = pd.concat([l_0, A_comp], axis=1,
                             keys=[self.matrix_refractive_index.columns[0], "Absorbance, a.u."])
        self.f_abs = pd.concat([l_0, f_abs(eps_m, e_1, e_2, l_0, C_tot, t)], axis=1,
                               keys=[self.matrix_refractive_index.columns[0], "f_abs"])
        self.f_sca = pd.concat([l_0, f_sca(eps_m, e_1, e_2, l_0, C_tot, t)], axis=1,
                               keys=[self.matrix_refractive_index.columns[0], "f_sca"])

        sigma_sca = 1e4 * f_sca(eps_m, e_1, e_2, l_0, C_tot, t) * I_1(6, mean, sd) / (0.4343 * C_tot * t)
        sigma_abs = 1e4 * f_abs(eps_m, e_1, e_2, l_0, C_tot, t) * I_1(3, mean, sd) / (0.4343 * C_tot * t)
        self.cross_sections = pd.concat([l_0, sigma_abs, sigma_sca], axis=1,
                               keys=[self.matrix_refractive_index.columns[0], "Cross section for absorption, cm2",
                                     "Cross section for scattering, cm2"])

    def T(self, R_corr: bool=True):

        """
        Calculates the transmittance of the Matrix-Nanoparticle composite with the preset thickness.

        Parameters:
            R_corr:
                Reflectance correction. When True takes reflection losses into account. Uses Fresnel equation
                approximation for 90 deg incidence and transition between composite and air.

        Returns:
            T:
                Transmittance of the composite. pd.DataFrame format.
        """

        Abs = self.Abs.iloc[:,1]
        # calculates the transmittance of the sample composite
        transmittance = None

        if R_corr == True:
            # evaluate corrected transmittance
            transmittance = (1 - 0.01*self.mreflectance.iloc[:,1]) * 10**(-Abs)
        else:
            # evaluate without correction
            transmittance = 10**(-(Abs))

        transmittance = 100*transmittance
        transmittance = pd.concat([self.matrix_refractive_index.iloc[:,0], transmittance], axis=1,
                                  keys=[self.matrix_refractive_index.columns[0], "Transmittance, %"])

        return transmittance

    def plot_results(self, feature: str = "T", low: float = None, high: float = None):
        """
        Plots the specified feature with respect to the wavelength.
        """

        dane = None

        if feature == "cross_sections":
            dane = self.cross_sections
            xlabel = dane.columns[0]
            title = "Cross-sections for absorption and scattering"
            ylabel_1 = dane.columns[1]
            ylabel_2 = dane.columns[2]
            x_data = dane.iloc[:,0]
            y_1 = dane.iloc[:,1]
            y_2 = dane.iloc[:,2]

            lower_xlimit = dane.iloc[0,0]
            higher_xlimit = dane.iloc[-1,0]
            if low != None:
                lower_xlimit = low
            if high != None:
                higher_xlimit = high

            color_1 = "red"
            color_2 = "blue"

            fig = plt.figure(figsize=(10,5.625))

            plt.subplot(1,2,1)
            plt.plot(x_data, y_1, color=color_1)
            st = plt.suptitle(title, size=19)
            plt.ylabel(ylabel_1, color=color_1, size=19)
            plt.xlabel(xlabel, size=19)
            plt.yticks(size=13)
            plt.xticks(size=13)
            plt.xlim(lower_xlimit,higher_xlimit)
            plt.tick_params(axis='y', labelcolor=color_1)

            plt.subplot(1,2,2)
            plt.plot(x_data, y_2, color=color_2)
            plt.ylabel(ylabel_2, color=color_2, size=19)
            plt.xlabel(xlabel, size=19)
            plt.yticks(size=13)
            plt.xticks(size=13)
            plt.xlim(lower_xlimit,higher_xlimit)
            plt.tick_params(axis='y', labelcolor=color_2)

            plt.tight_layout()
            st.set_y(0.95)
            fig.subplots_adjust(top = 0.85)
        else:
            if feature == "T":
                dane = self.Abs
                title = "Transmittance without reflectance correction"
                xlabel = dane.columns[0]
                ylabel = "Transmittance, %"
                x_data = dane.iloc[:,0]
                y_data = 100 * 10**(-dane.iloc[:,1])

            elif feature == "T_corr":
                dane = self.Abs
                title = "Transmittance with reflectance correction"
                xlabel = dane.columns[0]
                ylabel = "Transmittance, %"
                x_data = dane.iloc[:,0]
                y_data = 100 * (1 - 0.01*self.mreflectance.iloc[:,1]) * 10**(-dane.iloc[:,1])

            elif feature == "Abs":
                dane = self.Abs
                title = "Absorbance"
                xlabel = dane.columns[0]
                ylabel = "Absorbance, a.u."
                x_data = dane.iloc[:,0]
                y_data = dane.iloc[:,1]

            elif feature == "R":
                dane = self.mreflectance
                title = "Reflectance"
                xlabel = dane.columns[0]
                ylabel = "Reflectance, %"
                x_data = dane.iloc[:,0]
                y_data = dane.iloc[:,1]


            lower_xlimit = dane.iloc[0,0]
            higher_xlimit = dane.iloc[-1,0]
            if low != None:
                lower_xlimit = low
            if high != None:
                higher_xlimit = high

            fig = plt.figure(figsize=(10,5.625))

            color = "red"
            plt.plot(x_data, y_data, color=color)
            st = plt.suptitle(title, size=19)
            plt.ylabel(ylabel, color=color, size=19)
            plt.xlabel(xlabel, size=19)
            plt.yticks(size=13)
            plt.xticks(size=13)
            plt.xlim(lower_xlimit,higher_xlimit)
            plt.tick_params(axis='y', labelcolor=color)

            plt.tight_layout()
            st.set_y(0.95)
            fig.subplots_adjust(top = 0.85)


    def align_datasets(self):

        """
        Takes the combined Wavelength sets of Matrix and Nanoparticle (outer join)
        and interpolates the missing data using linear approximation of neighbouring points.

        WARNING!
        Can produce errors when large difference in datasets is present.
        """

        # add x subsets of matrix and nanoparticle
        mat_per = self.matrix.permittivity.iloc[:,0]
        nano_per = self.nanoparticle.material.permittivity.iloc[:,0]
        aligned_x_subset = pd.concat([mat_per, nano_per]).sort_values().reset_index(drop=True)

        # interpolate data based on new x subset
        interpolated_matrix_RI = pd.merge(aligned_x_subset, self.matrix.refractive_index,
                                          how="outer").interpolate()

        interpolated_nanomaterial_RI = pd.merge(aligned_x_subset,
                                                self.nanoparticle.material.refractive_index,
                                                how="outer").interpolate()

        # Generate interpolated matrix and nanoparticle instances
        nano_name = self.nanoparticle.name
        nano_conc = self.nanoparticle.concentration
        nano_size = self.nanoparticle.size
        nano_stdev = self.nanoparticle.stdev
        inter_mat_name = "Interpolated" + self.matrix.name
        inter_nanomat_name = "Interpolated" + self.nanoparticle.material.name
        inter_nanopart_name = "Interpolated" + nano_name

        interpolated_matrix = OpticalMaterial(inter_mat_name, interpolated_matrix_RI)
        interpolated_nano_material = OpticalMaterial(inter_nanomat_name, interpolated_nanomaterial_RI)
        interpolated_nanoparticle = Nanoparticle(inter_nanopart_name, interpolated_nano_material,
                                                 size = nano_size, stdev = nano_stdev,
                                                 concentration = nano_conc)

        # initialize a copy SpectralMeasurement
        aligned_measurement = SpectralMeasurement("Aligned data measurement", matrix=interpolated_matrix,
                                                  nanoparticle=interpolated_nanoparticle)

        return aligned_measurement
