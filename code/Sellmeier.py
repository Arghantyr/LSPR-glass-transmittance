class Sellmeier:

    """
    Class for fitting Sellmeier or Forouhi-Bloomer model to refractive index experimental data.

    Available methods:

    model:
        input: x, model
            x: int, float, np.array, pd.Series
                Wavelengths for which the refractive index and extinction coefficient will be computed.
            model: str = "Sellmeier"
                Model to be used.
                Available:
                - "Sellmeier" - a 3rd order Sellmeier formula
                - "Forouhi-Bloomer" - Forouhi-Bloomer model for glass
        output: pd.DataFrame with columns "Wavelength", "n", "k"

    fit_data:
        input: model
            model: str = "Sellmeier"
                Model to be used for fitting.
                Available:
                - "Sellmeier" - a 3rd order Sellmeier formula
                - "Forouhi-Bloomer" - Forouhi-Bloomer model for glass

        output: Constants for the specified model fitting
            - Sellmeier constants for n: B_i, C_i (i=1...3) and for k: b_i (i=1...3)
            - Forouhi-Bloomer contants: n_inf, Eg, A, B, C
    """
    
    # Get the important physical constants
    h = constants.physical_constants["Planck constant"][0]
    c = constants.physical_constants["speed of light in vacuum"][0]
    Jev = constants.physical_constants["joule-electron volt relationship"][0]


    def __init__(self, data):
        """
        data:
            A table of shape (a,b), where b=3;
            1st column - wavelength, 2nd column - real refractive index, 3rd column - imaginary refractive index.
            Available formats: pd.DataFrame, np.array, .csv sourcefile
        """

        self.constants = {}

        # READ the input data
        typ = type(data)
        d = None

        if typ == pd.DataFrame:
            d = data
        elif typ == np.ndarray:
            d = pd.Dataframe(data, columns=["Wavelength", "n", "k"])
        elif ((typ == str)&(data[-4:] == ".csv")):
            d = pd.read_csv(data)

        d.columns = ["Wavelength", "n", "k"]

        self.data = d

        # Add "Energy, eV" column - 1e6 scaling factor for assumed wavelength in microns
        self.data["Wavelength, eV"] = 1e6 * h*c*Jev / self.data["Wavelength"]

        # Initiate an empty dict with tables of constants

        # Initiate Sellmeier constants table with 0's
        terms = []

        for j in range(1,4):
            symb = ("B", "C", "b")
            for k in symb:
                terms.append(k+str(j))

        self.constants["Sellmeier"] = pd.Series([0 for _ in range(6)], name="Value").reindex(terms)


        # Initiate Forouhi-Bloomer constants table with 0's
        terms = ["n_inf", "Eg", "A", "B", "C"]

        self.constants["Forouhi-Bloomer"] = pd.Series([0 for _ in range(4)], name="Value").reindex(terms)


    def model(self, x, model: str = "Sellmeier"):

        if model == "Sellmeier":
            """
            A 3rd order Sellmeier function.

            Parameters:
                x: int, float, list, numpy.ndarray, pandas.Series, pandas.DataFrame

            Returns:
                Refractive index as fitted with the 3rd order Sellmeier function:
                - pd.DataFrame

            """
            B1 = self.constants["Sellmeier"]["B1"]
            C1 = self.constants["Sellmeier"]["C1"]
            B2 = self.constants["Sellmeier"]["B2"]
            C2 = self.constants["Sellmeier"]["C2"]
            B3 = self.constants["Sellmeier"]["B3"]
            C3 = self.constants["Sellmeier"]["C3"]
            b1 = self.constants["Sellmeier"]["b1"]
            b2 = self.constants["Sellmeier"]["b2"]
            b3 = self.constants["Sellmeier"]["b3"]

            # READ the input data
            typ = type(x)

            if typ == pd.DataFrame:
                x = x.iloc[:,0]
            elif typ == np.ndarray:
                x = pd.DataFrame(x).iloc[:,0]
            elif typ == list:
                x = pd.Series(x, name="Wavelength")
            elif ((typ == int)|(typ == float)):
                x = pd.DataFrame(x, columns="Wavelength")

            Sellm_n = (1 + (B1 * x**2)/(x**2 - C1) + (B2 * x**2)/(x**2 - C2) + (B3 * x**2)/(x**2 - C3))**0.5
            Sellm_k = (Sellm_n*((b1*x)+(b2/x)+(b3/(x**3))))**(-1)

            return pd.concat([x, Sellm_n, Sellm_k], axis=1, keys=["Wavelength", "n", "k"])

        elif model == "Forouhi-Bloomer":
            """
            Fit Forouhi-Bloomer model based on 5 parameters:

            n_inf     -   refractive index at infinity energy
            Eg        -   energy band gap
            A, B, C   -   energy band structure-dependent constants
            """

            n_inf = self.constants["Forouhi-Bloomer"]["n_inf"]
            Eg = self.constants["Forouhi-Bloomer"]["Eg"]
            A = self.constants["Forouhi-Bloomer"]["A"]
            B = self.constants["Forouhi-Bloomer"]["B"]
            C = self.constants["Forouhi-Bloomer"]["C"]

            # READ the input data
            typ = type(x)

            if typ == pd.DataFrame:
                x = x.iloc[:,0]
            elif typ == np.ndarray:
                x = pd.DataFrame(x).iloc[:,0]
            elif typ == list:
                x = pd.Series(x, name="Wavelength")
            elif ((typ == int)|(typ == float)):
                x = pd.DataFrame(x, columns="Wavelength")

            # Transform wavelengths into energies.
            # WARNING! Assumption is made that wavelength is in um - scaling factor of 1e6
            En = 1e6 * h*c*Jev / x

            # Initiate Forouhi-Bloomer complex constants
            Q = 0.5 * (4*C - B**2)**0.5
            B_0 = (A / Q) * (-0.5*B**2 + Eg*B - Eg**2 + C)
            C_0 = (A / Q) * (0.5*B * (Eg**2 + C) - 2*Eg*C)


            FH_n = n_inf + (B_0*En + C_0) / (En**2 - B*En - C)
            FH_k = (A * (En - Eg)**2) / (En**2 - B*En + C)

            return pd.concat([x, FH_n, FH_k], axis=1, keys=["Wavelength", "n", "k"])

        else:
            raise Exception("Invalid model. Expected 'Sellmeier' or 'Forouhi-Bloomer'." )

    def fit_data(self, model: str = "Sellmeier"):
        """
        Fit source data according to the specified model.

        Parameters:
        -----------

        model: str = "Sellmeier"
            Model to be used t=for fitting n, k source data. Accepts two values: "Sellmeier or "Forouhi-Bloomer".
            Defaults to 3rd order Sellmeier for n and k.

        Returns:
        -----------

        Fitting constants for the specified model:
            Sellmeier:
                B_i, C_i   1=(1...3) - refractive index
                b_i        i=(1...3) - extinction coefficient

            Forouhi-Bloomer:
                n_inf     -   refractive index at infinite energy
                Eg        -   energy band gap of the glass
                A, B, C   -   band structure dependent constants
        """

        dane = self.data

        # Preprocess the data
        dane_n = dane.drop(columns="k").dropna().reset_index(drop=True)
        dane_k = dane.drop(columns="n").dropna().reset_index(drop=True)

        if model == "Sellmeier":

            # Here define the Sellmeier_n function
            def Sellmeier_n(x, B1, C1, B2, C2, B3, C3):
                """
                A 3rd order Sellmeier function.

                Parameters:
                    x: int, float, np.array, pd.Series
                    B,C = constants of the Sellmeier equation for real part of the refractive index

                Returns:
                The output of the 3rd order Sellmeier function:
                        - a pd.DataFrame

                """

                Sellm_n = 1 + (B1 * x**2)/(x**2 - C1) + (B2 * x**2)/(x**2 - C2) + (B3 * x**2)/(x**2 - C3)

                return Sellm_n

            # Run optimization on the glass data
            """
            The optimizing parameters are based on Sellmeier constants for existing glasses.
            p0 initializing values are set on the order of the constants for fused silica.
            """

            popt_n, pcov_n = optimize.curve_fit(Sellmeier_n, dane_n.iloc[:,0], dane_n.iloc[:,1]**2,
                                                bounds = (0, [2., 1., 10., 0.01, 0.1, 200.]),
                                                p0 = [0.5, 0.5, 0.7, 0.005, 0.01, 100])


            # fill the fitted constants (popt) to the table of constants
            self.constants["Sellmeier"]["B1"] = popt_n[0]
            self.constants["Sellmeier"]["C1"] = popt_n[1]
            self.constants["Sellmeier"]["B2"] = popt_n[2]
            self.constants["Sellmeier"]["C2"] = popt_n[3]
            self.constants["Sellmeier"]["B3"] = popt_n[4]
            self.constants["Sellmeier"]["C3"] = popt_n[5]

            def Sellmeier_k(x, b1, b2, b3, B1 = popt_n[0], C1 = popt_n[1], B2 = popt_n[2], C2 = popt_n[3],
                            B3 = popt_n[4], C3 = popt_n[5]):

                """
                Fits the 3rd order Sellmeier formula for the extinction part of the refractive index.
                """

                n = (1 + (B1 * x**2)/(x**2 - C1) + (B2 * x**2)/(x**2 - C2) + (B3 * x**2)/(x**2 - C3))**0.5
                Sellm_k = (n*((b1*x)+(b2/x)+(b3/(x**3))))**(-1)

                return Sellm_k

            _Sellmeier_k = lambda x, b1, b2, b3: Sellmeier_k(x, b1, b2, b3)

            popt_k, pcov_k = optimize.curve_fit(_Sellmeier_k ,dane_k.iloc[:,0], dane_k.iloc[:,1], p0 = [0.5, 0.5, 0.5])

            # fill the fitted constants (popt) to the table of constants
            self.constants["Sellmeier"]["b1"] = popt_k[0]
            self.constants["Sellmeier"]["b2"] = popt_k[1]
            self.constants["Sellmeier"]["b3"] = popt_k[2]

            return self.constants["Sellmeier"]

        elif model == "Forouhi-Bloomer":

            # Define Forouhi-Bloomer formula for modelling k
            def FH_k(x, Eg, A, B, C):

                k = (A * (x - Eg)**2) / (x**2 - B*x + C)

                return k

            # Optimize the fitting parameters
            FH_k_parameters, FHk_cov = optimize.curve_fit(FH_k, dane_k.iloc[:,2], dane_k.iloc[:,1],
                                                          bounds=(0,[10, 2, 75, 300]),
                                                          p0=[2.5, 0.02, 10, 25])

            Eg = FH_k_parameters[0]
            A = FH_k_parameters[1]
            B = FH_k_parameters[2]
            C = FH_k_parameters[3]

            self.constants["Forouhi-Bloomer"]["Eg"] = Eg
            self.constants["Forouhi-Bloomer"]["A"] = A
            self.constants["Forouhi-Bloomer"]["B"] = B
            self.constants["Forouhi-Bloomer"]["C"] = C

            # Define the Forouhi-Bloomer formula for modelling n
            def FH_n(x, n_inf, Eg = Eg, A = A, B = B, C = C):

                Q = 0.5 * (4*C - B**2)**0.5
                B_0 = (A / Q) * (-0.5*B**2 + Eg*B - Eg**2 + C)
                C_0 = (A / Q) * (0.5*B * (Eg**2 + C) - 2*Eg*C)


                n = n_inf + (B_0*x + C_0) / (x**2 - B*x - C)

                return n

            # Optimize the fitting parameters
            FH_n_parameters, FHn_cov = optimize.curve_fit(lambda x, n_inf: FH_n(x, n_inf, Eg = Eg, A = A, B = B, C = C),
                                                          dane_n.iloc[:,2], dane_n.iloc[:,1])

            self.constants["Forouhi-Bloomer"]["n_inf"] = FH_n_parameters[0]

            return self.constants["Forouhi-Bloomer"]

        else:
            raise Exception("Invalid fitting model. Must be 'Sellmeier' or 'Forouhi-Bloomer'.")

    def plot_fit(self, model: str="Sellmeier"):

        """
        Plots the fitting results with experimental data as points and fitted model as a solid line.

        Parameters:
        -----------

        model: str="Sellmeier"
            Model to be used for fitting and plotted.
            Accepted values: "Sellmeier", "Forouhi-Bloomer"
            Defaults to "Sellmeier".

        Returns:
        -----------

        A plot with fitted model and experimental data, where:
        ___   fitted model
        -o-   experimental data

        """

        plt.figure(figsize = (12,7))

        x_min = self.data.iloc[0,0]
        x_max = self.data.iloc[-1,0]
        ymin_k = self.data["k"].min() * 1.1
        ymax_k = self.data["k"].max() * 1.1
        x = np.linspace(x_min, x_max, 1000)
        fit = self.model(x, model=model)

        # plot refractive index
        plt.subplot(1,2,1)
        experimental_n = plt.scatter(self.data["Wavelength"], self.data["n"], alpha=0.5, s=20)
        plt.title("{} fit for refractive index".format(model), size=13)
        plt.plot(fit["Wavelength"], fit["n"], c="black")
        plt.xlabel("Wavelength", size=19)
        plt.ylabel("Refractive index, a.u.", size=19)
        plt.xticks(size=13)
        plt.yticks(size=13)

        # plot the extinction coefficient
        plt.subplot(1,2,2)

        experimental_k = plt.scatter (self.data["Wavelength"], self.data["k"], alpha=0.5, s=20)
        plt.title("{} fit for extinction coefficient".format(model), size=13)
        plt.plot(fit["Wavelength"], fit["k"], c="black")
        plt.ylim(ymin_k, ymax_k)
        plt.xlabel("Wavelength", size=19)
        plt.ylabel("Extinction coefficient, a.u.", size=19)
        plt.xticks(size=13)
        plt.yticks(size=13)
        plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))

        plt.tight_layout()
