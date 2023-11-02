class OpticalMaterial:
    """
    A class for any material to be considered for its optical properties.

    Takes a tabular source file (.csv or pd.DataFrame) of form: Wavelength, n, k.
    Stores the refractive index and electric permittivity.
    """

    def __init__(self, name, source):

        """
        Initialize basic properties of the optical material based on the source.

        Parameters:
            name: str
                Name given to the optical material, e.g. "Silver_99" for optical data on 99% purity silver.
            source: .csv file, pd.DataFrame
                Source of the refractive index data. Expected a 3 column table: Wavelength, n, k.

        """
        
        import pandas as pd
        import numpy as np
        
        # Handle the name
        if isinstance(name, str):
            if name.endswith(".csv"):
                raise NameError("Invalid type for the 'name'. Must be 'str'.")
            else:
                self.name = name
        else:
            raise NameError("Invalid type for the 'name'. Must be 'str'.")


        # Handle the source
        if isinstance(source, str):
            if source.endswith('.csv'):
                self.refractive_index = pd.read_csv(source)
                self.refractive_index.columns = ["Wavelength", "n", "k"]
            else:
                raise Exception("Invalid sourcefile. Only '.csv' or 'pd.DataFrame' type accepted.")
        elif isinstance(source, pd.DataFrame):
            self.refractive_index = source
            self.refractive_index.columns = ["Wavelength", "n", "k"]
        else:
            raise Exception("Invalid sourcefile. Only '.csv' or 'pd.DataFrame' type accepted.")

        p = self.refractive_index
        self.n = p.drop(columns=['k'])
        self.k = p.drop(columns=['n'])

        self.permittivity = pd.concat([p.iloc[:,0], p.iloc[:,1]**2 - p.iloc[:,2]**2,
                                       2*p.iloc[:,1] * p.iloc[:,2]],
                                      axis=1,
                                      keys=[p.columns[0],
                                            "Real permittivity, a.u.",
                                            "Imaginary permittivity, a.u."])
        self.real_permittivity = self.permittivity.drop(columns=["Imaginary permittivity, a.u."])
        self.imaginary_permittivity = self.permittivity.drop(columns=["Real permittivity, a.u."])

    def even_xspacing(self, optical_property, left_bound: float = None,
                      right_bound: float = None, spacing: float = None):
        """
        Takes optical_property, e.g. refractive_index and returns evenly spaced data with defined spacing.
        Using spacing smaller than average spacing in the original data with throw a warning.
        """

        # identify a typical spectral spacing
        p = optical_property
        ds = p.iloc[:,0]
        ds.name = p.columns[0]
                          
        # Modify the bounds
        if left_bound == None:
            minds = min(ds)
        else:
            minds = left_bound

        if right_bound == None:
            maxds = max(ds)
        else:
            maxds = right_bound
        lds = len(ds)
        spacings = [ds.iloc[k+1]-ds.iloc[k] for k in range(lds-1)]
        avg_spacing = sum(spacings) / len(spacings)
        x_new = pd.Series(np.arange(minds, maxds, spacing),
                          name="Wavelength")
        x_additional = pd.merge(ds, x_new, how="outer").sort_values(by="Wavelength").reset_index(drop=True)

        # Warn if new spacing is smaller than old spacing
        len_spacing = len(str(spacing))
        if spacing < avg_spacing:
            print(f"""Warning: the specified value of spacing {spacing} is smaller than the average
            value of the spacing from the dataset {avg_spacing}. Interpolation may cause errors,
            depending on the difference between the new spacing and the average old spacing""")


        # Preprocess the data
        # add values separated by new spacing, with NaN as missing values if the value is not there
        new_x = pd.DataFrame(x_additional,
                             columns=['Wavelength'])
        new_empty = pd.DataFrame(columns=list(p.columns[-2:])).fillna(np.nan)
        reg_points = pd.concat([new_x, new_empty], axis=0, keys=list(p.columns))

        # Dataset with additional data to be interpolated
        extended = pd.concat([p, reg_points])

        # sort the new dataframe and drop duplicates
        extended = extended.drop_duplicates(subset='Wavelength')
        extended = extended.sort_values(by='Wavelength')

        # reset the index and drop the old one
        extended = extended.reset_index(drop=True)

        # interpolate the NaNs with dataframe.interpolate() - linear interpolation
        extended = extended.interpolate()

        # extract the values separated by 1nm
        regspace = extended[extended.iloc[:,0].isin(x_new)]

        # reset index
        regspace = regspace.reset_index(drop=True)

        return regspace

    def __repr__(self):
        return f"{self.name}"
