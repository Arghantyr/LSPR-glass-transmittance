class Nanoparticle:

    """
    Typical properties of nanoparticles:
    size, concentration, size distribution.
    """

    def __init__(self, name: str=None, material: OpticalMaterial=None, size: float=10,
                 stdev: float=1, concentration: float=1e18):

        """
        Initializes the Nanoparticle object.

        Parameters:
            name: str
                The name of the Nanoparticle object.
            material: OpticalMaterial
                The material the nanoparticles will be made of.
            size:
                Mean size of the nanoparticle population. [nm]
            stdev:
                Standard deviation of the nanoparticle size distribution. [nm]
            concentration:
                Total concentration of nanoparticles. [1/cm3]

        """

        self.name = name
        self.material = material
        self.size = size
        self.stdev = stdev
        self.concentration = concentration

    def modify(self, feature: str, new):

        if feature == "material":
            self.material = new
        elif feature == "size":
            self.size = new
        elif feature == "concentration":
            self.concentration = new
        elif feature == "stdev":
            self.stdev = new
        elif feature == "name":
            self.name = new
        else:
            raise AttributeError("Invalid attribute in Nanoparticle class.")

    def __repr__(self):
        return "{} nanoparticles with size {}±{} nm and concentration {:.2e} cm^-3. Size distribution: Normal.".format(
            self.material, self.size, self.stdev, self.concentration)
