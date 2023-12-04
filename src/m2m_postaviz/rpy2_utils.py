import rpy2.robjects as ro
import rpy2.robjects.packages as rpackages
from rpy2.robjects import default_converter
from rpy2.robjects import numpy2ri
from rpy2.robjects import pandas2ri


class Rconverter:
    numpy_converter = default_converter + numpy2ri.converter
    pandas_converter = default_converter + pandas2ri.converter
    r = ro.r

    def __init__(self, *packages: list) -> None:
        self.packages = packages

    def pcoa(self, pandas_df, nb_factor: int):
        metadata = pandas_df.iloc[:, :nb_factor]
        pandas_df = pandas_df.drop(pandas_df.iloc[:, :nb_factor], axis=1)

        vegan = rpackages.importr("vegan")
        with self.pandas_converter.context():
            r_df = self.pandas_converter.py2rpy(pandas_df)

        nrow = len(pandas_df.index) - 1

        dist_matrix = self.r.vegdist(r_df, "jaccard", binary=True)

        pcoa = self.r.cmdscale(dist_matrix, k=nrow, eig=True)

        asdata = self.r("as.data.frame")
        scores = self.r("scores")

        pcoa_df = asdata(scores(pcoa))

        with default_converter.context():
            pcoa_pandas_dataframe = pandas2ri.rpy2py_dataframe(pcoa_df)

        for factor in metadata:
            pcoa_pandas_dataframe.insert(0, factor, metadata[factor].values)

        return pcoa_pandas_dataframe
