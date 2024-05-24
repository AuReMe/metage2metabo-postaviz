import pandas as pd

from m2m_postaviz import data_utils



def test_abundance():
    mock_cscope1 = pd.DataFrame(
        data=[
            [1, 0, 0, 0, 1, 1, 0],
            [0, 0, 1, 0, 1, 1, 0],
            [1, 0, 0, 1, 1, 1, 1],
            [0, 1, 0, 0, 0, 0, 0],
        ],
        index=["bin1", "bin3", "spec437", "bin1030"],
        columns=["CPD1", "CPD2", "CPD3", "CPD4", "CPD5", "CPD6", "CPD7"],
    )
    mock_cscope1.index.name = "smplID"
    mock_cscope2 = pd.DataFrame(
        data=[
            [1, 0, 0, 0, 1, 1, 0],
            [0, 0, 1, 0, 1, 1, 0],
            [1, 0, 0, 1, 1, 1, 1],
            [0, 1, 0, 0, 0, 0, 0],
        ],
        index=["bin1", "spec88", "bin1030", "bin502"],
        columns=["CPD1", "CPD2", "CPD3", "CPD4", "CPD5", "CPD6", "CPD7"],
    )
    mock_cscope2.index.name = "smplID"

    mock_abundance_df = pd.DataFrame(
        data=[
            [15, 25, 2, 5],
            [1, 30, 12, 2],
            [8, 0, 0, 1],
            [2, 1, 18, 0],
            [8, 2, 2, 5],
            [0, 10, 12, 2],
            [7, 2, 6, 1],
            [0, 1, 18, 0],
        ],
        index=["bin1", "spec88", "bin1030", "bin502", "bin3", "spec437", "bin999", "spec999"],
        columns=["mock1", "mock2", "mock3", "mock4"],
    )

    sample_mock = dict()
    sample_mock["mock1"] = {}
    sample_mock["mock1"]["cscope"] = mock_cscope1
    sample_mock["mock2"] = {}
    sample_mock["mock2"]["cscope"] = mock_cscope2

    normalised_mock_ab = mock_abundance_df.apply(lambda x: x / x.sum(), axis=0)
    expected_results = sample_mock["mock1"]["cscope"].T.dot(
        normalised_mock_ab.loc[normalised_mock_ab.index.isin(sample_mock["mock1"]["cscope"].index)]["mock1"]
    )
    expected_results2 = sample_mock["mock2"]["cscope"].T.dot(
        normalised_mock_ab.loc[normalised_mock_ab.index.isin(sample_mock["mock2"]["cscope"].index)]["mock2"]
    )

    expected_results.rename("mock1", inplace=True)
    expected_results2.rename("mock2", inplace=True)
    expected_df = pd.DataFrame([expected_results, expected_results2])

    norm_observed_results, non_norm_observed_results = data_utils.relative_abundance_calc(mock_abundance_df, sample_mock)

    assert expected_df.equals(
        norm_observed_results
    ), "Expected abundance dataframe from unit_test_abundance() and abundance dataframe from tested function are not equals."
    assert expected_df.iloc[1, 1] == norm_observed_results.iloc[1, 1], "dataframe.iloc on [1,1] coordinate did not match."
