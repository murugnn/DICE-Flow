from pathlib import Path
import pandas as pd
import pytest

FIXTURES = Path(__file__).parent / "fixtures"
SEED = 20260507


@pytest.fixture(scope="session")
def tiny_graph_path():
    return FIXTURES / "tiny_graph.tsv"


@pytest.fixture(scope="session")
def tiny_expr_normal():
    return pd.read_csv(FIXTURES / "tiny_expr_normal.tsv", sep="\t", index_col=0)


@pytest.fixture(scope="session")
def tiny_expr_tumor():
    return pd.read_csv(FIXTURES / "tiny_expr_tumor.tsv", sep="\t", index_col=0)


@pytest.fixture(autouse=True)
def _seed():
    import random
    import numpy as np
    random.seed(SEED)
    np.random.seed(SEED)
