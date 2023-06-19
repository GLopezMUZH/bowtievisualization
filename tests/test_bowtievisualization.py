# %%
import pytest
from bowtievisualization.bowTieVisualization import (
    getBowTieNetworkValues,
    BowTieNetworkValues,
)


@pytest.fixture
def example_data():
    # Example input data for testing
    return {
        "IN": [1, 2, 3, 4],
        "SCC": [5, 6, 7],
        "OUT": [8, 9],
        "TENDRIL": [10],
        "OTHER": [11, 12, 13],
    }


# %%
def test_getBowTieNetworkValues(example_data):
    # Test case for BowTieNetworkValues class

    # Create example data
    nrNodesTubes = 10
    nrNodesTendrilsIn = 5
    nr_nodes_in = 2
    nrNodesSCC = 3
    nrNodesOut = 4
    nrNodesTendrilsOut = 6
    nr_nodes_OCC = 1
    connectedComponentsSizes = [1, 2, 3]

    # Create an instance of BowTieNetworkValues
    network_values = BowTieNetworkValues(
        nrNodesTubes,
        nrNodesTendrilsIn,
        nr_nodes_in,
        nrNodesSCC,
        nrNodesOut,
        nrNodesTendrilsOut,
        nr_nodes_OCC,
        connectedComponentsSizes,
    )

    # Test the getter functions
    assert network_values.get_nrNodesTubes() == nrNodesTubes
    assert network_values.get_nrNodesTendrilsIn() == nrNodesTendrilsIn
    assert network_values.get_nr_nodes_in() == nr_nodes_in
    assert network_values.get_nrNodesSCC() == nrNodesSCC
    assert network_values.get_nrNodesOut() == nrNodesOut
    assert network_values.get_nrNodesTendrilsOut() == nrNodesTendrilsOut
    assert network_values.get_nr_nodes_OCC() == nr_nodes_OCC
    assert network_values.get_connectedComponentsSizes() == connectedComponentsSizes


# %%
