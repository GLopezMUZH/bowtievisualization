import pytest
from bowtievisualization.bowTieVisualization import *


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


def test_get_node_type(example_data):
    # Test cases for get_node_type function
    assert get_node_type(1, example_data) == "IN"
    assert get_node_type(6, example_data) == "SCC"
    assert get_node_type(9, example_data) == "OUT"
    assert get_node_type(10, example_data) == "TENDRIL"
    assert get_node_type(13, example_data) == "OTHER"
    assert get_node_type(14, example_data) == "OTHER"  # Test case for node not found


def test_get_node_category(example_data):
    # Test cases for get_node_category function
    assert get_node_category(1, example_data) == "SOURCE"
    assert get_node_category(6, example_data) == "SCC"
    assert get_node_category(9, example_data) == "SINK"
    assert get_node_category(10, example_data) == "TENDRIL"
    assert get_node_category(13, example_data) == "OTHER"
    assert (
        get_node_category(14, example_data) == "OTHER"
    )  # Test case for node not found


def test_get_incoming_nodes(example_data):
    # Test cases for get_incoming_nodes function
    assert get_incoming_nodes(1, example_data) == []
    assert get_incoming_nodes(6, example_data) == [1, 2, 3, 4]
    assert get_incoming_nodes(9, example_data) == [5, 6, 7]
    assert get_incoming_nodes(10, example_data) == [8, 9]
    assert get_incoming_nodes(13, example_data) == [10]
    assert get_incoming_nodes(14, example_data) == []  # Test case for node not found
