import pytest
import json
import rasterio
import numpy.testing as ntest


def raster(filename):
    '''
    Helper function to load a raster from a filename into a numpy array.
    '''
    out_data = rasterio.open(filename)
    return out_data.read(1)


@pytest.fixture
def params():
    '''
    Fixture to build parameter sets based on the paths.json file.
    '''
    with open('fixtures/paths.json') as f:
        fixtures = json.loads(f.read())
    params = []

    for fix in fixtures:
        params.append(pytest.param(raster(fixtures[fix]['result']),
                                   raster(fixtures[fix]['expected']),
                                   id=fix))
    return params


class TestingLSD():
    @pytest.mark.parametrize('result,expected', params())
    def test_basic_metrics(self, result, expected):
        ntest.assert_allclose(result, expected, rtol=1e-03)
