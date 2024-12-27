import time
import pytest
import numpy as np
from smallestcircle import load_population_data, smallest_circle, plot_circle

def test_load_population_data():
    """
    Test the load_population_data function with a dataset.
    """
    file_path = "C:/Users/Moein/Desktop/dataset/gpw_v4_population_count_rev11_2pt5_min.nc"

    country_code='IRn'
    year=2020
    bounds=None#[10,50,-10,50]

    population_data, latitudes, longitudes  = load_population_data(file_path, year=year, country_code=country_code,bounds=bounds)
    assert len(population_data) > 0, "Populations array should not be empty"
    assert len(latitudes) > 0, "Latitudes array should not be empty"
    assert len(longitudes) > 0, "Longitudes array should not be empty"


    start_time = time.time()

    # Run the smallest_circle function
    best_center, best_radius = smallest_circle(
        population_data, latitudes, longitudes, target_population_ratio=0.5, plot=False, details=True
    )

    print("--- %s seconds ---" % (time.time() - start_time))

    # Assert results
    assert isinstance(best_center, tuple), "Center should be a tuple of (latitude, longitude)"
    assert len(best_center) == 2, "Center should have exactly two coordinates"
    assert isinstance(best_radius, float), "Radius should be a float"
    assert best_radius > 0, "Radius should be positive"


    # Run the smallest_circle function with plot
    best_center, best_radius = smallest_circle(
        population_data, latitudes, longitudes, target_population_ratio=0.5, plot=True, details=False
    )


if __name__ == "__main__":
    test_load_population_data()