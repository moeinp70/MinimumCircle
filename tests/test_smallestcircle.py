from src.smallestcircle import load_population_data, smallest_circle

def test_load_population_data():
    # Replace "path_to_dataset.nc" with a real file path
    file_path = "path_to_dataset.nc"
    populations, latitudes, longitudes, _ = load_population_data(file_path, country_code=380)

    assert len(populations) > 0, "Populations should not be empty"

def test_smallest_circle():
    # Mock data
    populations = [100, 200, 300, 400, 500]
    latitudes = [10, 20, 30, 40, 50]
    longitudes = [10, 20, 30, 40, 50]

    best_center, best_radius = smallest_circle(populations, latitudes, longitudes, target_population_ratio=0.5)

    assert best_center is not None, "Center should not be None"
    assert best_radius > 0, "Radius should be greater than zero"
