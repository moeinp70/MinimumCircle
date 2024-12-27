# SmallestCircle Library

The `SmallestCircle` library is designed to compute the smallest enclosing circle that covers a specified fraction of a population based on input geospatial population data.

---

## Features

- Load and process population data from NetCDF files or zipped archives.
- Support for filtering data by country or geographical bounds.
- Calculate the smallest circle covering a given fraction of the total population.
- Visualize the results with population density and circle overlay.

---

## Installation

### Requirements

The library requires Python 3.7 or higher. Install the required dependencies using:

```bash
pip install -r requirements.txt
```

### Install the Library

To install the library, run:

```bash
python setup.py install
```

---

## Data Requirements

The library operates on geospatial population data in NetCDF format. You can download the required dataset from:

- [NASA EarthData GPWv4 Catalog](https://earthdata.nasa.gov/data/catalog/sedac-ciesin-sedac-gpwv4-popcount-r11-4.11)
- Direct link: [GPWv4 Population Count NetCDF File](https://sedac.ciesin.columbia.edu/downloads/data/gpw-v4/gpw-v4-population-count-rev11/gpw-v4-population-count-rev11_totpop_2pt5_min_nc.zip)

Once downloaded, provide the path to the dataset when using the library.

---

## Usage

### 1. Load Population Data

The `load_population_data` function allows loading and preprocessing of population data.

```python
from smallestcircle import load_population_data

file_path = "path_to_population_data.zip"  # Path to the NetCDF or zip file
data, latitudes, longitudes, population_array = load_population_data(file_path, year=2020, country_code="ITA")
```

### 2. Compute the Smallest Circle

Use the `smallest_circle` function to compute the smallest circle that covers a specified fraction of the population.

```python
from smallestcircle import smallest_circle

center, radius = smallest_circle(data, population_array, latitudes, longitudes, target_population_ratio=0.5)
print(f"Center: {center}, Radius: {radius} km")
```

### 3. Visualize Results

Visualize the population density with the calculated circle overlay.

```python
from smallestcircle import plot_circle

plot_circle(population_array, center, radius, latitudes, longitudes)
```

---

## Testing

Run the test suite using:

```bash
pytest tests/
```

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## Contributing

Contributions are welcome! Please fork the repository and submit a pull request with your changes.

---

## Contact

For any inquiries or issues, please contact:

**Author:** Moein Zadeh  
**Email:** seyed.peyghambar@mail.polimi.it

