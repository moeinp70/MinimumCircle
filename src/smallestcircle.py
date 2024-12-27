import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.patches import Circle
import xarray as xr

def load_population_data(file_path, country_code=None, bounds=None):
    """
    Load and preprocess population data from a NetCDF file.

    Parameters:
        file_path (str): Path to the NetCDF file.
        country_code (int, optional): Numeric country code to filter data.
        bounds (list, optional): [lower_lat, upper_lat, lower_lon, upper_lon] to filter data by geographical bounds.

    Returns:
        tuple:
            populations (np.ndarray): Population values of valid cells.
            latitudes (np.ndarray): Latitude coordinates of valid cells.
            longitudes (np.ndarray): Longitude coordinates of valid cells.
            population_data (np.ndarray): Full 2D population data grid.
    """
    data = xr.open_dataset(file_path)

    # Load population data
    population_data = data['Population Count, v4.11 (2000, 2005, 2010, 2015, 2020): 2.5 arc-minutes'].isel(raster=4).values

    if country_code:
        # Filter by country code
        country_identifier = data['Population Count, v4.11 (2000, 2005, 2010, 2015, 2020): 2.5 arc-minutes'].isel(raster=10).values
        population_data = np.where((country_identifier == country_code) & (population_data >= 0), population_data, np.nan)

    if bounds:
        # Filter by geographical bounds
        lower_lat, upper_lat, lower_lon, upper_lon = bounds
        latitudes = data['latitude'].values
        longitudes = data['longitude'].values

        lat_mask = (latitudes >= lower_lat) & (latitudes <= upper_lat)
        lon_mask = (longitudes >= lower_lon) & (longitudes <= upper_lon)

        population_data = np.where(
            lat_mask[:, None] & lon_mask[None, :],
            population_data,
            np.nan
        )

    # Replace invalid data
    population_data = np.where(population_data < 0, np.nan, population_data)

    latitudes = data['latitude'].values
    longitudes = data['longitude'].values
    lat_grid, lon_grid = np.meshgrid(latitudes, longitudes, indexing="ij")

    # Extract valid data
    valid_mask = ~np.isnan(population_data) & (population_data > 0)
    populations = population_data[valid_mask]
    latitudes = lat_grid[valid_mask]
    longitudes = lon_grid[valid_mask]

    return populations, latitudes, longitudes, population_data



def plot_circle(population_data, best_center, best_radius, latitudes, longitudes):
    """
    Plot population density with the minimal circle overlay.

    Parameters:
        population_data (np.ndarray): Full 2D population data grid.
        best_center (tuple): Center of the circle (latitude, longitude).
        best_radius (float): Radius of the circle in kilometers.
        latitudes (np.ndarray): Latitude grid.
        longitudes (np.ndarray): Longitude grid.
    """
    plt.figure(figsize=(14, 8))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_global()

    ax.set_extent([longitudes.min(), longitudes.max(), latitudes.min(), latitudes.max()], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.BORDERS, linestyle=":")
    ax.add_feature(cfeature.COASTLINE)

    population_log = np.log1p(population_data)
    plt.imshow(population_log, transform=ccrs.PlateCarree(), extent=[-180, 180, -90, 90], cmap="viridis", origin="upper")

    circle_radius_deg = best_radius / 111
    circle = Circle((best_center[1], best_center[0]), circle_radius_deg, color="red", fill=False,
                    transform=ccrs.PlateCarree(), linewidth=2, linestyle="--", label="Valeriepieris Circle")
    ax.add_patch(circle)

    plt.plot(best_center[1], best_center[0], marker="o", color="blue", markersize=8, transform=ccrs.PlateCarree(),
             label="Center Point")

    plt.colorbar(label="Log(Population Density)")
    plt.title("Population Density with Valeriepieris Circle Covering 50% of the Population")
    plt.legend(handles=[circle], loc="lower left")
    plt.show()



def select_data(input_data, data_bounds, target_bounds=None):
    """
    Filter input data based on data bounds and optionally target bounds.

    Parameters:
        input_data (np.ndarray): 2D array of input data.
        data_bounds (list): [lower_lat, upper_lat, lower_lon, upper_lon].
        target_bounds (list, optional): [lower_lat, upper_lat, lower_lon, upper_lon].

    Returns:
        tuple:
            latitudes (np.ndarray): Array of latitude grid points.
            longitudes (np.ndarray): Array of longitude grid points.
            filtered_data (np.ndarray): Subset of input data.
            updated_bounds (list): [lower_lat, upper_lat, lower_lon, upper_lon].
    """
    rows, cols = input_data.shape
    lower_lat, upper_lat, lower_lon, upper_lon = data_bounds
    lat_step = (upper_lat - lower_lat) / rows
    lon_step = (upper_lon - lower_lon) / cols

    # Generate latitude and longitude arrays
    latitudes = np.array([lower_lat + (i + 0.5) * lat_step for i in range(rows)])
    longitudes = np.array([lower_lon + (j + 0.5) * lon_step for j in range(cols)])

    filtered_data = input_data  # Initially, use the full input data

    # Apply target bounds filtering if provided
    if target_bounds:
        lower_lat_t, upper_lat_t, lower_lon_t, upper_lon_t = target_bounds
        lat_mask = (latitudes >= lower_lat_t) & (latitudes <= upper_lat_t)
        lon_mask = (longitudes >= lower_lon_t) & (longitudes <= upper_lon_t)

        lat_indices = np.where(lat_mask)[0]
        lon_indices = np.where(lon_mask)[0]

        filtered_data = filtered_data[np.ix_(lat_indices, lon_indices)]
        latitudes = latitudes[lat_indices]
        longitudes = longitudes[lon_indices]

        updated_bounds = [
            latitudes[0], latitudes[-1],
            longitudes[0], longitudes[-1]
        ]
    else:
        updated_bounds = [lower_lat, upper_lat, lower_lon, upper_lon]

    return latitudes, longitudes, filtered_data, updated_bounds



def precompute_trig(lower_lat, upper_lat, lower_lon, upper_lon, rows, cols, centered=True):
    """
    Precompute trigonometric values for latitudes and longitudes.

    Parameters:
        lower_lat, upper_lat, lower_lon, upper_lon (float): Bounds of the grid.
        rows, cols (int): Number of grid points along each dimension.
        centered (bool): Whether to center grid points.

    Returns:
        tuple:
            latitudes (np.ndarray): Latitude grid points.
            longitudes (np.ndarray): Longitude grid points.
            sin_latitudes (np.ndarray): Sine of latitudes.
            cos_latitudes (np.ndarray): Cosine of latitudes.
    """
    latitudes = np.linspace(lower_lat, upper_lat, rows)
    longitudes = np.linspace(lower_lon, upper_lon, cols)

    if centered:
        latitudes = latitudes + (latitudes[1] - latitudes[0]) / 2
        longitudes = longitudes + (longitudes[1] - longitudes[0]) / 2

    sin_latitudes = np.sin(np.radians(latitudes))
    cos_latitudes = np.cos(np.radians(latitudes))

    return latitudes, longitudes, sin_latitudes, cos_latitudes


def compute_distances(center_lat, center_lon, latitudes, longitudes, sin_lat=None, cos_lat=None):
    """
    Compute great-circle distances between a center point and all other points.

    Optimized to avoid redundant trigonometric calculations.

    Parameters:
        center_lat, center_lon (float): Latitude and longitude of the center point.
        latitudes, longitudes (np.ndarray): Arrays of latitudes and longitudes of grid points.
        sin_lat, cos_lat (np.ndarray): Precomputed sine and cosine of latitudes.

    Returns:
        np.ndarray: Array of distances in kilometers.
    """
    sin_center_lat, cos_center_lat = np.sin(np.radians(center_lat)), np.cos(np.radians(center_lat))
    delta_lat = np.radians(latitudes - center_lat)
    delta_lon = np.radians(longitudes - center_lon)

    if sin_lat is None or cos_lat is None:
        sin_lat, cos_lat = np.sin(np.radians(latitudes)), np.cos(np.radians(latitudes))

    a = np.sin(delta_lat / 2) ** 2 + cos_center_lat * cos_lat * np.sin(delta_lon / 2) ** 2
    distances = 2 * 6371 * np.arcsin(np.sqrt(a))  # Earth's radius in km
    return distances


def expand_circle(distances, populations, target_population, sorted_indices):
    """
    Expand a circle incrementally by adding points based on sorted distances.

    Optimized to avoid redundant radius calculations.

    Parameters:
        distances (np.ndarray): Array of distances from the center.
        populations (np.ndarray): Array of population values at each grid point.
        target_population (float): Population target to cover.
        sorted_indices (np.ndarray): Indices sorted by distance.

    Returns:
        tuple:
            cumulative_population (float): Total population within the circle.
            radius (float): Radius of the circle in kilometers.
    """
    cumulative_population = 0.0

    for idx in sorted_indices:
        cumulative_population += populations[idx]
        if cumulative_population >= target_population * (1 - 0.01):
            radius = distances[idx]
            return cumulative_population, radius

    # In case the loop completes without meeting the target
    return cumulative_population, distances[sorted_indices[-1]]


def expand_circle222(distances, populations, target_population, sorted_indices):
    """
    Expand a circle incrementally by adding points based on sorted distances.

    Parameters:
        distances (np.ndarray): Array of distances from the center.
        populations (np.ndarray): Array of population values at each grid point.
        target_population (float): Population target to cover.
        sorted_indices (np.ndarray): Indices sorted by distance.

    Returns:
        tuple:
            cumulative_population (float): Total population within the circle.
            radius (float): Radius of the circle in kilometers.
    """
    cumulative_population = 0.0
    radius = 0.0

    for idx in sorted_indices:
        cumulative_population += populations[idx]
        radius = distances[idx]
        if cumulative_population >= target_population * (1 - 0.01):
            break

    return cumulative_population, radius




def smallest_circle(populations, latitudes, longitudes, target_population_ratio=0.50, tolerance=0.01, max_candidates=500, details=True, plot=True):
    """
    Find the smallest circle that covers the target population fraction.

    Optimized for efficiency.

    Parameters:
        populations (np.ndarray): Population values of valid cells.
        latitudes (np.ndarray): Latitude grid points.
        longitudes (np.ndarray): Longitude grid points.
        target_population_ratio (float): Target fraction of total population.
        tolerance (float): Precision threshold for radius refinement in km.
        max_candidates (int): Number of top-populated points to use as candidate centers.

    Returns:
        tuple:
            best_center (tuple): Center of the circle (latitude, longitude).
            best_radius (float): Radius of the circle in kilometers.
    """
    total_population = np.sum(populations)
    target_population = target_population_ratio * total_population

    # Efficiently select top candidate centers
    sorted_indices = np.argpartition(populations, -max_candidates)[-max_candidates:]
    sorted_indices = sorted_indices[np.argsort(populations[sorted_indices])[::-1]]
    candidate_centers = [(latitudes[i], longitudes[i]) for i in sorted_indices]

    best_radius = float('inf')
    best_center = None

    # Precompute trigonometric values for efficiency
    sin_lat, cos_lat = np.sin(np.radians(latitudes)), np.cos(np.radians(latitudes))

    for center_lat, center_lon in candidate_centers:
        # Compute distances from the current center
        distances = compute_distances(center_lat, center_lon, latitudes, longitudes, sin_lat, cos_lat)
        sorted_distances_indices = np.argsort(distances)

        # Expand the circle until the target population is covered
        cumulative_population, radius = expand_circle(distances, populations, target_population, sorted_distances_indices)

        # Update the best circle if the radius is smaller
        if cumulative_population >= target_population * (1 - tolerance) and radius < best_radius:
            best_radius = radius
            best_center = (center_lat, center_lon)

    # Final validation
    if best_radius == float('inf'):
        raise ValueError("Unable to find a circle meeting the population fraction criteria.")
    if details:
        print(f"Optimized Valeriepieris Circle center: {best_center}")
        print(f"Minimum radius to cover at least {cumulative_population / total_population * 100:.2f}% of population: {best_radius:.2f} km")
        print(f"Population within final circle: {cumulative_population:.2f}")
        print(f"Expected population (50% of country): {target_population:.2f}")
    
    if plot:
        # Plot the result
        plot_circle(populations, best_center, best_radius, latitudes, longitudes)

    return best_center, best_radius
