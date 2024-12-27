import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import xarray as xr
import zipfile
import os

def load_population_data(file_path, year= 2020, country_code=None, bounds=None):
    """
    Load and preprocess population data from a NetCDF file or a zip archive containing a NetCDF file.

    Parameters:
        file_path (str): Path to the zip file or NetCDF file.
        country_code (int, optional): Numeric country code to filter data.
        year (int): The year for which to extract population data.


    Returns:
        tuple:
            populations (np.ndarray): Flattened array of valid population values.
            latitudes (np.ndarray): Corresponding latitudes for valid data.
            longitudes (np.ndarray): Corresponding longitudes for valid data.
            population_data (np.ndarray): 2D array of population data.
    """

    COUNTRY_MAPPING = {
        'afg': 4, '4': 4, 'afghanistan': 4, 'alb': 8, '8': 8, 'albania': 8, 'dza': 12, '12': 12, 'algeria': 12, 'asm': 16, '16': 16, 'american samoa': 16, 'and': 20, '20': 20, 'andorra': 20, 'ago': 24, '24': 24, 'angola': 24, 'atg': 28, '28': 28, 'antigua and barbuda': 28, 'aze': 31, '31': 31, 'azerbaijan': 31, 'arg': 32, '32': 32, 'argentina': 32, 'aus': 36, '36': 36, 'australia': 36, 'aut': 40, '40': 40, 'austria': 40, 'bhs': 44, '44': 44, 'bahamas': 44, 'bhr': 48, '48': 48, 'bahrain': 48, 'bgd': 50, '50': 50, 'bangladesh': 50, 'arm': 51, '51': 51, 'armenia': 51, 'brb': 52, '52': 52, 'barbados': 52, 'bel': 56, '56': 56, 'belgium': 56, 'bmu': 60, '60': 60, 'bermuda': 60, 'btn': 64, '64': 64, 'bhutan': 64, 'bol': 68, '68': 68, 'bolivia (plurinational state of)': 68, 'bih': 70, '70': 70, 'bosnia and herzegovina': 70, 'bwa': 72, '72': 72, 'botswana': 72, 'bra': 76, '76': 76, 'brazil': 76, 'blz': 84, '84': 84, 'belize': 84, 'slb': 90, '90': 90, 'solomon islands': 90, 'vgb': 92, '92': 92, 'british virgin islands': 92, 'brn': 96, '96': 96, 
        'brunei darussalam': 96, 'bgr': 100, '100': 100, 'bulgaria': 100, 'mmr': 104, '104': 104, 'myanmar': 104, 'bdi': 108, '108': 108, 'burundi': 108, 'blr': 112, '112': 112, 'belarus': 112, 'khm': 116, '116': 116, 'cambodia': 116, 'cmr': 120, '120': 120, 'cameroon': 120, 'can': 124, '124': 124, 'canada': 124, 'cpv': 132, '132': 132, 'cape verde': 132, 'cym': 136, '136': 136, 'cayman islands': 136, 'caf': 140, '140': 140, 'central african republic': 140, 'lka': 144, '144': 144, 'sri lanka': 144, 'tcd': 148, '148': 148, 'chad': 148, 'chl': 152, '152': 152, 'chile': 152, 'chn': 156, '156': 156, 'china': 156, 'twn': 158, '158': 158, 'taiwan': 158, 'col': 170, '170': 170, 'colombia': 170, 'com': 
        174, '174': 174, 'comoros': 174, 'myt': 175, '175': 175, 'mayotte': 175, 'cog': 178, '178': 178, 'congo': 178, 'cod': 180, '180': 180, 'democratic republic of the congo': 180, 'cok': 184, '184': 184, 'cook islands': 184, 'cri': 188, '188': 188, 'costa rica': 188, 'hrv': 191, '191': 191, 'croatia': 191, 'cub': 192, '192': 192, 'cuba': 192, 'cyp': 
        196, '196': 196, 'cyprus': 196, 'cze': 203, '203': 203, 'czech republic': 203, 'ben': 204, '204': 204, 'benin': 204, 'dnk': 208, '208': 208, 'denmark': 208, 'dma': 212, '212': 212, 'dominica': 212, 'dom': 214, '214': 214, 'dominican republic': 214, 'ecu': 218, '218': 218, 'ecuador': 218, 'slv': 222, '222': 222, 'el salvador': 222, 'gnq': 226, '226': 226, 'equatorial guinea': 226, 'eth': 231, '231': 231, 'ethiopia': 231, 'eri': 232, '232': 232, 'eritrea': 232, 'est': 233, '233': 233, 'estonia': 233, 'fro': 234, '234': 234, 'faeroe islands': 234, 'flk': 238, '238': 238, 'falkland islands (malvinas)': 238, 'fji': 242, '242': 242, 'fiji': 242, 'fin': 246, '246': 246, 'finland': 246, 'ala': 0, '0': 0, 'aland islands': 0, 'fra': 250, '250': 250, 'france': 250, 'guf': 254, '254': 254, 'french guiana': 254, 'pyf': 258, '258': 258, 'french polynesia': 258, 'dji': 262, '262': 262, 'djibouti': 262, 'gab': 266, '266': 266, 'gabon': 266, 'geo': 268, '268': 268, 'georgia': 268, 'gmb': 270, '270': 270, 'gambia': 270, 'pse': 275, '275': 275, 'state of palestine': 275, 'deu': 276, '276': 276, 'germany': 276, 'gha': 288, '288': 288, 'ghana': 288, 'gib': 292, '292': 292, 'gibraltar': 292, 'kir': 296, '296': 296, 'kiribati': 296, 'grc': 300, '300': 300, 'greece': 300, 'grl': 304, '304': 304, 'greenland': 304, 'grd': 308, '308': 308, 'grenada': 308, 'glp': 312, '312': 312, 'guadeloupe': 312, 'gum': 316, '316': 316, 'guam': 316, 'gtm': 320, '320': 320, 'guatemala': 320, 'gin': 324, '324': 324, 'guinea': 324, 'guy': 328, '328': 328, 'guyana': 328, 'hti': 332, '332': 332, 'haiti': 332, 'vat': 336, '336': 336, 'holy see': 336, 'hnd': 340, '340': 340, 'honduras': 340, 'hkg': 344, '344': 344, 'china hong kong special administrative region': 344, 'hun': 348, '348': 348, 'hungary': 348, 'isl': 352, '352': 352, 'iceland': 352, 'ind': 356, '356': 356, 'india': 356, 'idn': 360, '360': 360, 'indonesia': 360, 'irn': 364, '364': 364, 'iran (islamic republic of)': 364, 'irq': 368, '368': 368, 'iraq': 368, 'irl': 372, '372': 372, 'ireland': 372, 'isr': 376, '376': 376, 'israel': 376, 'ita': 380, '380': 380, 'italy': 380, 'civ': 384, '384': 384, "cote d'ivoire": 384, 'jam': 388, '388': 388, 'jamaica': 388, 'jpn': 392, '392': 392, 'japan': 392, 'kaz': 398, '398': 398, 'kazakhstan': 398, 'jor': 400, '400': 400, 'jordan': 400, 'ken': 404, '404': 404, 'kenya': 404, 'prk': 408, '408': 408, "democratic people's republic of korea": 408, 
        'kor': 410, '410': 410, 'republic of korea': 410, 'kwt': 414, '414': 414, 'kuwait': 414, 'kgz': 417, '417': 417, 'kyrgyzstan': 417, 'lao': 418, '418': 418, "lao people's democratic republic": 418, 'lbn': 422, '422': 422, 'lebanon': 422, 'lso': 426, '426': 426, 'lesotho': 426, 'lva': 428, '428': 428, 'latvia': 428, 'lbr': 430, '430': 430, 'liberia': 430, 'lby': 434, '434': 434, 'libya': 434, 'lie': 438, '438': 438, 'liechtenstein': 438, 'ltu': 440, '440': 440, 'lithuania': 440, 'lux': 442, '442': 442, 'luxembourg': 442, 'mac': 446, '446': 446, 'china macao special administrative region': 446, 'mdg': 450, '450': 450, 'madagascar': 450, 'mwi': 454, '454': 454, 'malawi': 454, 'mys': 458, '458': 458, 'malaysia': 458, 'mdv': 462, '462': 462, 'maldives': 462, 'mli': 466, '466': 466, 'mali': 466, 'mlt': 470, '470': 470, 'malta': 470, 'mtq': 474, '474': 474, 'martinique': 474, 'mrt': 478, '478': 478, 'mauritania': 478, 'mus': 480, '480': 480, 'mauritius': 480, 'mex': 484, '484': 484, 'mexico': 484, 'mco': 492, '492': 492, 'monaco': 492, 
        'mng': 496, '496': 496, 'mongolia': 496, 'mda': 498, '498': 498, 'republic of moldova': 498, 'mne': 499, '499': 499, 'montenegro': 499, 'msr': 500, '500': 500, 'montserrat': 
        500, 'mar': 504, '504': 504, 'morocco': 504, 'moz': 508, '508': 508, 'mozambique': 508, 'omn': 512, '512': 512, 'oman': 512, 'nam': 516, '516': 516, 'namibia': 516, 'nru': 520, '520': 520, 'nauru': 520, 'npl': 524, '524': 524, 'nepal': 524, 'nld': 528, '528': 528, 'netherlands': 528, 'cuw': 531, '531': 531, 'curacao': 531, 'abw': 533, '533': 533, 'aruba': 533, 'sxm': 534, '534': 534, 'sint maarten (dutch part)': 534, 'bes': 535, '535': 535, 'bonaire saint eustatius and saba': 535, 'ncl': 540, '540': 540, 'new caledonia': 540, 'vut': 548, '548': 548, 'vanuatu': 548, 'nzl': 554, '554': 554, 'new zealand': 554, 'nic': 558, '558': 558, 'nicaragua': 558, 'ner': 562, '562': 562, 'niger': 562, 
        'nga': 566, '566': 566, 'nigeria': 566, 'niu': 570, '570': 570, 'niue': 570, 'nfk': 0, 'norfolk island': 0, 'nor': 578, '578': 578, 'norway': 578, 'mnp': 580, '580': 580, 'northern mariana islands': 580, 'fsm': 583, '583': 583, 'micronesia (federated states of)': 583, 'mhl': 584, '584': 584, 'marshall islands': 584, 'plw': 585, '585': 585, 'palau': 585, 'pak': 586, '586': 586, 'pakistan': 586, 'pan': 591, '591': 591, 'panama': 591, 'png': 598, '598': 598, 'papua new guinea': 598, 'pry': 600, '600': 600, 'paraguay': 600, 'per': 604, '604': 604, 'peru': 604, 'phl': 608, '608': 608, 'philippines': 608, 'pcn': 0, 'pitcairn': 0, 'pol': 616, '616': 616, 'poland': 616, 'prt': 620, '620': 620, 'portugal': 620, 'gnb': 624, '624': 624, 'guinea-bissau': 624, 'tls': 626, '626': 626, 'timor-leste': 626, 'pri': 630, '630': 630, 'puerto rico': 630, 'qat': 634, '634': 634, 
        'qatar': 634, 'reu': 638, '638': 638, 'reunion': 638, 'rou': 642, '642': 642, 'romania': 642, 'rus': 643, '643': 643, 'russian federation': 643, 'rwa': 646, '646': 646, 'rwanda': 646, 'blm': 0, 'saint-barthelemy': 0, 'shn': 654, '654': 654, 'saint helena': 654, 'kna': 659, '659': 659, 'saint kitts and nevis': 659, 'aia': 660, '660': 660, 'anguilla': 660, 'lca': 662, '662': 662, 'saint lucia': 662, 'maf': 0, 'saint-martin (french part)': 0, 'spm': 666, '666': 666, 'saint pierre and miquelon': 666, 'vct': 670, '670': 670, 'saint vincent and the grenadines': 670, 'smr': 674, '674': 674, 'san marino': 674, 'stp': 678, '678': 678, 'sao tome and principe': 678, 'sau': 682, '682': 682, 'saudi arabia': 682, 'sen': 686, '686': 686, 'senegal': 686, 'srb': 688, '688': 688, 'serbia': 688, 'syc': 690, '690': 690, 'seychelles': 690, 'sle': 694, '694': 694, 'sierra leone': 694, 'sgp': 702, '702': 702, 'singapore': 702, 'svk': 703, '703': 703, 'slovakia': 703, 'vnm': 704, '704': 704, 'viet nam': 704, 'svn': 705, '705': 705, 'slovenia': 705, 'som': 706, '706': 706, 'somalia': 706, 'zaf': 710, '710': 710, 'south africa': 710, 'zwe': 716, '716': 716, 'zimbabwe': 716, 'esp': 724, '724': 724, 'spain': 724, 'ssd': 728, '728': 728, 'south sudan': 728, 'sdn': 0, 'sudan': 0, 'esh': 732, '732': 732, 'western sahara': 732, 'sur': 740, '740': 740, 'suriname': 740, 'sjm': 0, 'svalbard and jan mayen islands': 0, 'swz': 748, '748': 748, 'swaziland': 748, 'swe': 752, '752': 752, 'sweden': 752, 'che': 756, '756': 756, 'switzerland': 756, 'syr': 760, '760': 760, 'syrian arab republic': 760, 'tjk': 762, '762': 762, 'tajikistan': 762, 'tha': 764, '764': 764, 'thailand': 764, 'tgo': 768, '768': 768, 'togo': 768, 'tkl': 772, '772': 772, 'tokelau': 
        772, 'ton': 776, '776': 776, 'tonga': 776, 'tto': 780, '780': 780, 'trinidad and tobago': 780, 'are': 784, '784': 784, 'united arab emirates': 784, 'tun': 788, '788': 788, 'tunisia': 788, 'tur': 792, '792': 792, 'turkey': 792, 'tkm': 795, '795': 795, 'turkmenistan': 795, 'tca': 796, '796': 796, 'turks and caicos islands': 796, 'tuv': 798, '798': 
        798, 'tuvalu': 798, 'uga': 800, '800': 800, 'uganda': 800, 'ukr': 804, '804': 804, 'ukraine': 804, 'mkd': 807, '807': 807, 'the former yugoslav republic of macedonia': 807, 'egy': 818, '818': 818, 'egypt': 818, 'gbr': 826, '826': 826, 'united kingdom of great britain and northern ireland': 826, 'ggy': 0, 'guernsey': 0, 'jey': 0, 'jersey': 0, 'imn': 833, '833': 833, 'isle of man': 833, 'tza': 834, '834': 834, 'united republic of tanzania': 834, 'usa': 840, '840': 840, 'united states of america': 840, 'vir': 850, '850': 850, 'united states virgin islands': 850, 'bfa': 854, '854': 854, 'burkina faso': 854, 'ury': 858, '858': 858, 'uruguay': 858, 'uzb': 860, '860': 860, 'uzbekistan': 860, 'ven': 862, '862': 862, 'venezuela (bolivarian republic of)': 862, 'wlf': 876, '876': 876, 'wallis and futuna islands': 876, 'wsm': 882, '882': 882, 'western samoa': 882, 'yem': 887, '887': 887, 'yemen': 887, 'zmb': 894, '894': 894, 'zambia': 894, 'atf': 0, 'french southern territories': 0, 'bvt': 0, 'bouvet island': 0, 'hmd': 0, 'heard island and mcdonald islands': 0, 'iot': 0, 'british indian ocean territory': 0, 'sgs': 0, 'south georgia and the south sandwich islands': 0, 'spr': 0, 'spratly islands': 0, 'umi': 0, 'united states minor outlying islands': 0, 'kos': 0, 'kosovo': 0}


        # Map year to idyear
    if year < 2000:
        raise ValueError("Year must be 2000 or later.")
    elif 2000 <= year <= 2004:
        idyear = 0
    elif 2005 <= year <= 2009:
        idyear = 1
    elif 2010 <= year <= 2014:
        idyear = 2
    elif 2015 <= year <= 2019:
        idyear = 3
    elif year >= 2020:
        idyear = 4
    else:
        raise ValueError("Invalid year provided.")
    

    if file_path.endswith('.zip'):
        # If the file is a zip archive, extract the NetCDF file
        with zipfile.ZipFile(file_path, 'r') as z:
            nc_files = [f for f in z.namelist() if f.endswith('.nc')]
            if not nc_files:
                raise ValueError("No NetCDF files found in the zip archive.")
            nc_file_name = nc_files[0]
            
            # Extract the NetCDF file to a temporary location
            temp_dir = os.path.join(os.path.dirname(file_path), "temp_extracted")
            os.makedirs(temp_dir, exist_ok=True)
            extracted_nc_path = z.extract(nc_file_name, path=temp_dir)
            file_to_load = extracted_nc_path
    elif file_path.endswith('.nc'):
        # If the file is a NetCDF file, use it directly
        file_to_load = file_path
    else:
        raise ValueError("Unsupported file format. Please provide a .zip or .nc file.")

    # Load the NetCDF file using xarray
    data = xr.open_dataset(file_to_load)


    #data = xr.open_dataset(file_path)

    # Load population data
    population_data = data['Population Count, v4.11 (2000, 2005, 2010, 2015, 2020): 2.5 arc-minutes'].isel(raster=idyear).values

    if country_code:
        key = str(country_code).strip().lower()
        if key in COUNTRY_MAPPING:
            country_code = COUNTRY_MAPPING[key]
        else:
            raise ValueError(f"Invalid country code '{country_code}'. Please provide a valid ISOCODE, UNSDCODE, or NAME.")

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

    return population_data, latitudes, longitudes





def plot_circle(population_data, best_center, best_radius, latitudes, longitudes):
    """
    Visualize the population data with the calculated circle overlay (no Cartopy).

    Parameters:
        population_data (np.ndarray): 2D array of population data.
        best_center (tuple): Latitude and longitude of the circle's center.
        best_radius (float): Radius of the circle in kilometers.
        latitudes (np.ndarray): Latitude grid points.
        longitudes (np.ndarray): Longitude grid points.
    """

    # Convert radius from kilometers to degrees
    circle_radius_deg = best_radius / 111  # Approximation for Earth's surface

    # Set up the figure and axes
    plt.figure(figsize=(14, 8))
    ax = plt.gca()
    # Set the extent of the plot (equivalent to Cartopy's set_extent)
    ax.set_xlim([longitudes.min(), longitudes.max()])
    ax.set_ylim([latitudes.min(), latitudes.max()])

    # Plot the population data as a heatmap
    population_log = np.log1p(population_data)  # Avoid log of zero
    extent = [-180, 180, -90, 90]
    im = ax.imshow(
        population_log,
        extent=extent,
        origin="upper",
        cmap="viridis",
        interpolation="nearest",
    )

    # Add a circle for the calculated smallest circle
    circle = Circle(
        (best_center[1], best_center[0]),  # Longitude, Latitude
        circle_radius_deg,
        color="red",
        fill=False,
        linewidth=2,
        linestyle="--",
        label="Smallest Circle",
    )
    ax.add_patch(circle)

    # Mark the center point
    ax.plot(
        best_center[1],
        best_center[0],
        marker="o",
        color="blue",
        markersize=8,
        label="Center Point",
    )

    # Add labels, title, and legend
    plt.colorbar(im, label="Log(Population Density)")
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.title("Population Density with Calculated Smallest Circle")
    plt.legend(loc="lower left")
    plt.grid(True)
    plt.show()


'''

def plot_circle_cartopy(population_data, best_center, best_radius, latitudes, longitudes):
    """
    Visualize the population data with the calculated circle overlay.

    Parameters:
        population_data (np.ndarray): 2D array of population data.
        best_center (tuple): Latitude and longitude of the circle's center.
        best_radius (float): Radius of the circle in kilometers.
        latitudes (np.ndarray): Latitude grid points.
        longitudes (np.ndarray): Longitude grid points.
    """

    plt.figure(figsize=(14, 8))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_global()

    ax.set_extent([longitudes.min(), longitudes.max(), latitudes.min(), latitudes.max()], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.BORDERS, linestyle=":")
    ax.add_feature(cfeature.COASTLINE)

    # Use the full 2D population_data for visualization
    population_log = np.log1p(population_data)  # Avoid log of zero
    plt.imshow(
        population_log,
        transform=ccrs.PlateCarree(),
        extent=[-180, 180, -90, 90],
        cmap="viridis",
        origin="upper",
    )

    # Convert radius from kilometers to degrees
    circle_radius_deg = best_radius / 111
    circle = Circle(
        (best_center[1], best_center[0]),  # Longitude, Latitude
        circle_radius_deg,
        color="red",
        fill=False,
        transform=ccrs.PlateCarree(),
        linewidth=2,
        linestyle="--",
        label="Smallest Circle",
    )
    ax.add_patch(circle)

    plt.plot(
        best_center[1],
        best_center[0],
        marker="o",
        color="blue",
        markersize=8,
        transform=ccrs.PlateCarree(),
        label="Center Point",
    )

    plt.colorbar(label="Log(Population Density)")
    plt.title("Population Density with Calculated Smallest Circle")
    plt.legend(loc="lower left")
    plt.show()

'''

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



def smallest_circle(population_data, latitudes, longitudes, target_population_ratio=0.50, tolerance=0.01, max_candidates=100, details=True, plot=True):
    """
    Find the smallest circle that covers the target population fraction.

    Optimized for efficiency.
    """


        # Replace invalid data
    population_data = np.where(population_data < 0, np.nan, population_data)


    lat_grid, lon_grid = np.meshgrid(latitudes, longitudes, indexing="ij")

    # Extract valid data
    valid_mask = ~np.isnan(population_data) & (population_data > 0)
    populations = population_data[valid_mask]
    latitudes = lat_grid[valid_mask]
    longitudes = lon_grid[valid_mask]


    
    # Ensure inputs are numpy arrays
    populations = np.array(populations)
    latitudes = np.array(latitudes)
    longitudes = np.array(longitudes)

    total_population = np.sum(populations)
    target_population = target_population_ratio * total_population

    # Ensure max_candidates does not exceed the number of available data points
    num_candidates = min(max_candidates, len(populations))

    # Select top candidate centers based on population
    sorted_indices = np.argpartition(populations, -num_candidates)[-num_candidates:]
    sorted_indices = sorted_indices[np.argsort(populations[sorted_indices])[::-1]]
    candidate_centers = [(latitudes[i], longitudes[i]) for i in sorted_indices]

    best_radius = float('inf')
    best_center = None
    cumulative_population=0.0
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
        print(
            f"Minimum radius to cover at least {cumulative_population / total_population * 100:.2f}% of population: {best_radius:.2f} km")
        print(f"Population within final circle: {cumulative_population:.2f}")
        print(f"Expected population (50% of country): {target_population:.2f}")
        
    if plot:
        # Plot the result
        plot_circle(population_data, best_center, best_radius, latitudes, longitudes)
    return best_center, best_radius
