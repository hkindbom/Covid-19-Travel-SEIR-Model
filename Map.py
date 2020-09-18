# https://stackoverflow.com/questions/60147431/how-to-put-a-label-on-a-country-with-python-cartopy
import matplotlib.pyplot as plt
import cartopy
import cartopy.io.shapereader as shpreader
import cartopy.crs as ccrs

def plot_map(country_data, compartment):
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.add_feature(cartopy.feature.LAND)
    ax.add_feature(cartopy.feature.OCEAN)
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.95)
    ax.add_feature(cartopy.feature.RIVERS)

    ax.set_extent([-150, 60, -25, 60])

    shpfilename = shpreader.natural_earth(resolution='110m',
                                          category='cultural',
                                          name='admin_0_countries')

    reader = shpreader.Reader(shpfilename)
    countries = reader.records()

    for country in countries:
        country_name = country.attributes['SOVEREIGNT']

        if country_name in country_data.keys():
            r = min(1, 2 * 1 * (country_data[country_name]))
            g = min(1, 2 * 1 * (1 - country_data[country_name]))
            b = 0
            rgb_color = (r, g, b)
            g = ax.add_geometries(country.geometry, ccrs.PlateCarree(), facecolor=rgb_color, label=country_name)
            x = country.geometry.centroid.x
            y = country.geometry.centroid.y
            ax.text(x, y, compartment + " = " + str(country_data[country_name]), color='black', size=15, ha='center', va='center', transform=ccrs.PlateCarree())
        else:
            ax.add_geometries(country.geometry, ccrs.PlateCarree(), facecolor=(1, 1, 1), label = country_name)

    plt.rcParams["figure.figsize"] = (50,50)
    plt.show()

plot_map({"Sweden": 0.8, "Norway": 0.5, "Denmark": 0.3, "Finland": 0.1}, "I")

