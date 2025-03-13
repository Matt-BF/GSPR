import requests
import time
import sys

# Initialize a session for improved performance
session = requests.Session()


# Function to handle rate limited query
def rate_limited_query(url):
    try:
        response = session.get(url)
        response.raise_for_status()  # This will raise an HTTPError for bad responses
        return response.json()
    except requests.exceptions.HTTPError as e:
        print(f"HTTPError for URL {url}: {e}")
    except requests.exceptions.RequestException as e:
        print(f"RequestException for URL {url}: {e}")
    except UnboundLocalError:
        print(f"Lat or long missing")
    # Return None for any error
    #return e


# Rate limit settings
queries_per_minute = 5
delay_seconds = 60 / queries_per_minute

# Open input file and append to output file if it exists; create if not
with open("temp_df_for_map.tsv") as f, open("soil_grid_results.txt", "a") as j:
    next(f)  # Skip the header line
    for line in f:
        try:
            oid, latitude, longitude = line.strip().split("\t")
            print(f"working on {oid}")
        except ValueError:
            print("No lat long!")
            continue
        url = f"https://rest.isric.org/soilgrids/v2.0/properties/query?lon={longitude}&lat={latitude}&property=bdod&property=cec&property=cfvo&property=clay&property=nitrogen&property=ocd&property=ocs&property=phh2o&property=sand&property=silt&property=soc&depth=0-5cm&value=mean"

        result = rate_limited_query(url)

        if result:
            try:
                properties = result["properties"]['layers']
                prop_results = {}
                for prop in properties:
                    prop_results[f"{prop['name']} ({prop['unit_measure']['mapped_units']})"] = prop['depths'][0]['values']['mean']
                
                j.write(f"{oid}\t{str(prop_results)}\n")
                j.flush()
            except Exception as e:
                print(oid,e)
        else:
            print(f"Failed to fetch data for {oid}")
            j.write(f"{oid}\terror\n")

        time.sleep(delay_seconds)  # Pause to comply with rate limiting
