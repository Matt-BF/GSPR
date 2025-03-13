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
    # Return None for any error
    return e


# Rate limit settings
queries_per_minute = 5
delay_seconds = 60 / queries_per_minute

# Load already queried IDs
queried_ids = set()
try:
    with open("soil_grid_res_2.txt", "r") as file:
        for line in file:
            oid = line.split("\t")[0]
            queried_ids.add(oid)
except FileNotFoundError:
    # If the file does not exist, that's fine; we'll create it
    pass

# Open input file and append to output file if it exists; create if not
with open("temp_df_for_map_errors.csv") as f, open("soil_grid_res_2.txt", "a") as j:
    next(f)  # Skip the header line
    for line in f:
        oid, _, latitude, longitude = line.strip().split(",")[:4]

        # Skip this OID if we've already queried it
        if oid in queried_ids:
            continue

        url = f"https://rest.isric.org/soilgrids/v2.0/classification/query?lon={float(longitude)}&lat={float(latitude)}&number_classes=1"

        result = rate_limited_query(url)

        if result:
            try:
                wrb_class_name = result["wrb_class_name"]
                j.write(f"{oid}\t{wrb_class_name}\n")
                queried_ids.add(oid)  # Update the set of queried IDs
            except KeyError as e:
                print(f"{oid} KeyError: No 'wrb_class_name' found in response.")
        else:
            print(f"Failed to fetch data for {oid}")
            j.write(f"{oid}\terror\n")

        time.sleep(delay_seconds)  # Pause to comply with rate limiting
