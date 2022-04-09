# How to obtain CMIP data

### 1. Get an account at [Climate Explorer](http://climexp.knmi.nl/)

This is fairly self-explanatory, just need to sign up so you can download large datasets.

### 2. Pick from the righthand side what data collection to look at

For me this is either "Monthly CMIP3+ scenario runs" or "Monthly CMIP5 scenario runs".

### 3. Pick Model and variable

There should be a big table of models and variables, for example, scroll until you find "CSIRO Mk 3.5" in the model column, and I used the "20c3m" (historical) experiment, and then clicking the button for the `tas` variable, which means "near-surface air temperature", and is more or less equivalent to sea surface temperature (SST) over the ocean, which is what I want to analyze.

### 4. Click "Select Field" button at the very top

### 5. Selections

Now, we need to make a few selections on this screen. First, I inputted the latitude and longitude ranges I wanted to use:

* Latitude: -50N to 50N
* Longitude: 120E to 300E

Then, select "subset of the field" in order to get all the data in those ranges. I also clicked under "Considering" the option "only sea points", but it didn't seem to make a difference here. I also left the units as Kelvin, but it won't matter since the anomalies are what we care about for now.

### 6. Click "Make Time Series"

### 6.5 (Optional) Click the compute ensemble mean button to save yourself the work later (I chose not to do this)

### 7. Computing Anomalies

On the same screen, scroll down to the "Compute Anomalies" section, input your time range (1871-2020 for the whole historical experiment) and then click the "Generate anomaly field" button.

### 8. Download!

At the very bottom, you should now be able to download the netcdf files for each ensemble member.