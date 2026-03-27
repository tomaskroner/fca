# FCA.py -- Floating Catchment Area

User manual\
March 21, 2026\
Author: Tomáš Kröner

------------------------------------------------------------------------

# FCA.py

Python script for calculating spatial accessibility using Floating
Catchment Area (FCA) methods.

## Overview

`FCA.py` is a Python script designed to compute **spatial accessibility
of services** using methods from the **Floating Catchment Area (FCA)**
family.

The script currently supports the following methods:

-   **2SFCA** -- Two-Step Floating Catchment Area
-   **E2SFCA** -- Enhanced Two-Step Floating Catchment Area
-   **M2SFCA** -- Modified Two-Step Floating Catchment Area
-   **3SFCA** -- Three-Step Floating Catchment Area
-   **E3SFCA** -- Enhanced Three-Step Floating Catchment Area

The output of the script is a **CSV table containing accessibility
values (SPAI -- Spatial Accessibility Index)** for each location.

------------------------------------------------------------------------

# Requirements

## Software

The script requires:

-   **ArcGIS Pro** (for the `arcpy` library)
-   **Python 3.x** (included in the ArcGIS Pro environment)

## Python Libraries

The script uses the following Python libraries:

-   `arcpy`
-   `pandas`
-   `numpy`

These libraries are usually already available in the **ArcGIS Pro Python
environment**.

------------------------------------------------------------------------

# Input Data

The script requires **three input datasets**.

## 1. OD Cost Matrix

Output from the **OD Cost Matrix tool** (Network Analyst).

To create an OD Cost Matrix, the following are required:

-   **Network Analyst extension**
-   A working **Network Dataset**

The OD Cost Matrix table must contain these attributes:

  Field           Description
  --------------- --------------------------------------------
  OriginID        ID of the origin location
  DestinationID   ID of the destination location
  Travel time     Travel time between origin and destination

Example:

  OriginID   Name                         Total_time
  ---------- ---------------------------- ------------
  1          Location 1 - Destination 5   12
  1          Location 1 - Destination 8   18
  2          Location 2 - Destination 5   6.5

------------------------------------------------------------------------

## 2. Origins Layer (Population)

A **point layer representing locations containing population**.

Required attributes:

  -----------------------------------------------------------------------
  Variable          Default Field               Description
  ----------------- --------------------------- -------------------------
  `SourceIDField`   `pointid`                   Unique location ID (e.g.,
                                                administrative unit code)

  `PopField`        `grid_code`                 Population count

  `JoinIDField`     `OriginID_Point`            Join ID linking to OD
                                                matrix (must equal
                                                `OBJECTID`, must be
                                                created manually)
  -----------------------------------------------------------------------

------------------------------------------------------------------------

## 3. Destination Layer (Facilities)

Layer representing **service facilities** (e.g., hospitals).

Required attributes:

  Variable               Default Field   Description
  ---------------------- --------------- -----------------------------
  `destinationIdField`   `OBJECTID`      Facility ID (do not modify)
  `supplyField`          `capacity`      Facility capacity

------------------------------------------------------------------------

# Data Path Configuration

This version of the script uses **absolute file paths**, which must be
defined at the beginning of the script.

Example:

``` python
ODmatrix = r"C:\Users\JohnDoe\Documents\ArcGIS\Projects\Project1\Project1.gdb\ODCostMatrixSolver1234xyz\ODLines1234xyz"

origins_layer = r"C:\Users\JohnDoe\Documents\ArcGIS\Projects\Project1\Project1.gdb\PopulationPoints"

destination_layer = r"C:\Users\JohnDoe\Documents\ArcGIS\Projects\Project1\Project1.gdb\Hospitals"

output_csv = r"C:\Users\JohnDoe\Desktop\Results.csv"
```

Description:

  Variable              Meaning
  --------------------- -------------------------
  `ODmatrix`            OD Cost Matrix table
  `origins_layer`       Population layer
  `destination_layer`   Facilities layer
  `output_csv`          Path to output CSV file

------------------------------------------------------------------------

# Calculation Parameters

The script uses the following parameters with default values:

``` python
d0 = 30
decay_f = 15
```

  Parameter   Description
  ----------- -----------------------------------
  `d0`        Maximum travel time (minutes)
  `decay_f`   Gaussian decay function parameter

Example:

``` python
d0 = 30
```

This means the script considers only facilities reachable within **30
minutes travel time**.

------------------------------------------------------------------------

# Calculation Workflow

The script performs the following steps.

## 1. Data Loading

The script loads:

-   OD matrix
-   population data
-   facility capacities

------------------------------------------------------------------------

## 2. Distance Filtering

OD pairs with **travel time greater than `d0`** are removed.

------------------------------------------------------------------------

## 3. Weight Calculation

By default, the script calculates weights using a **Gaussian decay
function**.

A **linear function** is also defined in the code but requires minor
modifications to be used.

The **2SFCA method does not use weights**, therefore it does not apply
any decay function.

------------------------------------------------------------------------

## 4. Supply-to-Demand Ratio (Rj)

For each facility, the **supply-to-demand ratio** is calculated:

    Rj = supply / population

The exact calculation depends on the **specific FCA method** used.

------------------------------------------------------------------------

## 5. Accessibility Index (SPAI)

For each origin location, the **Spatial Accessibility Index (SPAI)** is
calculated as:

    SPAI = sum of accessible facility ratios

------------------------------------------------------------------------

# Script Output

The output is a **CSV file**.

Example:

  ------------------------------------------------------------------------------
  pointid   SPAI_2sfca   SPAI_e2sfca   SPAI_m2sfca   SPAI_3sfca    SPAI_e3sfca
  --------- ------------ ------------- ------------- ------------- -------------
  101       0.34         0.29          0.31          0.27          0.25

  102       0.21         0.19          0.20          0.18          0.17
  ------------------------------------------------------------------------------

Each column corresponds to **one FCA method**.

Note:\
Resulting values are often very small. In practice, it may be useful to
**multiply all columns by a coefficient (e.g., 1,000,000)** to improve
interpretability.

------------------------------------------------------------------------

# Running the Script

The script can be executed directly in **VS Code**.

After completion, the script prints:

    Calculation complete! Results saved to ...

------------------------------------------------------------------------

# Script Structure

The script consists of three main parts.

## 1. Data Loading Functions

-   `load_od_matrix()`
-   `load_population_data()`
-   `load_facility_capacities()`

------------------------------------------------------------------------

## 2. Supply--Demand Ratio Calculations

-   `compute_rj_2sfca()`
-   `compute_rj_e2sfca()`
-   `compute_rj_m2sfca()`
-   `compute_rj_3sfca()`
-   `compute_rj_e3sfca()`

------------------------------------------------------------------------

## 3. Spatial Accessibility Calculation (SPAI)

-   `compute_accessibility_2sfca()`
-   `compute_accessibility_e2sfca()`
-   `compute_accessibility_m2sfca()`
-   `compute_accessibility_3sfca()`
-   `compute_accessibility_e3sfca()`

------------------------------------------------------------------------

# Common Errors

## ID Mismatch

The following fields must match:

    OriginID
    OriginID_Point

------------------------------------------------------------------------

## Incorrect Destination Name

The `Name` field must contain the format:

    Location X - Destination Y

because the script extracts the destination ID from this string.

If the destination layer already contains a field called **`name`**, it
must be **deleted before running the "Import Destinations" step** when
creating the OD Cost Matrix.

------------------------------------------------------------------------

## Missing Capacity Field

The destination layer must contain the field:

    capacity

------------------------------------------------------------------------

# Recommendations and Limitations

1.  Always carefully check that **absolute paths are correct and contain
    no typos**.

2.  The script is **not optimized for extremely large datasets**.

Problems may occur when the **OD Cost Matrix contains more than \~100
million rows**.

Possible solutions:

-   Divide the study area into **smaller regions** (recommended)
-   Compute **only one FCA method per run**
-   Modify the script (e.g., **use lookup instead of join** in certain
    cases)
