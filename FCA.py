# The following script computes various FCA methods: 2SFCA, E2SFCA, M2SFCA, 3SFCA, E3SFCA; and writes the results to a CSV file
# Import libraries
import arcpy
import pandas as pd
import numpy as np

# Paths to input files (Absolute paths)
ODmatrix = r"C:\Users\sholl\Documents\ArcGIS\Packages\E2SFCA_skore_c76f9b\commondata\e2sfca_skore.gdb\ODCostMatrixSolver17gr90w\ODLines1lbcnyk"
origins_layer = r"C:\Users\sholl\Desktop\UP\BC_3\BAK1_BAK2_HOTOVO\Zdrav_zarizeni_body.gdb\Obce_def_singleV2"
destination_layer = r"C:\Users\sholl\Desktop\UP\BC_3\BAK1_BAK2_HOTOVO\Zdrav_zarizeni_body.gdb\NemocniceV3"
output_csv = r"C:\Users\sholl\Desktop\UP\MGR_2\PRODA\fca_results.csv"
# ODmatrix = r"C:\Users\sholl\Documents\ArcGIS\Packages\E2SFCA_skore_c76f9b\commondata\e2sfca_skore.gdb\ODCostMatrixSolverotluf8\ODLines114ri4w"
# destination_layer = r"C:\Users\sholl\Desktop\UP\BC_3\BAK1_BAK2_HOTOVO\Zdrav_zarizeni_body.gdb\Gynekologie"

# Define field names (Global Variables)
OD_SourceIDField = "OriginID"
OD_DestinationNameField = "Name"
OD_TotalTimeField = "Total_Time"
SourceIDField = "kod"
JoinIDField = "OriginID_Obec"
PopField = "Poc_obyv_S"
destinationIdField = "OBJECTID"
supplyField = "PocetLuzek"

# Set Parameters (Global Variables)
d0 = 35     # (in minutes) 
passband_ratio = 0.25
d_passband = d0 * passband_ratio
power_f = 10
decay_f = 15

############################### DEFINE FUNCTIONS ###############################

def gauss_weight(time, threshold, b):

    """
    Computes the Gaussian distance decay weight

    Args:
        time (float): The travel time or (d_ij) between origin i and destination j
        threshold (float): The maximum allowed travel time (d_rel_i)
        b (float): The decay parameter controlling the steepness of the curve

    Returns:
        float: The weight value based on the Gaussian decay function
                Returns 0 if time exceeds the threshold
    """

    return np.exp(- (time / b) ** 2) if time <= threshold else 0

def linear_weight(time, threshold):

    """
    Computes a linear distance decay weight

    The weight decrease is linear as travel time increases, reaching zero when the threshold is exceeded

    Args:
        time (float): The travel time or (d_ij) between origin i and destination j
        threshold (float): The maximum allowed travel time (d0)

    Returns:
        float: The weight value, linearly scaled between 1 (at time = 0) and 0 (at time = threshold)
               Returns 0 if time exceeds the threshold
    """
    
    return max((threshold - time) / threshold, 0)

def no_ddecay_weight(time, threshold):

    """
    This function returns a weight of 1 if the travel time is less than or equal to the specified threshold, 
    and 0 if the travel time exceeds the threshold. 

    Args:
        time (float): The travel time or distance between the origin and destination.
        threshold (float): The maximum allowed time or distance for the weight to be considered as 1.

    Returns:
        int: 1 if the travel time is less than or equal to the threshold, otherwise 0.
    """
    
    return 1 if time <= threshold else 0

def load_od_matrix():

    """
    Loads the Origin-Destination (OD) cost matrix from the specified ArcGIS Pro table

    This function reads the OD matrix using an ArcPy SearchCursor and creates a DataFrame 
    containing travel times between origin and destination pairs.

    Returns:
        pd.DataFrame: A DataFrame with columns:
            - OriginID (origin location ID)
            - Name (destination name, typically contains a unique ID or label)
            - Total_Time (travel time or cost between origin and destination)
            - DestinationID (unique identifier for the destination, from the OD Matrix)
    """
        
    print("Loading OD matrix...")
    df_od = [row for row in arcpy.da.SearchCursor(ODmatrix, [OD_SourceIDField, OD_DestinationNameField, OD_TotalTimeField])]
    df_od = pd.DataFrame(df_od, columns=[OD_SourceIDField, OD_DestinationNameField, OD_TotalTimeField])
    df_od['DestinationID'] = df_od[OD_DestinationNameField].astype(str).str.split().str[-1]
    return df_od

def load_population_data():

    """
    Loads population data from the origins layer.

    This function reads information about origin locations (e.g., municipalities)
    using an ArcPy SearchCursor. It returns a DataFrame that includes each location's 
    unique identifier, population size, and join field for matching with the OD matrix.

    Returns:
        pd.DataFrame: A DataFrame with columns:
            - kod (unique identifier for each origin, e.g., municipality code)
            - Poc_obyv_S (population count at each origin location)
            - OriginID_Obec (join field used to match origins to OD matrix)
    """

    print("Loading population data...")
    df_pop = [row for row in arcpy.da.SearchCursor(origins_layer, [SourceIDField, PopField, JoinIDField])]
    return pd.DataFrame(df_pop, columns=[SourceIDField, PopField, JoinIDField])

def load_facility_capacities():

    """
    Loads facility (destination) capacities from the destination layer

    This function reads data about healthcare facilities or service providers,
    including their unique IDs and supply capacity values (e.g., number of physicians, beds, etc.).
    It converts the OBJECTID to a string and sets it as the index to enable joining with OD matrix data.

    Returns:
        pd.DataFrame: A DataFrame indexed by the facility ID as a string ('id_str') and containing:
            - OBJECTID (unique identifier of the facility)
            - kapacita (supply capacity of the facility)
    """
        
    print("Loading facility capacities...")
    df_dest = [row for row in arcpy.da.SearchCursor(destination_layer, [destinationIdField, supplyField])]
    df_dest = pd.DataFrame(df_dest, columns=[destinationIdField, supplyField])
    df_dest['id_str'] = df_dest[destinationIdField].astype(str)
    return df_dest.set_index('id_str')

def filter_od_matrix(df_od, threshold):

    """
    Filters the OD matrix to remove rows where Total_Time exceeds the threshold (d0).
    
    Args:
        df_od (pd.DataFrame): Origin-Destination DataFrame containing:
            - 'Total_Time': travel time from origin to destination
        threshold (float): maximum allowed travel time (d0 in minutes)
    
    Returns:
        pd.DataFrame: Filtered DataFrame containing only rows with Total_Time <= threshold
    """
    
    print(f"Filtering OD matrix with threshold {threshold} minutes...")
    filtered_df = df_od[df_od[OD_TotalTimeField] <= threshold].copy()
    print(f"Remaining OD pairs after filtering: {len(filtered_df)}")
    return filtered_df

def complete_spai_output(df_spai, df_pop, spai_column):

    """
    Ensures all origins (municipalities) are represented in the SPAI output,
    even if their accessibility score is zero (i.e., no facilities within threshold).
    
    Args:
        df_spai (pd.DataFrame): Calculated SPAI values. Must include 'kod' and [spai_column].
        df_pop (pd.DataFrame): Population DataFrame. Must include 'kod'.
        spai_column (str): The name of the column containing SPAI scores.
        
    Returns:
        pd.DataFrame: Final merged DataFrame with one row per municipality.
    """

    all_origins = df_pop[[SourceIDField]].drop_duplicates()
    df_final = all_origins.merge(df_spai, on=SourceIDField, how='left')
    df_final[spai_column] = df_final[spai_column].fillna(0)
    return df_final

def compute_selection_weights(df_od):

    """
    Computes the selection weight Gij for each origin-destination pair based on a distance decay function.

    Formula:
        Gij = f(d_ij) / sum_j( f(d_ij) ) where j within catchment threshold

    Args:
        df_od (pd.DataFrame): DataFrame containing:
            - 'OriginID': population site i
            - 'DestinationID': facility j
            - 'Total_Time': travel time from i to j

    Returns:
        pd.DataFrame: DataFrame with:
            - 'OriginID'
            - 'DestinationID'
            - 'Gij': selection weight
    """

    print("Computing selection weights Gij...")

    # Compute Tij = f(d_ij) for all OD pairs
    df_od['Tij'] = df_od[OD_TotalTimeField].apply(lambda x: gauss_weight(x, d0, decay_f))
    # Tij is the same as f(d_ij) in the context of the function
    # Normalize f(d_ij) within each OriginID
    df_od['Gij'] = df_od.groupby(OD_SourceIDField)['Tij'].transform(lambda x: x / x.sum()).fillna(0)
    return df_od

def compute_huff_probability(df_od, df_dest):
    # !note: this function should be rewritten using dataframes instead of dictionaries (performance? and code consistency!)

    """
    Computes Huff probabilities (Huff_ij) for each origin-destination pair

    The Huff model estimates the probability that a population at origin i will choose
    a facility j based on:
      - the supply capacity S_j of facility j
      - the travel impedance (distance decay) from i to j using a Gaussian function

    The formula used:
        Huff_ij = (S_j * f(d_ij)) / sum(S_j * f(d_ij))
    Where:
        - f(d_ij): Gaussian weight based on travel time from i to j
        - d0: maximum catchment threshold (in minutes)
        - decay_f: Gaussian decay parameter

    Only facility locations j within d0 minutes of origin i are considered.

    Args:
        df_od (pd.DataFrame): DataFrame with origin-destination pairs containing:
            - OriginID (origin location identifier)
            - DestinationID (destination identifier)
            - Total_Time (travel time between origin and destination)
        df_dest (pd.DataFrame): DataFrame with facility capacities containing:
            - OBJECTID (facility identifier)
            - kapacita (supply capacity of each facility)

    Returns:
        pd.DataFrame: DataFrame with columns:
            - OriginID
            - DestinationID
            - Huff_ij: Probability of population at i choosing facility j
    """

    print("Computing Huff probabilities...")

# Create a dictionary of {(OriginID, DestinationID): travel_time} from the OD matrix
    distances = {
        (row[OD_SourceIDField], row['DestinationID']): row[OD_TotalTimeField]
        for _, row in df_od.iterrows()
    }

    # Convert the facility capacities into a dictionary {DestinationID: capacity}
    supply_capacities = df_dest[supplyField].to_dict()

    # Prepare an empty dictionary to hold the computed Huff probabilities
    huff_scores = {}

    # Get a list of unique OriginIDs from the OD matrix
    origins = set(i for i, j in distances.keys())

    # Loop through each origin i
    for i in origins:
        # Calculate the denominator for Huff_ij: sum of (S_j * f(d_ij)) for all reachable destinations j
        denominator = sum(
            supply_capacities.get(j, 0) * gauss_weight(distances[i, j], d0, decay_f)
            for j in supply_capacities.keys()
            if (i, j) in distances and distances[i, j] <= d0  # only include reachable destinations
        )

        # Loop through destinations j for the same origin i
        for j in supply_capacities.keys():
            if (i, j) in distances and distances[i, j] <= d0:
                # Calculate the numerator: S_j * f(d_ij)
                numerator = supply_capacities[j] * gauss_weight(distances[i, j], d0, decay_f)

                # Final Huff_ij = numerator / denominator, if denominator is not zero
                huff_scores[(i, j)] = numerator / denominator if denominator > 0 else 0

    # Convert the dictionary into a DataFrame with three columns
    huff_df = pd.DataFrame(
        [(i, j, prob) for (i, j), prob in huff_scores.items()],
        columns=[OD_SourceIDField, 'DestinationID', 'Huff_ij']
    )

    return huff_df

def compute_rj_2sfca(df_od, df_pop, df_dest):

    """
    Computes Rj values for the 2-Step Floating Catchment Area (2SFCA) method.

    The formula implemented is:

        Rj = Sj / SUM_over_i{d_ij <= dmax} [ Pi ]

    Where:
        - Rj: accessibility index for facility j
        - Sj: supply or capacity at facility j
        - Pi: population at location i
        - dij: travel time or distance between population i and facility j
        - dmax: catchment threshold (set globally as d0)

    Args:
        df_od (pd.DataFrame): Origin-Destination matrix.
            Required columns: ['OriginID', 'DestinationID', 'Total_Time']
        df_pop (pd.DataFrame): Population data.
            Required columns: ['kod', 'Poc_obyv_S', 'OriginID_Obec']
        df_dest (pd.DataFrame): Facility supply data, indexed by DestinationID.
            Required columns: ['kapacita']

    Returns:
        pd.DataFrame: DataFrame with:
            - 'DestinationID': facility ID
            - 'Rj': supply-to-demand ratio for each facility
    """

    print("Computing Rj for 2SFCA...")

    # Apply distance decay function
    df_od["f_dij"] = df_od[OD_TotalTimeField].apply(lambda x: no_ddecay_weight(x, d0))

    # Merge with population data to get Pi
    df_merged = df_od.merge(df_pop, left_on=OD_SourceIDField, right_on=JoinIDField)

    # Compute denominator: SUM( Pi * f(d_ij) )
    df_merged["Pi_f_dij"] = df_merged[PopField] * df_merged["f_dij"]
    df_denom = df_merged.groupby("DestinationID")["Pi_f_dij"].sum().reset_index(name="Denom")

    # Merge with facility capacities S_j
    df_rj_2sfca = df_denom.merge(df_dest, left_on="DestinationID", right_index=True, how="left")

    # Compute Rj = ( S_j ) / SUM( Pi )
    df_rj_2sfca["Rj"] = df_rj_2sfca[supplyField] / df_rj_2sfca["Denom"]
    return df_rj_2sfca[["DestinationID", "Rj"]]

def compute_rj_e2sfca(df_od, df_pop, df_dest):

    """
    Computes Rj values for the Enhanced 2-Step Floating Catchment Area (E2SFCA) method.

    The formula implemented is:

        Rj = Sj / SUM_over_i{d_ij <= dmax} [ Pi * f(d_ij) ]

    Where:
        - Rj: accessibility index for facility j
        - Sj: supply or capacity at facility j
        - Pi: population at location i
        - dij: travel time or distance between population i and facility j
        - f(dij): distance decay function (e.g., linear decay)
        - dmax: catchment threshold (set globally as d0)

    Args:
        df_od (pd.DataFrame): Origin-Destination matrix.
            Required columns: ['OriginID', 'DestinationID', 'Total_Time']
        df_pop (pd.DataFrame): Population data.
            Required columns: ['kod', 'Poc_obyv_S', 'OriginID_Obec']
        df_dest (pd.DataFrame): Facility supply data, indexed by DestinationID.
            Required columns: ['kapacita']

    Returns:
        pd.DataFrame: DataFrame with:
            - 'DestinationID': facility ID
            - 'Rj': supply-to-demand ratio for each facility
    """

    print("Computing Rj for E2SFCA...")

    # Apply distance decay function
    df_od["f_dij"] = df_od[OD_TotalTimeField].apply(lambda x: linear_weight(x, d0))

    # Merge with population data to get Pi
    df_merged = df_od.merge(df_pop, left_on=OD_SourceIDField, right_on=JoinIDField)

    # Compute denominator: SUM( Pi * f(d_ij) )
    df_merged["Pi_f_dij"] = df_merged[PopField] * df_merged["f_dij"]
    df_denom = df_merged.groupby("DestinationID")["Pi_f_dij"].sum().reset_index(name="Denom")

    # Merge with facility capacities S_j
    df_rj_e2sfca = df_denom.merge(df_dest, left_on="DestinationID", right_index=True, how="left")

    # Compute Rj = ( S_j ) / SUM( Pi * f(d_ij) )
    df_rj_e2sfca["Rj"] = df_rj_e2sfca[supplyField] / df_rj_e2sfca["Denom"]
    return df_rj_e2sfca[["DestinationID", "Rj"]]

def compute_rj_m2sfca(df_od, df_pop, df_dest):

    """
    Computes Rj values for the Modified 2-Step Floating Catchment Area (M2SFCA) method using:

        R_ij = (S_j * f(d_ij)) / SUM_over_i(P_i * f(d_ij))

    Where:
        - S_j: supply (capacity) of facility j
        - P_i: population at location i
        - d_ij: travel time between population location i and facility j
        - f(d_ij): distance decay function (e.g., linear decay)
        - R_ij: intermediate accessibility contribution from facility j to location i

    This function returns the sum of R_ij values per facility j, resulting in the final Rj values.

    Args:
        df_od (pd.DataFrame): Origin-Destination matrix with travel times.
            Required columns: ['OriginID', 'DestinationID', 'Total_Time']
        df_pop (pd.DataFrame): Population data.
            Required columns: [SourceIDField, PopField, JoinIDField]
        df_dest (pd.DataFrame): Facility capacity data, indexed by destination ID (id_str).
            Required columns: [supplyField]

    Returns:
        pd.DataFrame: DataFrame with:
            - 'DestinationID': destination/facility identifier
            - 'Rij': computed supply-to-demand ratio for each facility
    """

    print("Computing Rj for M2SFCA...")

    # Apply distance decay function
    df_od["f_dij"] = df_od[OD_TotalTimeField].apply(lambda x: linear_weight(x, d0))

    # Merge with population data to get Pi
    df_merged = df_od.merge(df_pop, left_on=OD_SourceIDField, right_on=JoinIDField)

    # Compute denominator: SUM( Pi * f(d_ij) )
    df_merged["Pi_f_dij"] = df_merged[PopField] * df_merged["f_dij"]
    df_denom = df_merged.groupby("DestinationID")["Pi_f_dij"].sum().reset_index(name="Denom")

    # Merge with facility capacities S_j
    df_rj_m2sfca = df_denom.merge(df_dest, left_on="DestinationID", right_index=True, how="left")

    # Compute Rj = ( S_j * f(d_ij) ) / SUM( Pi * f(d_ij) )
    df_od = df_od.merge(df_rj_m2sfca[["DestinationID", supplyField, "Denom"]], on="DestinationID", how="left")
    df_od["Numerator"] = df_od[supplyField] * df_od["f_dij"]
    df_rj_m2sfca["Rij"] = (df_od["Numerator"]) / df_rj_m2sfca["Denom"]
    return df_rj_m2sfca[["DestinationID", "Rij"]]

def compute_rj_3sfca(df_od, df_pop, df_dest):
    
    """
    Computes Rj for the 3-Step Floating Catchment Area (3SFCA) method using:

        Rj = Sj / SUM_i ( Gij * Pi * f(dij) )

    Where:
        - Sj is the supply capacity at facility j
        - Gij is the selection weight between i and j
        - Pi is the population at location i
        - f(dij) is a linear distance decay function

    Args:
        df_od (pd.DataFrame): Origin-Destination matrix with:
            - 'OriginID'
            - 'DestinationID'
            - 'Total_Time'
            - 'Gij' (selection weight)
        df_pop (pd.DataFrame): Population data with:
            - 'OriginID_Obec'
            - 'Poc_obyv_S'
        df_dest (pd.DataFrame): Facility capacities, indexed by destination ID ('id_str'), with:
            - 'kapacita'

    Returns:
        pd.DataFrame: DataFrame with:
            - 'DestinationID'
            - 'Rj' (computed accessibility score)
    """

    print("Computing Rj for 3SFCA...")

    # Apply distance decay function f(dij)
    df_od["f_dij"] = df_od[OD_TotalTimeField].apply(lambda x: linear_weight(x, d0))

    # Merge with population (to get Pi)
    df_merged = df_od.merge(df_pop, left_on=OD_SourceIDField, right_on=JoinIDField)

    # Compute weighted demand: Gij * Pi * f(dij)
    df_merged["weighted_demand"] = df_merged["Gij"] * df_merged[PopField] * df_merged["f_dij"]

    # Aggregate weighted demand per destination
    df_denom = df_merged.groupby("DestinationID")["weighted_demand"].sum().reset_index(name="Denom")

    # Merge with facility capacities (Sj)
    df_rj_3sfca = df_denom.merge(df_dest, left_on="DestinationID", right_index=True, how="left")

    # Calculate Rj = Sj / Denom
    df_rj_3sfca["Rj"] = df_rj_3sfca[supplyField] / df_rj_3sfca["Denom"]
    df_rj_3sfca["Rj"] = df_rj_3sfca["Rj"].replace([np.inf, -np.inf], 0).fillna(0)

    return df_rj_3sfca[["DestinationID", "Rj"]]

def compute_rj_e3sfca(df_huff, df_od, df_pop, df_dest):

    """
    Computes Rj values for each destination for the Enhanced 3-Step Floating Catchment Area (E3SFCA) method using:

        Rj = S_j / sum(Huff_ij * P_i * f(d_ij))

    The denominator is calculated over all i (origins) within the catchment area of j (i.e., where d_ij <= d0)

    Args:
        df_huff (pd.DataFrame): DataFrame with Huff probabilities for each (i, j) pair.
            Required columns: ['OriginID', 'DestinationID', 'Huff_ij']
        df_od (pd.DataFrame): OD matrix with travel times.
            Required columns: ['OriginID', 'DestinationID', 'Total_Time']
        df_pop (pd.DataFrame): Population data.
            Required columns: [SourceIDField, PopField, JoinIDField]
        df_dest (pd.DataFrame): Facility capacity data (indexed by destination ID).
            Required index: 'id_str'
            Required column: supplyField

    Returns:
        pd.DataFrame: DataFrame with columns:
            - 'DestinationID': ID of each facility
            - 'Rj': computed availability ratio for each facility
    """

    print("Computing Rj for E3SFCA...")

    # Merge Huff probabilities with OD matrix (for distances)
    df_merged = df_huff.merge(df_od, on=[OD_SourceIDField, 'DestinationID'])

    # Merge with population data
    df_merged = df_merged.merge(df_pop, left_on=OD_SourceIDField, right_on=JoinIDField)

    # Apply distance decay function
    df_merged['f_dij'] = df_merged[OD_TotalTimeField].apply(lambda x: linear_weight(x, d0))

    # Compute denominator: SUM(Huff_ij * P_i * f(d_ij))
    df_merged['denominator'] = df_merged['Huff_ij'] * df_merged[PopField] * df_merged['f_dij']
    denominators = df_merged.groupby('DestinationID')['denominator'].sum()

    # Compute Rj for each facility
    df_rj_e3sfca = df_dest.copy()
    df_rj_e3sfca['Rj'] = df_rj_e3sfca.index.map(lambda j: df_rj_e3sfca.loc[j, supplyField] / denominators[j] if j in denominators and denominators[j] > 0 else 0)

    # Reset index for output
    df_rj_e3sfca = df_rj_e3sfca.reset_index()[['id_str', 'Rj']].rename(columns={'id_str': 'DestinationID'})

    return df_rj_e3sfca[["DestinationID", "Rj"]]

def compute_accessibility_2sfca(df_od, df_rj_2sfca):

    """
    Computes SPAi for 2SFCA:

        SPAi = SUM( Rj )

    Parameters:
        - df_od: DataFrame containing distances (columns: OriginID, DestinationID, Total_Time)
        - df_rj_2sfca: DataFrame containing Rj values (columns: DestinationID, Rj)

    Returns:
        - DataFrame with SPAi values.
    """

    print("Computing SPAI for 2SFCA...")

    # Apply distance decay function f(d_ij)
    df_od['f_dij'] = df_od[OD_TotalTimeField].apply(lambda x: no_ddecay_weight(x, d0))

    # Merge with Rj values
    df_merged = df_od.merge(df_rj_2sfca, on="DestinationID", how="left")

    # Compute SPA_i = SUM( R_j * f(d_ij) ) for each origin
    df_spai_2sfca = df_merged.groupby(OD_SourceIDField).apply(lambda group: (group["Rj"] * group["f_dij"]).sum()).reset_index(name="SPAI_2sfca")

    # Merge with population data to get SourceIDField
    df_spai_2sfca = df_spai_2sfca.merge(df_pop[[JoinIDField, SourceIDField]], left_on=OD_SourceIDField, right_on=JoinIDField, how='left')

    # Reorder columns for clarity
    df_spai_2sfca = df_spai_2sfca[[SourceIDField, 'SPAI_2sfca']]

    return df_spai_2sfca

def compute_accessibility_e2sfca(df_od, df_rj_e2sfca):

    """
    Computes SPAi for E2SFCA:

        SPAi = SUM(Rj * f(d_ij))

    Parameters:
        - df_od: DataFrame containing distances (columns: OriginID, DestinationID, Total_Time)
        - df_rj_e2sfca: DataFrame containing Rj values (columns: DestinationID, Rj)

    Returns:
        - DataFrame with SPAi values.
    """

    print("Computing SPAI for E2SFCA...")

    # Apply distance decay function f(d_ij)
    df_od['f_dij'] = df_od[OD_TotalTimeField].apply(lambda x: linear_weight(x, d0))

    # Merge with Rj values
    df_merged = df_od.merge(df_rj_e2sfca, on="DestinationID", how="left")

    # Compute SPA_i = SUM( R_j * f(d_ij) ) for each origin
    df_spai_e2sfca = df_merged.groupby(OD_SourceIDField).apply(lambda group: (group["Rj"] * group["f_dij"]).sum()).reset_index(name="SPAI_e2sfca")

    # Merge with population data to get SourceIDField
    df_spai_e2sfca = df_spai_e2sfca.merge(df_pop[[JoinIDField, SourceIDField]], left_on=OD_SourceIDField, right_on=JoinIDField, how='left')

    # Reorder columns for clarity
    df_spai_e2sfca = df_spai_e2sfca[[SourceIDField, 'SPAI_e2sfca']]

    return df_spai_e2sfca

def compute_accessibility_m2sfca(df_od, df_rj_m2sfca):

    """
    Computes SPAi for M2SFCA:

        SPAi = SUM(Rij * f(d_ij))

    Parameters:
        - df_od: DataFrame containing distances (columns: OriginID, DestinationID, Total_Time)
        - df_rj_m2sfca: DataFrame containing Rj values (columns: DestinationID, Rj)

    Returns:
        - DataFrame with SPAi values.
    """

    print("Computing SPAI for M2SFCA...")

    # Apply distance decay function f(d_ij)
    df_od['f_dij'] = df_od[OD_TotalTimeField].apply(lambda x: linear_weight(x, d0))

    # Merge with Rj values
    df_merged = df_od.merge(df_rj_m2sfca, on="DestinationID", how="left")

    # Compute SPA_i = SUM( R_j * f(d_ij) ) for each origin
    df_spai_m2sfca = df_merged.groupby(OD_SourceIDField).apply(lambda group: (group["Rij"] * group["f_dij"]).sum()).reset_index(name="SPAI_m2sfca")

    # Merge with population data to get SourceIDField
    df_spai_m2sfca = df_spai_m2sfca.merge(df_pop[[JoinIDField, SourceIDField]], left_on=OD_SourceIDField, right_on=JoinIDField, how='left')

    # Reorder columns for clarity
    df_spai_m2sfca = df_spai_m2sfca[[SourceIDField, 'SPAI_m2sfca']]

    return df_spai_m2sfca

def compute_accessibility_3sfca(df_od, df_rj_3sfca):

    """
    Computes SPAi values for the 3-Step Floating Catchment Area (3SFCA) method.

    Formula:
        SPAi = sum_j ( Gij * Rj * f(dij) )

    Args:
        df_od (pd.DataFrame): Origin-Destination matrix with travel times and Gij weights.
            Required columns: ['OriginID', 'DestinationID', 'Total_Time', 'Gij']
        df_rj (pd.DataFrame): Rj values per facility.
            Required columns: ['DestinationID', 'Rj']

    Returns:
        pd.DataFrame: SPAi accessibility scores per OriginID
            Columns: ['OriginID', 'SPAi']
    """

    print("Computing SPAi for 3SFCA...")

    # Compute f(dij)
    df_od["f_dij"] = df_od[OD_TotalTimeField].apply(lambda x: linear_weight(x, d0))

    # Merge Rj values onto OD table
    df_merged = df_od.merge(df_rj_3sfca, on="DestinationID", how="left")

    # Calculate SPAi = Gij * Rj * f(dij)
    df_merged["SPAi_term"] = df_merged["Gij"] * df_merged["Rj"] * df_merged["f_dij"]

    # Aggregate for each OriginID
    df_spai_3sfca = df_merged.groupby(OD_SourceIDField)["SPAi_term"].sum().reset_index(name="SPAI_3sfca")

    # Merge with population data to get SourceIDField
    df_spai_3sfca = df_spai_3sfca.merge(df_pop[[JoinIDField, SourceIDField]], left_on=OD_SourceIDField, right_on=JoinIDField, how='left')

    # Reorder columns for clarity
    df_spai_3sfca = df_spai_3sfca[[SourceIDField, 'SPAI_3sfca']]

    return df_spai_3sfca

def compute_accessibility_e3sfca(df_huff, df_od, df_rj_e3sfca):
    
    """
    Computes SPAi for each location i based on:

        SPAi = SUM(Huff_ij * R_j * f(d_ij))

    Parameters:
        - df_huff: DataFrame containing Huff probabilities (columns: OriginID, DestinationID, Huff_ij)
        - df_od: DataFrame containing distances (columns: OriginID, DestinationID, Total_Time)
        - df_rj_e3sfca: DataFrame containing Rj values (columns: DestinationID, Rj)

    Returns:
        - DataFrame with SPAi values for each location i.
    """

    print("Computing SPAI for E3SFCA...")

    # Merge Huff probabilities with OD matrix (to get distances)
    df_merged = df_huff.merge(df_od, on=[OD_SourceIDField, 'DestinationID'])

    # Apply distance decay function
    df_merged['f_dij'] = df_merged[OD_TotalTimeField].apply(lambda x: linear_weight(x, d0))

    # Merge with Rj values
    df_merged = df_merged.merge(df_rj_e3sfca, on='DestinationID')

    # Compute SPAi components
    df_merged['SPAi_component'] = df_merged['Huff_ij'] * df_merged['Rj'] * df_merged['f_dij']

    # Sum across all destinations for each origin
    df_spai_e3sfca = df_merged.groupby(OD_SourceIDField)['SPAi_component'].sum().reset_index()
    df_spai_e3sfca.rename(columns={'SPAi_component': 'SPAI_e3sfca'}, inplace=True)

    # Merge with population data to get SourceIDField
    df_spai_e3sfca = df_spai_e3sfca.merge(df_pop[[JoinIDField, SourceIDField]], left_on=OD_SourceIDField, right_on=JoinIDField, how='left')

    # Reorder columns for clarity
    df_spai_e3sfca = df_spai_e3sfca[[SourceIDField, 'SPAI_e3sfca']]

    return df_spai_e3sfca

############################### EXECUTE ###############################

# Load Data
df_od = load_od_matrix()
df_od = filter_od_matrix(df_od, d0)
df_pop = load_population_data()
df_dest = load_facility_capacities()

# Compute Huff Probability and Gij
df_huff = compute_huff_probability(df_od, df_dest)
df_od = compute_selection_weights(df_od)

# Compute Rj scores
df_rj_2sfca = compute_rj_2sfca(df_od, df_pop, df_dest)
df_rj_e2sfca = compute_rj_e2sfca(df_od, df_pop, df_dest)
df_rj_m2sfca = compute_rj_m2sfca(df_od, df_pop, df_dest)
df_rj_3sfca = compute_rj_3sfca(df_od, df_pop, df_dest)
df_rj_e3sfca = compute_rj_e3sfca(df_huff, df_od, df_pop, df_dest)

# Compute accessibility scores
df_spai_2sfca = compute_accessibility_2sfca(df_od, df_rj_2sfca)
df_spai_e2sfca = compute_accessibility_e2sfca(df_od, df_rj_e2sfca)
df_spai_m2sfca = compute_accessibility_m2sfca(df_od, df_rj_m2sfca)
df_spai_3sfca = compute_accessibility_3sfca(df_od, df_rj_3sfca)
df_spai_e3sfca = compute_accessibility_e3sfca(df_huff, df_od, df_rj_e3sfca)

# Complete output (fill with zeros) for municipalities with zero SPAI
df_spai_2sfca = complete_spai_output(df_spai_2sfca, df_pop, spai_column="SPAI_2sfca")
df_spai_e2sfca = complete_spai_output(df_spai_e2sfca, df_pop, spai_column="SPAI_e2sfca")
df_spai_m2sfca = complete_spai_output(df_spai_m2sfca, df_pop, spai_column="SPAI_m2sfca")
df_spai_3sfca = complete_spai_output(df_spai_3sfca, df_pop, spai_column="SPAI_3sfca")
df_spai_e3sfca = complete_spai_output(df_spai_e3sfca, df_pop, spai_column="SPAI_e3sfca")

# Merge results into a single DataFrame
df_spai_merged = pd.merge(df_spai_2sfca, df_spai_e2sfca, on=SourceIDField, how='outer')
df_spai_merged = pd.merge(df_spai_merged, df_spai_m2sfca, on=SourceIDField, how='outer')
df_spai_merged = pd.merge(df_spai_merged, df_spai_3sfca, on=SourceIDField, how='outer')
df_spai_merged = pd.merge(df_spai_merged, df_spai_e3sfca, on=SourceIDField, how='outer')

# Export to CSV
df_spai_merged.to_csv(output_csv, index=False)
print(f"Calculation complete! Results saved to {output_csv}")