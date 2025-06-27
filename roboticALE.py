import pandas as pd
import os
import io
from datetime import datetime
import numpy as np
from minio import Minio
from minio.error import S3Error
import re
from utilities import read_local_file, read_minio_file, list_local_files, list_minio_files, get_well_name
import numpy as np
from amiga.libs.growth import GrowthPlate
from amiga.libs.model import GrowthModel

def extract_from_robotic_ALE(
    path_to_data: str = None,
    minio_config: dict = None, 
    minio_path_to_data: str = None,
    exp_meta: dict = None,
    fname_pattern = r'(?P<experiment>\w+)_(?P<timestamp>\d+)_(?P<uniqueID>\w+)_(?P<series>\w+)_(?P<transfer>\d+)_(?P<timepoint>\d+).txt'
):
    """
    Extracts data from robotic ALE experiments.
    
    Args:
        path_to_data: Path to local directory containing plate reader files
        minio_config: Dictionary containing MinIO configuration:
            {
                'endpoint': 'minio-server:port',
                'access_key': 'your-access-key',
                'secret_key': 'your-secret-key',
                'bucket': 'your-bucket-name'
            }
        minio_path_to_data: Path prefix in MinIO bucket containing plate reader files
        
    Returns:
        pd.DataFrame: Raw data containing optical density measurements
    """

    # Initialize data df
    data = pd.DataFrame()
    
    # Define file pattern for plate reader files
    fname_pattern = re.compile(fname_pattern)

    # Get list of files based on source
    if path_to_data is not None:
        files = list_local_files(path_to_data, pattern=fname_pattern)
    elif minio_config is not None and minio_path_to_data is not None:
        files = list_minio_files(minio_config, prefix=minio_path_to_data, pattern=fname_pattern)
    else:
        raise ValueError("Either path_to_data or both minio_config and minio_path_to_data must be provided")

    if len(files) == 0:
        print('No files matching specified file name pattern.')
        return None
    
    # Read info from plate reader file names and file content into df
    for f in files:
        # Initialize row in dataframe
        data_row = {}
        match = fname_pattern.match(os.path.basename(f))
        
        # Parse info contained in plate reader file name
        data_row['experiment'] = exp_meta['experiment_id'] if exp_meta else str(match.group('experiment')) 
        data_row['file_ID'] = str(match.group('uniqueID'))
        data_row['timestamp'] = int(match.group('timestamp'))
        data_row['series'] = str(match.group('series'))
        data_row['plate_index'] = int(match.group('transfer'))
        # t_transfer indicates which plate cols were most recently innoculated.
        # data_row['t_transfer'] = int(match.group('transfer'))
        data_row['transfer'] = int(match.group('transfer'))
   
        # Read plate reader file
        if path_to_data is not None:
            plate_data = read_local_file(f)
        else:
            plate_data = read_minio_file(minio_config, f)

        # Process plate data
        for row in range(8):
            for col in range(12):
                data_row['row'] = row
                data_row['column'] = col
                data_row['od'] = plate_data.iloc[row, col]
                data = pd.concat([data, pd.Series(data_row).to_frame().T])

    data.reset_index(inplace=True, drop=True)

    # Translate row and column numbers to well names
    data['well'] = data.apply(
        lambda x: get_well_name(x['row'], x['column']), axis=1
    )
    
    # Translate timestamp into isoformat time
    data['datetime'] = data['timestamp'].apply(
        lambda x: datetime.fromtimestamp(x).isoformat()
    )
    
    # Appending metadata
    data['filename'] = path_to_data if path_to_data else minio_path_to_data
    data['measurement_type'] = 'OD'
    data['culture_container'] = 'plate_well'
    data['plate_type'] = exp_meta['plate_type'] if exp_meta else'96_shallow'  # Could be parameterized
    data['start_date'] = exp_meta['start_date'] if exp_meta else pd.NA
    # data['exp_index'] = self.exp_index
    # data['operation_id'] = self.op_id

    return data


def map_metadata(df, path_to_metadata: str = None, minio_config: str = None, minio_path_to_metadata: str = None):
    """
    Maps metadata in layout files to a df containing OD measurements from a robotic ALE experiment.
    
    Args:
        path_to_data: Path to local directory containing plate reader files
        minio_config: Dictionary containing MinIO configuration:
            {
                'endpoint': 'minio-server:port',
                'access_key': 'your-access-key',
                'secret_key': 'your-secret-key',
                'bucket': 'your-bucket-name'
            }
        minio_path_to_data: Path prefix in MinIO bucket containing plate layout files
        
    Returns:
        pd.DataFrame: Optical density measurement data and associated metadata.
    """

    # Load layout files
    
    layout_path = path_to_metadata if path_to_metadata else minio_path_to_metadata
    if path_to_metadata is not None:
        strain_layout = read_local_file(os.path.join(layout_path, 'strain_layout.csv'))
        rep_layout = read_local_file(os.path.join(layout_path, 'replicate_layout.csv'))
        gc_layout = read_local_file(os.path.join(layout_path, 'growth_condition_layout.csv'))
        transfer_layout = read_local_file(os.path.join(layout_path, 'transfer_layout.csv'))
    else:
        strain_layout = read_minio_file(minio_config, f"{layout_path}/strain_layout.csv")
        rep_layout = read_minio_file(minio_config, f"{layout_path}/replicate_layout.csv")
        gc_layout = read_minio_file(minio_config, f"{layout_path}/growth_condition_layout.csv")
        transfer_layout = read_minio_file(minio_config, f"{layout_path}/transfer_layout.csv")
        
    # Apply layouts to df
    df['layout_filename'] = layout_path
    df['strain'] = df.apply(
        lambda x: strain_layout.iloc[x['row'], x['column']], axis=1
    ).astype("Int64")
    df['replicate'] = df.apply(
        lambda x: rep_layout.iloc[x['row'], x['column']], axis=1
    ).astype("Int64")
    df['gc'] = df.apply(
        lambda x: gc_layout.iloc[x['row'], x['column']], axis=1
    ).astype("Int64")
    df['l_transfer'] = df.apply(
        lambda x: transfer_layout.iloc[x['row'], x['column']], axis=1
    ).astype("Int64")

    return df


def get_sample_names(df):
    sample_names = 'E:' + df['experiment'] + \
                  '.P:' + df['series'] + '-' + df['plate_index'].astype(str) + \
                  '.W:' + df['well'].astype(str) + \
                  '.S:' + df['strain'].astype(str) + \
                  '.C:' + df['gc'].astype(str) + \
                  '.R:' + df['replicate'].astype(str) + \
                  '.T:' + df['transfer'].astype(str)

    # 'Sample' can only exist where there is media
    df['sample_name'] = np.where(~pd.isna(df['gc']), sample_names, pd.NA)

    return df


def get_plate_names(df):
    df['plate_name'] = (
        'E:' + df['experiment'] + '.P:' + df['series'] + '-' + df['plate_index'].astype(str)
    )
    return df
    

def get_parent_samples(df, parent_frame = None):

    if 'parent_sample' in df.columns:
        print('Parent samples were already assigned.')
        return df

    if not 'sample_name' in df.columns:
        print('Must call get_sample_names() first.')
        return df
        
    parent_sample_names = 'E:' + df['experiment'] + \
                          '.P:' + df['series'] + '-' + (df['plate_index']-1).astype(str) + \
                          '.W:' + df['well'].astype(str) + \
                          '.S:' + df['strain'].astype(str) + \
                          '.C:' + df['gc'].astype(str) + \
                          '.R:' + df['replicate'].astype(str) + \
                          '.T:' + (df['transfer'] - 1).astype(str)

    # Parent sample can only exist where there is an inoculated sample, ie media and strain
    df['parent_sample'] = np.where(~pd.isna(df['sample_name']) & ~pd.isna(df['strain']) & (df['plate_index'] != 1), parent_sample_names, pd.NA)

    # If a dataframe with parent sample names for the first transfer was provided, assign the provided names.
    if not parent_frame is None:
        parent_sample_names_y = df.merge(parent_frame, on=['plate_index', 'well'], how='left')["parent_sample_y"]
        df['parent_sample'] = np.where(~pd.isna(df['sample_name']) & ~pd.isna(df['strain']) & (df['plate_index'] == 1), parent_sample_names_y, df['parent_sample'])
        
    return df


def compute_background(df): #### NOT SURE IF THIS IS CORRECT!!!
    
    # Calculate a background value for each plate reader measurement
    # (based on the wells that only contain media)
    
    df['background'] = pd.NA
    df['background'] = df.groupby(
        ['experiment', 'series', 'plate_index', 'timestamp']
        )['od'].transform(
        lambda x: x[df.loc[x.index, 'strain'].isna()].mean()
        )
    
    return df


def compute_inoculation(df):

    # Compute innoculation timestamp based on the oldest timestamp
    # Compute timepoint in hours for each timestamp
    df['innoculation_timestamp'] = pd.NA
    df.loc[
        ~(pd.isna(df['transfer'])),
        'innoculation_timestamp'
        ] = df.groupby(
            ['series','plate_index', 'transfer']
            )['datetime'].transform("min")

    def calc_timepoint(time_0, time):
        
        if not pd.isna(time_0):
            timepoint = (
                datetime.fromisoformat(time) - datetime.fromisoformat(time_0)
            ).total_seconds()/3600
    
        else:
            timepoint = pd.NA
    
        return timepoint
    
    df['timepoint'] =  df.apply(
        lambda x: calc_timepoint(x['innoculation_timestamp'], x['datetime']),
        axis=1
    )

    return df


def create_plates(df):

    # req_cols are the required fields for the db
    req_cols = ['plate_name','experiment', 'plate_type', 'plate_index', 'layout_filename']
    
    if not set(req_cols) <= set(df.columns):
        print('df is missing required col(s):', set(req_cols) - set(df.columns))
        return None
    
    plates = df[req_cols].drop_duplicates().reset_index(drop=True)
    plates.rename(
        columns={'plate_name': 'id', 'experiment': 'experiment_id'}, inplace=True
        )
    
    return plates


def create_samples(df):

    # req_cols are the required fields for the db
    req_cols = ['sample_name', 'experiment', 'plate_name', 'well', 'transfer',
              'gc', 'strain', 'innoculation_timestamp', 'replicate', 'parent_sample', 'culture_container']
    
    if not set(req_cols) <= set(df.columns):
        print('df is missing required col(s):', set(req_cols) - set(df.columns))
        return None
    
    samples = df[req_cols].loc[
        ~pd.isna(df['sample_name'])
        ].drop_duplicates(
            (req_cols)
                     ).reset_index(drop=True)
    samples.rename(
        columns={
            'sample_name': 'name', 'experiment': 'experiment_id', 
            'plate_name': 'plate', 'transfer': 'passage',
            'gc': 'growth_condition_id', 'strain': 'strain_id',
            'parent_sample': 'parent_sample_name'
            }, inplace=True)
    
    return samples


def create_inoc_procedures(df, inoc_protocol, lab_id=None, contact_id=None): 
# it would be nice if could give experiment procedure instead of the inoc protocol, lab, contact
    
    # req_cols are the required fields for creating operation id and
    # populating other fields in the db
    req_cols = ['sample_name', 'innoculation_timestamp']
    
    if not set(req_cols) <= set(df.columns):
        print('df is missing required col(s):', set(req_cols) - set(df.columns))
        return None
    
    inoc_procedures = df[req_cols].loc[
        ~pd.isna(df['sample_name'])
        ].drop_duplicates(
            (req_cols)
                     ).reset_index(drop=True)
    
    inoc_procedures['id'] = inoc_procedures['sample_name'] + '_inoculation'
    inoc_procedures['protocol_id'] = inoc_protocol
    inoc_procedures['lab_id'] = lab_id
    inoc_procedures['contact_id'] = contact_id

    inoc_procedures.rename(
        columns={'innoculation_timestamp': 'timestamp'},
        inplace=True
    )
    inoc_procedures.drop('sample_name', axis=1, inplace=True)
    
    return inoc_procedures

def create_od_procedures(df, od_protocol, lab_id=None, contact_id=None):
    # req_cols are the required fields for creating operation id and
    # populating other fields in the db
    req_cols = ['sample_name', 'innoculation_timestamp']
    
    if not set(req_cols) <= set(df.columns):
        print('df is missing required col(s):', set(req_cols) - set(df.columns))
        return None
    
    od_procedures = df[req_cols].loc[
        ~pd.isna(df['sample_name'])
        ].drop_duplicates(
            (req_cols)
                     ).reset_index(drop=True)
    
    od_procedures['id'] = od_procedures['sample_name'] + '_od'
    od_procedures['parent_operation'] = od_procedures['sample_name'] + '_inoculation'
    od_procedures['protocol_id'] = od_protocol
    od_procedures['lab_id'] = lab_id
    od_procedures['contact_id'] = contact_id

    od_procedures.rename(
        columns={'innoculation_timestamp': 'timestamp'},
        inplace=True
    )
    od_procedures.drop('sample_name', axis=1, inplace=True)

    return od_procedures

def create_measurements(df):
    req_cols = ['sample_name', 'filename', 'measurement_type']

    if not set(req_cols) <= set(df.columns):
        print('df is missing required col(s):', set(req_cols) - set(df.columns))
        return None

    measurements = df[req_cols].loc[
        ~pd.isna(df['sample_name'])
        ].drop_duplicates(
            (req_cols)
                     ).reset_index(drop=True)
    
    measurements['operation_id'] = measurements['sample_name'] + '_od'
    measurements.rename(
        columns={'sample_name': 'sample_id', 'measurement_type': 'type'},
        inplace=True
    )

    return measurements
        

def create_od_measurements(df):
    req_cols = ['sample_name', 'datetime', 'timepoint', 'od', 'background']

    if not set(req_cols) <= set(df.columns):
        print('df is missing required col(s):', set(req_cols) - set(df.columns))
        return None

    od_measurements = df.loc[~pd.isna(df['sample_name'])][req_cols].reset_index(drop=True)
    od_measurements['operation_id'] = od_measurements['sample_name'] + '_od'
    od_measurements.drop('sample_name', axis=1, inplace=True)
    
    return od_measurements


def get_amiga_metrics(df, subtract_background = False):
    # Given timepoint-od data for a set of samples (operation_ids)
    # calculate growth rate, doubling time, lag time

    # Check if df has required columns
    req_cols = ['operation_id', 'timepoint', 'od']

    if not set(req_cols) <= set(df.columns):
        print('df is missing required col(s):', set(req_cols) - set(df.columns))
        return None

    value = 'od'

    # Calculate background-subtracted values if requested
    if subtract_background:

        if not 'background' in df.columns:
            print('df is missing column "background".')
            return None
            
        df['od_background_subtracted'] = df['od'] - df['background']
        value = 'od_background_subtracted'

    # Reformat df for df GrowthPlate object
    pivot_data = pd.pivot(df, 
                  index = 'timepoint',
                  values = value,
                  columns = 'operation_id'
                 ).reset_index(
                 ).rename(columns=({'timepoint':'time'})
                 ).sort_values('time'
                 ).astype(float)
    
    mapping = pd.DataFrame(data = {
        'm_id':df['operation_id'].unique(),
        'plate_ID':['1'] * df['operation_id'].nunique()
    }
                          ).set_index('m_id')
    
    plate = GrowthPlate(data=pivot_data,key=mapping)
    plate.computeBasicSummary()
    plate.raiseData()
    plate.logData()
    plate.subtractBaseline(True,poly=False)

    # Create df growth curves for each sample, collect into df
    metrics = plate.key[['OD_Max']]
    metrics.insert(0, 'growth_rate', [pd.NA]*len(metrics))
    metrics.insert(0, 'doubling_time', [pd.NA]*len(metrics))
    metrics.insert(0, 'lag_time', [pd.NA]*len(metrics))
    
    for sample in plate.data.columns:
        # fit curve
        thisSample = pd.concat([plate.time, plate.data[sample]], axis=1).dropna()
        thisSample.columns = ['Time', 'OD']
        gm = GrowthModel(thisSample,ARD=True)
        gm.fit()
        curve = gm.run()
        metrics.loc[sample, 'growth_rate'] = curve.gr
        metrics.loc[sample, 'doubling_time'] = curve.td * 60 # convert to minutes
        metrics.loc[sample, 'lag_time'] = curve.lagC * 60 # convert to minutes

    metrics.reset_index(inplace=True)
            
    return metrics
    

def create_growth_measurements(amiga):
    # Massage data so it will fit into database
    # i.e., replace infinity values, replace very large values with NA,
    # and round to 3 dec places
    
    df = amiga.rename(columns={
        'OD_Max': 'max_od',
        'm_id': 'operation_id'
    })
    df['lag_time'] = df['lag_time'].replace([np.inf, -np.inf], 10000)
    df['doubling_time'] = df['doubling_time'].replace([np.inf, -np.inf], 10000)
    df['max_od'] = df['max_od'].apply(lambda x: round(x, 3))
    df['growth_rate'] = df['growth_rate'].apply(lambda x: round(x, 3))
    
    crazy_lag_time = (df['lag_time'] >= 10000) | (df['lag_time'] < 0)
    df['lag_time'] = np.where(
        crazy_lag_time,
        [pd.NA]*len(df),
        df['lag_time'].apply(round)
    )
    
    crazy_doubling_time = (df['doubling_time'] >= 10000) | (df['doubling_time'] < 0)
    df['doubling_time'] = np.where(
        crazy_doubling_time,
        [pd.NA]*len(df),
        df['doubling_time'].apply(round)
    )
    
    od_max_near_0 = df['max_od'] < 0.0001
    df['max_od'] = np.where(
        od_max_near_0,
        [pd.NA]*len(df),
        df['max_od']
    )
    growth_rate_near_0 = df['growth_rate'] < 0.0001
    df['growth_rate'] = np.where(
        growth_rate_near_0,
        [pd.NA]*len(df),
        df['growth_rate']
    )

    return df














                                         
        