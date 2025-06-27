from sqlalchemy import create_engine
import pandas as pd
import numpy as np
import minio
import io
import re
from datetime import datetime
from datetime import timedelta
from amiga.libs.growth import GrowthPlate
from amiga.libs.model import GrowthModel
from utilities import read_mio_csv

def etl(engine, mio):
    '''
    Extracts and transforms plate reader data from MinIO storage,
    and stores it in MySQL database.

    Args:

    Returns:
    
    '''
    
    
    # Functions
    
    
    def get_plate_reader_filenames(
        minio_bucket, minio_path_to_plate_reader_files, regexpattern
    ):
        # Return list of file names that match the expected file name pattern.
        objects = mio.list_objects(
            minio_bucket, prefix=minio_path_to_plate_reader_files
        )
        filepaths = [obj.object_name for obj in objects]
        filenames = [fn.split('/')[1] for fn in filepaths]
        plate_reader_files = [fn for fn in filenames if re.match(regexpattern, fn)]
        
        return plate_reader_files
    
    
    def get_transfers(batch_dict, batch, plate, transfer):
        # Which transfer since beginning of experiment?
        return batch_dict[batch] + (plate-1)*3 + (transfer)
    
    
    def get_plates(batch_dict, batch, plate):
        # How many plates in total?
        return batch_dict[batch] + plate
    
    
    def get_parent_plate(plate, column):
        # Which plate was the parent sample on?
        if column in (1, 4, 7):
            parent_plate = plate-1
        else:
            parent_plate = plate
        return parent_plate
    
    
    def get_parent_col(col):
        # Which column was the parent sample in?
        cols = (1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 10, 11)
        parent_cols = (3, 1, 2, 6, 4, 5, 9, 7, 8, 0, 0, 0)
        parent_col_dict = dict(zip(cols, parent_cols))
        return parent_col_dict[col]
    
    
    def get_well_name(row, col):
        # Return well name based on plate rows/cols.
        well_name = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'][row] + str(col+1)
        return well_name
    
    
    # Hardcoded variables. 
    # Later, these things should be supplied when the experimental team 
    # registers the experiment through the web UI.
    minio_bucket_name = 'synbio'
    path_to_plate_reader_files = 'ALE1b_OD_data/'
    path_to_layout_files = 'plate_layouts/'
    fname_pattern = re.compile(
        r'(?P<experiment>\w+)_(?P<timestamp>\d+)_(?P<uniqueID>\w+)_(?P<plate>\d+)_(?P<transfer>[1-3])_(?P<timepoint>\d+).txt'
        )
    experiment_id = path_to_plate_reader_files.split('_')[0]
    plate_type = '96_shallow' # could also be deep well plate
    start_date = '2025-04-04'
    exp_index = 1
    exp_type = 'autoALE'
    description = ''
    protocol_id = 'mock_1b_protocol'
    lab_id = 1
    contact_id = 1
    operation_id = f"{experiment_id}_operation"
    measurement_type = 'growth'
    plate_reader_filenames = get_plate_reader_filenames(
        minio_bucket_name, path_to_plate_reader_files, fname_pattern
        )
    transfer_layout = read_mio_csv(
        mio, minio_bucket_name, path_to_layout_files + 'transfer_layout.csv'
        )
    rep_layout = read_mio_csv(
        mio, minio_bucket_name, path_to_layout_files + 'replicate_layout.csv'
        )
    strain_layout = read_mio_csv(
        mio, minio_bucket_name, path_to_layout_files + 'strain_layout.csv'
        )
    gc_layout = read_mio_csv(
        mio, minio_bucket_name, path_to_layout_files + 'growth_condition_layout.csv'
        )
    
    
    # Upload this experiment and its operation (i.e., procedure) to the db.
    operation_dict = {
        'id': [operation_id],
        'protocol_id': [protocol_id],
        'lab_id': [lab_id],
        'contact_id': [contact_id],
        'timestamp': [start_date]
    }
    exp_operation_df = pd.DataFrame.from_dict(operation_dict)
    exp_operation_df.to_sql('operation', engine, index=False, if_exists='append')
    
    exp_dict = {
        'id': [experiment_id],
        'type': [exp_type],
        'start_date': [start_date],
        'index': [exp_index],
        'description': [description],
        'operation_id': [operation_id]
    }
    new_exp_df = pd.DataFrame.from_dict(exp_dict)
    new_exp_df.to_sql('experiment', engine, index=False, if_exists='append')
    
    # Initialize data df
    data = pd.DataFrame()
    
    # Read info from plate reader file names and file content into df
    for f in plate_reader_filenames:
    
        try:
    
            # Initialize row in dataframe
            data_row = {}
            match = fname_pattern.match(f)
    
            # Parse info contained in plate reader file name
            data_row['experiment'] = str(match.group('experiment'))
            data_row['file_ID'] = str(match.group('uniqueID'))
            data_row['timestamp'] = int(match.group('timestamp'))
            data_row['plate_index'] = int(match.group('plate'))
            # t_transfer indicates which plate cols were most recently innoculated.
            data_row['t_transfer'] = int(match.group('transfer'))
    
            # Read plate reader files into dataframe
            response = mio.get_object(
                minio_bucket_name, path_to_plate_reader_files+f
                )
            csv_data = response.data
            plate_data = pd.read_csv(io.BytesIO(csv_data), header=None)
            for row in range(8):
                for col in range(12):
                    data_row['row'] = row
                    data_row['column'] = col
                    data_row['OD'] = plate_data.iloc[row, col]
                    data = pd.concat([data, pd.Series(data_row).to_frame().T])
                    
        except Exception as e:
            print(f"Error: {e}")
            
        finally: 
            response.close()
            response.release_conn()
    
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
    data['filename'] = path_to_plate_reader_files
    data['measurement_type'] = measurement_type
    data['experiment'] = experiment_id
    data['plate_type'] = plate_type
    data['start_date'] = start_date
    data['exp_index'] = exp_index
    data['operation_id'] = operation_id
    data['layout_filename'] = path_to_layout_files 
    # Given location on plate (row, col) and layout files, get the strain, 
    # growth condition, replicate number, and transfer_l for each well.
    # (transfer_l indicates transfer based on plate location.)
    data['strain'] = data.apply(
        lambda x: strain_layout.iloc[x['row'], x['column']], axis=1
        ).astype("Int64")
    data['replicate'] = data.apply(
        lambda x: rep_layout.iloc[x['row'], x['column']], axis=1
        ).astype("Int64")
    data['gc'] = data.apply(
        lambda x: gc_layout.iloc[x['row'], x['column']], axis=1
        ).astype("Int64")
    data['l_transfer'] = data.apply(
        lambda x: transfer_layout.iloc[x['row'], x['column']], axis=1
        ).astype("Int64")
    
    # Determine innoculation status
    # Explanation: All wells on a plate are read at each timepoint, but only some
    # of the wells will be inoculated at any given timepoint. Only wells that have
    # media in them (i.e., growth condition not NA) are considered *samples*.
    # Samples that do not have a strain designation are negative controls.
    # Depending on the location on the plate, samples can have 11, 22, or 33
    # timepoints. The calculation below determines which transfer to assign a
    # given OD reading to.
    
    # Calculate actual transfer number for a given reading
    data['l-t_transfer'] = data.apply(
        lambda x: x['l_transfer']-x['t_transfer'], axis=1
        )
    data['transfer'] = np.where(
        (data['l-t_transfer']<=0), data['l_transfer'], 1
        )
    
    # Compute total number of plates, transfers for this experiment
    # There are 3 transfers per plate, x plates per batch, and y batches per
    # experiment. A batch is defined as a continuous run of the robot without
    # any interruption. All measurements in a batch will have the same file_ID.
    # Each batch starts with plate=1, transfer=1, timepoint=1. To accurately
    # calculate the cummulative number of transfers since the beginning of the
    # experiment for each data point, we need to calculate the number of transfers
    # for the individual batches and their cumulative number (in their correct 
    # order).
    
    batches = (
        data.groupby('file_ID')['datetime']
        .agg("min")
        .reset_index()
        .sort_values('datetime')['file_ID']
        .to_list()
    )
    batch_n_transfers = [0]  # transfers for each of the batches.
    # Start with 0, because the first batch has no prior transfers.
    batch_n_plates = [0]  # number of plates for each of the batches.
    # Start with 0, because first batch has no prior plates.
    
    for batch in batches:
        # Loop through each batch and calculate the number of transfers
        batch_data = data.loc[
            (data['file_ID'] == batch) & (~pd.isna(data['transfer']))
        ]
        # From last measurement in batch, total number of plates and transfers in
        # this batch
        max_plate = batch_data.sort_values(['plate_index', 'transfer']).iloc[-1]
        batch_n_plates.append(max_plate['plate_index'])
        batch_transfers = (max_plate['plate_index']-1)*3 + (max_plate['transfer'])
        batch_n_transfers.append(batch_transfers)
    
    # cumulative transfers (# of transfers since the start of the experiment)
    batch_cumsum_t = np.cumsum(np.array(batch_n_transfers[:-1])) 
    batch_dict_t = dict(zip(batches, batch_cumsum_t))
    # cumulative plates (# of plates since the start of the experiment)
    batch_cumsum_p = np.cumsum(np.array(batch_n_plates[:-1])) 
    batch_dict_p = dict(zip(batches, batch_cumsum_p))
    
    data['cum_transfer'] = np.where(
        data['transfer']==0, np.array([0]*len(data)), 
        data.apply(
            lambda x: get_transfers(
                batch_dict_t, x['file_ID'], x['plate_index'], x['transfer']
                ), axis=1)
                ) 
    data['cum_plate'] = data.apply(
        lambda x: get_plates(
            batch_dict_p, x['file_ID'], x['plate_index']
            ), axis=1
            )
    
    # data['strain'] is NA (i.e., negative control) if not yet inoculated
    data['strain'] = np.where(
        (data['l-t_transfer']>0), 
        np.array([pd.NA]*len(data)), 
        data['strain']
        )
    data['strain'] = data['strain'].astype("Int64")
    data['replicate'] = np.where(
        pd.isna(data['strain']), 
        np.array([pd.NA]*len(data)), 
        data['replicate']
        )
    data['replicate'] = data['replicate'].astype("Int64")
    
    # Calculate a background value for each plate reader measurement
    # (based on the wells that only contain media)
    data['background'] = pd.NA
    data['background'] = data.groupby(
        ['experiment', 'plate_index', 'timestamp']
        )['OD'].transform(
        lambda x: x[data.loc[x.index, 'strain'].isna()].mean()
        )
    
    # Compute innoculation timestamp based on the oldest timestamp
    data['innoculation_timestamp'] = pd.NA
    data.loc[
        (~(pd.isna(data['cum_transfer']))) & (~(pd.isna(data['strain']))),
        'innoculation_timestamp'
        ] = data.groupby(
            ['cum_plate', 'cum_transfer']
            )['datetime'].transform("min")
    
    # Assign parent samples
    # Only innoculated samples (i.e., not neg. controls) can have parent samples.
    # Only samples after passage 1 have parent samples (passage 1 parents will 
    # have to be manually assigned)
    data['parent_plate'] = pd.NA
    data.loc[
        ~pd.isna(data['strain']) & (data['cum_transfer'] != 1),
        'parent_plate'] = data.apply(
            lambda x: get_parent_plate(
                x['cum_plate'], x['column']
                ), axis=1)
    data['parent_well'] = pd.NA
    data.loc[
        ~pd.isna(data['strain']) & (data['cum_transfer'] != 1),
        'parent_well'] = data.apply(
            lambda x: get_well_name(x['row'], get_parent_col(x['column'])), axis=1)
    
    # Assign plate, well, and sample names
    data = data.convert_dtypes()
    data['plate_name'] = (
        'E:' + data['experiment'] + '.P:' + data['cum_plate'].astype(str)
    )
    data['well_name'] = data['plate_name'] + '.W:' + data['well'].astype(str)
    data['sample_name']= pd.NA
    data.loc[(~(pd.isna(data['gc']))), 'sample_name'] = (
        data['well_name'] + '.S:' + data['strain'].astype(str) + '.C:' +
        data['gc'].astype(str) + '.R:' + data['replicate'].astype(str) +
        '.T:' + data['cum_transfer'].astype(str)
    )
    
    # get parent_id
    data['parent_id'] = pd.NA
    data.loc[
        (~(pd.isna(data['strain']))) &
        (~(pd.isna(data['gc']))) &
        (data['cum_transfer'] != 1),
        'parent_id'] = (
            'E:' + data['experiment'] + '.P:' + data['parent_plate'].astype(str) + 
            '.W:' + data['parent_well'].astype(str) + '.S:' + data['strain'].astype(str) + 
            '.C:' + data['gc'].astype(str) + '.R:' + data['replicate'].astype(str) + 
            '.T:' + (data['cum_transfer']-1).astype(str)
        )
    # data = data.convert_dtypes()
    
    
    # Reformat data for database upload
    # 1) Plates
    plates = data[
        ['plate_name', 'experiment', 'plate_type', 'plate_index', 'layout_filename']
        ].drop_duplicates().reset_index(drop=True)
    plates.rename(
        columns={'plate_name': 'id', 'experiment': 'experiment_id'}, inplace=True
        )
    plates.to_sql('plate', engine, index=False, if_exists='append')
    
    # 2) Samples and associated measurements
    # Each sample has a measurement of type 'growth' to which all od_measurements map.
    sample_meas = data.loc[
        ~pd.isna(data['sample_name']) & (~(pd.isna(data['gc'])))
        ].drop_duplicates(
            (['sample_name', 'experiment', 'plate_name', 'well', 'cum_transfer',
              'gc', 'strain', 'innoculation_timestamp', 'replicate', 'parent_id'])
                     ).sort_values('innoculation_timestamp')
    samples = sample_meas[
        (['sample_name', 'experiment', 'plate_name', 'well', 'cum_transfer','gc', 
          'strain', 'innoculation_timestamp', 'replicate', 'parent_id'])
          ].copy()
    samples.rename(
        columns={
            'sample_name': 'name', 'experiment': 'experiment_id', 
            'plate_name': 'plate', 'cum_transfer': 'passage',
            'gc': 'growth_condition_id', 'strain': 'strain_id',
            'parent_id': 'parent_sample_name'
            }, inplace=True)
    samples.to_sql('sample', engine, index=False, if_exists='append')
    measurements = sample_meas[
        ['sample_name', 'operation_id', 'measurement_type', 'filename']
        ].copy()
    measurements.rename(
        columns={'sample_name': 'sample_id', 'measurement_type': 'type'},
        inplace=True
        )
    measurements.to_sql('measurement', engine, index=False, if_exists='append')
    
    # 3) OD_measurements
    # The index of the db table `measurements` autoincrements, so the measurement
    # ids to which the od_measurements of the current dataset belong are not known
    # ahead of time. Therefore, need to download all measurements whith sample ids
    # from this dataset to get their ids, then merge with OD readings.
    sample_names = tuple(sample_meas['sample_name'])
    meas_from_db = pd.read_sql(
        f"SELECT `id`, `sample_id` FROM `measurement` WHERE `sample_id` IN {sample_names};",
        engine
        ).rename(columns={'id': 'measurement_id'})
    od_meas = meas_from_db.merge(
        data, left_on='sample_id', right_on='sample_name', how='inner'
        )[['plate_name', 'measurement_id', 'datetime', 'OD', 'background']].rename(
            columns={'OD': 'od'}
            )
    od_meas.drop('plate_name', axis=1).to_sql('od_measurement', engine, index=False, if_exists='append')

    
    # 4) Generate growth curves with AMiGA code and store in db   
    # Calculate background value, calculate timepoints in hours.
    
    # First drop duplicate measurement.
    # (There are always two timepoints that are the same in the ALE1a dataset
    # because measurement taken before and after transfer received same
    # timepoint.)
    od_meas = od_meas.drop_duplicates(['measurement_id', 'datetime']).copy()
    od_meas['od_sub_background'] = od_meas['od'] - od_meas['background']
    
    # For each sample (measurement_id), calculate timepoints from datetime in hours
    od_meas['datetime'] = od_meas['datetime'].apply(lambda x: datetime.fromisoformat(x))
    baseline_timepoints = od_meas.groupby('measurement_id')['datetime'].transform("min")
    od_meas['timepoint'] = od_meas['datetime'] - baseline_timepoints
    od_meas['timepoint'] = od_meas['timepoint'].apply(lambda x: x.total_seconds()/3600)
    
    # Collection data frame
    metrics_all = pd.DataFrame()
    
    # For every plate_name, create a GrowthPlate
    
    value = 'od' #'OD_sub_background' # or
    
    for plate_name, plate_data in od_meas.groupby('plate_name'):
    
        pivot_data = pd.pivot(plate_data, 
                              index = 'timepoint',
                              values = value,
                              columns = 'measurement_id'
                             ).reset_index(
                             ).rename(columns=({'timepoint':'time'})
                             ).sort_values('time'
                             ).astype(float)
        mapping = pd.DataFrame(data = {
            'm_id':plate_data['measurement_id'].unique(),
            'plate_ID':[plate_name]*len(plate_data['measurement_id'].unique())
        }
                              ).set_index('m_id')
    
        plate = GrowthPlate(data=pivot_data,key=mapping)
        plate.computeBasicSummary()
        plate.raiseData()
        plate.logData()
        plate.subtractBaseline(True,poly=False)
        
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
    
        metrics_all = pd.concat([metrics_all, metrics])
    
    
    # Massage data so it will fit into database
    # i.e., replace infinity values, replace very large values with NA,
    # and round to 3 dec places, 
    
    metrics_all = metrics_all.reset_index().rename(columns={
        'm_id':'measurement_id',
        'OD_Max': 'max_od'}
                                                  ).copy()
    metrics_all['lag_time'] = metrics_all['lag_time'].replace([np.inf, -np.inf], 10000)
    metrics_all['doubling_time'] = metrics_all['doubling_time'].replace([np.inf, -np.inf], 10000)
    metrics_all['max_od'] = metrics_all['max_od'].apply(lambda x: round(x, 3))
    metrics_all['growth_rate'] = metrics_all['growth_rate'].apply(lambda x: round(x, 3))
    
    crazy_lag_time = (metrics_all['lag_time'] >= 10000) | (metrics_all['lag_time'] < 0)
    metrics_all['lag_time'] = np.where(
        crazy_lag_time,
        [pd.NA]*len(metrics_all),
        metrics_all['lag_time'].apply(round)
    )
    
    crazy_doubling_time = (metrics_all['doubling_time'] >= 10000) | (metrics_all['doubling_time'] < 0)
    metrics_all['doubling_time'] = np.where(
        crazy_doubling_time,
        [pd.NA]*len(metrics_all),
        metrics_all['doubling_time'].apply(round)
    )
    
    od_max_near_0 = metrics_all['max_od'] < 0.0001
    metrics_all['max_od'] = np.where(
        od_max_near_0,
        [pd.NA]*len(metrics_all),
        metrics_all['max_od']
    )
    growth_rate_near_0 = metrics_all['growth_rate'] < 0.0001
    metrics_all['growth_rate'] = np.where(
        growth_rate_near_0,
        [pd.NA]*len(metrics_all),
        metrics_all['growth_rate']
    )
    
    metrics_all.to_sql('growth_measurement', engine, index=False, if_exists='append')


def query_OD(engine, experiment_id, strain_id):
    '''
    Queries db for all samples from specified experiment and returns all 
    associated od_measurements as pandas DataFrame. Plots all samples of
    specified strain_id.

    Args:
        experiment_id (str): Experiment id. Must be in db.
        strain_id (str): Strain id. Must be in db.
    Returns:
        DataFrame    
    '''
    

    # # Hardcode experiment id. To be changed later
    # experiment_id = 'ALE1b'

    # Check validity of passed arguments
    db_experiments = pd.read_sql(
        "SELECT experiment.id FROM experiment", engine
    )['id'].to_list()

    db_strains = pd.read_sql(
        "SELECT strain.id FROM strain", engine
    )['id'].to_list()

    if (
        experiment_id not in db_experiments
    ) or (
        strain_id not in db_strains
    ):
        print(f"Check if experiment and strain are registered in the db.")
        
        return None
        
    # Hardcoded query. This can be made more flexible later.
    query = """
    SELECT 
        experiment.id,
        sample.name, sample.passage, sample.plate,
        sample.strain_id, strain.long_name,
        sample.growth_condition_id, growth_condition.long_name AS gc_name, growth_condition.carbon_source,
        measurement.type,
        od_measurement.datetime, od_measurement.timepoint, od_measurement.od, od_measurement.background
    FROM 
        experiment
        INNER JOIN sample ON sample.experiment_id = experiment.id
        INNER JOIN measurement ON measurement.sample_id = sample.name
        # INNER JOIN od_measurement ON od_measurement.measurement_id = measurement.id
        INNER JOIN od_measurement ON od_measurement.operation_id = measurement.operation_id
        INNER JOIN strain ON strain.id = sample.strain_id
        INNER JOIN growth_condition ON growth_condition.id = sample.growth_condition_id
        
    WHERE 
        (experiment.id=%(experiment)s) AND (sample.strain_id=%(strain)s)
    """
    
    selection = pd.read_sql(
        query, engine, params={'experiment': experiment_id, 'strain': str(strain_id)}
    ).rename(
        columns={'id': 'experiment_id',
                 'name': 'sample_name',
                 'type': 'measurement_type',
                 'long_name':'strain_name'}
    )

    return selection


def query_all_OD(engine, strain_id):
    '''
    Queries db for all samples from specified strain and returns all 
    associated od_measurements as pandas DataFrame. Plots all samples of
    specified strain_id.

    Args:
        experiment_id (str): Experiment id. Must be in db.
        strain_id (str): Strain id. Must be in db.
    Returns:
        DataFrame    
    '''
    

    # # Hardcode experiment id. To be changed later
    # experiment_id = 'ALE1b'

    db_strains = pd.read_sql(
        "SELECT strain.id FROM strain", engine
    )['id'].to_list()

    if strain_id not in db_strains:
        print(f"Check if strain is registered in the db.")
        
        return None
        
    # Hardcoded query. This can be made more flexible later.
    query = """
    SELECT 
        sample.experiment_id, sample.name, sample.replicate, sample.passage, 
        sample.strain_id, strain.long_name,
        sample.growth_condition_id, growth_condition.carbon_source,
        measurement.type,
        od_measurement.timepoint, od_measurement.od, od_measurement.background
    FROM 
        sample
        INNER JOIN measurement ON measurement.sample_id = sample.name
        INNER JOIN od_measurement ON od_measurement.measurement_id = measurement.id
        INNER JOIN strain ON strain.id = sample.strain_id
        INNER JOIN growth_condition ON growth_condition.id = sample.growth_condition_id
        
    WHERE 
        (sample.strain_id=%(strain)s)
    """
    
    selection = pd.read_sql(
        query, engine, params={'strain': str(strain_id)}
    ).rename(
        columns={'experiment_id': 'experiment_id',
                 'name': 'sample_name',
                 'type': 'measurement_type',
                 'long_name':'strain_name'}
    )

    return selection
    
def query_growth_rate(engine, experiment_id, strain_id):
    '''
    Queries db for all samples from specified experiment and returns all 
    associated growth_measurements as pandas DataFrame.

    Args:
        experiment_id (str): Experiment id. Must be in db.
        strain_id (str): Strain id. Must be in db.
    Returns:
        DataFrame    
    '''


    # # Hardcode experiment id. To be changed later
    # experiment_id = 'ALE1b'

    # Check validity of passed arguments
    db_experiments = pd.read_sql(
        "SELECT experiment.id FROM experiment", engine
    )['id'].to_list()

    db_strains = pd.read_sql(
        "SELECT strain.id FROM strain", engine
    )['id'].to_list()

    if (
        experiment_id not in db_experiments
    ) or (
        strain_id not in db_strains
    ):
        print(f"Check if experiment and strain are registered in the db.")
        
        return None
        
    # Hardcoded query. This can be made more flexible later.
    query = """
    SELECT
        experiment.id,
        sample.name, sample.passage, sample.plate,
        sample.strain_id, strain.long_name,
        sample.growth_condition_id, growth_condition.long_name AS gc_name, growth_condition.carbon_source,
        measurement.type,
        growth_measurement.growth_rate, growth_measurement.doubling_time, growth_measurement.max_od
    FROM 
        experiment
        INNER JOIN sample ON sample.experiment_id = experiment.id
        INNER JOIN measurement ON measurement.sample_id = sample.name
        # INNER JOIN growth_measurement ON growth_measurement.measurement_id = measurement.id
        INNER JOIN growth_measurement ON growth_measurement.operation_id = measurement.operation_id
        INNER JOIN strain ON strain.id = sample.strain_id
        INNER JOIN growth_condition ON growth_condition.id = sample.growth_condition_id
        
    WHERE 
        (experiment.id=%(experiment)s) AND (sample.strain_id=%(strain)s)
    """
    
    selection = pd.read_sql(
        query, engine, params={'experiment': experiment_id, 'strain': str(strain_id)}
    ).rename(
        columns={'id': 'experiment_id',
                 'name': 'sample_name',
                 'type': 'measurement_type',
                 'long_name':'strain_name'}
    )

    return selection


def growth_exp_ETL(engine, mio, operation, experiment, minio_folder):

    minio_bucket_name = 'synbio'
    path_to_data = minio_folder + 'hour_OD.csv'
    path_to_metadata = minio_folder + 'metadata.csv'

    # Read metadata
    response = mio.get_object(
        minio_bucket_name, path_to_metadata
        )
    csv_data = response.data
    mapping = pd.read_csv(io.BytesIO(csv_data))
    
    # Read data
    response = mio.get_object(
        minio_bucket_name, path_to_data
        )
    csv_data = response.data
    raw = pd.read_csv(io.BytesIO(csv_data))
    
    data_cols = ['time'] + mapping['sample_number'].to_list()
    
    data = raw.copy()
    data.columns = data_cols
    
    # Prepare data tables for db upload
    mapping['experiment_id'] = experiment
    mapping['operation_id'] = operation
    mapping['measurement_type'] = 'growth'
    mapping['filename'] = path_to_data
    mapping['sample_name'] = ('E:' + experiment + 
                              '.S:' + mapping['strain_id'].astype(str) + 
                              '.C:' +  mapping['growth_condition_id'].astype(str) + 
                              '.R:' + mapping['replicate'].astype(str))
    
    samples = mapping[['sample_name', 'experiment_id', 'growth_condition_id', 'strain_id', 'replicate']]
    samples.rename(columns={'sample_name':'name'}, inplace=True)
    samples.to_sql('sample', engine, index=False, if_exists='append')
    
    measurements = mapping[['sample_name', 'operation_id', 'measurement_type', 'filename']]
    measurements.rename(columns={'sample_name':'sample_id', 'measurement_type':'type'}, inplace=True)
    measurements.to_sql('measurement', engine, index=False, if_exists='append')
    
    od = pd.melt(data, id_vars=['time'], var_name='sample_number', value_name='od')
    
    sample_names = tuple(samples['name'])
    
    meas_from_db = pd.read_sql(
        f"SELECT `id`, `sample_id` FROM `measurement` WHERE `sample_id` IN {sample_names};",
        engine
        ).rename(columns={'id': 'measurement_id'})
    
    od_meas = od.merge(mapping[['sample_name', 'sample_number']], on='sample_number', how='inner'
                        ).merge(meas_from_db, left_on='sample_name', right_on='sample_id', how='inner'
                               ).drop(['sample_name','sample_id', 'sample_number'], axis=1
                                     ).rename(columns={'time':'timepoint'})
    
    od_meas.to_sql('od_measurement', engine, index=False, if_exists='append')
    
    
    mapping.set_index('sample_number', inplace=True)
    
    plate = GrowthPlate(data=data,key=mapping)
    plate.computeBasicSummary()
    plate.raiseData()
    plate.logData()
    plate.subtractBaseline(True,poly=False) 
    
    metrics = plate.key[['sample_name','OD_Max']]
    
    for sample in plate.data.columns:
        # fit curve
        thisSample = pd.concat([plate.time, plate.data[sample]], axis=1).dropna() #### There were some missing values, so had to drop those, otherwise curve fitting won't work
        thisSample.columns = ['Time', 'OD']
        gm = GrowthModel(thisSample,ARD=True)
        gm.fit()
        curve = gm.run()
        metrics.loc[sample, 'growth_rate'] = curve.gr
        metrics.loc[sample, 'doubling_time'] = curve.td * 60
        metrics.loc[sample, 'lag_time'] = curve.lagC * 60
    
    # Massage data so it will fit into database
    # i.e., replace infinity values, replace very large values with NA,
    # and round to 3 dec places, 
    metrics = metrics.merge(meas_from_db, left_on='sample_name',right_on='sample_id', how='inner'
                           ).rename(columns={'OD_Max': 'max_od'}
                                   ).drop(['sample_name', 'sample_id'], axis=1)
    
    metrics['lag_time'] = metrics['lag_time'].replace([np.inf, -np.inf], 10000)
    metrics['doubling_time'] = metrics['doubling_time'].replace([np.inf, -np.inf], 10000)
    metrics['max_od'] = metrics['max_od'].apply(lambda x: round(x, 3))
    metrics['growth_rate'] = metrics['growth_rate'].apply(lambda x: round(x, 3))
    
    crazy_lag_time = (metrics['lag_time'] >= 10000) | (metrics['lag_time'] < 0)
    metrics['lag_time'] = np.where(
        crazy_lag_time,
        [pd.NA]*len(metrics),
        metrics['lag_time'].apply(round)
    )
    
    crazy_doubling_time = (metrics['doubling_time'] >= 10000) | (metrics['doubling_time'] < 0)
    metrics['doubling_time'] = np.where(
        crazy_doubling_time,
        [pd.NA]*len(metrics),
        metrics['doubling_time'].apply(round)
    )
    
    od_max_near_0 = metrics['max_od'] < 0.0001
    metrics['max_od'] = np.where(
        od_max_near_0,
        [pd.NA]*len(metrics),
        metrics['max_od']
    )
    growth_rate_near_0 = metrics['growth_rate'] < 0.0001
    metrics['growth_rate'] = np.where(
        growth_rate_near_0,
        [pd.NA]*len(metrics),
        metrics['growth_rate']
    )
    
    metrics.to_sql('growth_measurement', engine, index=False, if_exists='append')



def flask_data_from_ALE_samples_ETL(engine, mio, operation, experiment, minio_folder):

    minio_bucket_name = 'synbio'
    path_to_data = minio_folder + 'hour_OD.csv'
    path_to_metadata = minio_folder + 'metadata.csv'

    # Read metadata
    response = mio.get_object(
        minio_bucket_name, path_to_metadata
        )
    csv_data = response.data
    mapping = pd.read_csv(io.BytesIO(csv_data))
    
    # Read data
    response = mio.get_object(
        minio_bucket_name, path_to_data
        )
    csv_data = response.data
    raw = pd.read_csv(io.BytesIO(csv_data))
    
    data_cols = ['time'] + mapping['sample_number'].to_list()
    
    data = raw.copy()
    data.columns = data_cols
    
    # Prepare data tables for db upload
    mapping['experiment_id'] = experiment
    mapping['operation_id'] = operation
    mapping['measurement_type'] = 'growth'
    mapping['filename'] = path_to_data
    mapping['parent_plate_name'] = 'E:' + experiment + '.' + 'P:' + mapping['parent_plate'].astype(str)
    mapping['sample_name'] = ('E:' + experiment +
                              '.PP:' + mapping['parent_plate'].astype(str) +
                              '.PW:' + mapping['parent_well'].astype(str) +
                              '.S:' + mapping['strain_id'].astype(str) + 
                              '.C:' +  mapping['growth_condition_id'].astype(str) + 
                              '.R:' + mapping['replicate'].astype(str))
    
    samples = mapping[['sample_name', 'experiment_id', 'parent_plate_name', 'parent_well', 'growth_condition_id', 'strain_id', 'replicate', 'parent_sample_name']]
    samples.rename(columns={'sample_name':'name', 'parent_plate_name':'plate', 'parent_well':'well'}, inplace=True)

    samples.to_sql('sample', engine, index=False, if_exists='append')
    
    measurements = mapping[['sample_name', 'operation_id', 'measurement_type', 'filename']]
    measurements.rename(columns={'sample_name':'sample_id', 'measurement_type':'type'}, inplace=True)
    measurements.to_sql('measurement', engine, index=False, if_exists='append')
    
    od = pd.melt(data, id_vars=['time'], var_name='sample_number', value_name='od')
    
    sample_names = tuple(samples['name'])
    
    meas_from_db = pd.read_sql(
        f"SELECT `id`, `sample_id` FROM `measurement` WHERE `sample_id` IN {sample_names};",
        engine
        ).rename(columns={'id': 'measurement_id'})
    
    od_meas = od.merge(mapping[['sample_name', 'sample_number']], on='sample_number', how='inner'
                        ).merge(meas_from_db, left_on='sample_name', right_on='sample_id', how='inner'
                               ).drop(['sample_name','sample_id', 'sample_number'], axis=1
                                     ).rename(columns={'time':'timepoint'})
    
    od_meas.to_sql('od_measurement', engine, index=False, if_exists='append')
    
    
    mapping.set_index('sample_number', inplace=True)
    
    plate = GrowthPlate(data=data,key=mapping)
    plate.computeBasicSummary()
    plate.raiseData()
    plate.logData()
    plate.subtractBaseline(True,poly=False) 
    
    metrics = plate.key[['sample_name','OD_Max']]
    
    for sample in plate.data.columns:
        # fit curve
        thisSample = pd.concat([plate.time, plate.data[sample]], axis=1).dropna() #### There were some missing values, so had to drop those, otherwise curve fitting won't work
        thisSample.columns = ['Time', 'OD']
        gm = GrowthModel(thisSample,ARD=True)
        gm.fit()
        curve = gm.run()
        metrics.loc[sample, 'growth_rate'] = curve.gr
        metrics.loc[sample, 'doubling_time'] = curve.td * 60
        metrics.loc[sample, 'lag_time'] = curve.lagC * 60
    
    # Massage data so it will fit into database
    # i.e., replace infinity values, replace very large values with NA,
    # and round to 3 dec places, 
    metrics = metrics.merge(meas_from_db, left_on='sample_name',right_on='sample_id', how='inner'
                           ).rename(columns={'OD_Max': 'max_od'}
                                   ).drop(['sample_name', 'sample_id'], axis=1)
    
    metrics['lag_time'] = metrics['lag_time'].replace([np.inf, -np.inf], 10000)
    metrics['doubling_time'] = metrics['doubling_time'].replace([np.inf, -np.inf], 10000)
    metrics['max_od'] = metrics['max_od'].apply(lambda x: round(x, 3))
    metrics['growth_rate'] = metrics['growth_rate'].apply(lambda x: round(x, 3))
    
    crazy_lag_time = (metrics['lag_time'] >= 10000) | (metrics['lag_time'] < 0)
    metrics['lag_time'] = np.where(
        crazy_lag_time,
        [pd.NA]*len(metrics),
        metrics['lag_time'].apply(round)
    )
    
    crazy_doubling_time = (metrics['doubling_time'] >= 10000) | (metrics['doubling_time'] < 0)
    metrics['doubling_time'] = np.where(
        crazy_doubling_time,
        [pd.NA]*len(metrics),
        metrics['doubling_time'].apply(round)
    )
    
    od_max_near_0 = metrics['max_od'] < 0.0001
    metrics['max_od'] = np.where(
        od_max_near_0,
        [pd.NA]*len(metrics),
        metrics['max_od']
    )
    growth_rate_near_0 = metrics['growth_rate'] < 0.0001
    metrics['growth_rate'] = np.where(
        growth_rate_near_0,
        [pd.NA]*len(metrics),
        metrics['growth_rate']
    )
    
    metrics.to_sql('growth_measurement', engine, index=False, if_exists='append')


def etl_ALE1a(engine, mio):
    '''
    Extracts and transforms plate reader data from MinIO storage,
    and stores it in MySQL database.

    Args:

    Returns:
    
    '''   
    
    # Functions
    
    
    def get_plate_reader_filenames(
        minio_bucket, minio_path_to_plate_reader_files, regexpattern
    ):
        # Return list of file names that match the expected file name pattern.
        objects = mio.list_objects(
            minio_bucket, prefix=minio_path_to_plate_reader_files
        )
        filepaths = [obj.object_name for obj in objects]
        filenames = [fn.split('/')[1] for fn in filepaths]
        plate_reader_files = [fn for fn in filenames if re.match(regexpattern, fn)]
        
        return plate_reader_files
    
    
    def get_transfers(batch_dict, batch, plate, transfer):
        # Which transfer since beginning of experiment?
        return batch_dict[batch] + (plate-1)*3 + (transfer)
    
    
    def get_plates(batch_dict, batch, plate):
        # How many plates in total?
        return batch_dict[batch] + plate
    
    
    def get_parent_plate(plate, column):
        # Which plate was the parent sample on?
        if column in (1, 4, 7):
            parent_plate = plate-1
        else:
            parent_plate = plate
        return parent_plate
    
    
    def get_parent_col(col):
        # Which column was the parent sample in?
        cols = (1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 10, 11)
        parent_cols = (3, 1, 2, 6, 4, 5, 9, 7, 8, 0, 0, 0)
        parent_col_dict = dict(zip(cols, parent_cols))
        return parent_col_dict[col]
    
    
    def get_well_name(row, col):
        # Return well name based on plate rows/cols.
        well_name = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'][row] + str(col+1)
        return well_name
    
    
    # Hardcoded variables. 
    # Later, these things should be supplied when the experimental team 
    # registers the experiment through the web UI.
    minio_bucket_name = 'synbio'
    path_to_plate_reader_files = 'ALE1a_OD_data/'
    path_to_layout_files = 'plate_layouts_ALE1a/'
    fname_pattern = re.compile(
        r'(?P<uniqueID>\w+)_(?P<plate>\d+)_(?P<transfer>[1-3])_(?P<timepoint>\d+).txt'
        )
    experiment_id = path_to_plate_reader_files.split('_')[0]
    plate_type = '96_shallow' # could also be deep well plate
    start_date = '2025-03-01'
    exp_index = 1
    exp_type = 'autoALE'
    description = ''
    protocol_id = 'mock_1a_protocol'
    lab_id = 1
    contact_id = 1
    operation_id = f"{experiment_id}_operation"
    measurement_type = 'growth'
    plate_reader_filenames = get_plate_reader_filenames(
        minio_bucket_name, path_to_plate_reader_files, fname_pattern
        )
    transfer_layout = read_mio_csv(
        mio, minio_bucket_name, path_to_layout_files + 'transfer_layout.csv'
        )
    rep_layout = read_mio_csv(
        mio, minio_bucket_name, path_to_layout_files + 'replicate_layout.csv'
        )
    strain_layout = read_mio_csv(
        mio, minio_bucket_name, path_to_layout_files + 'strain_layout.csv'
        )
    gc_layout = read_mio_csv(
        mio, minio_bucket_name, path_to_layout_files + 'growth_condition_layout.csv'
        )
    
    # Upload this experiment and its operation (i.e., procedure) to the db.
    operation_dict = {
        'id': [operation_id],
        'protocol_id': [protocol_id],
        'lab_id': [lab_id],
        'contact_id': [contact_id],
        'timestamp': [start_date]
    }
    exp_operation_df = pd.DataFrame.from_dict(operation_dict)
    exp_operation_df.to_sql('operation', engine, index=False, if_exists='append')
    
    exp_dict = {
        'id': [experiment_id],
        'type': [exp_type],
        'start_date': [start_date],
        'index': [exp_index],
        'description': [description],
        'operation_id': [operation_id]
    }
    new_exp_df = pd.DataFrame.from_dict(exp_dict)
    new_exp_df.to_sql('experiment', engine, index=False, if_exists='append')
    
    
    # Initialize data df
    data = pd.DataFrame()
    
    # Read info from plate reader file names and file content into df
    for f in plate_reader_filenames:
    
        try:
    
            # Initialize row in dataframe
            data_row = {}
            match = fname_pattern.match(f)
    
            # Parse info contained in plate reader file name
            data_row['experiment'] = 'ALE1a'  #str(match.group('experiment'))
            data_row['file_ID'] = str(match.group('uniqueID'))
            data_row['timestamp'] = int(match.group('timepoint'))
            data_row['plate_index'] = int(match.group('plate'))
            # t_transfer indicates which plate cols were most recently innoculated.
            data_row['t_transfer'] = int(match.group('transfer'))
    
            # Read plate reader files into dataframe
            response = mio.get_object(
                minio_bucket_name, path_to_plate_reader_files+f
                )
            csv_data = response.data
            plate_data = pd.read_csv(io.BytesIO(csv_data), header=None)
            for row in range(8):
                for col in range(12):
                    data_row['row'] = row
                    data_row['column'] = col
                    data_row['OD'] = plate_data.iloc[row, col]
                    data = pd.concat([data, pd.Series(data_row).to_frame().T])
                    
        except Exception as e:
            print(f"Error: {e}")
            
        finally: 
            response.close()
            response.release_conn()
    
    data.reset_index(inplace=True, drop=True)
    
    # Translate row and column numbers to well names
    data['well'] = data.apply(
        lambda x: get_well_name(x['row'], x['column']), axis=1
        )
    # # Translate timestamp into isoformat time
    # data['datetime'] = data['timestamp'].apply(
    #     lambda x: datetime.fromtimestamp(x).isoformat()
    #     )
    
    # Appending metadata
    data['filename'] = path_to_plate_reader_files
    data['measurement_type'] = measurement_type
    data['experiment'] = experiment_id
    data['plate_type'] = plate_type
    data['start_date'] = start_date
    data['exp_index'] = exp_index
    data['operation_id'] = operation_id
    data['layout_filename'] = path_to_layout_files 
    # Given location on plate (row, col) and layout files, get the strain, 
    # growth condition, replicate number, and transfer_l for each well.
    # (transfer_l indicates transfer based on plate location.)
    data['strain'] = data.apply(
        lambda x: strain_layout.iloc[x['row'], x['column']], axis=1
        ).astype("Int64")
    data['replicate'] = data.apply(
        lambda x: rep_layout.iloc[x['row'], x['column']], axis=1
        ).astype("Int64")
    data['gc'] = data.apply(
        lambda x: gc_layout.iloc[x['row'], x['column']], axis=1
        ).astype("Int64")
    data['l_transfer'] = data.apply(
        lambda x: transfer_layout.iloc[x['row'], x['column']], axis=1
        ).astype("Int64")
    
    # Determine innoculation status
    # Explanation: All wells on a plate are read at each timepoint, but only some
    # of the wells will be inoculated at any given timepoint. Only wells that have
    # media in them (i.e., growth condition not NA) are considered *samples*.
    # Samples that do not have a strain designation are negative controls.
    # Depending on the location on the plate, samples can have 11, 22, or 33
    # timepoints. The calculation below determines which transfer to assign a
    # given OD reading to.
    
    # Calculate actual transfer number for a given reading
    data['l-t_transfer'] = data.apply(
        lambda x: x['l_transfer']-x['t_transfer'], axis=1
        )
    data['transfer'] = np.where(
        (data['l-t_transfer']<=0), data['l_transfer'], 1
        )   
    
    # This experiment only had two transfers on the 5th plate.
    # This breaks the code that computes the parent sample later.
    # For now, I AM DELETING ALL DATA FROM PLATE 5.
    # There needs to be a solution for that!!!
    
    data = data.loc[~(
        (data['file_ID']=='01JNKYBNV34M1FFSDQQ8898C91') &
        (data['plate_index']==5)
    )]
    
    
    # Compute total number of plates, transfers for this experiment
    # There are 3 transfers per plate, x plates per batch, and y batches per
    # experiment. A batch is defined as a continuous run of the robot without
    # any interruption. All measurements in a batch will have the same file_ID.
    # Each batch starts with plate=1, transfer=1, timepoint=1. To accurately
    # calculate the cummulative number of transfers since the beginning of the
    # experiment for each data point, we need to calculate the number of transfers
    # for the individual batches and their cumulative number (in their correct 
    # order).
    
    batches = ['01JNKYBNV34M1FFSDQQ8898C91', '01JPAMF2MM1SE25WH7KQS31ZS4']
    
    batch_n_transfers = [0]  # transfers for each of the batches.
    # Start with 0, because the first batch has no prior transfers.
    batch_n_plates = [0]  # number of plates for each of the batches.
    # Start with 0, because first batch has no prior plates.
    
    for batch in batches:
        # Loop through each batch and calculate the number of transfers
        batch_data = data.loc[
            (data['file_ID'] == batch) & (~pd.isna(data['transfer']))
        ]
        # From last measurement in batch, total number of plates and transfers in
        # this batch
        max_plate = batch_data.sort_values(['plate_index', 'transfer']).iloc[-1]
        batch_n_plates.append(max_plate['plate_index'])
        batch_transfers = (max_plate['plate_index']-1)*3 + (max_plate['transfer'])
        batch_n_transfers.append(batch_transfers)
    
    # cumulative transfers (# of transfers since the start of the experiment)
    batch_cumsum_t = np.cumsum(np.array(batch_n_transfers[:-1])) 
    batch_dict_t = dict(zip(batches, batch_cumsum_t))
    # cumulative plates (# of plates since the start of the experiment)
    batch_cumsum_p = np.cumsum(np.array(batch_n_plates[:-1])) 
    batch_dict_p = dict(zip(batches, batch_cumsum_p))
    
    data['cum_transfer'] = np.where( ################### It's never 0!!!!!!!!!!!!!!!!!!!!!!! what's the point?
        data['transfer']==0, np.array([0]*len(data)), 
        data.apply(
            lambda x: get_transfers(
                batch_dict_t, x['file_ID'], x['plate_index'], x['transfer']
                ), axis=1)
                ) 
    data['cum_plate'] = data.apply(
        lambda x: get_plates(
            batch_dict_p, x['file_ID'], x['plate_index']
            ), axis=1
            )
    
    
    # Translate timestamp into datetime
    # data['datetime'] = data['timestamp'].apply(
    #     lambda x: datetime.fromtimestamp(x).isoformat()
    #     )
    
    # data['total_h'] = pd.NA
    data['total_h'] = np.where(
        ((data['l-t_transfer']<0) &  ~(data['l-t_transfer'].isna())), 
        (data['cum_transfer'] - 1)*10 + data['l-t_transfer']*(-11) + (data['timestamp'] - 1), 
        (data['cum_transfer'] - 1)*10 + (data['timestamp'] - 1)
    )
    
    # data['total_h'] = (data['cum_transfer'] - 1)*10 + (data['timestamp'] - 1)
    start_datetime = datetime.fromisoformat('2025-03-01T09:00:00')  # MOCK START DATE!!
    data['datetime'] = start_datetime + data['total_h'].apply(lambda x: timedelta(hours=x))
    
    
    # data['strain'] is NA (i.e., negative control) if not yet inoculated
    data['strain'] = np.where(
        (data['l-t_transfer']>0), 
        np.array([pd.NA]*len(data)), 
        data['strain']
        )
    data['strain'] = data['strain'].astype("Int64")
    data['replicate'] = np.where(
        pd.isna(data['strain']), 
        np.array([pd.NA]*len(data)), 
        data['replicate']
        )
    data['replicate'] = data['replicate'].astype("Int64")
    
    # Calculate a background value for each plate reader measurement
    # (based on the wells that only contain media)
    data['background'] = pd.NA
    data['background'] = data.groupby(
        ['experiment', 'plate_index', 'timestamp']
        )['OD'].transform(
        lambda x: x[data.loc[x.index, 'strain'].isna()].mean()
        )
    
    # Compute innoculation timestamp based on the oldest timestamp
    data['innoculation_timestamp'] = pd.NA
    data.loc[
        (~(pd.isna(data['cum_transfer']))) & (~(pd.isna(data['strain']))),
        'innoculation_timestamp'
        ] = data.groupby(
            ['cum_plate', 'cum_transfer']
            )['datetime'].transform("min")
    
    # Assign parent samples
    # Only innoculated samples (i.e., not neg. controls) can have parent samples.
    # Only samples after passage 1 have parent samples (passage 1 parents will 
    # have to be manually assigned)
    data['parent_plate'] = pd.NA
    data.loc[
        ~(pd.isna(data['strain'])) & (data['cum_transfer'] != 1),
        'parent_plate'] = data.apply(
            lambda x: get_parent_plate(
                x['cum_plate'], x['column']
                ), axis=1)
    data['parent_well'] = pd.NA
    data.loc[
        ~pd.isna(data['strain']) & (data['cum_transfer'] != 1),
        'parent_well'] = data.apply(
            lambda x: get_well_name(x['row'], get_parent_col(x['column'])), axis=1)
    
    # Assign plate, well, and sample names
    data = data.convert_dtypes()
    data['plate_name'] = (
        'E:' + data['experiment'] + '.P:' + data['cum_plate'].astype(str)
    )
    data['well_name'] = data['plate_name'] + '.W:' + data['well'].astype(str)
    data['sample_name']= pd.NA
    data.loc[(~(pd.isna(data['gc']))), 'sample_name'] = (
        data['well_name'] + '.S:' + data['strain'].astype(str) + '.C:' +
        data['gc'].astype(str) + '.R:' + data['replicate'].astype(str) +
        '.T:' + data['cum_transfer'].astype(str)
    )
    
    # get parent_id
    data['parent_id'] = pd.NA
    data.loc[
        (~(pd.isna(data['strain']))) &
        (~(pd.isna(data['gc']))) &
        (data['cum_transfer'] != 1),
        'parent_id'] = (
            'E:' + data['experiment'] + '.P:' + data['parent_plate'].astype(str) + 
            '.W:' + data['parent_well'].astype(str) + '.S:' + data['strain'].astype(str) + 
            '.C:' + data['gc'].astype(str) + '.R:' + data['replicate'].astype(str) + 
            '.T:' + (data['cum_transfer']-1).astype(str)
        )
    data = data.convert_dtypes()
    
    # Reformat data for database upload
    # 1) Plates
    plates = data[
        ['plate_name', 'experiment', 'plate_type', 'plate_index', 'layout_filename']
        ].drop_duplicates().reset_index(drop=True)
    plates.rename(
        columns={'plate_name': 'id', 'experiment': 'experiment_id'}, inplace=True
        )
    plates.to_sql('plate', engine, index=False, if_exists='append')
    
    # 2) Samples and associated measurements
    # Each sample has a measurement of type 'growth' to which all od_measurements map.
    sample_meas = data.loc[
        ~pd.isna(data['sample_name']) & (~(pd.isna(data['gc'])))
        ].drop_duplicates(
            (['sample_name', 'experiment', 'plate_name', 'well', 'cum_transfer',
              'gc', 'strain', 'innoculation_timestamp', 'replicate', 'parent_id'])
                     ).sort_values('innoculation_timestamp')
    samples = sample_meas[
        (['sample_name', 'experiment', 'plate_name', 'well', 'cum_transfer','gc', 
          'strain', 'innoculation_timestamp', 'replicate', 'parent_id'])
          ].copy()
    samples.rename(
        columns={
            'sample_name': 'name', 'experiment': 'experiment_id', 
            'plate_name': 'plate', 'cum_transfer': 'passage',
            'gc': 'growth_condition_id', 'strain': 'strain_id',
            'parent_id': 'parent_sample_name'
            }, inplace=True)
    samples.to_sql('sample', engine, index=False, if_exists='append')
    measurements = sample_meas[
        ['sample_name', 'operation_id', 'measurement_type', 'filename']
        ].copy()
    measurements.rename(
        columns={'sample_name': 'sample_id', 'measurement_type': 'type'},
        inplace=True
        )
    measurements.to_sql('measurement', engine, index=False, if_exists='append')
    
    # 3) OD_measurements
    # The index of the db table `measurements` autoincrements, so the measurement
    # ids to which the od_measurements of the current dataset belong are not known
    # ahead of time. Therefore, need to download all measurements whith sample ids
    # from this dataset to get their ids, then merge with OD readings.
    sample_names = tuple(sample_meas['sample_name'])
    meas_from_db = pd.read_sql(
        f"SELECT `id`, `sample_id` FROM `measurement` WHERE `sample_id` IN {sample_names};",
        engine
        ).rename(columns={'id': 'measurement_id'})
    od_meas = meas_from_db.merge(
        data, left_on='sample_id', right_on='sample_name', how='inner'
        )[['plate_name', 'measurement_id', 'datetime', 'OD', 'background']].rename(
            columns={'OD': 'od'}
            )
    od_meas.drop('plate_name', axis=1).to_sql('od_measurement', engine, index=False, if_exists='append')

    # 4) Generate growth curves with AMiGA code and store in db   
    # Calculate background value, calculate timepoints in hours.
    
    # First drop duplicate measurement.
    # (There are always two timepoints that are the same in the ALE1a dataset
    # because measurement taken before and after transfer received same
    # timepoint.)
    od_meas = od_meas.drop_duplicates(['measurement_id', 'datetime']).copy()
    od_meas['od_sub_background'] = od_meas['od'] - od_meas['background']
    
    # For each sample (measurement_id), calculate timepoints from datetime in hours
    # od_meas['datetime'] = od_meas['datetime'].apply(lambda x: datetime.fromisoformat(x))
    baseline_timepoints = od_meas.groupby('measurement_id')['datetime'].transform("min")
    od_meas['timepoint'] = od_meas['datetime'] - baseline_timepoints
    od_meas['timepoint'] = od_meas['timepoint'].apply(lambda x: x.total_seconds()/3600)
    
    # Collection data frame
    metrics_all = pd.DataFrame()
    
    # For every plate_name, create a GrowthPlate
    
    value = 'od' #'OD_sub_background' # or
    
    for plate_name, plate_data in od_meas.groupby('plate_name'):
    
        pivot_data = pd.pivot(plate_data, 
                              index = 'timepoint',
                              values = value,
                              columns = 'measurement_id'
                             ).reset_index(
                             ).rename(columns=({'timepoint':'time'})
                             ).sort_values('time'
                             ).astype(float)
        mapping = pd.DataFrame(data = {
            'm_id':plate_data['measurement_id'].unique(),
            'plate_ID':[plate_name]*len(plate_data['measurement_id'].unique())
        }
                              ).set_index('m_id')
    
        plate = GrowthPlate(data=pivot_data,key=mapping)
        plate.computeBasicSummary()
        plate.raiseData()
        plate.logData()
        plate.subtractBaseline(True,poly=False)
        
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
    
        metrics_all = pd.concat([metrics_all, metrics])
    
    
    # Massage data so it will fit into database
    # i.e., replace infinity values, replace very large values with NA,
    # and round to 3 dec places, 
    
    metrics_all = metrics_all.reset_index().rename(columns={
        'm_id':'measurement_id',
        'OD_Max': 'max_od'}
                                                  ).copy()
    metrics_all['lag_time'] = metrics_all['lag_time'].replace([np.inf, -np.inf], 10000)
    metrics_all['doubling_time'] = metrics_all['doubling_time'].replace([np.inf, -np.inf], 10000)
    metrics_all['max_od'] = metrics_all['max_od'].apply(lambda x: round(x, 3))
    metrics_all['growth_rate'] = metrics_all['growth_rate'].apply(lambda x: round(x, 3))
    
    crazy_lag_time = (metrics_all['lag_time'] >= 10000) | (metrics_all['lag_time'] < 0)
    metrics_all['lag_time'] = np.where(
        crazy_lag_time,
        [pd.NA]*len(metrics_all),
        metrics_all['lag_time'].apply(round)
    )
    
    crazy_doubling_time = (metrics_all['doubling_time'] >= 10000) | (metrics_all['doubling_time'] < 0)
    metrics_all['doubling_time'] = np.where(
        crazy_doubling_time,
        [pd.NA]*len(metrics_all),
        metrics_all['doubling_time'].apply(round)
    )
    
    od_max_near_0 = metrics_all['max_od'] < 0.0001
    metrics_all['max_od'] = np.where(
        od_max_near_0,
        [pd.NA]*len(metrics_all),
        metrics_all['max_od']
    )
    growth_rate_near_0 = metrics_all['growth_rate'] < 0.0001
    metrics_all['growth_rate'] = np.where(
        growth_rate_near_0,
        [pd.NA]*len(metrics_all),
        metrics_all['growth_rate']
    )
    
    metrics_all.to_sql('growth_measurement', engine, index=False, if_exists='append')