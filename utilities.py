import pandas as pd
from minio import Minio
import io
from sqlalchemy.sql import text
import os
from minio.error import S3Error
from typing import Dict, Any, Optional, List
import re


def read_mio_csv(mio, minio_bucket_name, file_path):

    # Read mio csv.
    try:
        response = mio.get_object(minio_bucket_name, file_path)
        csv_data = response.data
        df = pd.read_csv(io.BytesIO(csv_data))
    except Exception as e:
        print(f"Error: {e}")
        df = None
    finally:
        response.close()
        response.release_conn()
        return df


def del_from_table_where(engine, table, field, value):

    statement = f"DELETE FROM `{table}` WHERE `{field}` = '{value}';"

    with engine.connect() as con:
        con.execute(text(statement))
        con.commit()

    print(statement)


def select_from_table_where_in(engine, table, field, value_list):
        
    query = f"""
    SELECT * FROM `{table}` WHERE `{field}` IN %(value_list)s;
    """

    selection = pd.read_sql(
        query, engine, params={'value_list': value_list}
    )

    return selection


def read_local_file(file_path: str) -> pd.DataFrame:
    """
    Read and process data from a local file.
    
    Args:
        file_path: Path to the local file
        
    Returns:
        pd.DataFrame: Processed data
        
    Raises:
        FileNotFoundError: If the file doesn't exist
        ValueError: If the file type is not supported
    """
    
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Local file not found: {file_path}")
        
    # Determine file type and read accordingly
    if file_path.endswith('.csv'):
        df = pd.read_csv(file_path, header=None)
    elif file_path.endswith('.txt'):
        df = pd.read_csv(file_path, header=None)
    # elif file_path.endswith('.xlsx'):
    #     df = pd.read_excel(file_path)
    else:
        raise ValueError(f"Unsupported file type: {file_path}")
        
    return df


def read_minio_file(minio_client, minio_bucket, minio_path: str) -> pd.DataFrame:

    try:
        
        # Get the object from MinIO
        response = minio_client.get_object(
            bucket_name= minio_bucket,
            object_name=minio_path
        )
        
        # Read the data into a DataFrame
        if minio_path.endswith('.csv') or minio_path.endswith('.txt'):
            df = pd.read_csv(io.BytesIO(response.read()), header=None)
        elif minio_path.endswith('.xlsx'):
            df = pd.read_excel(io.BytesIO(response.read()), header=None)
        else:
            raise ValueError(f"Unsupported file type: {minio_path}")
        
        return df
        
    except S3Error as e:
        raise S3Error(f"Error accessing MinIO: {str(e)}")
    except Exception as e:
        raise Exception(f"Error reading from MinIO: {str(e)}")
    finally:
        if 'response' in locals():
            response.close()
            response.release_conn()


def list_local_files(directory_path: str, pattern: str = None) -> List[str]:
    """
    List files in a local directory, optionally filtered by a pattern.
    
    Args:
        directory_path: Path to the directory
        pattern: Optional regex pattern to filter filenames
        
    Returns:
        List[str]: List of file paths
        
    Raises:
        FileNotFoundError: If the directory doesn't exist
    """
    if not os.path.exists(directory_path):
        raise FileNotFoundError(f"Directory not found: {directory_path}")
        
    files = []

    for file in os.listdir(directory_path):
        if pattern is None or re.match(pattern, file.split('/')[-1]):
            files.append(os.path.join(directory_path, file))

    return files


def list_minio_files(minio_client, minio_bucket, prefix: str = "", pattern: str = None) -> List[str]:

    try:
       
        # List objects in the bucket
        objects = minio_client.list_objects(
            bucket_name= minio_bucket,
            prefix=prefix,
            recursive=True
        )
        
        # Filter objects by pattern if provided
        files = []
        for obj in objects:
            if pattern is None or re.match(pattern, obj.object_name.split('/')[-1]):
                files.append(obj.object_name)
                
        return files
        
    except S3Error as e:
        raise S3Error(f"Error accessing MinIO: {str(e)}")
    except Exception as e:
        raise Exception(f"Error listing MinIO files: {str(e)}")


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
        sample.name, sample.passage, 
        sample.strain_id, strain.long_name,
        sample.growth_condition_id, growth_condition.carbon_source,
        measurement.type,
        od_measurement.datetime, od_measurement.od, od_measurement.background
    FROM 
        experiment
        INNER JOIN sample ON sample.experiment_id = experiment.id
        INNER JOIN measurement ON measurement.sample_id = sample.name
        INNER JOIN od_measurement ON od_measurement.measurement_id = measurement.id
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
        sample.name, sample.passage, 
        sample.strain_id, strain.long_name,
        sample.growth_condition_id, growth_condition.carbon_source,
        measurement.type,
        growth_measurement.growth_rate, growth_measurement.doubling_time, growth_measurement.max_od
    FROM 
        experiment
        INNER JOIN sample ON sample.experiment_id = experiment.id
        INNER JOIN measurement ON measurement.sample_id = sample.name
        INNER JOIN growth_measurement ON growth_measurement.measurement_id = measurement.id
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


def get_table_col_names(engine, table):
    query = """
    SELECT COLUMN_NAME 
    FROM INFORMATION_SCHEMA.COLUMNS 
    WHERE TABLE_NAME = %(table)s AND TABLE_SCHEMA = 'anl_synbio';
    """
    selection = pd.read_sql(
    query, engine, params={'table': table}
    )
    return selection


def get_well_name(row, col):
    # Return well name based on plate rows/cols.
    well_name = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'][row] + str(col+1)
    return well_name


def query_samples_by_name(engine, experiment_id, sample_names):
    '''
    Queries db for specified samples from specified experiment and returns
    media and strain info.  
    '''


    # Check validity of passed arguments
    db_experiments = pd.read_sql(
        "SELECT experiment.id FROM experiment", engine
    )['id'].to_list()

    if (experiment_id not in db_experiments):
        print(f"Check if experiment and strain are registered in the db.")
        
        return None
        
    # Hardcoded query. This can be made more flexible later.
    query = """
    SELECT 
        experiment.id,
        sample.name, sample.parent_sample_name,
        sample.plate, sample.well,
        sample.strain_id, strain.long_name,
        sample.growth_condition_id, growth_condition.carbon_source

    FROM 
        experiment
        INNER JOIN sample ON sample.experiment_id = experiment.id
        INNER JOIN strain ON strain.id = sample.strain_id
        INNER JOIN growth_condition ON growth_condition.id = sample.growth_condition_id
        
    WHERE 
        (experiment.id=%(experiment)s) AND (sample.name IN %(sample_names)s)
    """
    
    selection = pd.read_sql(
        query, engine, params={'experiment': experiment_id, 'sample_names': tuple(sample_names)}
    ).rename(
        columns={'id': 'experiment_id',
                 'name': 'sample_name',
                 'long_name':'strain_name'}
    )

    return selection
