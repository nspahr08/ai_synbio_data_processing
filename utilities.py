import pandas as pd
import minio
import io
from sqlalchemy.sql import text

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


def del_from_table_where(engine, table, pk, value):

    statement = f"DELETE FROM `{table}` WHERE `{pk}` = '{value}';"
        
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