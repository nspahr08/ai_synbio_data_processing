import pandas as pd


def register_lab(engine, name, location):
    lab_dict = {
        'name': [name],
        'location': [location]
    }
    lab_df = pd.DataFrame.from_dict(lab_dict)
    lab_df.to_sql('lab', engine, index=False, if_exists='append')


def register_people(engine, first_name, last_name, email, lab_id):
    people_dict = {
        'first_name': [first_name],
        'last_name': [last_name],
        'email': [email],
        'lab_id': [lab_id]
    }
    people_df = pd.DataFrame.from_dict(people_dict)
    people_df.to_sql('people', engine, index=False, if_exists='append')


def register_protocol(engine, protocol_name, description, filename):
    protocol_dict = {
        'id': [protocol_name],
        'description': [description],
        'filename': [filename]
    }
    protocol_df = pd.DataFrame.from_dict(protocol_dict)
    protocol_df.to_sql('protocol', engine, index=False, if_exists='append')
    return protocol_name


def register_strain(engine, long_name, short_name, culture = None, parent_strain=None):
    strain_dict = {
        'long_name': [long_name],
        'short_name': [short_name],
        'culture': [culture],
        'parent_strain_id': [parent_strain]
    }
    strain_df = pd.DataFrame.from_dict(strain_dict)
    strain_df.to_sql('strain', engine, index=False, if_exists='append')


def register_growth_condition(engine, 
                              long_name,
                              short_name,
                              carbon_source,
                              temperature=None,
                              agitation_speed=0,
                              minimal_media="Ellen's media",
                              nitrogen_source=None,
                              carbon_concentration=20, 
                              nitrogen_concentration=None,
                              antibiotics=None,
                              antibiotic_concentration=None,
                              filename=None
                              ):
    growth_condition_dict = {
        'long_name': [long_name],
        'short_name': [short_name],
        'temperature': [temperature],
        'agitation_speed': [agitation_speed],
        'minimal_media': [minimal_media],
        'carbon_source': [carbon_source],
        'nitrogen_source': [nitrogen_source],
        'carbon_concentration': [carbon_concentration],
        'nitrogen_concentration': [nitrogen_concentration],
        'antibiotics': [antibiotics],
        'ab_concentration': [antibiotic_concentration],
        'filename': [filename]
    }
    growth_condition_df = pd.DataFrame.from_dict(growth_condition_dict)
    growth_condition_df.to_sql('growth_condition', engine, index=False, if_exists='append')
    

def register_operation(engine, op_id, protocol_id, lab_id, contact_id, timestamp):
    operation_dict = {
        'id': [op_id],
        'protocol_id': [protocol_id],
        'lab_id': [lab_id],
        'contact_id': [contact_id],
        'timestamp': [timestamp]
    }
    operation_df = pd.DataFrame.from_dict(operation_dict)
    operation_df.to_sql('operation', engine, index=False, if_exists='append')
    
    return op_id


def register_experiment(engine, experiment_id, exp_type, start_date, exp_index, description, op_id):
    exp_dict = {
        'id': [experiment_id],
        'type': [exp_type],
        'start_date': [start_date],
        'index': [exp_index],
        'description': [description],
        'operation_id': [op_id]
    }
    new_exp_df = pd.DataFrame.from_dict(exp_dict)
    new_exp_df.to_sql('experiment', engine, index=False, if_exists='append')

def register_sample(engine, name, experiment_id, plate, well, growth_condition_id, strain_id, replicate, passage, parent_sample, innoculation_timestamp):
    sample_dict = {
        'name': [name],
        'experiment_id': [experiment_id],
        'plate': [plate],
        'well': [well],
        'growth_condition_id': [growth_condition_id],
        'strain_id': [strain_id],
        'replicate': [replicate],
        'passage': [passage],
        'parent_sample': [parent_sample],
        'innoculation_timestamp': [innoculation_timestamp]
    }
    new_sample_df = pd.DataFrame.from_dict(sample_dict)
    new_sample_df.to_sql('sample', engine, index=False, if_exists='append')


def register_sample(engine, name, experiment_id, plate, well, growth_condition_id, strain_id, replicate, passage, parent_sample, innoculation_timestamp):
    sample_dict = {
        'name': [name],
        'experiment_id': [experiment_id],
        'plate': [plate],
        'well': [well],
        'growth_condition_id': [growth_condition_id],
        'strain_id': [strain_id],
        'replicate': [replicate],
        'passage': [passage],
        'parent_sample': [parent_sample],
        'innoculation_timestamp': [innoculation_timestamp]
    }
    new_sample_df = pd.DataFrame.from_dict(sample_dict)
    new_sample_df.to_sql('sample', engine, index=False, if_exists='append')


def register_measurement(engine, sample_id, operation_id, measurement_type, filename):
    measurement_dict = {
        'sample_id': [sample_id],
        'operation_id': [operation_id],
        'type': [measurement_type],
        'filename': [filename]
    }
    new_measurement_df = pd.DataFrame.from_dict(measurement_dict)
    new_measurement_df.to_sql('measurement', engine, index=False, if_exists='append')


def register(engine, table, value_dict):
    df = pd.Series(value_dict).to_frame().T
    df.to_sql(table, engine, index=False, if_exists='append')
