import ssm.celestials as cel

a=cel.solar_system_model.config_JSON_load('data/SEM_config_snapshot_0.json')

a.config_obj_write("data/SEM_model.joblib")
b=cel.solar_system_model.config_obj_load("data/SEM_model.joblib")
print(b)