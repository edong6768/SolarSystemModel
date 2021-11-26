import ssm.celestials as cel
import matplotlib.pyplot as plt

a=cel.solar_system_model.config_JSON_load('data/SEM_config_snapshot_0.json')

a.config_obj_write("data/SEM_model.joblib")
#b=cel.solar_system_model.config_obj_load("data/SEM_model.joblib")
#print(b)

#a.timetravel("2021-10-01T00:00:00.0Z", 3600)
a.timetravel_animation(3600)
#a.plot()


# for i in range(20): 
#     #a.timetravel_diff(d=1)
#     print(a.search_eclipse(60))

a.eclps_judge()
