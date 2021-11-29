import ssm.celestials as cel
import matplotlib.pyplot as plt
from time import sleep

a=cel.solar_system_model.config_JSON_load('data/SEM_config_snapshot_0.json')


#a.config_obj_write("data/SEM_model.joblib")
#b=cel.solar_system_model.config_obj_load("data/SEM_model.joblib")
#print(b)

#a.timetravel("2021-10-01T00:00:00.0Z", 3600)
#a.timetravel_animation(3600)



a.plot()

plt.show()
# a.fig.show()
# for i in range(100):
#     a.plot()
#     a.fig.canvas.flush_events()
#     a.next_iter(3600*24)
#     #sleep(1)



# for i in range(3): 
#     a.timetravel_diff(d=1)
#     print(a.search_eclipse(3600))

# a.eclps_judge()
