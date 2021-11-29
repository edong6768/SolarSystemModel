# code to find moon's optimal perap angle
import ssm.celestials as cel
import math as m
import datetime
from dateutil.parser import parse
import numpy as np
import matplotlib.pyplot as plt


true_time=parse("2021-12-04T07:00:01.0Z")
curr_time=0
td=0

opt_time=0
opt_td=1000000000
opt_peri=0

center=-0.03371705485228767
count=10
precision=0.0001
pad=(m.pi/180)*precision*count
peri_total = list(np.linspace(center-pad, center-pad+1, 2))
td_total=list()

for peri in peri_total:
    a=cel.solar_system_model.config_JSON_load('data/SEM_config_snapshot_0.json', print_info=False)    # reset model

    a.elements[2].orb.periap_ang_rad=peri               # set this iteration's periap angle

    a.timetravel("2021-12-04T06:30:00.0Z", 3600)        # use timetravel to reduce computation time
    a.timetravel("2021-12-04T06:57:00.0Z", 1)        # use timetravel to reduce computation time
    a.timetravel("2021-12-04T06:59:00.0Z", 0.1)        # use timetravel to reduce computation time
    curr_time=a.search_eclipse(0.001, until="2021-12-04T07:01:00.0Z")                 # search eclipse

    td=(true_time-curr_time).total_seconds()            # differece from true date
    if abs(td)<abs(opt_td):                             # update value if it is closer to true date
        opt_peri=peri
        opt_time=curr_time
        opt_td=td
    if td==0: break
    td_total.append(td)            # add td to td_total(for graph)
# next step

print(" optimal peri : {}\n true time : {}\n opti time : {}\n time delta : {}".format(
    opt_peri, true_time.strftime('%Y-%m-%dT%H:%M:%S.0Z'), opt_time.strftime('%Y-%m-%dT%H:%M:%S.0Z'), opt_td))

print(peri_total, td_total)

# fig=plt.figure()
# ax = fig.add_subplot()
# ax.plot(peri_total, td_total)
# plt.show()