import ssm.celestials as cel
import matplotlib.pyplot as plt

a=cel.solar_system_model.config_JSON_load('data/SEM_config_snapshot_0.json', False)

fig=plt.figure()
ax=fig.add_subplot(projection='3d')

# 공전궤도상 천체들의 움직임
a.orbit_animation(fig, ax, 3600*24, 2)

# 천구상 천체들의 움직임
#a.cheongu_animation(fig, ax, 60)

# 천구상 천체들의 움직임
#ax2=fig.add_subplot()
#a.cheongu_animation2d(fig, ax2, 3600)
