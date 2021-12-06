import ssm.celestials as cel
import matplotlib.pyplot as plt

model=cel.solar_system_model.config_JSON_load('data/SEM_config_snapshot_0.json')

# 항성계의 각 천체와 공전 궤도 plot
fig1=plt.figure()
ax1=fig1.add_subplot(projection='3d')
model.plot(ax1, frame = 0, sphere = False)

# 관측자 천체의 북반구에서의 천구 현상 plot
fig2=plt.figure()
ax2=fig2.add_subplot(projection='3d')

fig3=plt.figure()
ax3=fig3.add_subplot()

for i in range(5):
    model.timetravel_diff(M=30)
    model.cheongu_plot(ax2)
    model.cheongu_plot2d(ax3)

# 각 plot을 보기 위해서 show
plt.show()
