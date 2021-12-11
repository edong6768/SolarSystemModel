import ssm.celestials as cel
import matplotlib.pyplot as plt

# 모델 load
model=cel.solar_system_model.config_JSON_load('data/SEM_config_snapshot_0.json', False)

# 해당 날짜에서의 천체 위치로 모델 재설정
model.timetravel("2021-10-01T03:00:00.0Z", 3600)

fig=plt.figure()
ax=fig.add_subplot()
# 일식 현상(부분) 가까운 미래에 순서대로 5개 검색
for i in range(5): 
    model.timetravel_diff(d=1) # 하루 뒤로 이동 후 다시 검색하기 위해
    print(model.search_eclipse(3600)) # 1시간 간격으로 검색
    #model.cheongu_plot2d(ax, 0)
    #plt.show()

# # 개기/금환 일식 가까운 과거에 역순으로 5개 검색
# for i in range(5):
#     model.timetravel_diff(d=-1) # 하루 앞으로 이동 후 다시 검색하기 위해
#     print(model.search_eclipse(-60, 1)) # 1분 간격으로 검색
