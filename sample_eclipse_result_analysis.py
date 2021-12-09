import ssm.celestials as cel
import matplotlib.pyplot as plt
import csv
import json
import datetime
from dateutil.parser import parse

# # 모델 load
# model=cel.solar_system_model.config_JSON_load('data/SEM_config_snapshot_0.json')

# # 해당 날짜에서의 천체 위치로 모델 재설정
# #model.timetravel("2021-10-01T03:00:00.0Z", 3600)

# # 일식 현상 가까운 미래에 순서대로 50개 실제 값과 비교
# diff=list()
# with open("./data/eclipse_time.csv", 'r') as f:
#     for line in csv.reader(f):
#         model.timetravel_diff(d=1) # 하루 뒤로 이동 후 다시 검색하기 위해
#         diff.append((model.search_eclipse(3600)-parse(line[0])).total_seconds())
#         print(diff[-1])

# with open("./data/eclipse_prediction_error.json", "w") as file_json:
#     json.dumps(diff, file_json)


# plt.plot([i for i in range(50)], diff)

diff=list()
with open("./data/eclipse_time.csv", 'r') as f:
    for line in csv.reader(f):
        diff.append((parse(line[1])-parse(line[0])).total_seconds()/3600/24)
        print(diff[-1])

plt.plot([i for i in range(50)], diff)
plt.show()