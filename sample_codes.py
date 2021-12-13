import ssm.celestials as cel
import matplotlib.pyplot as plt

def celestial_plots():
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
        model.cheongu_plot2d(ax3, 0)

    # 각 plot을 보기 위해서 show
    plt.show()


def search_eclipse():
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

def eclipse_result_analysis():
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


def periap_finding_test():
    # code to find moon's optimal perap angle
    import ssm.celestials as cel
    import math as m
    import datetime
    from dateutil.parser import parse
    import numpy as np
    import matplotlib.pyplot as plt


    true_time=parse("2021-12-04T05:29:11.3Z")
    curr_time=0
    td=0

    opt_time=0
    opt_td=1000000000
    opt_peri=0

    center=-0.11195463893540733
    count=4
    ang_precision=0.1
    t_precision=1
    pad=(m.pi/180)*ang_precision*count/2
    peri_total = list(np.linspace(center-pad, center+pad, count))
    td_total=list()
    #print(peri_total)
    for peri in peri_total:
        a=cel.solar_system_model.config_JSON_load('data/SEM_config_snapshot_0.json', print_info=False)    # reset model

        a.elements[2].orb.periap_ang_rad=peri               # set this iteration's periap angle
        a.elements[2].orb.curr_rad=-peri

        a.timetravel("2021-12-04T00:00:00.0Z", 3600)        # use timetravel to reduce computation time
        # a.timetravel("2021-12-04T05:29:00.0Z", 1)        # use timetravel to reduce computation time
        #a.timetravel("2021-12-04T06:59:00.0Z", 0.1)        # use timetravel to reduce computation time
        curr_time=a.search_eclipse(t_precision, until="2021-12-05T00:00:00.0Z")                 # search eclipse

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

    fig=plt.figure()
    ax = fig.add_subplot()
    ax.plot(np.array(peri_total)*180/m.pi, np.array(td_total)/3600)
    plt.show()


    # a=[-3.2385553898117, -2.9078614262759324, -2.5771674627401646, -2.246473499204397, -1.9157795356686294, -1.5850855721328616, -1.254391608597094, -0.9236976450613263, -0.5930036815255586, -0.26230971798979086, 0.06838424554597688, 0.39907820908174463, 0.7297721726175119, 1.0604661361532797, 1.3911600996890474, 1.7218540632248152, 2.052548026760583, 2.3832419902963506, 2.7139359538321184, 3.044629917367886]
    # b=[-11019.7, -39279.7, -64059.7, -82659.7, -92859.7, -93519.7, -84459.7, -66579.7, -41979.7, -13419.7, 15800.3, 42380.3, 63440.3, 77000.3, 81740.3, 77240.3, 64100.3, 43520.3, 17600.3, -11019.7]
    # plt.plot(np.array(a)*180/m.pi, np.array(b)/3600)
    # plt.plot((180, 180),(-20, 20))
    # plt.plot((-180, -180),(-20, 20))
    # plt.grid(True)
    # plt.show()