import celestials as cel
import matplotlib.pyplot as plt
    
def load_model():
    model=cel.solar_system_model.config_JSON_load('data/random_ssm_config_snapshot_0.json')
    #model=cel.solar_system_model.config_JSON_load('data/SEM_config_snapshot_0.json', print_info=False)
    return model


def animation_tool(model):
    fig=plt.figure()
    print("\n----------< Animation modes >----------\n",
    "--1. orbit animation\n",
    "--2. celestial sphere animation\n",
    "--3. 2d celestial sphere animation\n",
    "--4. local celectial sphere animation\n",
    "--5. 2d local celestial sphere animation\n",
    "---------------------------------------")
    mode = int(input(" Select mode : "))
    if mode==1:
        # 공전궤도상 천체들의 움직임
        print("1. orbit animation mode selected")
        ax=fig.add_subplot(projection='3d')
        model.orbit_animation(fig, ax, 3600*24)
    elif mode==2:
        # 천구상 천체들의 움직임
        print("2. celestial sphere animation mode selected")
        ax=fig.add_subplot(projection='3d')
        model.cheongu_animation(fig, ax, 3600*24)
    elif mode==3:
        # 평면에 투영한 천구상 천체들의 움직임
        print("3. 2d celestial sphere animation mode selected")
        ax2=fig.add_subplot()
        model.cheongu_animation2d(fig, ax2, 3600*24, 3, 1000)
    elif mode==4:
        # 현 관측자 위도/경도에서 천구상 천체들의 움직임
        print("4. local celestial sphere animation mode selected")
        ax=fig.add_subplot(projection='3d')
        model.local_cheongu_animation(fig, ax, 126, 37, 600)
    elif mode==5:
        # 현 관측자 위도/경도에서 평먄에 투영한 천구상 천체들의 움직임
        print("5. 2d local celestial sphere animation mode selected")
        ax2=fig.add_subplot()
        model.local_cheongu_animation2d(fig, ax2, 126, 37, 60, )
    else : print("Wrong number!")

    print("Program terminated")


def plot_tools(model):
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


def eclipse_search_tool(model):
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


model=load_model()
print("""-----------< Tool selection >-----------
--1. orbit animation
--2. celestial sphere animation
--3. 2d celestial sphere animation
--4. local celectial sphere animation
--5. 2d local celestial sphere animation
--6. Quit
-----------------------------------------""")

tools=[animation_tool, plot_tools, eclipse_search_tool]
while(True):
    modes=int(input(" Select mode : "))
    if modes not in list(range(6)): break
    tools[modes-1](model)
print("Program terminated.")
