import math as m
from matplotlib import artist
import numpy as np
import json
import joblib
import datetime
from dateutil.parser import parse
import ssm.coor_ref as cr
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from timezonefinder import TimezoneFinder 
import pytz


def set_data(plts, data):
    plts.set_data(data[:2])
    plts.set_3d_properties(data[2])



###############################################################
#                     solar system bodies                     #
###############################################################

class cel_body:
    G_kms = 6.67384 * 10**-20

    def __init__(self, name, mass_kg, diameter, color=None):
        #print("Loading {} body info...".format(name), end="\r")
        self.name = name
        self.mass_kg = mass_kg
        self.GM_kms = mass_kg * self.G_kms
        self.diameter_km = diameter
        self.color=color

    def body_data_dict(self):
        return {
            "name": self.name,
            "mass_kg": self.mass_kg,
            "diameter_km": self.diameter_km
        }
    
    def __str__(self):
        return ('\n\n{}\n- mass[kg] : {}\n- diameter[km] : {}\n'.format(self.name, self.mass_kg, self.diameter_km))


    ###################### plot ###########################
    def sphere_coors(self, disp=(0, 0, 0), ratio=1):
        # plot celestial body(sphere)
        u = np.linspace(0, np.pi, 30)
        v = np.linspace(0, 2 * np.pi, 30)
        r = self.diameter_km/2

        x = r*np.outer(np.sin(u), np.sin(v))*ratio+disp[0]
        y = r*np.outer(np.sin(u), np.cos(v))*ratio+disp[1]
        z = r*np.outer(np.cos(u), np.ones_like(v))*ratio+disp[2]
        return x, y, z


    def plot(self, ax, disp=(0, 0, 0), sphere=False):
        # plot celestial position(dot)
        if not sphere:
            if self.color : bp,=ax.plot(*disp, linestyle="", marker='o', color=self.color, label=self.name)
            else : bp,=ax.plot(*disp, linestyle="", marker='o', label=self.name)
        else:
            # plot celestial body(sphere)
            if self.color: bp=ax.plot_surface(*self.sphere_coors(disp), color=self.color, label=self.name)
            else : bp=ax.plot_surface(*self.sphere_coors(disp), label=self.name)
        return bp   # for animation
        


###############################################################
#                      celestial motion                       #
###############################################################
# 공전/자전과 같은 천체운동의 중요한 계수들과 위치 정보 관련 class


# 자전 관련 데이터
class rotate:
    def __init__(self, eqnx_lon_rad, incln_rad, ang_v_rad,rot_ang_rad, body, body_orbit=None):
        #print("Loading {} rotation info...".format(body.name), end="\r")
        self.eqnx_lon_rad = eqnx_lon_rad  # 춘분점 경도 : 기준면상 기준축과 적도승교점 사이 각도(자전축 방향)
        self.incln_rad = incln_rad  # 자전축 기울기

        self.ang_v_rad = ang_v_rad  # 자전각속도
        self.curr_rad = rot_ang_rad % (2 * m.pi)  # 현재 자전 위상각(중심기준은 적도좌표계)

        self.body=body  # 천체obj
        self.body_orbit=body_orbit    # 공전obj
        self.name=body.name+'_rot'  # 이름

        #print("Calculating {} rotation refference frame matrix...".format(body.name),  end="\r")
        # reffrm
        self.reffrm=cr.reffrms(eqnx_lon_rad, incln_rad, self.curr_rad, self.name, *body_orbit.base_coor.tuplify(0)) if body_orbit \
            else cr.reffrms(eqnx_lon_rad, incln_rad, self.curr_rad, self.name) 

    def set_curr_rad(self, new_curr_rad):
        self.curr_rad=new_curr_rad
        self.reset_reffrm()


    # t_diff step시간 후 각도변화(천체 자전)
    def rotational_increment_angle_rad(self, tdelta):
        self.curr_rad+=tdelta*self.ang_v_rad
        self.curr_rad%=2*m.pi

        self.reset_reffrm() # 기준틀 업데이트


    def reset_reffrm(self):
        if self.body_orbit:
            self.reffrm.reset_reffrm(self.eqnx_lon_rad, self.incln_rad, self.curr_rad, self.name, 
                                    *self.body_orbit.base_coor.tuplify(0))
        else: self.reffrm.reset_reffrm(self.eqnx_lon_rad, self.incln_rad, self.curr_rad, self.name) 


    # rotation data 딕셔너리 출력 
    def rot_data_dict(self):
        return {
            "eqnx_lon_rad": self.eqnx_lon_rad,
            "incln_rad": self.incln_rad,
            "ang_v_rad": self.ang_v_rad,
            "rot_ang_rad": self.curr_rad
        }


    def __str__(self):
        return ("rotation characteristics : " + "\n- rotation speed[deg/s] : " + str(self.ang_v_rad*180/m.pi)
                + "\n- orientational parameters[deg] : (" + str(self.eqnx_lon_rad*180/m.pi) + ", " + str(self.incln_rad*180/m.pi) + ")" 
                + "\n\nstatus : " + "\n- angle[deg] : " + str(self.curr_rad*180/m.pi) + "\n- refference frame :" + str(self.reffrm)).replace('\n', '\n  ')


# 공전 관련 데이터(케플러적 이체모델)
class orbit:
    def __init__(self, eccentricity, a_km, node_lon_rad, incln_rad, periap_ang_rad,
                 apsidal_prcs_speed, nodal_prcs_speed, anom_true_rad, body, center_body, center_body_orbit=None):
        #print("Loading {} orbit info...".format(body.name), end="\r")
        # 궤도 특징
        self.e = eccentricity  # 궤도 이심률
        self.a_km = a_km  # 긴반지름

        self.node_lon_rad = node_lon_rad  # 승교점 경도 : 기준면상 기준축과 궤도승교점 사이 각(궤도 진중심에서)
        self.incln_rad = incln_rad  # 궤도 경사 : 기준면에 대한 궤도기울기
        self.periap_ang_rad = periap_ang_rad  # 근점 편각 : 궤도 근점과 승교점 사이 각도(궤도 진중심에서)

        self.apsidal_prcs_speed=apsidal_prcs_speed  # 근점 편각 세차, +는 반시계
        self.nodal_prcs_speed=nodal_prcs_speed  # 승교점 경도 세차, +는 반시계

        # 대상 천체obj 및 궤도상 위치
        self.curr_rad = anom_true_rad  # 진근점 이각 : 중심천체가 있는 초점에서 근점 기준으로 현재 천체 각위치

        self.body=body  # 공전천체obj
        self.center_body=center_body    # 중심천체obj
        self.center_body_orbit=center_body_orbit    # 중심천체공전obj

        # 이름
        self.name=body.name+'_orb'

        #print("Calculating {} orbit refference frame matrix...".format(body.name), end="\r")
        # reffrm
        self.reffrm=cr.reffrms(node_lon_rad, incln_rad, periap_ang_rad, self.name, *center_body_orbit.base_coor.tuplify(0)) \
             if center_body_orbit else cr.reffrms(node_lon_rad, incln_rad, periap_ang_rad, self.name)
        #if center_body_orbit: print("skrtskrt", self.name, center_body_orbit.base_coor.tuplify(0))

        # coor3
        self.home_coor=cr.coor3('s', (self.dist_km, anom_true_rad, m.pi/2))
        self.home_coor.conv_coor_modeOS()
        self.base_coor=self.reffrm.base_conv(self.home_coor, 'sa')

    # 현 각위치 재설정
    def set_curr_rad(self, new_curr_rad):
        self.curr_rad = new_curr_rad
        self.reset_reffrm()

    # 입력 : 해당 천체의 공전궤도상 진근점이각(중심천체에서 바라본 근일점과의 각도)
    # 출력 : 중심천체와의 거리(km)
    @property
    def dist_km(self):
        return self.__conic(self.curr_rad)
    
    def __conic(self, ang):
        return self.a_km * (1 - self.e**2) / (1 +self.e * m.cos(ang))

    # 입력 : 해당 천체의 공전궤도상 진근점이각(중심천체에서 바라본 근일점과의 각도)
    # 출력 : 순간공전속도(km/s)
    @property
    def inst_speed_kms(self):
        mu_km=self.center_body.GM_kms
        return m.sqrt(mu_km*(2/self.dist_km-1/self.a_km))

    # t_diff step시간 후 각도변화(천체 공전궤도 및 위치 변화)
    def orbital_increment_angle_rad(self, tdelta):
        # 진근점 이각 변화 수치근사(euler's method)
        self.curr_rad+=tdelta*self.inst_speed_kms/self.dist_km
        self.curr_rad%=2*m.pi

        # 근점편각 변화
        self.periap_ang_rad+=tdelta*self.apsidal_prcs_speed
        self.periap_ang_rad%=2*m.pi
        
        # 승교점 경도 변화
        self.node_lon_rad+=tdelta*self.nodal_prcs_speed
        self.node_lon_rad%=2*m.pi

        # self.reset_reffrm() # 기준틀 업데이트
        # self.reset_coor()

    def reset_reffrm(self):
        if self.center_body_orbit:
            self.reffrm.reset_reffrm(self.node_lon_rad, self.incln_rad, self.periap_ang_rad, self.name,
                                    *self.center_body_orbit.base_coor.tuplify(0))
        else: self.reffrm.reset_reffrm(self.node_lon_rad, self.incln_rad, self.periap_ang_rad, self.name) 
    
    def reset_coor(self):
        self.home_coor.reset_coor('s', (self.dist_km, self.curr_rad, m.pi/2))
        self.home_coor.conv_coor_modeOS()
        self.base_coor=self.reffrm.base_conv(self.home_coor, 'sa')

    # 천체의 공전주기
    def orbital_period_s(self):
        mu_km=self.center_body.GM_kms
        return 2*m.pi*m.sqrt((self.a_km**3)/mu_km)

    # orbit data 딕셔너러 출력
    def orb_data_dict(self):
        return {
            "center_body": self.center_body.name,
            "e": self.e,
            "a_km": self.a_km,
            "node_lon_rad": self.node_lon_rad,
            "incln_rad": self.incln_rad,
            "periap_ang_rad": self.periap_ang_rad,
            "apsidal_prcs_speed": self.apsidal_prcs_speed,
            "nodal_prcs": self.nodal_prcs_speed,
            "anom_true_rad": self.curr_rad
        }

    def __str__(self):
        return ("\n\norbit characteristics : \n- eccentricity : " + str(self.e) + "\n- semi major axis[km] : " + str(self.a_km) 
        + "\n- orientational parameters[deg] : (" + str(self.node_lon_rad*180/m.pi) + ", " + str(self.incln_rad*180/m.pi) + ", " + str(self.periap_ang_rad*180/m.pi) + ")"
        + "\n- precession speeds[deg/s] : (apsidal:" + str(self.apsidal_prcs_speed*180/m.pi) + ", nodal:" + str(self.nodal_prcs_speed*180/m.pi) + ")"
        + "\n\nstatus :"
        + "\n- angle[deg] : " + str(self.curr_rad*180/m.pi) + "\n- coordinate info : " + str(self.base_coor)+'\n- referance frame : '+str(self.reffrm)).replace('\n', '\n  ')


    # makes list of X, Y, Z coordinates of the entire elliptic orbit
    def orbital_coors(self, ratio=1):
        th=np.linspace(0, 2* m.pi, 500)
        coor=[(self.__conic(ang)*ratio, ang, m.pi/2) for ang in list(th)]
        coor=cr.coor3('s', *coor)
        coor.conv_coor_modeOS()
        coor=self.reffrm.base_conv(coor, 'sa')
        X, Y, Z = tuple([list(i) for i in coor.vec])
        return X, Y, Z
    
    def plot(self, ax, ratio=1, sphere=False):
        # plot orbit(ellipse)
        X, Y, Z = self.orbital_coors(ratio)
        op, = ax.plot(X, Y, Z, color="black", linewidth=0.5)

        # 중심천체--근일점 선분
        cntr=self.center_body_orbit.base_coor.vec if self.center_body_orbit else [0, 0, 0]
        pl, = ax.plot(*np.array([[cntr[0], cntr[1], cntr[2]],[X[0], Y[0], Z[0]]]).transpose(), color="black", linewidth=0.5)

        # plot celestial body
        bp=self.body.plot(ax, (self.base_coor.vec-np.array(cntr))*ratio+np.array(cntr), sphere=sphere)

        return X, Y, Z, [pl, op, bp]
    
    def anim_update(self, plots, ratio=1):
        X, Y, Z = self.orbital_coors(ratio)
        cntr=self.center_body_orbit.base_coor.vec if self.center_body_orbit else [0, 0, 0]
        set_data(plots[0], np.array([[cntr[0], cntr[1], cntr[2]],[X[0], Y[0], Z[0]]]).transpose())

        set_data(plots[1], np.array([X, Y, Z]))

        set_data(plots[2], (self.base_coor.vec-np.array(cntr))*ratio+np.array(cntr))

        return X, Y, Z

###############################################################
#                     solar system model                      #
###############################################################
# 항성계 구성모델 & 시간경과에 따른 운동 모델구동 class
#   - 방향성 : 천체크기 데이터가 포함된 태양-지구-달(SEM) 모델 구축
#   - 우선목표 : 일식 예측

# 모델 특성
#   - 항성은 움직이지 않는다(자전도 X)
#   - 항성-행성 / 행성-위성 등 케플러궤도 독립적 이체운동 모델 (정해진 주기성 이용: 알려진 공전/자전주기)
#   - 공전궤도, 자전축등 천체운동의 상대적 위치 데이터 구성
#   - 천체는 구형으로 가정

# 초기 config data 구성방법론
#   - 일식 데이터 : SEM system ruler
#   - 일식 데이터(시각/날짜)  ->  날짜: 지구 공전궤도상 위치, 일식 : 달 궤도교점, 달 공전궤도상 위치
#   - 일식모델(기하적으로) 만들기

# 모델 작동 방식
#
# 각 iteration별 작업
# 1.

# class cel_evaluation_set:

class ssm_element:
    def __init__(self, elm_data_dict, cntr_elm=None):
        #print("\nLoading {}...".format(elm_data_dict['name']))
        self.name=elm_data_dict["name"]
        body = ["name", "mass_kg", "diameter_km"]
        if self.name=="Sun": color='red'
        elif self.name=="Earth": color='blue'
        else: color='grey'

        self.body=cel_body(*[elm_data_dict[i] for i in body], color=color)

        if self.body.name=='Earth':   # 지구의 초기위치 매개변수값의 결정은 json에서 주어진 근일점 시각에서부터 일식 시각까지 수치해석하여 구한다.
            #지구 공전
            orb=["e", "a_km", "node_lon_rad", "incln_rad", "periap_ang_rad"]
            self.orb=orbit(*(elm_data_dict["orbit"][i] for i in orb), 0, 0, 0, self.body, cntr_elm.body) # 근일점에서의 값들로 설정
            
            # 지구 자전: 근일점에 있을 때 춘분점(승교점)과 원점(위도0경도0, 근일점 지나는 시각으로 계산) 사이 이각으로 curr_rad 입력
            rot=["eqnx_lon_rad", "incln_rad", "ang_v_rad"]
            self.rot=rotate(*(elm_data_dict["rotate"][i] for i in rot), (self.orb.periap_ang_rad+m.pi)%(2*m.pi), self.body, self.orb)
            self.elm_timetravel(elm_data_dict["perihelion"], elm_data_dict["Eclipse"])  # 근일점으로 설정된 지구를 일식이 일어나는 지점으로 옮김

        else:   # 지구 외 태양, 달 또는 다른 행성/위성
            # 공전
            if 'orbit' in elm_data_dict:
                orb=["e", "a_km", "node_lon_rad", "incln_rad", "periap_ang_rad", "apsidal_prcs", "nodal_prcs", "anom_true_rad"]
                
                # 일식 때는 지구-달-태양이 모두 일직선에 있으므로, 지구에 대해 계산한 값을 달의 현재 승교점편각에 그대로 사용할 수 있다.
                if elm_data_dict["orbit"]["node_lon_rad"]=='None' : 
                    elm_data_dict["orbit"]["node_lon_rad"]=(cntr_elm.orb.curr_rad+m.pi)%(2*m.pi)

                # 일식이므로 현재 달의 진근점이각이 승교점에 있어야 한다.
                elm_data_dict["orbit"]["anom_true_rad"]=-elm_data_dict["orbit"]["periap_ang_rad"]
                    
                # 각 공전궤도요소의 세차 주기를 이용해 등속이란 가정 하에 각속도 계산(rad/s)
                prcs=["apsidal_prcs", "nodal_prcs"]
                elm_data_dict["orbit"][prcs[0]], elm_data_dict["orbit"][prcs[1]]=(2*m.pi/datetime.timedelta(days=elm_data_dict["orbit"][p]).total_seconds() for p in prcs)
                self.orb=orbit(*(elm_data_dict["orbit"][i] for i in orb), self.body, cntr_elm.body, cntr_elm.orb)
            else: self.orb=None

            #자전
            if 'rotate' in elm_data_dict:
                rot=["eqnx_lon_rad", "incln_rad", "ang_v_rad", "rot_ang_rad"]
                self.rot=rotate(*(elm_data_dict["rotate"][i] for i in rot), self.body, self.orbit)
            else: self.rot = None

    # 천체를 tdelta 후의 위치로 이동
    def next_iter(self, tdelta):
        if self.orb: self.orb.orbital_increment_angle_rad(tdelta)
        if self.rot: self.rot.rotational_increment_angle_rad(tdelta)

    # 천체에 할당된 위치와 공전/자전 좌표계 재설정
    def reset_reffrm_coor(self):
        if self.orb:
            self.orb.reset_reffrm()
            self.orb.reset_coor()
        if self.rot:
            self.rot.reset_reffrm()

    # 천체를 특정 시각에서의 위치로 옮김
    def elm_timetravel(self, curr_time, target_time, tdelta=None):
        tdiff=(parse(target_time) - parse(curr_time)).total_seconds()    # 현재와 목표 시간 사이 초
        if not tdiff: return
        if tdelta:  # 사용자가 tdelta를 입력했을 시
            tdelta=tdelta if tdiff>=0 else -tdelta     # 목표시간이 현재보다 과거일 시 tdelta에 - 붙임

            while True:
                self.next_iter(tdelta)
                tdiff_prev=tdiff
                tdiff-=tdelta
                if tdiff*tdiff_prev<=0: break
        else:
            for td in [3600*24*365.25, 3600*24, 3600, 60, 1]:    # 시간절약을 위해 1년 전까지는 1년단위, 하루전 까지는 하루단위, 1시간 전까지는 1시간단위, 1초단위로 나누어 계산
                if td>abs(tdiff): continue
                td=td if tdiff>=0 else -td
                while True:
                    self.next_iter(td)
                    tdiff_prev=tdiff
                    tdiff-=td
                    if (tdiff-td)*(tdiff_prev-td)<=0: break
        self.reset_reffrm_coor()

    # compile을 위한 각 reffrm 및 base_coor의 리스트화 : 미완성
    def list_reffrm_coor(self):
        return self.rot.reffrm+self.orb.reffrm+self.orb.base_coor

    # element data 딕셔너리 출력
    def element_data_dict(self):
        datas=self.body.body_data_dict()
        if self.orb: datas["orbit"]=self.orb.orb_data_dict()
        if self.rot: datas["rotate"]=self.rot.rot_data_dict()
        return datas 
            
    def __str__(self):
        return ("┌──"
                +str(self.body) + 
                "\n- orbit info : " + 
                str(self.orb) + 
                "\n- rotation info : " + 
                str(self.rot) + '\n').replace('\n', '\n│  ')+ \
                "\n└──"
                
        


class solar_system_model:
    def __init__(self, name, date_time):
        self.name = name
        self.date_time=date_time
        self.elements=[]
        self.elm_dict = dict()
        

###############file load/export & model_obj constructor###############

    @classmethod
    def config_JSON_load(cls, file_name, print_info=True):
        with open(file_name, "r") as file_json:
            # load json & make new obj
            #print("Loading json...\n", end="\r")
            config = json.load(file_json)

            #print("configuring {}...".format(config["model_name"]), end="\r")
            new_ssm=cls(config["model_name"], parse(config["elements"][1][config["Event"]]))

            # make solar system individual celestial element obj
            for elm_data in config['elements']: new_ssm.add_element(elm_data)
            
            if print_info : print("\nLoading Complete\n", new_ssm)
            return new_ssm


    def add_element(self, elm_data):
        if 'orbit' in elm_data:
            cntr=self.elm_dict[elm_data['orbit']['center_body']]
            elm=ssm_element(elm_data, cntr)
        else: elm=ssm_element(elm_data)
        self.elements.append(elm)
        self.elm_dict[elm.name]=elm


    def config_JSON_dump(self, file_name):
        with open(file_name, "w") as file_json:
            json.dumps(self.ssm_data_dict, file_json)

    @classmethod
    def config_obj_load(cls, file_name): return joblib.load(file_name)

    def config_obj_write(self, file_name): joblib.dump(self, file_name)

    # solar system model data 딕셔너리 출력
    def ssm_data_dict(self):
        data=self.elm_dict
        for k, v in data.values(): data[k]=v.element_data_dict()
        data["model_name"]=self.name
        data["datetime"]=self.date_time
        return data

#################### string & visualize #######################
    def __str__(self):
        rpr='\n\n   [[ model information ]]\n'.replace('\n', '\n          ') +"""
      ___           ___           ___     
     /\__\         /\__\         /\  \    
    /:/ _/_       /:/ _/_       |::\  \   
   /:/ /\  \     /:/ /\  \      |:|:\  \  
  /:/ /::\  \   /:/ /::\  \   __|:|\:\  \ 
 /:/_/:/\:\__\ /:/_/:/\:\__\ /::::|_\:\__\ 
 \:\/:/ /:/  / \:\/:/ /:/  / \:\~~\  \/__/
  \::/ /:/  /   \::/ /:/  /   \:\  \      
   \/_/:/  /     \/_/:/  /     \:\  \     
     /:/  /        /:/  /       \:\__\    
     \/__/         \/__/         \/__/          
        """.replace('\n', '\n    ')
        rpr+=('\nmodel name : ' + self.name + "\ndatetime : " + self.date_time.strftime('%Y-%m-%dT%H:%M:%S.0Z')).replace('\n', '\n          ') + "\n\n Celestial body info : \n"
        for elm in self.elements: rpr+=str(elm)+'\n\n'
        return rpr

    @classmethod
    def plot_setlim(cls, ax, X, Y, Z=None, ratio=1):
        if Z:
            # X, Y, Z 축 사이 비율 교정(없으면 서로 길이 비율이 달라져 일그러짐)
            X, Y, Z = map(np.array, (X, Y, Z))
            max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max() / 2.0 * ratio

            mid_x = (X.max()+X.min()) * 0.5
            mid_y = (Y.max()+Y.min()) * 0.5
            mid_z = (Z.max()+Z.min()) * 0.5
            ax.set_xlim(mid_x - max_range, mid_x + max_range)
            ax.set_ylim(mid_y - max_range, mid_y + max_range)
            ax.set_zlim(mid_z - max_range, mid_z + max_range)
        else:
            X, Y = map(np.array, (X, Y))
            max_range = np.array([X.max()-X.min(), Y.max()-Y.min()]).max() / 2.0 * ratio

            mid_x = (X.max()+X.min()) * 0.5
            mid_y = (Y.max()+Y.min()) * 0.5
            ax.set_xlim(mid_x - max_range, mid_x + max_range)
            ax.set_ylim(mid_y - max_range, mid_y + max_range)

    def plot(self, ax, frame=0, sphere=False):

        ax.set_title("Orbit", fontdict={'fontsize': 16,'fontweight': 'bold'}, loc='left', pad=20)
        ax.set_title(self.date_time.strftime('%Y-%m-%dT%H:%M:%S.0Z'), loc='right', pad=20)

        self.elements[0].body.plot(ax, sphere=sphere)

        body_plt=list()
        X, Y, Z =list(), list(), list()
        for n, e in enumerate(self.elements[1:]):
            Xe, Ye, Ze, bp = e.orb.plot(ax, ratio=(100 if n+1==2 else 1), sphere=sphere)
            if frame==n+1: X, Y, Z = Xe, Ye, Ze
            else : X, Y, Z = X+Xe, Y+Ye, Z+Ze
            body_plt.append(bp)
        
        self.plot_setlim(ax, X, Y, Z)
        ax.legend()
        
        return body_plt

    def __orb_anim_update(self, i, body_plts, tdelta, ax, frame=0):

        ax.set_title("Orbit", fontdict={'fontsize': 16,'fontweight': 'bold'}, loc='left', pad=20)
        ax.set_title(self.date_time.strftime('%Y-%m-%dT%H:%M:%S.0Z'), fontdict={'fontsize': 9}, loc='right', pad=20)

        self.next_iter(tdelta)
        X, Y, Z = list(), list(), list()
        for n, elm in enumerate(self.elements[1:]):
            Xe, Ye, Ze = elm.orb.anim_update(body_plts[n], ratio=(100 if n+1==2 else 1))
            if frame==n+1: X, Y, Z = Xe, Ye, Ze
            else: X, Y, Z = X+Xe, Y+Ye, Z+Ze

        self.plot_setlim(ax, X, Y, Z)
        


    def orbit_animation(self, fig, ax, tdelta, frame=0):
        body_plts=self.plot(ax, frame=frame)
        ani = FuncAnimation(fig, self.__orb_anim_update, 100, fargs=(body_plts, tdelta, ax, frame), interval=1, blit=False)
        plt.show()

    # 천구 상 특정 천체의 원모양 plot coordinates
    def __cg_elem_coor(self, e, cntr_ref, local_ref=None):
            ### 천구상 천체를 원으로 표시 ###
            # 천체의 시직경
            ang_d=self.__ang_d(e)

            # 원 생성(numpy 이용) : 일단 천정쪽에 생성
            coor=cr.coor3('s', *[(1, ang, ang_d) for ang in np.linspace(0, 2* m.pi, 500)])

            # 좌표변환을 통해 천정상에 있는 원을 원래 있을 위치로 이동
            if e.orb: disp=cntr_ref.base_conv(e.orb.base_coor, 'as')
            else : disp=cntr_ref.base_conv(cr.coor3('o', (0, 0, 0)), 'as')
            if local_ref: disp=local_ref.base_conv(disp, 'as')
            disp.normalize()
            disp.conv_coor_modeOS()
            
            return (cr.reffrms((disp.vec[1]+m.pi/2)%(2*m.pi), disp.vec[2])).base_conv(coor.conv_coor_modeOS(), 'sa')

    #천구 plot
    def cheongu_plot(self, ax):

        ax.set_title("Cheongu", fontdict={'fontsize': 16,'fontweight': 'bold'}, loc='left', pad=20)
        ax.set_title(self.date_time.strftime('%Y-%m-%dT%H:%M:%S.0Z'), fontdict={'fontsize': 9}, loc='right', pad=20)

        # 적도 원
        coor=cr.coor3('s', *[(1, ang, m.pi/2) for ang in np.linspace(0, 2* m.pi, 500)])
        ax.plot(*tuple([list(i) for i in coor.conv_coor_modeOS().vec]), color='midnightblue', linewidth=0.5)

        # 위도/경도 0,0 선
        ax.plot((0, 1), (0, 0), (0, 0), color='midnightblue', linewidth=0.5)

        # 관측자가 서 있는 천체의 reffrm obj
        cntr_ref=self.elements[1].rot.reffrm

        # 태양 외 천체 천구에 표시
        plts=list()
        for e in [self.elements[0]]+self.elements[2:]:
            coor=self.__cg_elem_coor(e, cntr_ref)
            cel_p, = ax.plot(*[i for i in coor.vec], color=e.body.color)
            plts.append(cel_p)

        # X, Y, Z 축 사이 비율 교정(없으면 서로 길이 비율이 달라져 일그러짐)
        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
        ax.set_zlim(-1, 1)     


        return plts 

    def __cg_anim_update(self, i, plts, tdelta, ax):
        self.next_iter(tdelta)

        ax.set_title("Cheongu", fontdict={'fontsize': 16,'fontweight': 'bold'}, loc='left', pad=20)
        ax.set_title(self.date_time.strftime('%Y-%m-%dT%H:%M:%S.0Z'), fontdict={'fontsize': 9}, loc='right', pad=20)

        # 관측자가 서 있는 천체의 reffrm obj
        cntr_ref=self.elements[1].rot.reffrm

        # 태양 외 천체 천구에 표시
        for n, e in enumerate([self.elements[0]]+self.elements[2:]):
            coor = self.__cg_elem_coor(e, cntr_ref)
            plts[n].set_data(coor.vec[0], coor.vec[1])
            plts[n].set_3d_properties(coor.vec[2])
        

    def cheongu_animation(self, fig, ax, tdelta):
        plts=self.cheongu_plot(ax)
        ani = FuncAnimation(fig, self.__cg_anim_update, 100, fargs=(plts, tdelta, ax), interval=1, blit=False)
        plt.show()


    # 천구 2D Plot
    def cheongu_plot2d(self, ax, focus=1):

        ax.set_title("Cheongu2D", fontdict={'fontsize': 16,'fontweight': 'bold'}, loc='left', pad=20)
        ax.set_title(self.date_time.strftime('%Y-%m-%dT%H:%M:%S.0Z'), fontdict={'fontsize': 9}, loc='right', pad=20)

        # 적도선 및 원점(위도/경도=0,0)
        ax.plot(np.linspace(-180, 180, 500), [0 for i in range(500)], color='midnightblue', linewidth=0.5)
        ax.scatter(0, 0, color='midnightblue', s=2)

        # 관측자가 서 있는 천체의 reffrm obj
        cntr_ref=self.elements[1].rot.reffrm

        # 태양 외 천체 천구에 표시
        focus = -1 if focus==1 else focus
        focus-= 1 if focus else 0
        plts=list()
        for n, e in enumerate([self.elements[0]]+self.elements[2:]):
            coor = self.__cg_elem_coor(e, cntr_ref).conv_coor_modeOS()
            cel_p, = ax.plot((-coor.vec[1]*180/m.pi+360)%360-180, 90-coor.vec[2]*180/m.pi, color=e.body.color)
            plts.append(cel_p)
            
            if n==focus: self.plot_setlim(ax, (-coor.vec[1]*180/m.pi+360)%360-180, 90-coor.vec[2]*180/m.pi, ratio=4)
        
        if focus==-2:
            ax.set_xlim(-180, 180)
            ax.set_ylim(-90, 90)


        return plts

    def __cg2d_anim_update(self, i, plts, tdelta, ax, focus=1):
        self.next_iter(tdelta)

        ax.set_title("Cheongu2D", fontdict={'fontsize': 16,'fontweight': 'bold'}, loc='left', pad=20)
        ax.set_title(self.date_time.strftime('%Y-%m-%dT%H:%M:%S.0Z'), fontdict={'fontsize': 9}, loc='right', pad=20)

        # 관측자가 서 있는 천체의 reffrm obj
        cntr_ref=self.elements[1].rot.reffrm

        # 태양 외 천체 천구에 표시
        focus = -1 if focus==1 else focus
        focus-= 1 if focus else 0
        for n, e in enumerate([self.elements[0]]+self.elements[2:]):
            coor=self.__cg_elem_coor(e, cntr_ref).conv_coor_modeOS()
            plts[n].set_data((-coor.vec[1]*180/m.pi+360)%360-180, 90-coor.vec[2]*180/m.pi)

            if n==focus: self.plot_setlim(ax, (-coor.vec[1]*180/m.pi+360)%360-180, 90-coor.vec[2]*180/m.pi, ratio=4)
        
        if focus==-2:
            ax.set_xlim(-180, 180)
            ax.set_ylim(-90, 90)
    

    def cheongu_animation2d(self, fig, ax, tdelta, focus=1):
        plts=self.cheongu_plot2d(ax, focus)
        ani = FuncAnimation(fig, self.__cg2d_anim_update, 100, fargs=(plts, tdelta, ax, focus), interval=1, blit=False)
        plt.show()


    def local_datetime_str(self, lon=0, lat=0):
        str_tz=TimezoneFinder().timezone_at(lng=lon, lat=lat)
        local_tz = pytz.timezone(str_tz)
        local_dt = self.date_time.replace(tzinfo=pytz.utc).astimezone(local_tz)
        return str_tz+" "+local_tz.normalize(local_dt).strftime('%Y-%m-%dT%H:%M:%S.0%z')
    
    #천구 plot
    def local_cheongu_plot(self, ax, lon, lat):
        ax.set_title("Local Cheongu", fontdict={'fontsize': 16,'fontweight': 'bold'}, loc='left', pad=20)
        ax.set_title(self.local_datetime_str(lon, lat), fontdict={'fontsize': 9}, loc='right', pad=20)

        # 지표면
        coor=cr.coor3('s', *[(1, ang, m.pi/2) for ang in np.linspace(0, 2* m.pi, 500)])
        ax.plot(*tuple([list(i) for i in coor.conv_coor_modeOS().vec]), color='midnightblue', linewidth=0.5)

        # 위도/경도 0,0 선
        ax.plot((0, 1), (0, 0), (0, 0), color='midnightblue', linewidth=0.5)

        # 관측자가 서 있는 천체의 reffrm obj
        cntr_ref=self.elements[1].rot.reffrm
        
        # 관측자가 천체 지표면상 서 있는 위치에서의 reffrm obj
        local_ref=self.localize_reffrm(lon, lat)

        # 적도 원
        lon, lat = lon*m.pi/180, (-lat*m.pi/180)+m.pi/2
        coor=cr.coor3('s', *[(1, ang, m.pi/2) for ang in np.linspace(0, 2* m.pi, 500)])
        coor=cr.reffrms(m.pi/2+lon, lat, m.pi/2).base_conv(coor.conv_coor_modeOS(), 'as')
        ax.plot(*tuple([list(i) for i in coor.vec]), color='midnightblue', linewidth=0.5)

        # 태양 외 천체 천구에 표시
        plts=list()
        for e in [self.elements[0]]+self.elements[2:]:
            coor=self.__cg_elem_coor(e, cntr_ref, local_ref)
            cel_p, = ax.plot(*[i for i in coor.vec], color=e.body.color)
            plts.append(cel_p)

        # X, Y, Z 축 사이 비율 교정(없으면 서로 길이 비율이 달라져 일그러짐)
        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
        ax.set_zlim(-1, 1)     

        return plts 

    def __lcg_anim_update(self, i, plts, tdelta, ax, lon, lat):
        self.next_iter(tdelta)

        ax.set_title("Local Cheongu", fontdict={'fontsize': 16,'fontweight': 'bold'}, loc='left', pad=20)
        ax.set_title(self.local_datetime_str(lon, lat), fontdict={'fontsize': 9}, loc='right', pad=20)

        # 관측자가 서 있는 천체의 reffrm obj
        cntr_ref=self.elements[1].rot.reffrm

        # 관측자가 천체 지표면상 서 있는 위치에서의 reffrm obj
        local_ref=self.localize_reffrm(lon, lat)

        # 태양 외 천체 천구에 표시
        for n, e in enumerate([self.elements[0]]+self.elements[2:]):
            coor=self.__cg_elem_coor(e, cntr_ref, local_ref)
            plts[n].set_data(coor.vec[0], coor.vec[1])
            plts[n].set_3d_properties(coor.vec[2])
        

    def local_cheongu_animation(self, fig, ax, lon, lat, tdelta):
        plts=self.local_cheongu_plot(ax, lon, lat)
        ani = FuncAnimation(fig, self.__lcg_anim_update, 100, fargs=(plts, tdelta, ax, lon, lat), interval=1, blit=False)
        plt.show()

    # 천구 2D Plot
    def local_cheongu_plot2d(self, ax, lon, lat, focus=1):
        ax.set_title("Local Cheongu2D", fontdict={'fontsize': 16,'fontweight': 'bold'}, loc='left', pad=20)
        ax.set_title(self.local_datetime_str(lon, lat), fontdict={'fontsize': 9}, loc='right', pad=20)

        # 적도선 및 원점(위도/경도=0,0)
        ax.plot(np.linspace(-180, 180, 500), [0 for i in range(500)], color='midnightblue', linewidth=0.5)
        ax.scatter(0, 0, color='midnightblue', s=2)

        # 관측자가 서 있는 천체의 reffrm obj
        cntr_ref=self.elements[1].rot.reffrm

        # 관측자가 천체 지표면상 서 있는 위치에서의 reffrm obj
        local_ref=self.localize_reffrm(lon, lat)

        # 태양 외 천체 천구에 표시
        focus = -1 if focus==1 else focus
        focus-= 1 if focus else 0
        plts=list()
        for n, e in enumerate([self.elements[0]]+self.elements[2:]):
            coor = self.__cg_elem_coor(e, cntr_ref, local_ref).conv_coor_modeOS()
            cel_p, = ax.plot((-coor.vec[1]*180/m.pi+360)%360-180, 90-coor.vec[2]*180/m.pi, color=e.body.color)
            plts.append(cel_p)
            
            if n==focus: self.plot_setlim(ax, (-coor.vec[1]*180/m.pi+360)%360-180, 90-coor.vec[2]*180/m.pi, ratio=4)
        
        if focus==-2:
            ax.set_xlim(-180, 180)
            ax.set_ylim(-90, 90)

        return plts

    def __lcg2d_anim_update(self, i, plts, tdelta, ax, lon, lat, focus=1):
        self.next_iter(tdelta)

        ax.set_title("Local Cheongu2D", fontdict={'fontsize': 16,'fontweight': 'bold'}, loc='left', pad=20)
        ax.set_title(self.local_datetime_str(lon, lat), fontdict={'fontsize': 9}, loc='right', pad=20)

        # 관측자가 서 있는 천체의 reffrm obj
        cntr_ref=self.elements[1].rot.reffrm

        # 관측자가 천체 지표면상 서 있는 위치에서의 reffrm obj
        local_ref=self.localize_reffrm(lon, lat)

        # 태양 외 천체 천구에 표시
        focus = -1 if focus==1 else focus
        focus-= 1 if focus else 0
        for n, e in enumerate([self.elements[0]]+self.elements[2:]):
            coor=self.__cg_elem_coor(e, cntr_ref, local_ref).conv_coor_modeOS()
            plts[n].set_data((-coor.vec[1]*180/m.pi+360)%360-180, 90-coor.vec[2]*180/m.pi)

            if n==focus: self.plot_setlim(ax, (-coor.vec[1]*180/m.pi+360)%360-180, 90-coor.vec[2]*180/m.pi, ratio=4)
        
        if focus==-2:
            ax.set_xlim(-180, 180)
            ax.set_ylim(-90, 90)


    def local_cheongu_animation2d(self, fig, ax, lon, lat, tdelta, focus=1):
        plts=self.local_cheongu_plot2d(ax, lon, lat, focus)
        ani = FuncAnimation(fig, self.__lcg2d_anim_update, 100, fargs=(plts, tdelta, ax, lon, lat, focus), interval=1, blit=False)
        plt.show()

    # def all_animation(self, fig, ax, ax2, tdelta, frame=0):
    #     body_plts=self.plot(ax, frame=frame)
    #     ani = FuncAnimation(fig, self.__orb_anim_update, 100, fargs=(body_plts, tdelta, ax, frame), interval=1, blit=False)

    #     plts=self.cheongu_plot(ax)
    #     ani = FuncAnimation(fig, self.__cg_anim_update, 100, fargs=(plts, tdelta), interval=1, blit=False)

    #     plts2=self.cheongu_plot2d(ax2)
    #     ani = FuncAnimation(fig, self.__cg2d_anim_update, 100, fargs=(plts2, tdelta), interval=1, blit=False)
        
    #     plt.show()


############### model compile tools ###################

    # def model_compile(self, cls, *snap_times):
    #     compiled = np.vstack(compiled, cls.__mono_iterational_compile())
    #     return solsys_compile_data(model_compile_code, cel_name, *snap_times, snap_count, compiled)

    # tdelta 뒤 위치 추정
    def next_iter(self, tdelta):
        for e in self.elements[1:]:
            e.next_iter(tdelta)
            e.reset_reffrm_coor()
            #print(e)
        self.date_time+=datetime.timedelta(seconds=tdelta)

    # 천체를 특정 시각에서의 위치로 옮김
    def timetravel(self, target_time, tdelta=None):
        tdiff=(parse(target_time) - self.date_time).total_seconds()    # 현재와 목표 시간 사이 초
        if not tdiff: return
        if tdelta:  # 사용자가 tdelta를 입력했을 시
            tdelta=tdelta if tdiff>=0 else -tdelta     # 목표시간이 현재보다 과거일 시 tdelta에 - 붙임
            while True:
                self.next_iter(tdelta)
                tdiff_prev=tdiff
                tdiff-=tdelta
                if tdiff*tdiff_prev<=0: break
                print(self.date_time, end="\r")
        else:
            for td in [3600*24*365.25, 3600*24, 3600, 60, 1]:    # 시간절약을 위해 1년 전까지는 1년단위, 하루전 까지는 하루단위, 1시간 전까지는 1시간단위, 1초단위로 나누어 계산
                if td>abs(tdiff): continue
                td=td if tdiff>=0 else -td
                while True:
                    self.next_iter(td)
                    tdiff_prev=tdiff
                    tdiff-=td
                    if (tdiff-td)*(tdiff_prev-td)<=0: break
                    print(self.date_time, end="\r")
        print("Time travel terminated : Arrived at ", self.date_time)

    def timetravel_diff(self, d=0, h=0, M=0, s=0, tdelta=None):
        diff=datetime.timedelta(days=d, hours=h, minutes=M, seconds=s)
        self.timetravel((self.date_time+diff).strftime('%Y-%m-%dT%H:%M:%S.0Z'), tdelta)

    
    def __ang_d(self, target):
        ref_c, target_c = self.elements[1], target
        if target_c.orb: dist = ref_c.orb.base_coor/ target_c.orb.base_coor
        else : dist = abs(ref_c.orb.base_coor)
        return m.atan(target_c.body.diameter_km / dist)

    # dist 거리에 있는 천체의 시직경
    def angular_diameter(self, target_cel):
        return self.__ang_d(self.elm_dict[target_cel])

    # 지구중심 -> 관측자가 천체 지표면상 서 있는 위치에서의 reffrm obj
    def localize_reffrm(self, lon, lat):
        lon, lat = lon*m.pi/180, (-lat*m.pi/180)+m.pi/2
        coor = cr.coor3('s', (self.elements[1].body.diameter_km, lon, lat))
        coor=coor.conv_coor_modeOS().tuplify(0)
        return cr.reffrms(m.pi/2+lon, lat, m.pi/2, "reffrm", *coor)

############## 일식분석툴 ###############

    # 현제 위치에서 거리관계로 계산한 일식경계이각, eclipse boundary elongation
    ## 지구에서 거리dist_s2e_km에 있는 달이 태양에서 오는 빛을 가리기 시작하는 지점(일식경계각)
    ## (지구중심과 태양중심을 이은 일직선 기준)(l_s2e)
    ## 입력 : 지구-달 거리, 태양-지구 거리
    ## 출력 : 지구중심에서 l_s2e을 기준으로 일식경계각
    def partial_eclps_elon_rad(self):
        rs, re, rm = (e.body.diameter_km/2 for e in self.elements)
        s2e_km=abs(self.elements[1].orb.base_coor)
        e2m_km=self.elements[1].orb.base_coor / self.elements[2].orb.base_coor
        return m.asin((rs-re)/ s2e_km) + m.asin((re+rm) / e2m_km)

    def full_eclps_elon_rad(self):
        rs, re, rm = (e.body.diameter_km/2 for e in self.elements)
        s2e_km=abs(self.elements[1].orb.base_coor)
        e2m_km=self.elements[1].orb.base_coor / self.elements[2].orb.base_coor
        return m.acos((rs+re)/ s2e_km) - m.acos((re+rm) / e2m_km)
    
    # 이각 : 지구에서 관측한 태양과 천체 사이 각도
    def elon_rad(self):
        er, mn= (e.orb.base_coor for e in self.elements[1:])
        return (mn-er)%(-er)

    # 현재 일식이 발생하고 있는지 확인
    def eclps_judge(self, mode=0): 
        if mode: return (self.elon_rad()<self.full_eclps_elon_rad())
        else: return (self.elon_rad()<self.partial_eclps_elon_rad())
    
    # 
    def search_eclipse(self, tdelta, mode=0, until=None):
        while(not self.eclps_judge(mode)):
            self.next_iter(tdelta)
            if (until and parse(until)<self.date_time): break
            print(self.date_time, end="\r")
        if self.eclps_judge(mode) : print("Search success : ", self.date_time.strftime('%Y-%m-%dT%H:%M:%S.0Z'))
        else : print("Search failed")
        return self.date_time

        




# class solsys_compile_data:
#     def __init__(self, model_compile_code, cel_name, t_ini, t_diff, t_fin,
#                  snap_count, compiled):
#         self.comp_code = model_compile_code
#         self.cel_name = cel_name  # tuple
#         self.snap_time = t_ini, t_diff, t_fin
#         self.snap_count = snap_count
#         self.compiled = compiled  # 2D numpy array
#         pass

#     def compdata_loadcsv(cls, file_name):
#         csv = np.loadtxt(file_name, delimiter=",")
#         celname = tuple(csv[0].tolist())
#         compiled = csv[1:, :]
#         with open(file_name, 'r') as f:
#             header = f.readline().replace(" ", "").replace("#", "").split(',')
#             comp_code, snap_time, snap_count = header[1], (
#                 header[2], float(header[3]), header[4]), header[5]
#         return cls(comp_code, celname, *snap_time, snap_count, compiled)

#     def compdata_writecsv(self, file_name):
#         data = np.vstack((np.array(self.cel_name), self.compiled))
#         headers = "{}, {}, {}, {}, {}".format(self.comp_code, *self.snap_time,
#                                               self.snap_count)
#         np.savetxt(file_name,
#                    data,
#                    fmt='%.18e',
#                    delimiter=',',
#                    newline='n',
#                    header=headers)

# ############### comp analysis tools ##################

#     def distance(self, cel1, cel2):
#         c1, c2 = self.cel_dict[cel1], self.cel_dict[cel2]


# ############## comp visual_sim tools #################

#     def run_compiled_data():
#         pass

