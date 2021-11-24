import math as m
import numpy as np
import json
import struct
import datetime
from dateutil.parser import parse
import coor_ref as cr


###############################################################
#                     solar system bodies                     #
###############################################################

class cel_body:
    G_kms = 6.67384 * 10**-20

    def __init__(self, name, mass_kg, diameter):
        print("Loading {} body info...".format(name))
        self.name = name
        self.mass_kg = mass_kg
        self.GM_kms = mass_kg * self.G_kms
        self.diameter_km = diameter
    
    def __str__(self):
        return ('\n\n{}\n- mass[kg] : {}\n- diameter[km] : {}\n'.format(self.name, self.mass_kg, self.diameter_km))
        


###############################################################
#                      celestial motion                       #
###############################################################
# 공전/자전과 같은 천체운동의 중요한 계수들과 위치 정보 관련 class


# 자전 관련 데이터
class rotate:
    def __init__(self, eqnx_lon_rad, incln_rad, ang_v_rad,rot_ang_rad, body, body_orbit=None):
        print("Loading {} rotation info...".format(body.name))
        self.eqnx_lon_rad = eqnx_lon_rad  # 춘분점 경도 : 기준면상 기준축과 적도승교점 사이 각도(자전축 방향)
        self.incln_rad = incln_rad  # 자전축 기울기

        self.ang_v_rad = ang_v_rad  # 자전각속도
        self.curr_rad = rot_ang_rad % (2 * m.pi)  # 현재 자전 위상각(중심기준은 적도좌표계)

        self.body=body  # 천체obj
        self.body_orbit=body_orbit    # 공전obj
        self.name=body.name+'_rot'  # 이름

        print("Calculating {} rotation refference frame matrix...".format(body.name))
        # reffrm
        self.reffrm=cr.reffrms(eqnx_lon_rad, incln_rad, self.curr_rad, self.name, *(body_orbit.home_coor.tuplify())) if body_orbit \
            else cr.reffrms(eqnx_lon_rad, incln_rad, self.curr_rad, self.name) 

    def set_curr_rad(self, new_curr_rad):
        self.curr_rad=new_curr_rad
        self.reset_reffrm()


    # t_diff step시간 후 각도변화(천체 자전)
    def rotational_increment_angle_rad(self, t_diff):
        self.curr_rad+=t_diff*self.ang_v_rad
        self.curr_rad%=2*m.pi

        self.reset_reffrm() # 기준틀 업데이트


    def reset_reffrm(self):
        if self.body_orbit:
            self.reffrm.reset_reffrm(self.eqnx_lon_rad, self.incln_rad, self.curr_rad, self.name, 
                                    self.body_orbit.dist_km, self.body_orbit.curr_rad, 0)
        else: self.reffrm.reset_reffrm(self.eqnx_lon_rad, self.incln_rad, self.curr_rad, self.name) 
    
    def __str__(self):
        return ("\n\n- refference frame :" + str(self.reffrm)).replace('\n', '\n  ')


# 공전 관련 데이터(케플러적 이체모델)
class orbit:
    def __init__(self, eccentricity, a_km, node_lon_rad, incln_rad, periap_ang_rad, anom_true_rad,
                 apsidal_prcs_speed, nodal_prcs_speed, body, center_body, center_body_orbit=None):
        print("Loading {} orbit info...".format(body.name))
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

        print("Calculating {} orbit refference frame matrix...".format(body.name))
        # reffrm
        self.reffrm=cr.reffrms(node_lon_rad, incln_rad, periap_ang_rad, self.name, *center_body_orbit.base_coor.tuplify()) \
             if center_body_orbit else cr.reffrms(node_lon_rad, incln_rad, periap_ang_rad, self.name)

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
        return self.a_km * (1 - self.e**2) / (1 +self.e * m.cos(self.curr_rad))

    # 입력 : 해당 천체의 공전궤도상 진근점이각(중심천체에서 바라본 근일점과의 각도)
    # 출력 : 순간공전속도(km/s)
    @property
    def inst_speed_kms(self):
        mu_km=self.center_body.GM_kms
        return m.sqrt(mu_km*(2/self.dist_km-1/self.a_km))

    # t_diff step시간 후 각도변화(천체 공전궤도 및 위치 변화)
    def orbital_increment_angle_rad(self, t_diff):
        # 진근점 이각 변화 수치근사(euler's method)
        self.curr_rad+=t_diff*self.inst_speed_kms/self.dist_km
        self.curr_rad%=2*m.pi

        # 근점편각 변화
        self.periap_ang_rad+=t_diff*self.apsidal_prcs_speed
        self.periap_ang_rad%=2*m.pi
        
        # 승교점 경도 변화
        self.node_lon_rad+=t_diff*self.nodal_prcs_speed
        self.node_lon_rad%=2*m.pi

        #self.reset_reffrm() # 기준틀 업데이트

    def reset_reffrm(self):
        if self.center_body_orbit:
            self.reffrm.reset_reffrm(self.node_lon_rad, self.incln_rad, self.periap_ang_rad, self.name,
                                    *tuple(self.center_body_orbit.base_coor.vec)) 
        else: self.reffrm.reset_reffrm(self.node_lon_rad, self.incln_rad, self.periap_ang_rad, self.name) 
    
    def reset_coor(self):
        self.home_coor.reset_coor('s', (self.dist_km, self.curr_rad, m.pi/2))
        self.home_coor.conv_coor_modeOS()
        self.base_coor=self.reffrm.base_conv(self.home_coor, 'sa')

    # 천체의 공전주기
    def orbital_period_s(self):
        mu_km=self.center_body.GM_kms
        return 2*m.pi*m.sqrt((self.a_km**3)/mu_km)

    def __repr__(self):
        return ("\n\n- coordinate info : " + str(self.home_coor)+'\n\n- referance frame : '+str(self.reffrm)).replace('\n', '\n  ')

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
        print("\nLoading {}...".format(elm_data_dict['name']))
        self.name=elm_data_dict["name"]
        body = ["name", "mass_kg", "diameter_km"]
        self.body=cel_body(*[elm_data_dict[i] for i in body])

        if self.body.name=='Earth':   # 지구의 초기위치 매개변수값의 결정은 json에서 주어진 근일점 시각에서부터 일식 시각까지 수치해석하여 구한다.
            #지구 공전
            orb=["e", "a_km", "node_lon_rad", "incln_rad", "periap_ang_rad"]
            self.orb=orbit(*(elm_data_dict["orbit"][i] for i in orb), 0, 0, 0, self.body, cntr_elm.body) # 근일점에서의 값들로 설정
            
            # 지구 자전: 근일점에 있을 때 춘분점(승교점)과 원점(위도0경도0, 근일점 지나는 시각으로 계산) 사이 이각으로 curr_rad 입력
            rot=["eqnx_lon_rad", "incln_rad", "ang_v_rad"]
            self.rot=rotate(*(elm_data_dict["rotate"][i] for i in rot), (self.orb.periap_ang_rad+m.pi)%(2*m.pi), self.body, self.orb)
            self.timetravel(elm_data_dict["perihelion"], elm_data_dict["Eclipse"])  # 근일점으로 설정된 지구를 일식이 일어나는 지점으로 옮김

        else:   # 지구 외 태양, 달 또는 다른 행성/위성
            # 공전
            if 'orbit' in elm_data_dict:
                orb=["e", "a_km", "node_lon_rad", "incln_rad", "periap_ang_rad", "apsidal_prcs", "nodal_prcs", "anom_true_rad"]
                
                # 일식 때는 지구-달-태양이 모두 일직선에 있으므로, 지구에 대해 계산한 값을 달의 현재 승교점편각에 그대로 사용할 수 있다.
                if elm_data_dict["orbit"]["node_lon_rad"]=='None' : elm_data_dict["orbit"]["node_lon_rad"]=0
                elm_data_dict["orbit"]["nodal_lon_rad"]=(cntr_elm.orb.curr_rad+m.pi)%(2*m.pi)
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
    def increment_motion(self, tdelta):
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
    def timetravel(self, curr_time, target_time, tdelta=None):
        tdiff=(parse(target_time) - parse(curr_time)).total_seconds()    # 현재와 목표 시간 사이 초
        if not tdiff: return
        if tdelta:  # 사용자가 tdelta를 입력했을 시
            tdelta=tdelta if tdiff>=0 else -tdelta     # 목표시간이 현재보다 과거일 시 tdelta에 - 붙임

            while True:
                self.increment_motion(tdelta)
                tdiff_prev=tdiff
                tdiff-=tdelta
                if tdiff*tdiff_prev<=0: break
        else:
            for td in [3600*24*365.25, 3600*24, 3600, 60, 1]:    # 시간절약을 위해 1년 전까지는 1년단위, 하루전 까지는 하루단위, 1시간 전까지는 1시간단위, 1초단위로 나누어 계산
                if td>abs(tdiff): continue
                td=td if tdiff>=0 else -td
                while True:
                    self.increment_motion(td)
                    tdiff_prev=tdiff
                    tdiff-=td
                    if (tdiff-td)*(tdiff_prev-td)<=0: break
        self.reset_reffrm_coor()

    # compile을 위한 각 reffrm 및 base_coor의 리스트화 : 미완성
    def list_reffrm_coor(self):
        return self.rot.reffrm+self.orb.reffrm+self.orb.base_coor

            
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
    def config_JSON_load(cls, file_name):
        with open(file_name, "r") as file_json:
            # load json & make new obj
            print("Loading json...\n")
            config = json.load(file_json)

            print("configuring {}...".format(config["model_name"]))
            new_ssm=cls(config["model_name"], parse(config["elements"][1][config["Event"]]))

            # make solar system individual celestial element obj
            for elm_data in config['elements']: new_ssm.add_element(elm_data)
            
            print("\nLoading Complete\n", new_ssm)

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
            json.dumps(file_json)

    def config_struct_load(file_name):
        pass

    def config_struct_write(file_name):
        pass

    def __str__(self):
        rpr='\n\n   [[ model information ]]\n'.replace('\n', '\n          ') +"""
      ___           ___           ___     
     /\__\         /\__\         /\  \    
    /:/ _/_       /:/ _/_       |::\  \   
   /:/ /\  \     /:/ /\  \      |:|:\  \  
  /:/ /::\  \   /:/ /::\  \   __|:|\:\  \ 
 /:/_/:/\:\__\ /:/_/:/\:\__\ /::::|_\:\__\\
 \:\/:/ /:/  / \:\/:/ /:/  / \:\~~\  \/__/
  \::/ /:/  /   \::/ /:/  /   \:\  \      
   \/_/:/  /     \/_/:/  /     \:\  \     
     /:/  /        /:/  /       \:\__\    
     \/__/         \/__/         \/__/          
        """.replace('\n', '\n    ')
        rpr+=('\nmodel name : ' + self.name + "\ndatetime : " + self.date_time.strftime('%Y-%m-%dT%H:%M:%S.0Z')).replace('\n', '\n          ') + "\n\n Celestial body info : \n"
        for elm in self.elements: rpr+=str(elm)+'\n\n'
        return rpr

############### model compile tools ###################

    def model_compile(self, cls, *snap_times):
        compiled = np.vstack(compiled, cls.__mono_iterational_compile())
        return solsys_compile_data(model_compile_code, cel_name, *snap_times, snap_count, compiled)

    # dist 거리에 있는 천체의 시직경
    def angular_diameter(self, cel1, cel2):
        dist = self.distance(cel1, cel2)
        return m.atan(cel2.diameter_km / dist)

    def ellipse_boundary_angle(dist_e2m_km, dist_s2e_km, sun_data, earth_data,
                               moon_data):
        rs = sun_data.diameter_km / 2
        re = earth_data.diameter_km / 2
        rm = moon_data.diameter_km / 2
        rse = rs - re
        rem = re + rm
        return m.asin(rse / dist_s2e_km) + m.asin(rem / dist_e2m_km)



class solsys_compile_data:
    def __init__(self, model_compile_code, cel_name, t_ini, t_diff, t_fin,
                 snap_count, compiled):
        self.comp_code = model_compile_code
        self.cel_name = cel_name  # tuple
        self.snap_time = t_ini, t_diff, t_fin
        self.snap_count = snap_count
        self.compiled = compiled  # 2D numpy array
        pass

    def compdata_loadcsv(cls, file_name):
        csv = np.loadtxt(file_name, delimiter=",")
        celname = tuple(csv[0].tolist())
        compiled = csv[1:, :]
        with open(file_name, 'r') as f:
            header = f.readline().replace(" ", "").replace("#", "").split(',')
            comp_code, snap_time, snap_count = header[1], (
                header[2], float(header[3]), header[4]), header[5]
        return cls(comp_code, celname, *snap_time, snap_count, compiled)

    def compdata_writecsv(self, file_name):
        data = np.vstack((np.array(self.cel_name), self.compiled))
        headers = "{}, {}, {}, {}, {}".format(self.comp_code, *self.snap_time,
                                              self.snap_count)
        np.savetxt(file_name,
                   data,
                   fmt='%.18e',
                   delimiter=',',
                   newline='n',
                   header=headers)

############### comp analysis tools ##################

    def distance(self, cel1, cel2):
        c1, c2 = self.cel_dict[cel1], self.cel_dict[cel2]


############## comp visual_sim tools #################

    def run_compiled_data():
        pass


class solar_system_analysis_prediction_tools:
    def __init__(self, solar_system_motion_model):
        self.system_model = solar_system_motion_model

    # 일식 분석툴

    # dist 거리에 있는 천체의 시직경
    def angular_diameter(dist, sph_cel):
        return m.atan(sph_cel.diameter_km / dist)

    # 지구에서 거리dist_s2e_km에 있는 달이 태양에서 오는 빛을 가리기 시작하는 지점(일식경계각)
    # (지구중심과 태양중심을 이은 일직선 기준)(l_s2e)
    # 입력 : 지구-달 거리, 태양-지구 거리
    # 출력 : 지구중심에서 l_s2e을 기준으로 일식경계각
    def ellipse_boundary_angle(dist_e2m_km, dist_s2e_km, sun_data, earth_data, moon_data):
        rs = sun_data.diameter_km / 2
        re = earth_data.diameter_km / 2
        rm = moon_data.diameter_km / 2
        rse = rs - re
        rem = re + rm
        return m.asin(rse / dist_s2e_km) + m.asin(rem / dist_e2m_km)
