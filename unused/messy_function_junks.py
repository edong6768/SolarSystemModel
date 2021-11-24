import math as m



# 옮김

# 입력 : 해당 천체의 공전궤도상 진근점이각(중심천체에서 바라본 근일점과의 각도)
# 출력 : 중심천체와의 거리(km)
def orbit_ang2dist_km(angle, cel_data):
  a=cel_data['long_r_km']
  e=cel_data['eccentr']
  return a*(1-e**2)/(1+e*m.cos(angle))

# 해당 천체의 공전주기
def orbital_period_s(cel_data):
  mu_km=cel_data['orbit_center']['GM_m3s-2']/(1000**3)
  a=cel_data['long_r_km']
  return 2*m.pi*m.sqrt((a**3)/mu_km)

# 입력 : 해당 천체의 공전궤도상 진근점이각(중심천체에서 바라본 근일점과의 각도)
# 출력 : 순간공전속도(km/s)
def inst_orbit_speed(angle, cel_data):
  r = orbit_ang2dist_km(angle, cel_data)
  mu_km=cel_data['orbit_center']['GM_m3s-2']/(1000**3)
  a=cel_data['long_r_km']
  return m.sqrt(mu_km*(2/r-1/a))



# 지구에서 거리dist_s2e_km에 있는 달이 태양에서 오는 빛을 가리기 시작하는 지점(일식경계각)
# (지구중심과 태양중심을 이은 일직선 기준)(l_s2e)
# 입력 : 지구-달 거리, 태양-지구 거리
# 출력 : 지구중심에서 l_s2e을 기준으로 일식경계각
def ellipse_boundary_angle(dist_e2m_km, dist_s2e_km, earth_data, moon_data, sun_data):
  re=earth_data['diameter_km']/2
  rm=moon_data['diameter_km']/2
  rs=sun_data['diameter_km']/2
  rse=rs-re
  rem=re+rm
  return m.asin(rse/dist_s2e_km)+m.asin(rem/dist_e2m_km)

# 해당 천체의 각운동량
def spec_ang_momentum(cel_data):
  mu_km=cel_data['orbit_center']['GM_m3s-2']/(1000**3)
  a=cel_data['long_r_km']
  e=cel_data['eccentr']
  return m.sqrt(a*mu_km*(1-e**2))



# 지구에서 거리dist_earth_km에 있는 달(천체)가 태양에서 오는 빛을 가리기 시작하는 지점
# l_s2e : 태양과 지구(혹은 행성)을 일직선으로 이은 선분
# 입력 : 지구(행성)에서 달(천체)까지 거리
# 출력 : l_s2e 수직방향 반지름 거리
def radius_eclipse(dist_earth_km, dist_froms_tocel_km, earth_data, moon_data, sun_data):
  re=earth_data['diameter_km']/2
  rm=moon_data['diameter_km']/2
  rs=sun_data['diameter_km']/2
  rse=rs-re
  d=m.sqrt(dist_froms_tocel_km**2-rse**2)
  return ((rm+re)*dist_froms_tocel_km+rse*dist_earth_km)/d




# dist 거리에 있는 천체의 시직경
def angular_diameter(dist, cel_data):
  return m.atan(cel_data['diameter_km']/dist)

# 공전주기 수치적으로 구하기
def orbital_period_num_s(t_diff, cel_data):
  a=0   # 근일점에서 각도
  t_total_sec=0 #공전주기
  while(a<2*m.pi):
    a=orbital_increment_angle_rad(a, t_diff, cel_data)
    t_total_sec+=t_diff
    #print(a/m.pi*180)
  return t_total_sec

# t_diff step시간 후 각도변화 수치근사
def orbital_increment_angle_rad(a, t_diff, cel_data):
    r=orbit_ang2dist_km(a, cel_data)
    a+=t_diff*inst_orbit_speed(a, cel_data)/r
    return a


def time_next_solar_ecc_sec(t_diff, earth_data, moon_data):
  a_earth_prev=0
  a_earth=0

  a_moon_prev=0
  a_moon=0
  
  a_moon_node=0
  t_total_sec=0
  b1, b2, b3, b4 = 0, 0, 0, 0


  ang_v_node=t_diff/moon_data['nodal_precession_day']/24/60/60
  while(not((b1 and b2) or (b3 and b4)) or ((t_total_sec/60/60/24)<365)):

    a_moon_node+=ang_v_node
    a_moon_node=divmod(a_moon_node, 2*m.pi)[1]

    a_moon_prev=a_moon
    a_moon=divmod(orbital_increment_angle_rad(a_moon, t_diff, moon_data),2*m.pi)[1]

    a_earth_prev=a_earth
    a_earth=divmod(orbital_increment_angle_rad(a_earth, t_diff, earth_data),2*m.pi)[1]
    
    t_total_sec+=t_diff
    #print(a_earth*180/m.pi)

    b1=(a_moon_prev<a_moon_node)and(a_moon_node<a_moon)
    b2=(a_earth_prev<a_moon_node)and(a_moon_node<a_earth)
    b3=(a_moon_prev<(divmod((a_moon_node+m.pi),2*m.pi)[1]))and((divmod((a_moon_node+m.pi),2*m.pi)[1])<a_moon)
    b4=(a_earth_prev<(divmod((a_moon_node+m.pi),2*m.pi)[1]))and((divmod((a_moon_node+m.pi),2*m.pi)[1])<a_earth)


    #print(b1, b2, b3, b4)
    # print((a_earth-a_moon_node)*180/m.pi,( a_moon-a_moon_node)*180/m.pi, t_total_sec/60/60/24)
    # file.write("{}, {}, {}\n".format((a_earth-a_moon_node)*180/m.pi,( a_moon-a_moon_node)*180/m.pi, t_total_sec/60/60/24))

    print((a_earth_prev-a_moon_node)*180/m.pi,( a_earth-a_moon_node)*180/m.pi, t_total_sec/60/60/24)
    file.write("{}, {}, {}\n".format((a_earth_prev-a_moon_node)*180/m.pi,( a_earth-a_moon_node)*180/m.pi, t_total_sec/60/60/24))
    #print(t_total_sec/60/60/24)
  return t_total_sec


sun_data={
          'diameter_km': 1392700,
          'GM_m3s-2': 1.32712440018*10**20
          }
earth_data={
            'eccentr': 0.016710219,
            'long_r_km': 149597887.5, 
            'sidereal_day_h':23.9344696,
            'diameter_km': 12742,
            'GM_m3s-2': 3.986004418*10**14,
            'orbit_center': sun_data
            }
moon_data={
           'eccentr': 0.0549,
           'long_r_km': 384400,
           'orbit_incline_deg': 5.145,
           'apisidal_precession_day': 3232.6054, 
           'nodal_precession_day': 6798.383,  # 달 공전궤도 세차주기
           'orbital_period_day': 27.2122,
           'diameter_km': 3474.8,
           'orbit_center': earth_data
           }


referance={
          'moon_longitude_asc_node':39
          }

init_config={1,2}

t_diff_s=3600
currtime=0




#print(orbital_period_num_s(t_diff_s, earth_data)/60/60/24)
#print(orbital_period_s(earth_data)/60/60/24)
#print(time_next_solar_ecc_sec(t_diff_s, earth_data, moon_data)/60/60/24)
#365.24861111111113
