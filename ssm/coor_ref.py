import math as m
import numpy as np
import struct



################### Custom Exception Class ####################

class SequenceSizeError(Exception):
  def __init__(self, valid_size, param_name, sequence):
    self.valid_size=valid_size
    self.param_name=param_name
    self.sequence=sequence

  def __str__(self):
    return "SequenceSizeError: {}, size of parameter {} should be {}".format(self.sequence, self.param_name, self.valid_size)

  @classmethod
  def len_test(cls, valid_size, param_name, coor):
    if len(coor)!=valid_size:
      raise cls(valid_size, param_name, coor)

### 좌표는 Numpy를 이용
### 다차원도 가능(여러 벡터의 나열))


class coor3:  # 좌표계상 수치쌍 class(기준틀 하나에 대해)

  # 입력 규칙
    # 1. coor_mode : 's'(구면)와 'a'(직교)만 가능
    # 2. *coors : 3가지 형태의 자료형만 인정한다.

    #      [type1] 길이 3인 기본 Sequence형
    #      [type2] 길이 3인 1D배열
    #      [type3] 높이(행길이)가 3이고 열길이 자유인 2D배열(shape: 3xn, ndim: 2)

    #    개수는 상관없으며, 위의 3개가 섞여있어도 된다

  # set of spatial points : 기본적으로 좌표들을 공간점에 대한 묘사로 나타내는 것이 목표이므로, 같은 좌표계에서 중복되는 좌표를 삭제
  # staticity 규칙 : 한번 생성된 coor3 obj는 내부 좌표값을 임의로 바꾸지 못한다(좌표계 변경 제외)

  __valid_coor_mode = ('s', 'o')

  def __init__(self, coor_mode='s', *coors):
    self.reset_coor(coor_mode, *coors)

  def reset_coor(self, coor_mode='s', *coors):
    self.vec=np.hstack(tuple(map(self.__typecast, coors)))  # 입력들을 각각 형변환하여 하나의 2D배열(3xn 행렬)로 저장
    self.coor_mode=self.__coor_mode_test(coor_mode)       # 직교'o' 또는 구면's'
    self.current=0
    self.__l=len(self)

  # 각각에 대해 길이가 맞는지 확인 후 [type]으로 형변환시켜 반환 
  @staticmethod
  def __typecast(coor):
    if type(coor)==np.ndarray:
      coor=coor[:,None] if np.ndim(coor)==1 else coor
    else:
      coor = np.array(coor)[:,None]
    SequenceSizeError.len_test(3, 'coor', coor) # 입력 좌표의 크기가 3이 아니면 오류
    return coor

  def tuplify(self, index=None): 
    lst=[tuple(v) for v in self.vec.transpose()]
    if index!=None: lst=lst[index]
    return lst

  def normalize(self):
    self.vec=self.vec/abs(self)  
  
  # 허용된 좌표계가 아닐 경우 오류
  def __coor_mode_test(self, coor_mode):
    if not (coor_mode in self.__valid_coor_mode) : raise Exception('Please enter valid coor_mode from {}'.format(str(self.__valid_coor_mode)))
    return coor_mode


  # self의 좌표계 변경
  def conv_coor_modeOS(self):
    conv_mode = {'s': ('so', 'o'), 'o': ('os', 'o')}[self.coor_mode]
    self.vec=self.__convOS_bidir(self.vec, conv_mode[0])
    self.coor_mode=conv_mode[1]
    return self
  

  # 같은 정보를 가진 새로운 coor3 객체 반환
  def copy_coor(self):
    return coor3(self.coor_mode, np.copy(self.vec))

  # XY평면 상 원점을 중심으로 하는 원의 coor obj 제작
  @classmethod
  def circle_coors(cls, radius=1, num=500, ang_d=m.pi/2):
    return cls('s', *[(radius, ang, ang_d) for ang in np.linspace(0, 2* m.pi, num)]).conv_coor_modeOS()

  # XY평면 상 원점을 중심으로 하는 원의 coor obj 제작
  @classmethod
  def conic_coors(cls, conic=None, ratio=1, e=None, a_km=None, num=500, ang_d=m.pi/2):
    if not conic : conic = lambda ang: a_km * (1 - e**2) / (1 + e * m.cos(ang))
    return cls('s', *[(conic(ang)*ratio, ang, ang_d) for ang in np.linspace(0, 2* m.pi, num)]).conv_coor_modeOS()

  ################ operation overload ################

  # add two vectors
  def __add__(self, other): # self + other
    return coor3('o', self.__vecO+other.__vecO)

  # sub two vectors
  def __sub__(self, other): # self - other
    return coor3('o', self.__vecO-other.__vecO)



  # vector norm
  def norm(self):
    return self.__normO(self.__vecO)
    
  # norm of each vector
  def __abs__(self):  # abs(self)
    return self.norm()


  # inner product
  def inner_pdt(self, other=None):
    other_vec= None if other==None else other.__vecO
    return self.__inner_pdtO(self.__vecO, other_vec)

  # dot product(returns scalar)
  def __mul__(self, other): # self * other
    return self.inner_pdt(other)


  # 두 점 사이 거리
  def dist(self, other):
    return (self-other).norm()

  # vector들을 divint배 줄이기
  def div(self, divint):
    new=self.copy_coor()
    new_mode=(new.coor_mode=='s')
    if new_mode : new.conv_coor_modeOS()
    new.vec=new.vec/divint
    if new_mode : new.conv_coor_modeOS()
    return new
  
  def __truediv__(self, other): # self / other
    if type(self)==type(other) : return self.dist(other)
    elif type(other)==type(int()) : return self.div(other)
      


  # self(선, 원점지나고 self.vec와 평행)에서 other(점)까지 거리
  def distFL(self, other): # self // other
    return abs(other.__vecO-(((self*other)/(self*self))*self.__vecO))

  def __floordiv__(self, other):
    return self.distFL(other)

  # self와 other 사이 라디안 각
  def ang_rad(self, other):
    return m.acos((self*other)/(abs(self)*abs(other)))

  def __mod__(self, other): # self % other
    return self.ang_rad(other)
  
  # negate self
  def __neg__(self):
    return coor3(self.coor_mode, -np.copy(self.vec))


  ################ Magic Methods ################

  # len() 함수를 사용하면 coor벡터가 몇개 있는지 출력(열 길이)
  def __len__(self):
    return len(self.vec.T)

  # indexing, sliceing 가능
  def __getitem__(self, index):
    if type(index)==int:
      if index < self.__l:
        return coor3(self.coor_mode, self.vec[:,index])
      else:
        raise IndexError
    elif type(index)==slice:
      start, stop, stride = index.indices(self.__l)
      if stop <= self.__l:
        return coor3(self.coor_mode, self.vec[:,start:stop:stride])
      else:
        raise IndexError
  
  # iterator
  def __next__(self):
    if self.current < self.__l:
      r = coor3(self.coor_mode, self.vec[:,self.current])   # 내부 np vec만 보낼까?, 어디에 for문 쓸지 나중에 보고
      self.current += 1
      return r
    else:
      raise StopIteration

  # for iteratability
  def __iter__(self):
    self.index = 0
    return self
  
  # vec 출력
  def __repr__(self):
    return str(self.tuplify())



  ################ for internal np operation ################
  # np broadcasting 덕분에 각 요소와 스칼라간의 사칙연산 가능
  # 행렬(벡터(좌표)의 나열)도 계산 가능

  # return orthogonal mode vector
  @property
  def __vecO(self):
    return self.vec if self.coor_mode=='o' else self.__convOS_bidir(self.vec, 'so')


  ############### 직교 <--> 구면 ###############
  # 구면의 r은 비례값이나 다름없으므로, 변환할 때는 r=1인 경우에 대해서 변환하도록 한다.
  # 변수명(좌표는 np 배열 객체로 표현) : vecS(구면좌표/3D), vecO(직교좌표/3D), vecO_1(r=1인 vecO/3D), vec_r(r/1D), vecA(경도, 위도/2D) 
  # 방식 : vecS <--> (vec_r, vecA) <--> (vec_r, vec_O1) <--> vecO
 
  v_sin, v_cos, v_asin, v_acos, v_atan, v_sqrt = tuple(map(np.vectorize, (m.sin, m.cos, m.asin, m.acos, m.atan, m.sqrt)))
  
  ###############################
  # (1) vecS <--> (vec_r, vecA)

  # vecS(r, theta , phi) -> vec_r(r), vecA(경도, 위도)
  @classmethod
  def __slice_r(cls, vecS):   # synomym to convSA
    return vecS[0,:], np.vstack((vecS[1,:], m.pi/2-vecS[2,:]))
  
  # vec_r(r), vecA(경도, 위도) -> vecS(r, phi, theta)
  @classmethod
  def __unslice_r(cls, vec_r, vecA):   # synomym to convAS
    return np.vstack((vec_r, vecA[0,:], m.pi/2-vecA[1,:]))
  

  ###############################
  # (2) vecO <--> (vec_r, vec_O1)
  
  @classmethod
  def __inner_pdtO(cls, vecO1, vecO2=None):
    if type(vecO2)!=np.ndarray: vecO2=vecO1
    vecnew=vecO1*vecO2
    return vecnew[0,:]+vecnew[1,:]+vecnew[2,:]

  # norm(r) of vecO : vecO(x, y, z) --> vec_r(r)
  @classmethod
  def __normO(cls, vecO):
    vecO_pow2=cls.__inner_pdtO(vecO)
    return cls.v_sqrt(vecO_pow2)

  # vecO(x, y, z) ---> vec_r(r), vecO_1(x, y, z)
  @classmethod
  def __normalizeO(cls, vecO):
    vec_r = cls.__normO(vecO)
    return vec_r, vecO/vec_r
  
  # vec_r(r), vecO_1(x, y, z) ---> vecO(x, y, z)
  @classmethod
  def __denormalizeO(cls, vec_r, vecO_1):
    return vec_r*vecO_1
  
  ###############################
  # (3) (vec_r, vecA) <--> (vec_r, vec_O1)
  @classmethod
  def __convOA_1_bidir(cls, vec_OA, conv_mode):
    if conv_mode=='oa':
      return np.vstack((cls.v_atan(vec_OA[1,:]/vec_OA[0,:])+(vec_OA[0,:]<0)*m.pi
      , cls.v_asin(vec_OA[2,:])))  # 위도시스템은 기준면에서 올라가므로 z에 asin
    elif conv_mode=='ao':
      return np.vstack((cls.v_cos(vec_OA[1,:])*cls.v_cos(vec_OA[0,:]), cls.v_cos(vec_OA[1,:])*cls.v_sin(vec_OA[0,:]), cls.v_sin(vec_OA[1,:])))

  ###############################
  # (4) final : vecO(구면) <--> vecS(직교))
  @classmethod
  def __convOS_bidir(cls, vec_OS, conv_mode):
    if conv_mode=='os' :
      vec_r, vecO_1 = cls.__normalizeO(vec_OS)
      vecA = cls.__convOA_1_bidir(vecO_1, 'oa')
      return cls.__unslice_r(vec_r, vecA)
    elif conv_mode=='so' : 
      vec_r, vecA = cls.__slice_r(vec_OS)
      vecO_1 = cls.__convOA_1_bidir(vecA, 'ao')
      return cls.__denormalizeO(vec_r, vecO_1)




class reffrms:   # 기준틀 데이터 class

  # 기저기준틀 base_reffrm 지정
  # 절대기준틀은 base_reffrm=None으로 표시(자신을 묘사하는 다른 기준틀이 없을 때))
  # 하부기준틀 : 자신을 기저기준틀로 사용하는 다른 기준틀(기준틀들의 관계를 상부node에서 search/return하기 위해 필요)

  def __init__(self, node_lon_rad=0, incln_rad=0, lon_ref=0, name='reffrm', x_dsp=0, y_dsp=0, z_dsp=0):
    self.reset_reffrm(node_lon_rad, incln_rad, lon_ref, name, x_dsp, y_dsp, z_dsp)


  def reset_reffrm(self, node_lon_rad=0, incln_rad=0, lon_ref=0, name='reffrm', x_dsp=0, y_dsp=0, z_dsp=0):
    self.convM_O=self.__total_conv_matmaker(node_lon_rad, incln_rad, lon_ref, x_dsp, y_dsp, z_dsp)
    self.name=name+('_hash' if name=='reffrm' else '') 

  # self<->base 좌표변환
  def base_conv(self, coor3_source, conv_mode):
    return self.__M_coor_conv(coor3_source, conv_mode, self.convM_O)

  # self->target_reffrm 좌표변환
  def coor_conv(self, coor3_source, target_reffrm):
    coorS_base=self.base_conv(coor3_source, 'sa')
    coorS_target=target_reffrm.base_conv(coorS_base, 'as')
    return coorS_target


################# Magic methods ##################

  def __repr__(self):
    return '\n\nname : {}\n - characteristic matrix : \n   {}\n'.format(self.name, str(self.convM_O['as']).replace('\n', '\n   ')).replace('\n', '\n    ')

  # 포함 관계 a in tree를 사용하기 위해
  def __contains__(self, other):
    if other==None: return True
    return True if hash(other) in [hash(i) for i in iter(self)] else False

  # 두 기준틀이 같은지를 파이썬이 판단하기 위해서 obj의 hash를 이용한다.
  # 여기서 coorM_O의 각 요소(실수, float)를 floating point binary string으로 바꾸어준 뒤 모두 concatinate하여 int()하여준다
  def __hash__(self):
    hsh=''
    for convL in list(self.convM_O['sa']):
      for co in convL:
        hsh+=self.binary(co)
    return int(hsh)

  # 인터넷에서 따온 것, 소수점을 floating point binary로 바꾸어준다
  def binary(self, num): return ''.join('{:0>8b}'.format(c) for c in struct.pack('!f', num))

################# 연산자 오버라이딩 ##################

  # 같은 기준틀 : 같은 기저기준틀을 가지며, 같은 특성행렬을 갖는 것
  def __eq__(self, other):
    if other==None:
      return False
    return hash(self)==hash(other)

######################### private ##########################
  # angle range correction(rad, list of angles), 각도(rad)의 크기를 2pi로 교정해주는 메소드
  def __angle_correction(self, *ang_rad):
    return (val%(2*m.pi) for val in ang_rad)

######################### 특성행렬 계산 ##########################

  # abs<->self 좌표계 회전+평행이동 변환 행렬 계산 메소드(4x4 행렬)
  def __total_conv_matmaker(self, *char_val):
    a=char_val
    eu1=np.array([[m.cos(a[0]), -m.sin(a[0]), 0, a[3]], [m.sin(a[0]), m.cos(a[0]), 0, a[4]], [0, 0, 1, a[5]], [0, 0, 0, 1]])
    eu2=np.array([[1, 0, 0, 0], [0, m.cos(a[1]), -m.sin(a[1]), 0], [0, m.sin(a[1]), m.cos(a[1]), 0], [0, 0, 0, 1]])
    eu3=np.array([[m.cos(a[2]), -m.sin(a[2]), 0, 0], [m.sin(a[2]), m.cos(a[2]), 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    M=eu1@eu2@eu3       # self -> abs
    M_inv=np.linalg.inv(M)  # abs -> self
    return {'sa' : M, 'as' : M_inv}

  # convM가 주어졌을 때 coor 변환
  @classmethod
  def __M_coor_conv(cls, coor_source, conv_mode, convM):
    coor_mode_conv = (coor_source.coor_mode=='s')
    if  coor_mode_conv : coor_source.conv_coor_modeOS()

    coorO_source=np.vstack([coor_source.vec, [1]*len(coor_source)])  # 1을 붙여 4차원 벡터로 만든다 (new)
    coorO_dest=convM[conv_mode]@coorO_source
    coor_dest=coorO_dest[0:3,:]                                      # 다시 끝에 1을 자른다(new)

    dest=coor3('o', coor_dest)
    if coor_mode_conv: dest.conv_coor_modeOS()
    return dest

