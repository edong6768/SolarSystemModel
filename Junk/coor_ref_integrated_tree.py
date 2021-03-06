import math as m
import numpy as np
import json
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


###############################################################
#         Referance frames for spherical coordinates          #
###############################################################

# 구면기준좌표계의 표현 : 절대기준틀에 대한 특정 기준틀(또는 기준좌표, Frame of referance)) 관계를 지정

# <관측 위치 이동을 고려한 좌표시스템>

# 0. 용어정리 (''표시가 있는 용어는 모두 여기서 정의된 용어들)
#    (1) '공간점'과 '수치쌍' : 좌표를 구하는 것은 공간상의 점을 숫자의 나열로 나타낸다는 것이라고 볼 수 있다.
#                            통상적으로 '좌표'라는 용어의 사용은 중의적이며 이를 명확히 구분하기 위해 다음과 같은 용어로 정리한다.
#                               - 공간점(spatial point) : 공간상의 점으로써의 좌표를 의미한다. 관측자의 관측위치나 방식에 상관 없다.
#                               - 수치쌍(numeric tuple) : 숫자쌍의 쌍으로써의 좌표를 의미한다.
#                            사용예시)
#                             1) 다른 숫자들로 표현되어도, 공간상 같은 지점을 나타낸다면, '같은 공간점을 다른 수치쌍으로 표현하였다'고 할 수 있다.     
#                             2) 같은 숫자들로 표현되어도, 공간상 다른 지점을 나타낸다면, '다른 공간점을 같은 수치쌍으로 표현하였다'고 할 수 있다.

#    (2) '좌표계'와 '기준틀' : (관측자의 위치/방향)와 (좌표방식)을 다음 두가지로 용어로 정리한다.
#                                - 좌표계(coordinate mode) : 좌표 방식을 나타낸다.(직교좌표 방식(o), 구면좌표 방식(s), 원통...)
#                                - 기준틀(referance frame) : 관측자 위치와 관측기준의 방향(orientation)을 나타낸다.
#                             좌표계와 기준틀에 따라, 공간점-수치쌍의 연결(대응)이 달라진다.

#    (3) '묘사' : '공간점'과 '수치쌍'이 연결된 것을 '묘사'라고 한다. 이는 '좌표계'와 '기준틀'을 특정하는 것과 동치이다.
#                 즉 공간상의 점을 숫자로 나타내기 위해서, 관찰자의 위치와 좌표를 표시하는 방식이 정해져야 한다.
#                     사용예시) 위와 같은 상황을 '공간점이 좌표계와 기준틀에 대하여 수치쌍으로 묘사되었다'라고 한다

#    (4) 좌표 변환 : 한 좌표계/기준틀에 대해 묘사된 공간점을 다른 좌표계/기준틀에 대해 묘사하는 것(동일한 공간점을 대상)

#                   위치의 측정방식만 다룰 때와 기준틀의 차이를 다룰 때의 변환을 문맥에 따라 '좌표계 변환'이라는 같은 용어로 사용한다.
#                   즉 일상적으로는 문맥에 따라 '좌표계 변환'을 2가지 의미로 모두 사용한다. 

#                   둘 사이 정확한 구분을 위해 다음과 같은 용어를 사용한다.
#
#                          변하는 요소 |    좌표계    |    기준틀      
#                        -------------|--------------|--------------
#                      - '좌표계 변경' |      O       |      X
#                      - '기준틀 변환' |      X       |      O
#            
# 1. 기준틀의 묘사
#    (0) 관측자와 '기준틀' : 공간상의 점을 묘사하는 관측자가 하나만 있다면 좌표계만 정하여 점 간의 관계를 나타낼 수 있다.
#                           하지만 관측자가 둘 이상 된다면, 관측자간의 관계를 알아야 각자가 '묘사'한 '수치쌍'이 나타내는 '공간점'을 연결할 수 있다.
#                           즉, 기준틀은 관측상황이 다른 둘 이상의 관측자가 존재할 때 필요하며, 둘 간의 관계를 알면 두 기준틀의 관계를 알 수 있다.
#    (1) '기준틀'의 도입 : 기준틀의 개념의 발생한 근간이 되는 질문은 '내가 묘사한 좌표가 저 사람이 묘사한 좌표와 같은 것을 나타낼까'이다. 
#                         즉 저 사람이 본 좌표가 내가 본 어떤 좌표와 같은지를 보는 것과 같으며, 그러면 '내와 저 사람이 어떤게 다르기 때문인가'에 대한
#                         물음이 생기게 된다. 이에 대한 대답을 함축하여 '기준틀이 다르다'라고 하는 것이다.
#                         여기서 볼 수 있듯이, 기준틀을 다룰 때의 관심사는, 내(또는 하나의) 기준틀에 대해 다른 기준틀을 어떻게 나타낼 수 있는가이다.

#    (2) 기저기준틀(base referance frame) : 한 기준틀이 수치적으로 묘사되는데 쓰이는 다른 기준틀
#    (3) 다른 기준틀 : 다른 기준틀에서 본 좌표가 내 기준틀에서는 어떤 좌표인지를 보기 위해서는, 다른 기준틀이 내 기준틀에 대해 어떻게 표현 될 수 있는지를
#                     보면 된다.
#                     앞서 기준틀을 [관측자의 위치]와 [관측기준의 방향]을 나타낸다고 하였다. 내 기준틀을 중심으로 다른 기준틀을 나타내는 것도
#                     이와 같은 2가지 카테고리로 나눌 수 있으며, 내 기준에 대한 다른 기준틀의 관계를 '관계파라미터(relation parameter, rp)'로 
#                     다음과 같이 수치적으로 나타낼 수 있다.
# 
#                          - 위치적 관계파라미터(positional - ,prp): 상대 관측자의 공간상 위치를 내 기준틀에서 좌표계로 나타낸 수치값(직교/구면 상관X)
#                          - 방향적 관계파라미터(orientational -, orp): 상대 관측자의 관측기준면 방향의 내 기준틀에 대한 오일러각(차후 설명)
#
#                     여기서 위치는 직관적으로 이해 되지만, 방향은 오일러각을 알아야 하므로 바로 이해되기는 어려울 수 있다. 여기서 오일러각을 설명하기는
#                     무리가 있으므로, 각각 변수에 다한 명칭과 통상적인 말로써의 뜻만 설명하겠다.

#    (4) 방향적 관계파라미터 : 용어를 알아듣기 덜 어렵게 하기 위해서는 '천체구면좌표계'적인 용어를 이용하면 쉬워진다. 오일러각을 안다면 다음 용어 순서대로
#                             제 1, 2, 3 오일러각이다.

#                                   잠시 천체구면좌표 관련 몇가지 정의를 소개하자면
#                                         - 구면좌표계에서의 위쪽의 정의 : 경도각 측정방향이 반시계방향으로 보이는 쪽이 위쪽이다
#                                         - 승교점 : 기준면이 겹치는 두 점 중 2개 교점 중 경도측정방향을 따라 기준좌표면 아래쪽에서 위쪽으로 올라오는 점
#         
#                             방향적 관계파라미터는 다음과 같다
#                                   1) 승교각(node_lon_rad) : 승교점과 공통기준좌표계 경도
#                                   2) 경도면 기울기(incln_rad)
#                                   3) 경도측정 기준점(lon_ref) : 승교점과의 경도차로 표현한다.
#
#                            위 각각 3가지, 합쳐서 6가지의 파라미터로 다른 기준틀의 정보를 표현할 수 있으며, 다른 기준틀에서 묘사된 좌표('수치쌍')을 
#                            모두 내 기준틀의 좌표로 변환할 수 있다. 둘 다 직교좌표계로 표현된 좌표일 경우(물론 구면이어도 되는지 여부는 모르지만), 
#                            위 6가지 매개변수를 적절히 조합하여 4x4 변환행렬로 표현할 수 있다.

#    (5) 표준기준틀 : 우리는 내 기준틀에서 다른 기준틀의 관계를 구할 수 있다는 것을 알았으며, 이는 내 기준틀 뿐만 아니라 다른 기준틀 사이의 관계들도 나타낼
#                    수 있음을 보인다. 그렇다고 우리는 내 기준틀와의 관계 이외의 다른 기준틀 사이의 관계들을 모두 지정해주지 않는다.
#                    위에서 우리는 다른기준계를 내 기준계에 대한 4x4 변환행렬로 나타낼 수 있음을 보았다. 여기서 행렬이란 것을 생각해보면, 행렬을 2번 취하면
#                    2번의 변환을 1번해 하여주는 새로운 4x4 변환행렬을 만들 수 있음을 생각해 볼 수 있다. 이는 즉 만약 두 개 다른 기준틀 사이 관계를 나타내고 싶을 때
#                    굳이 각자 기준틀을 중심으로 새로 관계를 계산할 필요가 없다는 것을 의미한다. 그렇게 한다면 모든 기준틀간의 관계는 기준틀이 늘어날 때 마다 n(n-1)/2
#                    개를 구하여야 한다. 이럴 필요 없이 우리는 내 기준틀과 다른 기준틀과의 관계만 나타낸다면(즉 n-1번만 계산하면), 다른 기준틀 사이의 관계도
#                    구할 수 있다는 것을 뜻한다.
#                    이렇게 내 기준계와 다른 기준계와의 관계만 있으면 된다고 했을 때, 여러 기준계 중에서 내 기준계는 특별한 위상을 지니게 된다.
#                    여기서 다른 기준틀의 표준이 된다는 의미에서 이때 내 기준계를 '표준기준틀'이라고 한다. 즉 다시 정리할 수 있다.

#                       - 표준기준틀 : 다른 기준틀의 관계를 기술하는 중심이 되는 기준계
          
#
#    (6) 내 기준틀 : 앞에서 보았다 싶듯이 기준계의 표현은 내 기준틀에 대한 다른 기준틀의 표현으로 이루어지므로, 여기서 내 기준틀을 논하는 것 자체는 의미론 적으로
#                   무의미 하다. 하지만 나의 위치가 달라질 때, 내 기준틀에 관한 묘사가 필요하게 된다. 이는 변하지 않는 특정 기준틀에 대해 내 기준틀을 묘사하면
#                   할 수 있다. 여기쯤 오면 나의 기준틀이라는 용어가 의미있는지에 대한 의문이 생긴다. 나의 기준틀이 변한다면 다른 기준틀로 내 기준틀이 움직이는
#                   양상을 표현해야 하기 때문이다. 즉 내 기준틀이 움직인다면 다른 기준틀의 절대적일 기준이 될 수 없다. 하지만 그렇다고 다른 기준틀들은 안 움직인다는
#                   보장 또한 없다. 그렇다고 움직이는 기준틀로 다른 움직이는 기준틀을 묘사하는 것은 계산적으로 너무 복잡해지며 비효율적이다.

#                   여기서 우리는 움직이지 않는 한 기준계를 중심으로 다른 기준계들을 묘사하고 싶어지며, 변하지 않는 가상의 기준틀의 개념이 필요하다.
#                   이 표준기준틀을 '절대기준틀'이라고 한다. 다시 정리하자면
#
#                      - 절대기준틀 : 변하지 않는 가상의 표준기준틀

#    (7) 절대기준틀의 역설 : 우리는 위에서 절대기준틀이라는 개념이 필요한 이유를 알았다. 하지만 여기서 우리는 한가지 역설에 도달한다.
#                             '기준틀의 위치는 다른 기준틀로 묘사되는데, 절대기준틀이 움직이지 않는 것은 어떻게 묘사할까?'
#                          이 말을 듣는다면 역설적으로 들린다. 하지만 이는 절대기준틀에 대한 잘못된 이해에서 오는 것이다.

#                          우리는 다른 기준틀을 수치적(관계파라미터)로 묘사하는 시도를 하는 시점부터, 그것의 중심이 되는 기준틀이 변하지 않는다는 신뢰를
#                          가정하고있다. 즉, 기준틀을 묘사하려는 시도에 이미 절대기준틀의 존재를 내포하고 있다.

#                          이러한 것들을 들으면서 자연스럽게 우리는 물리에서의 상대성이론과 대응시켜보게 된다. 상대성이론은 시간과 공간에 대한 절대성의 허물을
#                          깨어버렸다. 역설적이게도 상대성이론은 이 절대성의 역설을 깨기 위해 다른 절대성을 이용하고 있다. 바로 광속불변성이다. 상대성이론은
#                          광속불변이 물리적세계의 절대적 규칙이라 규정하였으며, 이를 통해 기준에 절대적이라고 여겼던 시간/공간을 광속불변에 대한 상대적인 관계를
#                          가진다는 것을 이끌어냈다. 또 다른 절대성을 도입하여 기존의 절대성을 가진다고 생각했던 시간/공간을 이에 종속시킨 것이다.
#                          위에서 내 기준틀이 움직인다는 것을 알고 나서 어떻게 해야 되냐는 의문과 같은 것은 겪은 것이다.

#                          여기서 물리는 절대적인 것을 과학적 방법론으로 사실이라고 결론내렸기 때문에 이에 대한 신뢰를 가지게 되었으며, 과거에 절대적이라고 믿고
#                          있던 것을 수정하였다. 이는 과거 사람들이 사실관계을 명확히 확인하지 않은 잘못이라고 할 수 있을까? 과거 사람들이 잘못 한것인가?

#                          잘못하였다고 볼 수 없다. 그저 여러 물리적인 '것'들 사이의 관계를 알아차렸을 때, 그 '관계'를 정리하다 보면 관계의 뿌리에 있는 절대적인
#                          것을 하나 지정하지 않고서는 불안을 느끼게 된다. 모든 알려진 관계내에서 이 관계들을 부정하지 않는 하나를 절대적으로 두는 것도 이러한
#                          사람이 생각하는 방식과 연관있다. 절대적인 것 없이는 사람은 그걸 이용할 동기가 없을 것이다.

#    (8) 태양계에서 절대기준틀 : 태양계를 묘사하는데 있어서 우리도 우리가 편하다고 느끼는 하나의 절대기준틀을 정하여야 한다. 물론 우리는 태양이 은하중심을 돌고 있다는
#                              것을 알고 있다. 하지만 여기서는 태양계만 생각할 것이기 때문에, 여기서는 '지구공전면'과 '춘분점'을 절대기준틀로 잡으려고 한다.





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
    self.vec=np.unique(np.hstack(tuple(map(self.__typecast, coors))), axis=1)   # 입력들을 각각 형변환하여 하나의 2D배열(3xn 행렬)로 저장, 중복벡터 삭제
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
  
  # 허용된 좌표계가 아닐 경우 오류
  def __coor_mode_test(self, coor_mode):
    if not (coor_mode in self.__valid_coor_mode) : raise Exception('Please enter valid coor_mode from {}'.format(str(self.__valid_coor_mode)))
    return coor_mode

  # self의 좌표계 변경
  def conv_coor_modeOS(self):
    conv_mode = {'s': ('so', 'o'), 'o': ('os', 'o')}[self.coor_mode]
    self.vec=self.__convOS_bidir(self.vec, conv_mode[0])
    self.coor_mode=conv_mode[1]
  

  # 같은 정보를 가진 새로운 coor3 객체 반환
  def copy_coor(self):
    return coor3(self.coor_mode, self.vec)


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
    return str(self.vec)


  ################ operation overload ################

  # create new coor3 obj containing vectors of both operands
  def __add__(self, other):
    if self.coor_mode!=other.coor_mode : 
      raise Exception('Cannot combine coor3s between different coor_modes')
    return coor3(self.coor_mode, self.vec, other.vec)

#  def __sub__(self, other):

  ################ for internal np operation ################
  # np broadcasting 덕분에 각 요소와 스칼라간의 사칙연산 가능
  # 행렬(벡터(좌표)의 나열)도 계산 가능


  ############### 직교 <--> 구면 ###############
  # 구면의 r은 비례값이나 다름없으므로, 변환할 때는 r=1인 경우에 대해서 변환하도록 한다.
  # 변수명(좌표는 np 배열 객체로 표현) : vecS(구면좌표/3D), vecO(직교좌표/3D), vecO_1(r=1인 vecO/3D), vec_r(r/1D), vecA(경도, 위도/2D) 
  # 방식 : vecS <--> (vec_r, vecA) <--> (vec_r, vec_O1) <--> vecO
 
  v_sin, v_cos, v_asin, v_acos, v_sqrt = tuple(map(np.vectorize, (m.sin, m.cos, m.asin, m.acos, m.sqrt)))
  
  ###############################
  # (1) vecS <--> (vec_r, vecA)

  # vecS(r, theta , phi) -> vec_r(r), vecA(경도, 위도)
  def __slice_r(cls, vecS):   # synomym to convSA
    return vecS[0,:], np.vstack((vecS[2,:], m.pi-vecS[1,:]))
  
  # vec_r(r), vecA(경도, 위도) -> vecS(r, theta , phi)
  def __unslice_r(self, vec_r, vecA):   # synomym to convAS
    return np.vstack((vec_r, m.pi-vecA[1,:], vecA[0,:]))
  

  ###############################
  # (2) vecO <--> (vec_r, vec_O1)
  
  # norm(r) of vecO : vecO(x, y, z) --> vec_r(r)
  def __normO(cls, vecO):
    vecO_pow2=np.power(vecO,2)
    vec_r = cls.v_sqrt(vecO_pow2[0,:]+vecO_pow2[1,:]+vecO_pow2[2,:])
    return vec_r

  # vecO(x, y, z) ---> vec_r(r), vecO_1(x, y, z)
  def __normalizeO(cls, vecO):
    vec_r = cls.__normO(vecO)
    return vec_r, vecO/vec_r
  
  # vec_r(r), vecO_1(x, y, z) ---> vecO(x, y, z)
  def __denormalizeO(cls, vec_r, vecO_1):
    return vec_r*vecO_1
  
  ###############################
  # (3) (vec_r, vecA) <--> (vec_r, vec_O1)
  def __convOA_1_bidir(cls, vec_OA, conv_mode):
    if conv_mode=='oa':
      return np.vstack((cls.v_atan(vec_OA[1,:]/vec_OA[0,:]), cls.v_asin(vec_OA[2,:])))  # 위도시스템은 기준면에서 올라가므로 z에 asin
    elif conv_mode=='ao':
      return np.vstack((cls.v_sin(vec_OA[1,:])*cls.v_cos(vec_OA[0,:]), cls.v_sin(vec_OA[1,:])*cls.v_sin(vec_OA[0,:]), cls.v_sin(vec_OA[1,:])))

  ###############################
  # (4) final : vecO(구면) <--> vecS(직교))
  def __convOS_bidir(cls, vec_OS, conv_mode):
    if conv_mode=='os' :
      vec_r, vecO_1 = cls.__normalizeO(vec_OS)
      vecA = cls.__convOA_1_bidir(vecO_1, 'oa')
      return cls.__unslice(vec_r, vecA)
    elif conv_mode=='so' : 
      vec_r, vecA = cls.__slice_r(vec_OS)
      vecO_1 = cls.__convOA_1_bidir(vecA, 'ao')
      return cls.__denormalizeO(vec_r, vecO_1)





class reffrm:   # 기준틀 데이터 class

  # 기저기준틀 base_reffrm 지정
  # 절대기준틀은 base_reffrm=None으로 표시(자신을 묘사하는 다른 기준틀이 없을 때))
  # 하부기준틀 : 자신을 기저기준틀로 사용하는 다른 기준틀(기준틀들의 관계를 상부node에서 search/return하기 위해 필요)
  #             self를 선언할 때 지정하지 않으며, 하위기준틀이 자신을 기저기준틀로 정할 때 추가된다.
  #             이는 자신을 하부기준틀로 지정하지 않은 채로 기저기준틀에서 지정할 경우 데이터구조의 무결성을 해칠 위험이 있기 때문이다.

  def __init__(self, r_dsp=0, lon_dsp=0, lat_dsp=0, node_lon_rad=0, incln_rad=0, lon_ref=0,  name='reffrm', base_reffrm=None):
    self.convM_O=self.__total_conv_matmaker(r_dsp, lon_dsp, lat_dsp, *self.__angle_correction(node_lon_rad, incln_rad, lon_ref))
    self.name=name+('_hash' if name=='reffrm' else '')
    self.set_base_reffrm(base_reffrm) #base_reffrm 지정, hash가 이름에 들어간다면 업데이트된 hash로 다시 이름짓는다
    self.sub_reffrm=set()  # 같은 기준틀에 대한 중복표현은 필요없으므로, set로 두어 자체적으로 중복관리를 하도록 한다.

    self.sub_size = 0      # 하위 기준틀 내부 tree 구조 전체의 reffrm 개수(이 또한 하위에서 self를 기저로 추가할때 결정된다)

    self.incident_sub_size = 0      # 인접 하위 기준틀 reffrm 개수(이 또한 하위에서 self를 기저로 추가할때 결정된다)
    # iterability를 위해 존재
    self.current = 0       # iterator 선언에 필요, tree 전체 구조 scan이므로 정지조건은 sub_size가 된다.
    self.iter_sub_reffrm=None   # 현재 sub_reffrm의 iter
    self.iter_curr_sub_reffrm=None   # 현재 search중인 sub_reffrm의 iter

  

  # 좌표계 방향변화에 맞추어 새로 설정
  def incremental_change_reffrm(self, r_dsp_diff=0, lon_dsp_diff=0, lat_dsp_diff=0, node_lon_rad_diff=0, incln_rad_diff=0, lon_ref_diff=0):
    trns_val_diff=[r_dsp_diff, lon_dsp_diff, lat_dsp_diff]
    rot_val_diff=[node_lon_rad_diff, incln_rad_diff, lon_ref_diff]
    new_val=[trns_val_diff[n]+self.char_trans_val[n] for n in range(3)]+[rot_val_diff[n]+self.char_rot_val[n] for n in range(3)]
    self.reset_reffrm(*new_val)

  # 자신의 tree node depth(깊이)와 절대기준틀 obj 반환
  def abs_search(self, n=-1):
    n+=1
    return self.base_reffrm.abs_search(n) if self.base_reffrm!=None else (n, self)

  # self<->base 좌표변환
  def base_conv(self, coorS_source, conv_mode):
    return self.__M_coor_conv(coorS_source, conv_mode, self.convM_O)

  # self-->abs 좌표변환
  def abs_conv(self, coorS_source, conv_mode):
    base_coor=self.__M_coor_conv(coorS_source, conv_mode, self.convM_O)
    base_frm= self.base_reffrm
    return base_frm.abs_conv(base_coor) if base_frm!=None else base_coor

  # self->target_reffrm 좌표변환
  def coor_conv(self, coor3_source, target_reffrm):
    coorS_abs=self.abs_conv(coor3_source, 'sa')
    coorS_target=target_reffrm.abs_conv(coorS_abs, 'as')
    return coorS_target

  # 여러 reffrm간 합성에 대한 coor 변환
  def composit_reffrm_coor_conv(cls, coor3_source, conv_mode, *reffrm):
    convM = [reff.convM_O for reff in reffrm]
    return cls.__composit_convM_coor_conv(cls, coor3_source, conv_mode, *convM)



########## reffrm간의 관계(tree 구조) 정보 구성 함수 ##########

  # 자신의 기저기준틀을 지정한다 (즉 이는 자신의 기저기준틀이 절대기준틀이 아니며, 이 또한 다른 기준틀을 기저로 가진다는 것을 의미한다)
  # 또한 지정된 기저기준틀의 하위기준틀 set에 자신을 추가한다.
  def set_base_reffrm(self, base_reffrm):
    if self==None:
      raise ValueError
    if base_reffrm!=None:
      if self in  base_reffrm.abs_search()[1]:
        raise Exception
      base_reffrm.__add_sub_reffrm(self)
      base_reffrm.__update_sub_size(self.sub_size) # 기저기준틀의 sub_size update
      base_reffrm.incident_sub_size+=1      # 기저기준틀의 incident sub_size update
    self.base_reffrm=base_reffrm
    self.name=(('reffrm_hash'+str(hash(self))) if self.name[0:11]=='reffrm_hash' else self.name)

  def __update_sub_size(self, subs_sub_size):
      self.sub_size+=subs_sub_size+1
      if self.base_reffrm:
        self.base_reffrm.__update_sub_size(subs_sub_size)

  def __add_sub_reffrm(self, sub_reffrm):
    self.sub_reffrm.add(sub_reffrm)



  # [[A, [ [B,[E, F]], [C,[G]] [D,[]] ] ]]--> A - {B, C, D} - {{E, F}, {G}, {}}
  @classmethod
  def l2t(cls, tree_list, base=None, n=0):
    for refs in tree_list:
      if len(refs)==2:
        if type(refs[1])==type(list()):
          refs[0].set_base_reffrm(base)
          cls.l2t(refs[1], refs[0], 1)
        else:
          refs.set_base_reffrm(base)
      else:
        refs.set_base_reffrm(base)
    if n==0 : return [roots[0] for roots in tree_list]


  # tree list가 각각의 이름으로 되어있을 때, '이름':reffrm dict와 함께 주어진 것을 tree구조 연결에 사용(재귀)
  # [['A', [ ['B',['E', 'F']], ['C',['G']] ['D',[]] ] ]]--> A - {B, C, D} - {{E, F}, {G}, {}}
  @classmethod
  def dict_l2t(cls, tree_name_list, reffrm_dict, base='', n=0):
    base=reffrm_dict[base] if base else None
    for r_name in tree_name_list:
      if len(r_name)==2:
        if type(r_name[1])==type(list()):
          reffrm_dict[r_name[0]].set_base_reffrm(base)
          cls.l2t(reffrm_dict, r_name[1], r_name[0], 1)
        else: reffrm_dict[r_name].set_base_reffrm(base)
      else: reffrm_dict[r_name].set_base_reffrm(base)
    if n==0 : return [reffrm_dict[roots[0]] for roots in tree_name_list]
  
  # A - {B, C, D} - {{E, F}, {G}, {}} --> [['A', [ ['B',['E', 'F']], ['C',['G']] ['D',[]] ] ](재귀)
  def t2l(self):
    if self.sub_reffrm: return [[self.name, [sub.t2l() for sub in self.sub_reffrm]]]
    else: return self.name

  # '이름' : reffrm 딕셔너리를 트리 전체 reffrm에 대하여 한다.  
  @property
  def dictify(self):
    frm_dict = dict()
    for frm in self:
      frm_dict[frm.name] = frm
    return frm_dict

  # reffrm 내부 트리 구조도
  @property
  def treeplot(self):
    branch_str = ['├──', '└──', '│  ', '   ' ]
    tree_str = ''
    for node in self:                           # in(contains)를 사용하여 트리 traversal
      if node==self: tree_str += '\n' + node.name # 가장 상단 node(root)는 그대로 name 출력
      else:                                 # 그 외 sub_node는 이름 앞 가지 문자 필요
        front_str = ''                          # 빈 str, 이름 앞의 가지 문자가 입력될 문자열
        curr_base_reffrm=node                   # 현재 node
        for depth in range(node.abs_search()[0]):   # 전체 조상노드 개수(node depth)
          curr_base_reffrm=curr_base_reffrm.base_reffrm   # 현재 노드의 조상(기저기준틀)
          idx=2 if depth else 0
          front_str=branch_str[(0 if (curr_base_reffrm.current < curr_base_reffrm.sub_size-1) else 1)+idx]+front_str
        tree_str += '\n'+ front_str + node.name
    return tree_str + '\n'

##################### magic methods ####################

  # 트리 구조도 출력
  def __str__(self):
    return self.treeplot

  # 트리 전체 내부 행렬 정보와 관계 
  def __repr__(self):
    # return '\nname : {}\nbase reffrm : {}\ncharacteristic matrix : \n{} \nsub reffrm : {}'.format(self.name, self.base_reffrm.name if (self.base_reffrm!=None) else 'Empty', self.convM_O['as'], [i.name for i in self.sub_reffrm] if self.sub_reffrm else 'Empty')
    return self.name

  # reffrm의 하부트리 크기 + 자신 = 전체트리크기(node수)
  def __len__(self):
    return self.sub_size + 1

# sub_reffrm의 하부 포함여부 (a in b, 트리 하부구조에 node가 존재하는지 여부 확인) 사용을 위한 iterable 선언


  def __next__(self):
    if self.current <= self.sub_size:
      try:
        if self.current==0:
          self.current+=1
          return self
        elif self.current==1: 
          self.iter_curr_sub_reffrm=iter(next(self.iter_sub_reffrm))
        self.current+=1
        return next(self.iter_curr_sub_reffrm)
      except StopIteration:
        self.iter_curr_sub_reffrm=iter(next(self.iter_sub_reffrm))
        return next(self.iter_curr_sub_reffrm)
    else:
      raise StopIteration

  # for iteratability
  def __iter__(self):
    self.current = 0
    self.iter_sub_reffrm=iter(self.sub_reffrm)
    return self

  # 포함 관계 a in tree를 사용하기 위해+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++>>>>>>>>>>>
  def __contains__(self, other):
    if other==None: return True
    return True if hash(other) in [hash(i) for i in iter(self)] else False

  # 두 기준틀이 같은지를 파이썬이 판단하기 위해서 obj의 hash를 이용한다.
  # 여기서 coorM_O의 각 요소(실수, float)를 floating point binary string으로 바꾸어준 뒤 모두 concatinate하여 int()하여준다
  def __hash__(self):
    convList=list(self.convM_O['sa'])
    if(self.base_reffrm!=None): convList+=list(self.base_reffrm.convM_O['sa'])
    hsh=''
    for convL in convList:
      for co in convL:
        hsh+=self.binary(co)
    return int(hsh)

  # 인터넷에서 따온 것, 소수점을 floating point binary로 바꾸어준다
  def binary(self, num):
      return ''.join('{:0>8b}'.format(c) for c in struct.pack('!f', num))

################# 연산자 오버라이딩 ##################

  # 같은 기준틀 : 같은 기저기준틀을 가지며, 같은 특성행렬을 갖는 것
  def __eq__(self, other):
    if other==None:
      return False
    elif self.base_reffrm==None or other.base_reffrm==None:
      return (hash(self)==hash(other))
    return (hash(self.base_reffrm)==hash(other.base_reffrm)) and (hash(self)==hash(other))

######################### private ##########################

  # 여러 coorM의 합성에 대한 coor 변환
  def __composit_convM_coor_conv(cls, coor3_source, conv_mode, *convM):
    tM=cls.__convM_composit(*convM)
    return cls.__M_coor_conv(cls, coor3_source, conv_mode, tM)

  # 여러 convM의 합성
  def __convM_composit(*convM):
      tM_inv=convM[0]['as']
      for i in convM[1:]:
        tM_inv@=i['as']
      tM = np.linalg.inv(tM_inv)
      return {'sa' : tM, 'as' : tM_inv}


  # angle range correction(rad, list of angles), 각도(rad)의 크기를 2pi로 교정해주는 메소드
  def __angle_correction(self, *ang_rad):
    return [val%(2*m.pi) for val in ang_rad]

######################### 특성행렬 계산 ##########################

  # abs<->self 좌표계 회전+평행이동 변환 행렬 계산 메소드(4x4 행렬)
  def __total_conv_matmaker(self, r_dsp, lon_dsp, lat_dsp, node_lon_rad, incln_rad, lon_ref):
      R=self.__rot_conv_matmaker((node_lon_rad, incln_rad, lon_ref))
      T=np.reshape(np.array((r_dsp, lon_dsp, lat_dsp)), (3, 1))
      M_inv=np.vstack([np.hstack([R, T]), np.array([0, 0, 0, 1])])   # abs -> self
      M=np.linalg.inv(M_inv)                                         # self -> abs
      
      return {'sa' : M, 'as' : M_inv}

  # abs<->self 좌표계 회전변환 행렬 계산 메소드
  def __rot_conv_matmaker(self, char_rot_val):
    if char_rot_val.count(0)==3:
      return np.eye(3)
    else:
      r = [np.eye(3)] * 3
      rot_M=lambda a : np.array([[1, 0, 0], [0, m.cos(a), m.sin(a)], [0, -m.sin(a), m.cos(a)]])
      a=char_rot_val
      for n in range(3):
        k, M = 0, rot_M(a[n])
        for raw in range(8):
          i, j = divmod(raw, 3)
          if i!=n and j!=n: 
            ki, kd = divmod(k, 2)
            r[n][i][j]=M[kd][ki]
            k+=1
      R_inv=r[0]@(np.transpose(r[1]))@r[2]  # abs -> self
      
      return R_inv

  # convM가 주어졌을 때 coor 변환
  def __M_coor_conv(cls, coor_source, conv_mode, convM):
    coor_mode_conv = (coor_source.coor_mode=='s')
    if  coor_mode_conv : coor_source.conv_coor_modeOS()
    coorO_source=np.vstack([coor_source.vec, [1]*len(coor_source)])  # 1을 붙여 4차원 벡터로 만든다 (new)
    coorO_dest=convM[conv_mode]@coorO_source
    coor_dest=coorO_dest[0:3,:]                                      # 다시 끝에 1을 자른다(new)
    
    return  coor3('o', coor_dest).conv_coor_modeOS() if coor_mode_conv else coor3('o', coor_dest)






# 좌표/좌표계 쌍
class coor_set:
  def __init__(self, coor , reffrm):
    self.coor=coor
    self.reffrm=reffrm

  # abs 좌표계의 coor_set
  def abs_conv_coorset(self, cls):
    return cls(self.reffrm.abs_conv(self.coor, 'sa'), reffrm())

  # 바뀐 좌표계의 coor_set
  def conv_coorset(self, cls, target_reffrm): 
    return cls(self.reffrm.coor_conv(target_reffrm, self.coor), target_reffrm)

  # 여러 reffrm간 합성에 대한 coor 변환
  def composit_reffrm_coor_conv(cls, coorS_source, conv_mode, *reffrm):
    convM = [reff.rot_convM_O for reff in reffrm]
    return cls.__composit_convM_coor_conv(cls, coorS_source, conv_mode, *convM)


##################### print str ####################

  def __str__(self):
    return 'coor_set info:\n  coor: {}\n  reference frame: {}'.format(self.coor, self.reffrm)

################# 연산자 오버라이딩 ##################

  # abs 좌표계의 coor_set (~coor_set)
  def __invert__(self):
    return self.abs_conv_coorset()

  # abs 좌표계상 같은 좌표인지 (cs1==cs2)
  def __eq__(self, other):
    a=self.abs_conv_coorset
    b=other.abs_conv_coorset
    return (a.coor==b.coor).all()

  # self좌표를 other의 좌표계 표현으로 변환 (cs1>>cs2)
  def __rshift__(self, other):
    return self.conv_coorset(other.reffrm)



