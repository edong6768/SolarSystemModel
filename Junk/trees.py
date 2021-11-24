class tree:

  def __init__(self, data, name='', parent=None, childmode=True):
    self.data=data  # any desirable object
    self.name=name
    self.size = 1
    self.set_parent(parent)
    self.children=set() if childmode else []

    self.current = 0       # iterator 선언에 필요, tree 전체 구조 scan이므로 정지조건은 sub_size가 된다.
    self.iter_children=None   # 현재 sub_children의 iter
    self.iter_curr_child=None   # 현재 search중인 child의 iter


  @property
  def depth(self): return self.__root_search()[0]

  @property
  def root(self): return self.__root_search()[1]


  # 자신의 기저기준틀을 지정한다 (즉 이는 자신의 기저기준틀이 절대기준틀이 아니며, 이 또한 다른 기준틀을 기저로 가진다는 것을 의미한다)
  # 또한 지정된 기저기준틀의 하위기준틀 set에 자신을 추가한다.
  def set_parent(self, parent):
    if self==None:
      raise ValueError
    if parent!=None:
      parent.__add_children(self)
      parent.__update_sub_size(self.size) # 기저기준틀의 sub_size update
    self.parent=parent

  def set_children(self, *children):
    for child in children:
      child.set_parent(self)

  # 자신의 tree node depth(깊이)와 절대기준틀 obj 반환
  def __root_search(self, n=-1):
    n+=1
    return self.parent.__root_search(n) if self.parent!=None else (n, self)

  def __update_sub_size(self, subs_size):
    self.size+=subs_size+1
    if self.parent: self.parent.__update_sub_size(subs_size)

  def __add_children(self, child):
    self.children.add(child)
  

  ##################### magic methods ####################

  # 트리 구조도 출력
  def __str__(self):
    return self.treeplot()

  # 트리 전체 내부 행렬 정보와 관계 
  def __repr__(self):
    # return '\nname : {}\nbase reffrm : {}\ncharacteristic matrix : \n{} \nsub reffrm : {}'.format(self.name, self.base_reffrm.name if (self.base_reffrm!=None) else 'Empty', self.convM_O['as'], [i.name for i in self.sub_reffrm] if self.sub_reffrm else 'Empty')
    return self.name

  # 전체트리크기(node수)
  def __len__(self):
    return self.size

  # sub_reffrm의 하부 포함여부 (a in b, 트리 하부구조에 node가 존재하는지 여부 확인) 사용을 위한 iterable 선언 
  def __next__(self):
    if self.current < self.size:
      try:
        if self.current==0:
          self.current+=1
          return self
        elif self.current==1: 
          self.iter_curr_child=iter(next(self.iter_children))
        self.current+=1
        return next(self.iter_curr_child)
      except StopIteration:
        self.iter_curr_child=iter(next(self.iter_children))
        return next(self.iter_curr_child)
    else:
      raise StopIteration

  # for iteratability
  def __iter__(self):
    self.current = 0
    self.iter_children=iter(self.children)
    return self

  # 포함 관계 a in tree를 사용하기 위해
  def __contains__(self, other):
    if other==None: return True
    return True if other.data in [i.data for i in iter(self)] else False

  # 같은 node = 같은 데이터 & 같은 부모
  def __eq__(self, other):
    if other==None:
      return False
    elif self.parent==None or other.parent==None:
      return self.data==other.data
    return (self.parent==other.parent) and (self.data==other.data)

  # 두 기준틀이 같은지를 파이썬이 판단하기 위해서 obj의 hash를 이용한다.
  # 여기서 coorM_O의 각 요소(실수, float)를 floating point binary string으로 바꾸어준 뒤 모두 concatinate하여 int()하여준다
  def __hash__(self):
    return hash(self.data) + (hash(self.parent.data) if self.parent!=None else 0)

  ########### Tree construction / export / print ###########

  
  # [[A, [ [B,[E, F]], [C,[G]] [D,[]] ] ]]--> A - {B, C, D} - {{E, F}, {G}, {}}
  @classmethod
  def l2t(cls, tree_list, base=None, n=0):
    for refs in tree_list:
      if len(refs)==2 and type(refs[1])==type(list()):
        refs[0].set_parent(base)
        cls.l2t(refs[1], refs[0], 1)
      else: refs.set_parent(base)
    if n==0 : return [roots[0] for roots in tree_list]

  # tree list가 각각의 이름으로 되어있을 때, '이름':reffrm dict와 함께 주어진 것을 tree구조 연결에 사용(재귀)
  # [['A', [ ['B',['E', 'F']], ['C',['G']] ['D',[]] ] ]]--> A - {B, C, D} - {{E, F}, {G}, {}}
  @classmethod
  def dict_l2t(cls, tree_name_list, node_dict, base=''):
    base=node_dict[base] if base else None
    if not base: tree_name_list=[tree_name_list]
    for r_name in tree_name_list:
      if len(r_name)==2 and type(r_name[1])==type(list()):
        node_dict[r_name[0]].set_parent(base)
        cls.dict_l2t(r_name[1], node_dict, r_name[0])
        if not base: return node_dict[r_name[0]]
      else: node_dict[r_name].set_parent(base)


  # A - {B, C, D} - {{E, F}, {G}, {}} --> [['A', [ ['B',['E', 'F']], ['C',['G']] ['D',[]] ] ](재귀)
  def t2l(self):
    if self.children: return [self.name, [sub.t2l() for sub in self.children]]
    else: return self.name

  # '이름' : node 딕셔너리를 트리 전체 node에 대하여 구성한다.
  @property
  def dictify(self):
    node_dict = dict()
    for node in self:
      node_dict[node.name] = node
    return node_dict

  # 트리 구조도
  def treeplot(self, front_str='', infront_str=''):
    branch_str = ['├──', '└──']
    space_str = { '': '', '├──': '│  ', '└──' :  '   '}
    tree_str = front_str + infront_str + self.name
    if self.children:
      l_child = len(self.children)
      front_str += space_str[infront_str]
      for i, child in enumerate(self.children):
        infront_str = branch_str[i==(l_child-1)]
        tree_str += '\n' + child.treeplot(front_str, infront_str)
    return tree_str


  # 리스트 형태 트리의 plot 출력
  @classmethod
  def treelist_plot(cls, tree_list, front_str='', infront_str=''):
    branch_str = ['├──', '└──']
    space_str = { '': '', '├──': '│  ', '└──' :  '   '}
    tree_str = front_str + infront_str + tree_list[0]
    if len(tree_list)==2 and type(tree_list[1])==type(list()):
      listsize = len(tree_list[1])
      front_str += space_str[infront_str]
      for i in range(listsize):
        infront_str = branch_str[i==(listsize-1)]
        tree_str += '\n' + cls.treeplot(tree_list[1], front_str, infront_str)
    return tree_str

