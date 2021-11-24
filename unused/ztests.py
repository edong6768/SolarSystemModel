#import celestials as c
import coor_ref as cr
import numpy as np
import celestials as cel
from trees import tree as tr



############################### tree test ###############################

a=tr(2, '2')
b=tr(3, '3')
c=tr(4, '4')
d=tr(5, '5')
e=tr(6, '6')

a.set_children(b, c)
b.set_children(d)
d.set_children(e)

print(a)
print(a.t2l())
print(a.dictify)

############################### coor3 test ###############################
b=np.array(((1, 2, 3), (3, 4, 5), (3, 4, 5))*3).T
dc=cr.coor3('o', b, (5, 6, 7), (3, 4, 5), (1, 2, 3))


acd=np.array([[2, 0, 0 ], [2, 0, 0], [1, 0, 0]])
#print(acd)
p=cr.coor3('o', acd)
#print(cc.coor3('o', acd))

#print(b.base_conv(p, 'sa'))

# full_coor=np.empty((0, 0), float)

# print(full_coor)
# j=np.array((1, 2, 3))

# k=np.column_stack([j, j])

# print(k)









############################### reffrm test ###############################

a=cr.reffrm(1, 2, 3, 4, 5, 6, 'A')
b=cr.reffrm(1, 2, 3, 4, 5, 7, 'B')
bed=cr.reffrm(1, 2, 3, 4, 5, 3, 'C')
deb=cr.reffrm(1, 3, 3, 4, 5, 3, 'D')
asy=cr.reffrm(1, 2, 3, 4,9, 10, 'E')
f=cr.reffrm(1, 2, 3, 4,4, 10, 'F')


tree_list = [[a, [ [b,[deb]], [bed, [[asy, [f]]]] ] ]]
k=cr.reffrm.l2t(tree_list)
# print(a)
# print(a.dictify)
# print(a.t2l())

#print(asdf[0], ' ,', asdf[1].name)
#print(a, b)

# a.reffrm_tree_plot()

# hello
#     goodbye
#     bed
#     hint
#현재는 goodbye까지 나옴








############################### celestials test ###############################


cel.solar_system_model.config_JSON_load('data/SEM_config_snapshot_0.json')