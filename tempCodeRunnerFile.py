fig=plt.figure()
ax=fig.add_subplot(projection='3d')
a.timetravel_animation(fig, ax, 3600*24)
plt.show()