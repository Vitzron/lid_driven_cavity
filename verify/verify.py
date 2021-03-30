import numpy as np
import matplotlib.pyplot as plt
import sys

X_data = np.loadtxt(sys.argv[1], skiprows=0)
Y_data = np.loadtxt(sys.argv[2], skiprows=0)
u_data = np.loadtxt(sys.argv[3], skiprows=0)
v_data = np.loadtxt(sys.argv[4], skiprows=0)
Re     = np.int(sys.argv[5])

if len(sys.argv) > 6:
    fname_fig = sys.argv[6]

Re_lst = [100, 400, 1000, 3200, 5000, 7500, 10000]
if Re not in Re_lst:
    print('WRONG Re number!')
    sys.exit()
col = Re_lst.index(Re) + 1
Ghia_y = np.loadtxt('Ghia1982_u.dat', skiprows=3, usecols=0)
Ghia_x = np.loadtxt('Ghia1982_v.dat', skiprows=3, usecols=0)
Ghia_u = np.loadtxt('Ghia1982_u.dat', skiprows=3, usecols=col)
Ghia_v = np.loadtxt('Ghia1982_v.dat', skiprows=3, usecols=col)

rows = v_data.shape[0]
cols = u_data.shape[1]

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
axes[0].set_ylim(0.0, 1.0)
axes[0].scatter(Ghia_u, Ghia_y, marker='.', linewidth=1.5, color='k')
axes[0].plot(u_data[:, int(cols / 2)], Y_data[:, 1], linewidth=1.0, color='k')
axes[0].set_xlabel(r'$u$', fontsize=15)
axes[0].set_ylabel(r'$y$', fontsize=15)
axes[0].set_title('$u$-velocity')
axes[1].set_xlim(0.0, 1.0)
axes[1].scatter(Ghia_x, Ghia_v, marker='.', linewidth=1.5, color='k')
axes[1].plot(X_data[1, :], v_data[int(rows / 2), :], linewidth=1.0, color='k')
axes[1].set_xlabel(r'$x$', fontsize=15)
axes[1].set_ylabel(r'$v$', fontsize=15)
axes[1].set_title('$v$-velocity')

def close_event():
    plt.close()

timer = fig.canvas.new_timer(interval=2000)
timer.add_callback(close_event)

if len(sys.argv) > 6:
    plt.savefig(fname_fig, dip=300)

timer.start()
plt.show()