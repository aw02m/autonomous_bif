import sys
import json
import numpy as np
from numpy import cos
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
# import poincare_section as psec
from scipy import linalg


def q(t, v, data):
    return v[data.dict['p_index']]-data.dict['p_place']


q.terminal = True
q.direction = -1

# Rossler
def func(t, x, data):
    return np.array([
        [-x[1][0]-x[2][0]],
        [x[0][0]+data.dict['params'][0]*x[1][0]],
        [data.dict['params'][1]*x[0][0] - data.dict['params'][2]*x[2][0]+x[0][0]*x[2][0]]]
    )

# BVP
# def func(t, x, data):
# 	return np.array([
# 		[-x[2][0] + data.dict['params'][0]*x[0][0] + np.tanh(data.dict['params'][1]*x[0][0])],
# 		[x[2][0] - data.dict['params'][2]*x[1][0]],
# 		[x[0][0] - x[1][0]]
# 	])


class DataStruct():
    def __init__(self):
        if len(sys.argv) != 2:
            print(f"Usage: python {sys.argv[0]} filename")
            sys.exit(0)
        fd = open(sys.argv[1], 'r')
        self.dict = json.load(fd)
        fd.close()
        self.param_ptr = 0
        self.ax = None
        self.fig = None


def keyin(event, s, data, tau, fixed):
    ptr = data.param_ptr
    if event.key == 'q':
        plt.close('all')
        print("quit")
        sys.exit()
    elif event.key == 'w':
        # data.dict['fixed'] = [[] if i == [] else i.ravel().tolist() for i in fixed]
        data.dict['fixed'] = list(fixed)
        data.dict['tau'] = tau
        jd = json.dumps(data.dict)
        # print(jd)
        with open("__ppout__.json", 'w') as fd:
            json.dump(data.dict, fd, indent=4)
        print("wrote a file")
    elif event.key == ' ' or event.key == 'e':
        plt.cla()
        initial_setup(data)
    elif event.key == 'f':
        plt.cla()
        initial_setup(data)
        data.visual_orbit = 1 - data.visual_orbit
    elif event.key == 's':
        for i in data.dict['params']:
            print(i, end=' ')
        print(s[0], s[1])
    elif event.key == 'p':
        data.param_ptr += 1
        if data.param_ptr >= len(data.dict['params']):
            data.param_ptr = 0
        print(f"changable parameter: {data.param_ptr}")
    elif event.key == 'up':
        ptr = data.param_ptr
        data.dict['params'][ptr] += data.dict['dparams'][ptr]
    elif event.key == 'down':
        ptr = data.param_ptr
        data.dict['params'][ptr] -= data.dict['dparams'][ptr]
    show_param(data)


def show_param(data):
    s = ""
    cnt = 0
    for key in data.dict['params']:
        s += " param{:d}: {:.5f}  ".format(cnt, key)
        cnt += 1
    plt.title(s, color='b')


def on_click(event, s0, lines, data):
    x = data.dict['axis'][0]
    y = data.dict['axis'][1]
    s0[x] = event.xdata
    s0[y] = event.ydata
    lines.set_data(s0[x], s0[y])
    plt.plot(s0[x], s0[y], 'o', markersize=2, color="red")
    print(s0[x], s0[y])
    initial_setup(data)
    show_param(data)


def on_close():
    running = False


def initial_setup(data):
    xr = data.dict['xrange']
    yr = data.dict['yrange']
    data.ax.set_xlim(xr[0], xr[1])
    data.ax.set_ylim(yr[0], yr[1])
    data.ax.set_xlabel('axis[0]')
    data.ax.set_ylabel('axis[1]')


def main():
    plt.rcParams['keymap.save'].remove('s')
    plt.rcParams['keymap.quit'].remove('q')
    data = DataStruct()

    data.fig = plt.figure(figsize=(10, 10))
    data.ax = data.fig.add_subplot(111)

    state0 = data.dict['x0']
    tick = data.dict['tick']

    initial_setup(data)
    data.visual_orbit = 1

    tau = 0
    fixed = []

    plt.connect('button_press_event',
                lambda event: on_click(event, state0, lines, data))
    plt.connect('key_press_event',
                lambda event: keyin(event, state0, data, tau, fixed))
    plt.ion()  # I/O non blocking

    duration = 1
    #tspan = np.arange(0, duration, tick)

    period_buf = 0.0

    running = True

    while running:
        state = solve_ivp(func, (0, duration), state0,
                          method='RK45', args=(data,),
                          events=q, max_step=tick, vectorized=True)
        # reset termination and direction settings
        # for i in range(len(poincare)):
        # 	poincare[i].terminal = True
        # 	poincare[i].direction = direction[i]

        if state.status == 0:
            period_buf += duration
        # and abs(np.linalg.norm(state0) - np.linalg.norm(state.y[:,-1])) > 0.001:
        elif state.status == 1:
            # hit_section = [x for x, y in enumerate(state.t_events) if y.size > 0][0]
            tau = period_buf + state.t_events[0][0]
            fixed = state.y_events[0][0]

            # poincare[hit_section[0]].terminal = False # avoid termination in next step
            # poincare[hit_section[0]].direction = -poincare[hit_section[0]].direction
            # if hit_section == 0:
            # 	fixed = state.y_events[0][0]

            sys.stdout.write("%s\n" % fixed)
            sys.stdout.write("%s" % tau)
            sys.stdout.flush()
            if data.visual_orbit == 1:
                lines, = plt.plot(state.y[data.dict['axis'][0], :], state.y[data.dict['axis'][1], :],
                                  linewidth=0.5, color=(0, 0, 0),
                                  ls="-")
            plt.plot(state.y_events[0][:, data.dict['axis'][0]],
                     state.y_events[0][:, data.dict['axis'][1]],
                     'o', markersize=2, color="red")
            state = solve_ivp(func, (0, tick), state.y[:, -1],
                              method='RK45', args=(data,),
                              max_step=tick, vectorized=True)
            period_buf = tick
        else:
            plt.close('all')
            print("Integration failed.")
            sys.exit()
        # if any((x.size > 0 for x in state.t_events)):
        # 	for i in state.y_events:
        # 		if len(i) != 0:
        # 			plt.plot(i[:, 0], i[:, 1],
        # 					 'o', markersize = 2, color="red")
        # 	hit_section = [x for x, y in enumerate(state.t_events) if y.size > 0]
        # 	period[hit_section[0]] = period_buf - (duration - state.t_events[hit_section[0]][0])
        # 	if len(hit_section) > 1:
        # 		print(hit_section)
        # 		for i in hit_section[1:]:
        # 			period[i] = state.t_events[i][0] - state.t_events[i-1][0]
        # 			period_buf = duration - state.t_events[i][0]
        if data.visual_orbit == 1:
            lines, = plt.plot(state.y[data.dict['axis'][0], :], state.y[data.dict['axis'][1], :],
                              linewidth=0.5, color=(0, 0, 0),
                              ls="-")
        state0 = state.y[:, -1]
        plt.pause(0.001)  # REQIRED


if __name__ == '__main__':
    main()
