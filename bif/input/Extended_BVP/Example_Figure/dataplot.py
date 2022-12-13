import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from PIL import Image
from matplotlib.patches import ArrowStyle

# データ読み取り用
def setdata(filename, xcol, ycol):
    x_list = []
    y_list = []
    fd = open(filename, "rt")  # specify appropriate data file here
    for line in fd:
        data = line[:-1].split(" ")
        x_list.append(float(data[xcol]))
        y_list.append(float(data[ycol]))
    return x_list, y_list


# 版組パラメタ
params = {
    "text.usetex": True,
    "text.latex.preamble": r"\usepackage{newtxtext,newtxmath}",
    "legend.fontsize": 12,
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "font.family": "serif",
    "grid.color": "k",
    "grid.linestyle": ":",
    "grid.linewidth": 0.5,
}
plt.rcParams.update(params)

# 図の大きさと比
fig = plt.figure(figsize=(5, 5))
ax = fig.add_subplot(111)
# ax.set_aspect('equal', 'datalim')

# plt.xticks(fontsize=16)
# plt.yticks(fontsize=16)

# 軸回り
x_limit = [0, 1.5]
y_limit = [0, 2]
x_tick = 0.5
y_tick = 0.5
ax.set_xlim(x_limit)
ax.set_ylim(y_limit)
plt.xticks(np.arange(x_limit[0], x_limit[1] + x_tick, x_tick), fontsize=16)
plt.yticks(np.arange(y_limit[0], y_limit[1] + y_tick, y_tick), fontsize=16)
ax.set_xlabel(r"$A \longrightarrow$", fontsize=16)
ax.set_ylabel(r"$B \longrightarrow $", fontsize=16)

# 格子
ax.grid(c="gainsboro", zorder=2)

# データ
# 周期解
(x_list, y_list) = setdata("NS1", 0, 1)
ax.plot(x_list, y_list, color="BLACK", linewidth=2)
(x_list, y_list) = setdata("NS2", 0, 1)
ax.plot(x_list, y_list, color="BLACK", linewidth=2)
(x_list, y_list) = setdata("5G1", 0, 1)
ax.plot(x_list, y_list, color="LIGHTGRAY", linewidth=0.75)
(x_list, y_list) = setdata("5G2", 0, 1)
ax.plot(x_list, y_list, color="LIGHTGRAY", linewidth=0.75)
(x_list, y_list) = setdata("5PD1", 0, 1)
ax.plot(x_list, y_list, color="LIGHTGRAY", linewidth=0.75)
(x_list, y_list) = setdata("5PD2", 0, 1)
ax.plot(x_list, y_list, color="LIGHTGRAY", linewidth=0.75)
(x_list, y_list) = setdata("1G1", 0, 1)
ax.plot(x_list, y_list, color="BLACK", linewidth=2)
(x_list, y_list) = setdata("1G2", 0, 1)
ax.plot(x_list, y_list, color="BLACK", linewidth=2)
(x_list, y_list) = setdata("2G1", 0, 1)
ax.plot(x_list, y_list, color="BLACK", linewidth=2)
(x_list, y_list) = setdata("2G2", 0, 1)
ax.plot(x_list, y_list, color="BLACK", linewidth=2)

# 平衡点
(x_list, y_list) = setdata("1EH1", 0, 1)
ax.plot(x_list, y_list, color="BLACK", linewidth=2)
(x_list, y_list) = setdata("1EH2", 0, 1)
ax.plot(x_list, y_list, color="BLACK", linewidth=2)
(x_list, y_list) = setdata("3EH1", 0, 1)
ax.plot(x_list, y_list, color="BLACK", linewidth=2)
(x_list, y_list) = setdata("3EH2", 0, 1)
ax.plot(x_list, y_list, color="BLACK", linewidth=2)
(x_list, y_list) = setdata("2EH1", 0, 1)
ax.plot(x_list, y_list, color="BLACK", linewidth=2)
(x_list, y_list) = setdata("2EH2", 0, 1)
ax.plot(x_list, y_list, color="BLACK", linewidth=2)
(x_list, y_list) = setdata("1EG1", 0, 1)
ax.plot(x_list, y_list, color="BLACK", linewidth=2)
(x_list, y_list) = setdata("1EG2", 0, 1)
ax.plot(x_list, y_list, color="BLACK", linewidth=2)

# テキスト
# ax.text(1.2, 0.035, r"Hopf Bifurcation", color="Black", fontsize="12", rotation=-49)
# ax.text(0.165, 0.15, r"$F_2$", color="Blue", fontsize="14", rotation=90)
# ax.text(1.44, -0.00475, r"$F_1$", color="Blue", fontsize="14", rotation=0)

# 矢印
# ax.annotate('', xy=(0.7, 0.07), xytext=(0.875, 0.15),
#                 arrowprops=dict(shrink=0, width=1, headwidth=8,
#                                 headlength=10, connectionstyle='arc3',
#                                 facecolor='gray', edgecolor='gray')
#                )
# ax.text(1.1, 0.024, r"Periodic Solution", color="gray", fontsize="14", rotation=0)
# ax.text(1.3925, 0.055, r"Stable Equilibrium", color="gray", fontsize="12", rotation=0)


# 背景画像
# xlim = ax.get_xlim()
# ylim = ax.get_ylim()
# im = Image.open("./inzhang3.png")
# ax.imshow(im, extent=[*xlim, *ylim], aspect='auto', resample=None)

# 出力
fig.subplots_adjust(left=0.175, right=0.94, bottom=0.135, top=0.975)
pdf = PdfPages("snapshot.pdf")
pdf.savefig(dpi=300)
pdf.close()
