import numpy as np

omega = 7.29211 * 10 ** (-5)
R = 6371302
mu = 398600.4415 * 10 ** 9
h = mu ** (1/3) * (3 / 46 / omega) ** (2/3) - R
v = np.sqrt(mu / (R + h))

def longitude(input):
    full_circle = np.floor(input / 360)
    answer = input - (full_circle * 360)
    return answer

n = 3  # число спутников на орбите
# вероятно, решение при 3-ёх спутниках не единственно, с фазой 0,8 работает, 2 спутника недостаточно
phase = 0.8  # расфазировка спутников

# рассчитываем долготу начального пролёта спутников через экватор
start = [i / n * phase * (R + h) / v * omega * 360 for i in range(n)]
ans = [np.array(start)]
# время прохождения каждым спутником экватора
T_start = [i / n * phase * 2 * np.pi * (R + h) / v / 60 for i in range(n)]
T = [np.array(T_start)]
semi_T = np.pi * (R + h) / v / 60

# 46 полувитка, расчитываем долготу прохождения восходящего и нисходящего узла
for i in range(46 * 2):
    new_points = [180 * (1 - omega * (R + h) / v) + ans[-1][i] for i in range(n)]
    ans.append(longitude(np.array(new_points)))
    T.append(T[-1] + semi_T)  # time, minutes


ans = np.array(ans)
# T = np.array(T)
# for i in range(ans.shape[0]):
#     ans[i][0] = np.mod(ans[i][0] + 180, 360) - 180
print(ans)
# print(T)
points = ans.reshape(-1) * 2 * np.pi * R / 360  # точки прохождения, км
points = np.sort(points) / 1000
diff = np.diff(points)
maxim = np.abs(diff).max()
print(max(maxim, np.abs(points[-1] - 40030.173)))  # максимальная разница по экватору между точками, км