import numpy as np

start = [0.0, 120.0, 240.0]
omega = 7.29211 * 10 ** (-5)
R = 6371302
mu = 398600.4415 * 10 ** 9
# оптимизировав орбиту одного спутника, получаем, что 3 - ёх хватит
h = mu ** (1/3) * (3 / 46 / omega) ** (2/3) - R
v = np.sqrt(mu / (R + h))

def longitude(input):
    full_circle = np.floor(input / 360)
    answer = input - (full_circle * 360)
    return answer

ans = [np.array(start)]
T = [0]

for i in range(46 * 2):
    new_point = 180 * (1 - omega * (R + h) / v) + ans[-1][0]
    ans.append(longitude(np.array([new_point, new_point + ans[0][1], new_point + ans[0][2]])))
    T.append(T[-1] + np.pi * (R + h) / v / 60)  # time, minutes

delta = 360 * omega * (R + h) / v  # запаздывание за виток в градусах
# в первом приближении период 1,5 часа, хотим полный период витков <= 100 часам,
# ближе всего 66 витков за 99 часов. delta хотим целое число * 360, ближайшее целое 4 (с delta_init 4,15)
# 66 витков много, орбита низкая; подбором 62 витка
# так как 4 - чётное, минимальный период 2 и мы неоптимально используем 1 спутник, нужно 3 полных оборота
h_opt = mu ** (1/3) * (3 / 46 / omega) ** (2/3) - R  # это то же h сверху, просто сверху подбиралось h
# print(delta)

ans = np.array(ans)
print("Времена прохождения восход/нисход узлов, минуты", T)
print("Высота орбиты:", h_opt)
points = ans.reshape(-1) * 2 * np.pi * R / 360  # точки прохождения, км
points = np.sort(points) / 1000
diff = np.diff(points)
maxim = np.abs(diff).max()
print("Максимальная разница между соседними точками, км:",
      max(maxim, np.abs(points[-1] - 40030.173)))

# нормализация долготы в связанной с Землёй неИСО, + восточная долгота, - западная долгота
for i in range(ans.shape[0]):
    for j in range(3):
        ans[i][j] = np.mod(ans[i][j] + 180, 360) - 180

print("Три столбика для 3-ёх спутников с прохождениями восходящего/нисходящего узлов в градусах:")
print(ans)