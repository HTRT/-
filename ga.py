import pandas as pd
import numpy as np
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import rc
import random

rc('font', family='YouYuan')  # 设置中文

GENERATIONS = 200  # 迭代次数
DNA_SIZE = 22  # DNA 链的长度
MUTATION_RATE = 0.1  # 基因突变概率
DOUBLE_CHILD_RATE = 0.5  # 二胎概率
dna_decode = []  # 种群DNA 解码列表
pop = []  # 种群DNA列表


# 种群生存环境(需要优化的函数)
def function_env_Ackley(x, y):
    if type(x) == list:
        x = np.array(x)
    if type(y) == list:
        y = np.array(y)
    a = 20
    b = 0.2
    c = 2 * np.pi
    sum1 = (x ** 2 + y ** 2)
    sum2 = np.cos(c * x) + np.cos(c * y)

    term1 = - a * np.exp(-b * np.sqrt(sum1 / 2))
    term2 = -np.exp(sum2 / 2)
    z = term1 + term2 + a + np.exp(1)
    return z


# 图像显示
def show_Ackley(fig, x1, x2, i, n):
    plt.clf()
    X, Y = np.mgrid[-40:40:35j, -40:40:35j]

    Z = function_env_Ackley(X, Y)  # Ackley函数图像
    Z1 = function_env_Ackley(x1, x2)  # 在图像上显示种群
    print(f"种群所在位置:{Z1}")

    # 图像显示
    ax = Axes3D(fig)
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='cool', alpha=0.4)  # alpha 设置透明度
    plt.suptitle(f'遗传算法：初始种群数量为{n}, 已进化到第{i + 1}代, ', c='blue')
    ax.scatter(x1, x2, Z1, color='red', alpha=1)
    # ax.plot3D(x1, x2, Z1, color='red', alpha=1)
    ax.set_xlabel('x1', c='blue')
    ax.set_ylabel('x2', c='blue')
    ax.set_zlabel('种群位置', c='blue')

    plt.pause(0.01)
    return Z1


# 初始化种群数量　
def initial_population(pop_num: int, DNA_SIZE):
    for i in range(pop_num):
        pop_dna = np.random.randint(0, 2, size=DNA_SIZE * 2)  # 随机生成种群个体的 DNA　序列
        pop.append(pop_dna)
    # print(pop)
    return pop  # 返回种群列表


# DNA解码
def decode_dna(dna):
    x1_t = 0
    x2_t = 0
    # list_xx = []

    dna1 = dna[0:DNA_SIZE]  # x1
    dna2 = dna[DNA_SIZE:]  # x2
    # 解码过程
    for i, data in enumerate(dna1):
        if data == 1:
            x1_t += 2 ** i
    # 解码过程
    for j, data in enumerate(dna2):
        if data == 1:
            x2_t += 2 ** j

    # 把 DNA　解码值约束到（-32, 32）
    x1 = x1_t * 2 / ((2 ** DNA_SIZE) - 1) - 1
    x2 = x2_t * 2 / ((2 ** DNA_SIZE) - 1) - 1

    x1 = x1 * 32
    x2 = x2 * 32
    x1 = round(x1, 4)  # 保留小数点后4位
    x2 = round(x2, 4)

    transDNA1.append(x1)
    transDNA2.append(x2)
    return transDNA1, transDNA2  # 返回个体 DNA 解码列表1, 解码列表2


# 基因突变
def mutation_dna(dna_2):
    if np.random.rand() < MUTATION_RATE:
        m = random.randint(1, 5)  # 基因突变的长度
        print("基因突变的编码个数：", m)
        m1 = random.sample(range(0, DNA_SIZE), m)  # 基因突变的位置
        for i in m1:
            if dna_2[i] == 1:
                dna_2[i] = 0
            else:
                dna_2[i] = 1
        print("子代发生突变！")
    return dna_2


# 基因交叉
def cross_dna(father, mother):
    father_arr = np.array(father[0:DNA_SIZE])
    mother_arr = np.array(mother[DNA_SIZE:])
    c1 = father[0:DNA_SIZE]  # 继承父亲前端基因
    c2 = np.bitwise_and(father_arr, mother_arr)  # 后端基因是父母基因的按位相与的结果
    c1 = list(c1)
    c2 = list(c2)
    child_dna = c1 + c2
    return child_dna


# 适应度
def get_fitness(x1, x2):
    z = function_env_Ackley(x1, x2)
    fitness = (z - np.min(z)) + 1e-4
    fitness = fitness / fitness.sum()
    fitness = -np.round(fitness, 4) + 1  # z值越大，适应度越小
    print(f"适应度:{fitness}")
    return fitness


# 配对
def make_pair(pop):
    father = []
    mother = []
    n = random.randrange(2, len(pop) + 1, 2)  # 配对个数
    idx = random.sample(range(0, len(pop)), n)  # 选择配对个体
    for i, data in enumerate(idx):
        if i % 2 == 0:
            father.append(pop[data])
        else:
            mother.append(pop[data])
    return father, mother


# 生成子代
def get_child(father, mother):
    child_list = []
    for fath, moth in zip(father, mother):
        child_dna = cross_dna(fath, moth)  # 生成子代
        child_dna = mutation_dna(child_dna)  # 子代有一定几率发生突变
        child_list.append(child_dna)
        if np.random.rand() > DOUBLE_CHILD_RATE:
            child_2_dna = cross_dna(moth, fath)
            child_list.append(child_2_dna)
            # print("已产生二胎")
    return child_list


# 选择
def select_pop(fitness, num, transDNA1, transDNA2, pop_num):
    sort_idx = np.argsort(-fitness)  # 按降序排列，输出原索引值
    sort_idx = list(sort_idx)  # 数组转列表
    sort_idx = sort_idx[:num]  # 取初始大小的种群存活，过滤掉适应度较小的个体

    select_dna1 = []
    select_dna2 = []
    select_pop_num = []
    for idx in sort_idx:
        select_dna1.append(transDNA1[idx])
        select_dna2.append(transDNA2[idx])
        select_pop_num.append(pop_num[idx])
    transDNA1.clear()
    transDNA2.clear()
    pop_num.clear()
    transDNA1 = select_dna1
    transDNA2 = select_dna2
    pop_num = select_pop_num

    return transDNA1, transDNA2, pop_num


# 把数据保存到csv
def save_csv(idx, dna1, dna2, fit):
    if (idx + 1) % 10 == 0:
        name = ['x1', "x2", "z", "fitness"]
        list_csv = zip(dna1, dna2, Z1, list(fit))
        test = pd.DataFrame(columns=name, data=list_csv)
        test.to_csv(f'G:/ackley/第{i + 1}代.csv', encoding="gbk")


# 随机生成颜色
def randomcolor():
    colArr = ['1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F']
    color = ""
    for i in range(6):
        color += colArr[random.randint(0, 14)]
    return "#" + color


# 生成种群位置随代数增加而变化的折线图
def graph_line(generation_num, z, population_num):
    fig = plt.figure()
    plt.title("种群位置随迭代次数变化图")
    plt.xlabel('迭代次数')
    plt.ylabel("种群位置")
    for x, y, n in zip(generation_num, z, population_num):
        plt.plot(x, y, color=randomcolor(), label=f'初始种群数量:{n}')
    plt.legend()


if __name__ == "__main__":
    # 生成种群
    # num = 50  # 初始种群数量
    num_list = [30, 50, 100, 150, 200]

    g_num = []
    z_m = []

    for num in num_list:
        pop_num = initial_population(num, DNA_SIZE)

        # 开启交互模式
        plt.ion()
        fig = plt.figure()

        generation_num_list = []
        z_mean = []
        for i in range(0, GENERATIONS):
            transDNA1 = []  # dna序列1
            transDNA2 = []  # dna序列2
            print(f"种群已繁衍到第{i + 1}代。")
            print(f"当前种群数量：{len(pop_num)}")

            father, mother = make_pair(pop_num)  # 生成父母
            child_list = get_child(father, mother)

            pop_num = pop_num + child_list  # 子代和父代（二进制）
            # pop_num_copy = pop_num.copy()
            # 种群dna解码
            for dna in pop_num:
                transDNA1, transDNA2 = decode_dna(dna)
            print(f"选择前的种群数量:{len(pop_num)}")
            # show_Ackley(transDNA1, transDNA2)

            # 适应度
            fitness = get_fitness(transDNA1, transDNA2)
            transDNA1, transDNA2, pop_num = select_pop(fitness, num, transDNA1, transDNA2, pop_num)
            # 选择

            # 选择后的种群
            print(f"选择后的种群{len(pop_num)}")

            print("解码列表x1", transDNA1)  # x1
            print("解码列表x2", transDNA2)  # x2
            Z1 = show_Ackley(fig, transDNA1, transDNA2, i, len(pop_num))
            print("\n")

            generation_num_list.append(i)
            z_mean.append(np.mean(Z1))

        g_num.append(generation_num_list)
        z_m.append(z_mean)

    graph_line(g_num, z_m, num_list)

    # 关闭交互模式
    plt.ioff()
    plt.show()