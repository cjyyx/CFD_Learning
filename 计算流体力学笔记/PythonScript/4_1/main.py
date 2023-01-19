# %%
from typing import List

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

# %%
# 设置参数
gif_name="E1"

f_E: float = 1
f_R_e: float = 5000.0
f_Time: float = 3200.0  # 运行的总时间
i_N: int = 20


f_delta_y: float = 1/i_N
f_delta_t: float = f_E*f_R_e*(f_delta_y*f_delta_y)

# u[n,j]的数据结构: list[np.array(N+1)]
array_u: List = []

# %%
# 初值条件
vec_temp = np.zeros(i_N+1)
array_u.append(vec_temp)

# %%


def u(n: int, j: int) -> float:
    """ 
    从array_u中读取[n,j]处的值
     """
    return array_u[n][j-1]


def K(n: int, j: int) -> float:
    """ 
    计算K_j^n
     """
    result = (1-f_E)*u(n, j)+(f_E/2)*(u(n, j+1)+u(n, j-1))
    return result


f_A = -f_E/2
f_B = 1+f_E

vec_temp1 = np.ones(i_N-2)
vec_temp2 = np.ones(i_N-1)

matrix_T = (
    f_A*np.diagflat(vec_temp1, 1) +
    f_B * np.diagflat(vec_temp2, 0) +
    f_A*np.diagflat(vec_temp1, -1)
)
matrix_S = np.linalg.inv(matrix_T)


def next_step():
    """ 
    计算一个时间步
     """
    i_n: int = len(array_u)-1
    vec_K = np.array(list(map(
        lambda j: K(i_n, j) if j != i_N else K(i_n, j)-f_A,
        list(range(2, i_N+1))
    )))

    vec_temp = np.dot(matrix_S, vec_K)
    vec_temp = np.insert(vec_temp, 0, 0)
    vec_temp = np.append(vec_temp, 1)
    array_u.append(vec_temp)


# %%

for i in tqdm(range(int(f_Time/f_delta_t))):
    next_step()

step_num = len(array_u)
print("共{}个时间步".format(step_num))

# %%

gif_time = 2
gif_frame_count = 20

choices = np.linspace(0, step_num, gif_frame_count).astype("int").tolist()

vec_y = np.linspace(0, 1, i_N+1)
fig = plt.figure()
ims = []
for choice in choices:
    if (choice > step_num-1):
        choice = choice-1
    im = plt.plot(array_u[choice], vec_y, color="blue")
    plt.xlabel("u")
    plt.ylabel("y")
    plt.title("y-u diagram")
    ims.append(im)

print("正在生成gif……")

ani = animation.ArtistAnimation(
    fig, ims, interval=1000*(gif_time/gif_frame_count), repeat_delay=1000)
ani.save("{}.gif".format(gif_name), writer='pillow')
