# %%
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib.animation as animation


# %%
# 设置参数
f_a = 1.0
f_l = 1.0  # 绳长1m
f_x_0 = 0.3
f_v_0 = 1.0
f_delta_x = 0.001
f_delta_t = 0.001

i_N = int(f_l / f_delta_x)
i_k = int(f_x_0 / f_delta_x)

# u[n,i]的数据结构: list[dict{n,u[i]}]
array_u = []

# %%

# n=0,1
dict_u = {
    "n": 0,
    "vec": np.zeros(i_N+2)
}
array_u.append(dict_u)

dict_u = {
    "n": 1,
    "vec": np.zeros(i_N+2)
}
dict_u["vec"][i_k] = (f_v_0*f_delta_t)/f_delta_x
array_u.append(dict_u)


# %%

def u(n: int, i: int) -> float:
    """ 
    从array_u中读取[n,i]处的值
     """
    return array_u[n]["vec"][i]


def next_vec() -> np.array:
    """ 
    返回new_vec
     """
    new_vec = np.zeros(i_N+2)
    n = len(array_u)-1

    # 这里应用边界条件
    new_vec[0] = 0.0
    new_vec[i_N+1] = 0.0

    for i in range(1, i_N+1):
        new_vec[i] = 2*u(n, i)-u(n-1, i)+((f_a*f_a*f_delta_t*f_delta_t) /
                                          (f_delta_x**2))*(u(n, i+1)-2*u(n, i)+u(n, i-1))

    return new_vec

# %%


for iii in tqdm(range(10000)):
    dict_u = {
        "n": len(array_u),
        "vec": next_vec()
    }
    array_u.append(dict_u)

# %%
# plt.ion()
# for index, dict_u in enumerate(array_u):
#     if index % 100 == 0:
#         plt.clf()
#         plt.plot(dict_u["vec"])
#         plt.pause(0.1)

# %%
fig = plt.figure()
ims = []
for index, dict_u in enumerate(array_u):
    if index % 100 == 0:
        im=plt.plot(dict_u["vec"],color="blue")
        ims.append(im)

ani = animation.ArtistAnimation(fig, ims, interval=100, repeat_delay=1000)
ani.save("2_1.gif",writer='pillow')       