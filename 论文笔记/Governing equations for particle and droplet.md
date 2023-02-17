# 粒子和液滴的控制方程

本文为文献 CHEN X, FENG Y, ZHONG W, et al. Numerical investigation of the interaction, transport and deposition of multicomponent roplets in a simple mouth-throat model[J]. Journal of Aerosol Science,2017,105: 108-127. 中Chapter 2.3 的内容。

## 1.粒子和液滴输运

由于空气中的微粒（aerosol）密度大于空气，旋转可以忽略不计，热迁移有限，因此可以根据**球形粒子假设**确定颗粒或液滴的运动轨迹。
$$
m_{d} \frac{d \vec{u}_{d}}{d t}=\frac{1}{8} \pi \rho d_{d}^{2} C_{D d}\left(\vec{u}-\vec{u}_{d}\right)\left|\vec{u}-\vec{u}_{d}\right|+m_{d} \vec{g} \tag{1}
$$
其中，$d_{d}$为液滴直径；$C_{D d}$为阻力系数（drag force coefficient），其定义如下
$$
C_{D d}=C_{D} / C_{c} \tag{2}
$$
其中$C_{D}$是Cunningham修正因子（Allen & Raabe, 1985）;$C_{c}$为阻力系数（drag coefficient），其定义如下（Morsi & Alexander,1972）
$$
C_{D}=a_{1}+\frac{a_{2}}{\operatorname{Re}_{p}}+\frac{a_{3}}{\operatorname{Re}_{p}^{2}} \tag{3}
$$
其中常数$a_{1},a_{2},a_{3}$由粒子雷诺数决定。

对于湍流流动，瞬时流体速度$u$由两部分组成
$$
u_{i}=\bar{u}_{i}+u_{i}^{\prime} \quad i=x, y, z \tag{4}
$$
其中$\bar{u}_{i}$是随时间变化的平均速度，$u_{i}^{\prime}$是表示随机方向涡流的波动分量。对于各向同性区域，波动分量$u_{i}^{\prime}$可以用(Gosman & Loannides, 1983)计算。
$$
u_{i}^{\prime}=\xi_{i} \sqrt{\frac{2}{3} k} \tag{5}
$$
其中，$\xi_{i}$为来自标准正态分布的随机数。然而，DNS结果表明，近壁区域的流体是各向异性的，流体的波动速度可以通过在不同方向上应用阻尼函数来模拟(Kim, Moin & Moser, 1987;Wang & James, 1999)，即
$$
\begin{array}{l}
u_{i}^{\prime} =f_{i} \xi_{i} \sqrt{\frac{2}{3} k} \\
f_{u} =1+0.285\left(y^{+}+6\right) \exp \left[-0.455\left(y^{+}+6\right)^{0.53}\right] \\
f_{v} =1-\exp \left(-0.02 y^{+}\right) \\
f_{w} =\sqrt{3-f_{u}^{2}-f_{v}^{2}}
\end{array} \tag{6}
$$
其中$f_{u}$为顺流方向函数，$f_{v}$垂直于最近的壁面，$f_{w}$垂直于$f_{u}$和$f_{v}$。方程式。以上式子适用于$y^{+}$值小于80的区域。

尽管这些函数是由管道流导出的，但它们对于更广泛的几何形状是相当准确的(Wang & James, 1999)。