# 粒子和液滴的控制方程

本文为文献 CHEN X, FENG Y, ZHONG W, et al. Numerical investigation of the interaction, transport and deposition of multicomponent roplets in a simple mouth-throat model[J]. Journal of Aerosol Science,2017,105: 108-127. 中Chapter 2.3 的内容。

## 1.粒子和液滴输运

> 球形粒子假设是指将一些物质或化学物质视为由无限小的、完全相同的、球形粒子组成的模型。在这个模型中，球形粒子的大小、形状和相互作用被认为是恒定的，这使得研究者可以将物质视为一系列简单的球形粒子，并对它们的运动和相互作用进行数学建模。——by ChatGPT

由于空气中的微粒（aerosol）密度大于空气，旋转可以忽略不计，热迁移有限，因此可以根据**球形粒子假设**确定颗粒或液滴的运动轨迹。
$$
m_{d} \frac{d \vec{u}_{d}}{d t}=\frac{1}{8} \pi \rho d_{d}^{2} C_{D d}\left(\vec{u}-\vec{u}_{d}\right)\left|\vec{u}-\vec{u}_{d}\right|+m_{d} \vec{g} \tag{1}
$$
其中，$d_{d}$为粒子直径；$C_{D d}$为阻力系数（drag force coefficient），其定义如下
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

## 2.多组分液滴-蒸汽相互作用


对于一个液滴或颗粒，其组分可以分为

1. 可溶且可蒸发组分：可蒸发组分通过液-气界面传质（应该指蒸发）直接导致液滴尺寸的变化
2. 可溶但不可蒸发的组分：可溶但不可蒸发的组分通常不会发生质量变化。然而，它们在液气界面上占据了部分表面积。因此，它们降低了液滴侧可蒸发组分的质量分数，从而影响液滴在平衡状态下的直径。
3. 不可溶且不可挥发组分：除对颗粒大小和热容量的贡献外，不可溶且不可挥发组分对颗粒的影响有限。

液滴(或颗粒)质量的变化可由各可蒸发组分的质量变化推导出来(Longest &
Kleinstreuer, 2005; Zhang et al., 2012; Feng, Castro, Kleinstreuer & Rostami, 2016)
$$
\frac{d m_{d}}{d t}=-\sum_{e=1}^{k} \int_{\text {surf }} n_{e} d A \approx-\sum_{e=1}^{k}\left(\bar{n}_{e} \cdot A\right)
$$

式中$n_{e}$为表面可蒸发组分$e$的平均质量通量，可以表示为
$$
\bar{n}_{e}=\frac{\rho_{g} S h \widetilde{D}_{e} C_{m}}{d_{d}} \ln \frac{1-Y_{e, \infty}}{1-Y_{e, surf}}
$$

其中$\rho_{g}$是周围气体的密度；$Sh$是Sherwood数 (Clift, Grace & Weber, 1978)，其定义为
$$
Sh=\sqrt[3]{1+\operatorname{Re}_{d} \cdot Sc} \cdot \max \left[1, Re_{d}^{0.077}\right]
$$
其中$Sc$是Schmidt数，定义为
$$
S c=\frac{\mu}{\rho D_{e}}
$$
其中$D_{e}$是组分$e$的质量扩散系数。

回到$\bar{n}_{e}$。$Y_{e, surf},Y_{e, \infty}$分别是，组分$e$在液滴表面和远离液滴的空气中的的质量分数；$C_{m}$是Fuchs-Knudsen数修正(Ferron, Kreyling & Haider, 1988)，定义为
$$
C_{m}=\frac{1+Kn}{1+\left(\frac{4}{3 \alpha_{m}}+0.377\right) Kn+\frac{4}{3 \alpha_{m}} Kn^{2}}
$$

其中$Kn$为Knudsen数
$$
Kn=\frac{2 \lambda }{d_{d}}
$$
$\lambda$为气体平均自由程。

$\alpha_{m}$是质量热调节系数(Fuchs & Sutugin, 1969)。

$Y_{e, surf}$可以用修正Raoult's定律获得 (Finlay, 2001; Longest & Xi, 2008; Zhang et al., 2012; Feng et al. 2016)
$$
Y_{e, \text { surf }}=\gamma_{e} x_{e} K_{e} \frac{P_{ve, sat}\left(T_{d}\right)}{\rho R_{e} T_{d}}
$$

其中$\gamma_{e}$为组分$e$的活度系数；$x_{e}$是组分$e$在液滴中的摩尔分数；$T_{d}$是液滴温度；$P_{ve, sat}\left(T_{d}\right)$为组分$e$在温度$T_{d}$下的饱和压力；$K_{e}$是对Kelvin效应的修正：与平面相比，高度弯曲的表面上的蒸汽浓度更高 (Finlay, 2001)。$K_{e}$定义如下
$$
K_{e}=e^{\frac{4 \sigma M_{e}}{R \rho_{d} d_{d} T_{d}}}
$$

其中$\sigma$是液滴的表面张力；$M_{e}$是组分$e$的摩尔质量；$R$是通用气体常数；$\rho_{d}$是液滴的密度。



