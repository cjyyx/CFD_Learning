# openfoam 修改 src 库

遇到一个问题，要**把 sprayFoam 求解器的蒸发模型修改为自定义蒸发模型**。

sprayFoam 求解器本身没有实现蒸发模型，而是调用 `$FOAM_SRC/lagrangian/intermediate/submodels/Reacting/PhaseChangeModel/` 中的模型。

## 添加蒸发模型库探究

这里先从 liquidEvaporationBoil 这个已有的蒸发模型入手

首先运行

```bash
cd $FOAM_SRC
cd ..
grep -r liquidEvaporationBoil .
```

可以看到结果中有

```
grep: ./platforms/linux64GccDPInt32Opt/lib/liblagrangianIntermediate.so: binary file matches
grep: ./platforms/linux64GccDPInt32Opt/lib/liblagrangianSpray.so: binary file matches
```

因此使用自定义的蒸发模型，需要重新编译 `liblagrangianIntermediate`、`liblagrangianSpray` 这两个库

因此做出以下尝试

## 复制并修改 src/lagrangian

```bash
mkdir -p $WM_PROJECT_USER_DIR/src
cp -r $FOAM_SRC/lagrangian $WM_PROJECT_USER_DIR/src
```

(1) 修改 `$WM_PROJECT_USER_DIR/src/lagrangian/intermediate/Make/files`

`LIB = $(FOAM_LIBBIN)/liblagrangianIntermediate` 改为

```
LIB = $(FOAM_USER_LIBBIN)/libmyLagrangianIntermediate
```

(2) 修改 `$WM_PROJECT_USER_DIR/src/lagrangian/spray/Make/files`

`LIB = $(FOAM_LIBBIN)/liblagrangianSpray` 改为

```
LIB = $(FOAM_USER_LIBBIN)/libmyLagrangianSpray
```

(3) 修改 `$WM_PROJECT_USER_DIR/src/lagrangian/spray/Make/options`

对 `EXE_INC`：

`-I$(LIB_SRC)/lagrangian/intermediate/lnInclude \` 改为

```
-I$(WM_PROJECT_USER_DIR)/src/lagrangian/intermediate/lnInclude \
```

对 `LIB_LIBS`：

在开头加上

```
-L$(FOAM_USER_LIBBIN) \
```

把 `-llagrangianIntermediate \` 修改为

```
-lmyLagrangianIntermediate \
```

(4) 由于测试的 case 调用的是 LiquidEvaporationBoil 模型，因此修改 `$WM_PROJECT_USER_DIR/src/lagrangian/intermediate/submodels/Reacting/PhaseChangeModel/LiquidEvaporationBoil/LiquidEvaporationBoil.C` 文件

首先 `#include <iostream>`，然后在构造器中加上

```cpp
std::cerr << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << std::endl;
std::cerr << "aaaaaaaaaaaaaaaaaaa Test !!!!!!!!!!!!!!" << std::endl;
std::cerr << "call LiquidEvaporationBoil Constructors" << std::endl;
std::cerr << "aaaaaaaaaaaaaaaaaaa Test !!!!!!!!!!!!!!" << std::endl;
std::cerr << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << std::endl;
```

这样可以证明修改蒸发模型有效

## 复制并修改 sparyFoam

```bash
cp -r $FOAM_SRC/../applications/solvers/lagrangian/sprayFoam <项目路径>
```

(1) 修改 `<项目路径>/sprayFoam/Make/files`

```
EXE = $(FOAM_USER_APPBIN)/mySprayFoam
```

(2) 修改 `<项目路径>/sprayFoam/Make/options`

对 `EXE_INC`：

```
-I$(LIB_SRC)/../applications/solvers/lagrangian/reactingParcelFoam \
```

```
-I$(WM_PROJECT_USER_DIR)/src/lagrangian/intermediate/lnInclude \
-I$(WM_PROJECT_USER_DIR)/src/lagrangian/spray/lnInclude \
```

对 `EXE_LIBS`：

在开头加上 `-L$(FOAM_USER_LIBBIN) \`

把 `-llagrangianIntermediate \` 修改为 `-lmyLagrangianIntermediate \`

把 `-llagrangianSpray \` 修改为 `-lmyLagrangianSpray \`

## 尝试编译并运行

编译

```bash
cd $WM_PROJECT_USER_DIR/src/lagrangian/intermediate
rm -r lnInclude
rm -r Make/linux*
wmake -j4 > wmake.log

cd $WM_PROJECT_USER_DIR/src/lagrangian/spray
rm -r lnInclude
rm -r Make/linux*
wmake -j4 > wmake.log

cd <项目路径>/sprayFoam
wmake -j4 > wmake.log
```

之后可以用 `$FOAM_USER_APPBIN/mySprayFoam` 运行使用 `liquidEvaporationBoil` 蒸发模型的算例，然后打印

```
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
aaaaaaaaaaaaaaaaaaa Test !!!!!!!!!!!!!!
call LiquidEvaporationBoil Constructors
aaaaaaaaaaaaaaaaaaa Test !!!!!!!!!!!!!!
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
```

但是同时看到有一堆报错

```
Duplicate entry dual in runtime table AveragingMethod
[stack trace]
=============
#1      platforms/linux64GccDPInt32Opt/lib/libmyLagrangianIntermediate.so(+0x4e564a) [0x7f14a396f64a]
#2      /lib64/ld-linux-x86-64.so.2(+0x647e) [0x7f14a5d3847e]
#3      /lib64/ld-linux-x86-64.so.2(+0x6568) [0x7f14a5d38568]
#4      /lib64/ld-linux-x86-64.so.2(+0x202ca) [0x7f14a5d522ca]
=============
Duplicate entry moment in runtime table AveragingMethod
[stack trace]
=============
#1      platforms/linux64GccDPInt32Opt/lib/libmyLagrangianIntermediate.so(+0x4e56dd) [0x7f14a396f6dd]
#2      /lib64/ld-linux-x86-64.so.2(+0x647e) [0x7f14a5d3847e]
#3      /lib64/ld-linux-x86-64.so.2(+0x6568) [0x7f14a5d38568]
#4      /lib64/ld-linux-x86-64.so.2(+0x202ca) [0x7f14a5d522ca]
=============
Duplicate entry moment in runtime table AveragingMethod
[stack trace]
=============
#1      platforms/linux64GccDPInt32Opt/lib/libmyLagrangianIntermediate.so(+0x4e576b) [0x7f14a396f76b]
#2      /lib64/ld-linux-x86-64.so.2(+0x647e) [0x7f14a5d3847e]
#3      /lib64/ld-linux-x86-64.so.2(+0x6568) [0x7f14a5d38568]
#4      /lib64/ld-linux-x86-64.so.2(+0x202ca) [0x7f14a5d522ca]
=============
```

这是 Duplicate Entry 问题，为解决该问题，运行

```bash
cd $WM_PROJECT_USER_DIR/src/lagrangian
grep -r llagrangianIntermediate .
```

结果为

```
./coalCombustion/Make/options:    -llagrangianIntermediate \
./turbulence/Make/options:    -llagrangianIntermediate \
```

这说明要重新编译 `coalCombustion`、`turbulence` (╬▔皿▔)╯

再运行 `grep -r llagrangianTurbulence .`，结果为

```
./coalCombustion/Make/options:    -llagrangianTurbulence \
./spray/Make/options:    -llagrangianTurbulence \
```

说明还要修改 `coalCombustion`、`spray`

于是进行以下操作

(1) 修改 `$WM_PROJECT_USER_DIR/src/lagrangian/turbulence/Make/files`

```
LIB = $(FOAM_USER_LIBBIN)/libmyLagrangianTurbulence
```

(2) 修改 `$WM_PROJECT_USER_DIR/src/lagrangian/turbulence/Make/options`

```
-I$(WM_PROJECT_USER_DIR)/src/lagrangian/intermediate/lnInclude \
```

```
-L$(FOAM_USER_LIBBIN) \
```

```
-lmyLagrangianIntermediate \
```

(3) 修改 `$WM_PROJECT_USER_DIR/src/lagrangian/spray/Make/options`

```
-I$(WM_PROJECT_USER_DIR)/src/lagrangian/turbulence/lnInclude \
```

```
-lmyLagrangianTurbulence \
```

(4) 修改 `<项目路径>/sprayFoam/Make/options`

```
-lmyLagrangianTurbulence \
```

(5) 编译

```bash
cd $WM_PROJECT_USER_DIR/src/lagrangian/intermediate
rm -r lnInclude
rm -r Make/linux*
wmake -j3 > wmake.log

cd $WM_PROJECT_USER_DIR/src/lagrangian/turbulence
rm -r lnInclude
rm -r Make/linux*
wmake -j3 > wmake.log

cd $WM_PROJECT_USER_DIR/src/lagrangian/spray
rm -r lnInclude
rm -r Make/linux*
wmake -j3 > wmake.log

cd <项目路径>/sprayFoam
wmake -j3 > wmake.log
```

## 创建自定义蒸发模型

首先要知道蒸发模型是如何被调用的，运行下面命令

```
cd $WM_PROJECT_USER_DIR/src/lagrangian/intermediate
grep -r LiquidEvaporationBoil.H .
```

结果为

```
./parcels/include/makeReactingParcelPhaseChangeModels.H:#include "LiquidEvaporationBoil.H"
```

这提示了如何添加自定义蒸发模型。

因此首先复制已有的蒸发模型至项目路径并重命名

```bash
cp -r $WM_PROJECT_USER_DIR/src/lagrangian/intermediate/submodels/Reacting/PhaseChangeModel/LiquidEvaporationBoil <项目路径>
cd <项目路径>
mv LiquidEvaporationBoil MyLiquidEvaporation

cd MyLiquidEvaporation
sed s/LiquidEvaporationBoil/MyLiquidEvaporation/g \
LiquidEvaporationBoil.C >MyLiquidEvaporation.C
sed s/LiquidEvaporationBoil/MyLiquidEvaporation/g \
LiquidEvaporationBoil.H >MyLiquidEvaporation.H
rm LiquidEvaporationBoil.C LiquidEvaporationBoil.H
```

也要修改 `MyLiquidEvaporation.H`

```
TypeName("myLiquidEvaporation");
```

对应的，要更新 case 文件，修改 `constant/sprayCloudProperties.subModels.phaseChangeModel` 为 `myLiquidEvaporation`，这应该与 `TypeName("myLiquidEvaporation")` 相对应。

接下来要创建软链接

```bash

ln -s <项目路径>/MyLiquidEvaporation/MyLiquidEvaporation.H $WM_PROJECT_USER_DIR/src/lagrangian/intermediate/submodels/Reacting/PhaseChangeModel/MyLiquidEvaporation.H
ln -s <项目路径>/MyLiquidEvaporation/MyLiquidEvaporation.C $WM_PROJECT_USER_DIR/src/lagrangian/intermediate/submodels/Reacting/PhaseChangeModel/MyLiquidEvaporation.C
``` 


然后在 `$WM_PROJECT_USER_DIR/src/lagrangian/intermediate/parcels/include/makeReactingParcelPhaseChangeModels.H` 中加上

```cpp
#include "MyLiquidEvaporation.H"
```

和

```cpp
makePhaseChangeModelType(MyLiquidEvaporation, CloudType);
```


## 4. 更新 case 文件

修改 `constant/sprayCloudProperties>subModels>phaseChangeModel` 为 `myLiquidEvaporation`，这应该与 `MyLiquidEvaporation.H` 中的 `TypeName("myLiquidEvaporation")` 相对应。


## 5. 尝试编译并运行

编译

```bash
cd $WM_PROJECT_USER_DIR/src/lagrangian/intermediate
rm -r lnInclude
rm -r Make/linux*
wmake > wmake.log

cd <项目路径>/sprayFoam
wmake > wmake.log
```

运行算例，结果报错

```
--> FOAM FATAL IO ERROR: (openfoam-2306)
Unknown phaseChangeModel type myLiquidEvaporation

Valid phaseChangeModel types :
4(liquidEvapFuchsKnudsen liquidEvaporation liquidEvaporationBoil none)


file: constant/sprayCloudProperties.subModels

    From static Foam::autoPtr<Foam::PhaseChangeModel<CloudType> > Foam::PhaseChangeModel<CloudType>::New(const Foam::dictionary&, CloudType&) [with CloudType = Foam::ReactingCloud<Foam::ThermoCloud<Foam::KinematicCloud<Foam::Cloud<Foam::SprayParcel<Foam::ReactingParcel<Foam::ThermoParcel<Foam::KinematicParcel<Foam::particle> > > > > > > >]
    in file /home/cjy/OpenFOAM/cjy-v2306/src/lagrangian/intermediate/lnInclude/PhaseChangeModelNew.C at line 47.

FOAM exiting
```

这里可以发现一个问题，就是 `Valid phaseChangeModel types` 中没有 `myLiquidEvaporation`。因此进行以下尝试：

删除 `$WM_PROJECT_USER_DIR/src/lagrangian` 中所有与 `liquidEvapFuchsKnudsen` 有关的内容，再进行编译和运行，结果发现 `Valid phaseChangeModel types` 中还有 `liquidEvapFuchsKnudsen`！！！

于是进行以下尝试

```bash
cd $FOAM_SRC
cd ..
grep -r liquidEvapFuchsKnudsen .
```

结果中有

```
grep: ./platforms/linux64GccDPInt32Opt/lib/liblagrangianSpray.so: binary file matches
```

**这说明不但要编译** `lagrangian/intermediate`，**还要编译** `lagrangian/spray`**！！！**

明白问题后，就继续以下操作

## 6. 修改并编译 src/lagrangian/spray

(1) 修改 `$WM_PROJECT_USER_DIR/src/lagrangian/spray/Make/files`

`LIB = $(FOAM_LIBBIN)/liblagrangianSpray` 改为

```
LIB = $(FOAM_USER_LIBBIN)/libmyLagrangianSpray
```

添加

```
<项目路径>/MyLiquidEvaporation/MyLiquidEvaporation.C
```

(2) 修改 `$WM_PROJECT_USER_DIR/src/lagrangian/spray/Make/options`

对 `EXE_INC`：

```
-I$(WM_PROJECT_USER_DIR)/src/lagrangian/intermediate/lnInclude \
```

添加

```
-I<项目路径>/MyLiquidEvaporation
```

对 `EXE_LIBS`：

在开头加上 `-L$(FOAM_USER_LIBBIN) \`

把 `-llagrangianIntermediate \` 修改为 `-lmyLagrangianIntermediate \`

(3) 编译

```bash
cd $WM_PROJECT_USER_DIR/src/lagrangian/spray
rm -r lnInclude
rm -r Make/linux*
wmake > wmake.log
```

(4) 最后修改 sprayFoam 的 options ，把 `-llagrangianSpray \` 修改为

```
-lmyLagrangianSpray \
```

---

然后编译，结果开始报与 `liquidEvapFuchsKnudsen` 有关的错误了……

于是把 `liquidEvapFuchsKnudsen` 加回来，再全部重新编译一遍

最终程序跑通，成功实现目标！！！