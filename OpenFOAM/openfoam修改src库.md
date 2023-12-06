#! https://zhuanlan.zhihu.com/p/670815002
# openfoam 修改 src 库

遇到一个问题，要**把 sprayFoam 求解器的蒸发模型修改为自定义蒸发模型**。

sprayFoam 求解器本身没有实现蒸发模型，而是调用 `$FOAM_SRC/lagrangian/intermediate/submodels/Reacting/PhaseChangeModel/` 中的模型。

## 蒸发模型库依赖探究

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

这样可以测试修改蒸发模型的有效性

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

之后运行，发现报错消失了

## 创建自定义蒸发模型

首先要知道蒸发模型是如何被调用的，运行下面命令

```bash
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

对应的，要更新 case 文件，修改 `constant/sprayCloudProperties.subModels.phaseChangeModel` 为 `myLiquidEvaporation`，这应该与 `TypeName("myLiquidEvaporation")` 相对应。同时添加

```
    myLiquidEvaporationCoeffs
    {
        enthalpyTransfer enthalpyDifference;

        activeLiquids    ( C7H16 );
    }
```

接下来要创建链接

```bash
mkdir -p $WM_PROJECT_USER_DIR/src/lagrangian/intermediate/submodels/Reacting/PhaseChangeModel/MyLiquidEvaporation
ln <项目路径>/MyLiquidEvaporation/MyLiquidEvaporation.H $WM_PROJECT_USER_DIR/src/lagrangian/intermediate/submodels/Reacting/PhaseChangeModel/MyLiquidEvaporation/MyLiquidEvaporation.H
ln <项目路径>/MyLiquidEvaporation/MyLiquidEvaporation.C $WM_PROJECT_USER_DIR/src/lagrangian/intermediate/submodels/Reacting/PhaseChangeModel/MyLiquidEvaporation/MyLiquidEvaporation.C
``` 


然后在 `$WM_PROJECT_USER_DIR/src/lagrangian/intermediate/parcels/include/makeReactingParcelPhaseChangeModels.H` 中加上

```cpp
#include "MyLiquidEvaporation.H"
```

```cpp
makePhaseChangeModelType(MyLiquidEvaporation, CloudType);
```

之后编译

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
```


为判断是否成功添加，运行以下命令

```bash
cd $WM_PROJECT_USER_DIR/src/lagrangian/intermediate
echo ---lagrangian/intermediate

echo MyLiquidEvaporation
grep -r -I MyLiquidEvaporation | wc -l
grep -r -I myLiquidEvaporation | wc -l

echo LiquidEvaporationBoil
grep -r -I LiquidEvaporationBoil | wc -l
grep -r -I liquidEvaporationBoil | wc -l

echo LiquidEvapFuchsKnudsen
grep -r -I LiquidEvapFuchsKnudsen | wc -l
grep -r -I liquidEvapFuchsKnudsen | wc -l

cd $WM_PROJECT_USER_DIR/src/lagrangian/spray
echo ---lagrangian/spray

echo MyLiquidEvaporation
grep -r -I MyLiquidEvaporation | wc -l
grep -r -I myLiquidEvaporation | wc -l

echo LiquidEvaporationBoil
grep -r -I LiquidEvaporationBoil | wc -l
grep -r -I liquidEvaporationBoil | wc -l

cd $FOAM_USER_LIBBIN
echo ---FOAM_USER_LIBBIN

echo MyLiquidEvaporation
nm -C libmyLagrangianIntermediate.so | grep MyLiquidEvaporation | wc -l
nm -C libmyLagrangianTurbulence.so | grep LiquidEvaporationBoil | wc -l
nm -C libmyLagrangianSpray.so | grep MyLiquidEvaporation | wc -l

echo LiquidEvaporationBoil
nm -C libmyLagrangianIntermediate.so | grep LiquidEvaporationBoil | wc -l
nm -C libmyLagrangianTurbulence.so | grep MyLiquidEvaporation | wc -l
nm -C libmyLagrangianSpray.so | grep LiquidEvaporationBoil | wc -l
```

如果 `MyLiquidEvaporation` 相关的符号数与 `LiquidEvaporationBoil` 相同，应该可以说明成功添加 `MyLiquidEvaporation`

之后再编译求解器，最后成功运行。至此，已满足把 sprayFoam 求解器的蒸发模型修改为自定义蒸发模型的要求。