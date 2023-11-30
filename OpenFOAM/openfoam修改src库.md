# openfoam 修改 src 库

遇到一个问题，要**把 sprayFoam 求解器的蒸发模型修改为自定义蒸发模型**。

sprayFoam 求解器本身没有实现蒸发模型，而是调用 `$FOAM_SRC/lagrangian/intermediate/submodels/Reacting/PhaseChangeModel/` 中的模型。

查看 sprayFoam 的 `Make/options` 文件，其中有

```
-I$(WM_PROJECT_USER_DIR)/src/lagrangian/intermediate/lnInclude
```

因此做出以下尝试

## 1. 复制并修改 src/lagrangian

```bash
mkdir -p $WM_PROJECT_USER_DIR/src
cp -r $FOAM_SRC/lagrangian $WM_PROJECT_USER_DIR/src
```

修改 `$WM_PROJECT_USER_DIR/src/lagrangian/intermediate/Make/files`

```
LIB = $(FOAM_USER_LIBBIN)/libmyLagrangianIntermediate
```

<!-- 由于测试的 case 调用的是 LiquidEvaporationBoil 模型，因此修改 `$LIB_SRC/src/lagrangian/intermediate/submodels/Reacting/PhaseChangeModel/LiquidEvaporationBoil.C` 文件

首先 `#include <iostream>`，然后在构造器中加上

```cpp
std::cerr << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << std::endl;
std::cerr << "call LiquidEvaporationBoil Constructors" << std::endl;
std::cerr << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << std::endl;
``` -->

## 2. 复制并修改 sparyFoam

```bash
cp -r $FOAM_SRC/../applications/solvers/lagrangian/sprayFoam <项目路径>
```

修改 `files`

```
EXE = $(FOAM_USER_APPBIN)/mySprayFoam
```

修改 `options`

对 `EXE_INC`：

```
-I$(LIB_SRC)/../applications/solvers/lagrangian/reactingParcelFoam \
```

```
-I$(WM_PROJECT_USER_DIR)/src/lagrangian/intermediate/lnInclude \
```

对 `EXE_LIBS`：

在开头加上 `-L$(FOAM_USER_LIBBIN) \`

把 `-llagrangianIntermediate \` 修改为 `-lmyLagrangianIntermediate \`


## 3. 创建自定义蒸发模型

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

在 `$WM_PROJECT_USER_DIR/src/lagrangian/intermediate/Make/files` 中添加

```
<项目路径>/MyLiquidEvaporation/MyLiquidEvaporation.C
```

之后要修改 `src/lagrangian/intermediate`

在 `$WM_PROJECT_USER_DIR/src/lagrangian/intermediate/Make/options` 中添加

```
-I<项目路径>/MyLiquidEvaporation
```

<!-- 创建软链接

```bash
ln -s <项目路径>/MyLiquidEvaporation/MyLiquidEvaporation.H $WM_PROJECT_USER_DIR/src/lagrangian/intermediate/lnInclude/MyLiquidEvaporation.H
ln -s <项目路径>/MyLiquidEvaporation/MyLiquidEvaporation.C $WM_PROJECT_USER_DIR/src/lagrangian/intermediate/lnInclude/MyLiquidEvaporation.C
ln -s <项目路径>/MyLiquidEvaporation/MyLiquidEvaporation.H $WM_PROJECT_USER_DIR/src/lagrangian/intermediate/submodels/Reacting/PhaseChangeModel/MyLiquidEvaporation.H
ln -s <项目路径>/MyLiquidEvaporation/MyLiquidEvaporation.C $WM_PROJECT_USER_DIR/src/lagrangian/intermediate/submodels/Reacting/PhaseChangeModel/MyLiquidEvaporation.C
``` -->

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