# openfoam修改src库

遇到一个问题，要把 sprayFoam 求解器的蒸发模型修改为自定义蒸发模型。

sprayFoam 求解器本身没有实现蒸发模型，而是调用 `$LIB_SRC/src/lagrangian/intermediate/submodels/Reacting/PhaseChangeModel/` 中的模型。

查看 sprayFoam 的 `Make/options` 文件，其中有

```
-I$(WM_PROJECT_USER_DIR)/src/lagrangian/intermediate/lnInclude
```

因此做出以下尝试

## 复制并修改 src/lagrangian

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

## 复制并修改 sparyFoam

```bash
cp -r $FOAM_SRC/../applications/solvers/lagrangian/sprayFoam <项目路径>
```

修改 files

```
EXE = $(FOAM_USER_APPBIN)/mySprayFoam
```

修改 options

```
-I$(LIB_SRC)/../applications/solvers/lagrangian/reactingParcelFoam \
```

```
-I$(WM_PROJECT_USER_DIR)/src/lagrangian/intermediate/lnInclude \
```

把 `-llagrangianIntermediate \` 修改为

```
-L$(FOAM_USER_LIBBIN)/libmyLagrangianIntermediate \
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

在 `$WM_PROJECT_USER_DIR/src/lagrangian/intermediate/Make/files` 中添加

```
<项目路径>/MyLiquidEvaporation/MyLiquidEvaporation.C
```

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

## 更新 case 文件

修改 `constant/sprayCloudProperties>subModels>phaseChangeModel` 为 `myLiquidEvaporation`，这与 `MyLiquidEvaporation.H` 中的 `TypeName("myLiquidEvaporation")` 相对应。


## 尝试编译

编译

```bash
cd $WM_PROJECT_USER_DIR/src/lagrangian/intermediate
rm -r lnInclude
rm -r Make/linux*
wmake > wmake.log

cd <项目路径>/sprayFoam
wmake > wmake.log
```
