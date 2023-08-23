## 验证安装成功

```bash
mkdir -p $FOAM_RUN
cp -r $FOAM_TUTORIALS/incompressible/icoFoam/cavity/cavity $FOAM_RUN
cd $FOAM_RUN/cavity/

blockMesh
icoFoam
paraFoam
```

或者

```bash
mkdir -p $FOAM_RUN
cp -r $FOAM_TUTORIALS/incompressible/icoFoam/cavity/cavity $FOAM_RUN
cd $FOAM_RUN

blockMesh -case cavity
icoFoam -case cavity
foamToVTK -case cavity -ascii
paraview --data=cavity/VTK/cavity.vtm.series
```

## 初始化工作区

```bash
rm -rf $FOAM_RUN
mkdir -p $FOAM_RUN
chmod -R +x $FOAM_RUN
```



## 复制示例、源码至工作区

```bash
mkdir -p $FOAM_RUN/openfoam
cp -r $FOAM_TUTORIALS $FOAM_RUN/openfoam
cp -r $FOAM_SRC $FOAM_RUN/openfoam
cp -r $FOAM_APP $FOAM_RUN/openfoam
```

