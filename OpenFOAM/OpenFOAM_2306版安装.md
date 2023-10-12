#! https://zhuanlan.zhihu.com/p/660892262
# OpenFOAM-v2306 安装

参考 https://develop.openfoam.com/Development/openfoam/-/wikis/precompiled/debian

```bash
curl https://dl.openfoam.com/add-debian-repo.sh | sudo bash
# 建议使用美国代理
sudo apt-get update && sudo apt-get upgrade
sudo apt-get install openfoam2306-default

echo "source /usr/lib/openfoam/openfoam2306/etc/bashrc" >> ~/.bashrc
source ~/.bashrc
```

**paraview 另行安装**，参考 https://zhuanlan.zhihu.com/p/660892102

验证安装成功

```bash
mkdir -p $FOAM_RUN
cp -r $FOAM_TUTORIALS/incompressible/icoFoam/cavity/cavity $FOAM_RUN
cd $FOAM_RUN/cavity/

blockMesh
icoFoam
paraFoam
```

初始化工作区，复制示例、源码至工作区

```bash
rm -rf $FOAM_RUN
mkdir -p $FOAM_RUN

mkdir -p $FOAM_RUN/openfoam
cp -r $FOAM_TUTORIALS $FOAM_RUN/openfoam
cp -r $FOAM_SRC $FOAM_RUN/openfoam
cp -r $FOAM_APP $FOAM_RUN/openfoam
```