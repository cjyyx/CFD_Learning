#! https://zhuanlan.zhihu.com/p/601418292
# 计算流体力学（六）OpenFOAM安装

OpenFOAM是用于处理cfd问题的开源软件，是在cfd学习过程中必然使用到的软件工具。这里介绍OpenFOAM的安装与简单使用。

## 安装Linux环境

OpenFOAM只能在Linux系统或Mac系统上运行。这里仅说Linux系统。如果是国产linux系统，比如Deepin、优麒麟等，可以直接进行下一步。

如果按装的是ubuntu等非国产linux系统，由于网络原因，这些系统内置软件包管理器不能更新或下载软件。因此需要换源操作。简单来说，就是让系统软件访问国内的服务器而不是国外的。具体方法可以参考[清华大学镜像站的教程](https://mirrors.tuna.tsinghua.edu.cn/help/ubuntu/)，直接在终端执行命令替换中的命令。

如果想要在Windows中运行OpenFOAM，有两种方案。一种是安装虚拟机软件，这里推荐VMware。另一种是安装WSL2，具体参考[微软官方教程](https://learn.microsoft.com/zh-cn/windows/wsl/)。

## 安装OpenFOAM

OpenFOAM有很多个版本，这里给出OpenFOAM v10的安装教程。当然，[官网也给出了教程](https://openfoam.org/download/10-ubuntu/)。但是与ubuntu一样，OpenFOAM的安装需要换源。这里使用[常恭大佬提供的镜像站](https://www.cfdem.cn/openfoam-mirror/)。

直接在终端输入以下命令。
```bash
sudo sh -c "wget -O - https://dl.openfoam.org/gpg.key | apt-key add -"
sudo add-apt-repository http://dl.cfdem.cn/ubuntu
sudo apt-get update
sudo apt-get -y install openfoam10
sudo apt install paraview
```
然后OpenFOAM v10就安装完成了。之后是添加环境变量。运行以下命令
```bash
cat /opt/openfoam10/etc/bashrc >> ~/.bashrc
source ~/.bashrc
```
之后运行`pisoFoam -help`，如果结果如下，初步说明安装完成了。
![](PasteImage/2023-01-28-14-59-17.png)

## 第一个算例

参考https://doc.cfd.direct/openfoam/user-guide-v10/cavity

将示例中的cavity算例复制到`$FOAM_RUN`目录中。
```bash
mkdir -p $FOAM_RUN
cp -r $FOAM_TUTORIALS/incompressible/icoFoam/cavity/cavity $FOAM_RUN
cd $FOAM_RUN/cavity/
```

运行如下命令
```bash
blockMesh
surfaceFeatureExtract
snappyHexMesh -overwrite
checkMesh
icoFoam
```

最后结果如下。这是计算结果。

![](PasteImage/2023-01-28-22-13-31.png)

<未完待续>

## 参考资料
1. https://zhuanlan.zhihu.com/p/85092318
2. https://cfd-china.com/category/6/openfoam
3. http://dyfluid.com/
4. https://www.cfdem.cn/openfoam-mirror/
5. https://doc.cfd.direct/openfoam/user-guide-v10/
6. https://www.bilibili.com/video/BV14S4y1c7NS

[目录](https://zhuanlan.zhihu.com/p/599909213)