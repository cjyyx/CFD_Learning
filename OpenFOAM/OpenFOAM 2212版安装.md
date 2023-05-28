参考 https://develop.openfoam.com/Development/openfoam/-/wikis/precompiled/debian

```bash
curl -s https://dl.openfoam.com/add-debian-repo.sh | sudo bash
wget -q -O - https://dl.openfoam.com/add-debian-repo.sh | sudo bash
sudo apt-get install openfoam2212-default

echo "source /usr/lib/openfoam/openfoam2212/etc/bashrc" >> ~/.bashrc
source ~/.bashrc

sudo apt install paraview
```

验证安装成功

```bash
mkdir -p $FOAM_RUN
cp -r $FOAM_TUTORIALS/incompressible/icoFoam/cavity/cavity $FOAM_RUN
cd $FOAM_RUN/cavity/

blockMesh
icoFoam
paraFoam
```

