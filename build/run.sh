#!/bin/bash

# 循环运行Geant4模拟，传入不同的i值
for i in {0..19}
do
  # 运行Geant4模拟
  ./muPos SingleEngMu.mac

  # 重命名产生的ROOT文件
  mv mu1GeV.root root_file/mu1GeV_$i.root
done
