#!/bin/bash

inputfiledir_const=/data/bond/yimengzhang/GEM/MuHitDM-11.20
inputfiledir_maxwell=/data/bond/yimengzhang/GEM/MuHitDM-11.20
inputfilename_const=("Mu_out_Constant(mDM=0.005GeV).dat"
"Mu_out_Constant(mDM=0.05GeV).dat"
"Mu_out_Constant(mDM=0.1GeV).dat"
"Mu_out_Constant(mDM=0.2GeV).dat"
"Mu_out_Constant(mDM=0.5GeV).dat"
"Mu_out_Constant(mDM=1.0GeV).dat"
"Mu_out_Constant(mDM=10.0GeV).dat"
"Mu_out_Constant(mDM=100.0GeV).dat"
)
inputfilename_maxwell=("Mu_out_Maxwell(mDM=0.005GeV).dat"
"Mu_out_Maxwell(mDM=0.05GeV).dat"
"Mu_out_Maxwell(mDM=0.1GeV).dat"
"Mu_out_Maxwell(mDM=0.2GeV).dat"
"Mu_out_Maxwell(mDM=0.5GeV).dat"
"Mu_out_Maxwell(mDM=1.0GeV).dat"
"Mu_out_Maxwell(mDM=10.0GeV).dat"
"Mu_out_Maxwell(mDM=100.0GeV).dat"
)

outputfiledir_const=./root_file/constant
outputfiledir_maxwell=./root_file/maxwell
outputfilename_const=("DMmuon_const_0p005.root"
"DMmuon_const_0p05.root"
"DMmuon_const_0p1.root"
"DMmuon_const_0p2.root"
"DMmuon_const_0p5.root"
"DMmuon_const_1.root"
"DMmuon_const_10.root"
"DMmuon_const_100.root"
)
outputfilename_maxwell=("DMmuon_maxwell_0p005.root"
"DMmuon_maxwell_0p05.root"
"DMmuon_maxwell_0p1.root"
"DMmuon_maxwell_0p2.root"
"DMmuon_maxwell_0p5.root"
"DMmuon_maxwell_1.root"
"DMmuon_maxwell_10.root"
"DMmuon_maxwell_100.root"
)

evtnumber_const=(3382521
3242793
2999700
2695788
2302523
2000241
1195145
912537
)

evtnumber_maxwell=(3381847
3241424
2999514
2693908
2300931
1998008
1194383
911307
)

mass_const=(0p005 0p05 0p1 0p2 0p5 1 10 100)
mass_maxwell=(0p005 0p05 0p1 0p2 0p5 1 10 100)

# loop constant model
for (( i=0; i<${#inputfilename_const[@]}; i++ ))
do

  cp DMmuon_template.mac DMmuon_const_${mass_const[$i]}.mac
  sed -i "s#inputfilename#${inputfiledir_const}/${inputfilename_const[$i]}#g" DMmuon_const_${mass_const[$i]}.mac
  sed -i "s#outputfilename#${outputfiledir_const}/${outputfilename_const[$i]}#g" DMmuon_const_${mass_const[$i]}.mac
  sed -i "s#evtnumber#${evtnumber_const[$i]}#g" DMmuon_const_${mass_const[$i]}.mac

  ./muPos DMmuon_const_${mass_const[$i]}.mac &

done

# loop maxwell model
for (( i=0; i<${#inputfilename_maxwell[@]}; i++ ))
do

  cp DMmuon_template.mac DMmuon_maxwell_${mass_maxwell[$i]}.mac
  sed -i "s#inputfilename#${inputfiledir_maxwell}/${inputfilename_maxwell[$i]}#g" DMmuon_maxwell_${mass_maxwell[$i]}.mac
  sed -i "s#outputfilename#${outputfiledir_maxwell}/${outputfilename_maxwell[$i]}#g" DMmuon_maxwell_${mass_maxwell[$i]}.mac
  sed -i "s#evtnumber#${evtnumber_maxwell[$i]}#g" DMmuon_maxwell_${mass_maxwell[$i]}.mac

  ./muPos DMmuon_maxwell_${mass_maxwell[$i]}.mac &

done
