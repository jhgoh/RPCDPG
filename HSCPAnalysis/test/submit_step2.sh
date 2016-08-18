#!/bin/bash

cmsDriver.py step2 --conditions auto:run2_mc -n -1 --era Phase2C1 --eventcontent FEVTDEBUGHLT --mc \
                   -s DIGI:pdigi_valid,L1,DIGI2RAW,HLT:@fake --datatier GEN-SIM-DIGI-RAW \
                   --customise SLHCUpgradeSimulations/Configuration/combinedCustoms.cust_2023tilted \
                   --geometry Extended2023D1 --filein file:step1.root --fileout step2.root \
                   --no_exec --python_filename step2_cfg.py

source /afs/cern.ch/project/eos/installation/cms/etc/setup.sh
INPUTDIR=/store/user/jhgoh/RPCUpgrade/20160816_1
for i in `eos ls $INPUTDIR | grep HSCP`; do
  eos ls $INPUTDIR/$i | sed -e "s;^;$INPUTDIR/$i/;g" > step2_$i.txt
  create-batch --fileList step2_$i.txt --jobName step2_$i --cfg step2_cfg.py --maxFiles 1 --transferDest $INPUTDIR/$i
done

