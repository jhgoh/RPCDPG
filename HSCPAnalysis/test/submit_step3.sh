#!/bin/bash

cmsDriver.py step3 --conditions auto:run2_mc -n 10 --era Phase2C1 --eventcontent FEVTDEBUGHLT,DQM \
                   --runUnscheduled -s RAW2DIGI,L1Reco,RECO,VALIDATION:@phase2Validation,DQM:@phase2 \
                   --datatier GEN-SIM-RECO,DQMIO \
                   --customise SLHCUpgradeSimulations/Configuration/combinedCustoms.cust_2023tilted \
                   --geometry Extended2023D1 --filein file:step2.root --fileout step3.root \
                   --no_exec --python_filename step3_cfg.py

source /afs/cern.ch/project/eos/installation/cms/etc/setup.sh
INPUTDIR=/store/user/jhgoh/RPCUpgrade/20160816_1
for i in `eos ls $INPUTDIR | grep HSCP`; do
  eos ls $INPUTDIR/$i | grep step2 | sed -e "s;^;$INPUTDIR/$i/;g" > step2_$i.txt
  create-batch --fileList step2_$i.txt --jobName step3_$i --cfg step3_dump_cfg.py --nJobs 1 --transferDest $INPUTDIR/$i
done

