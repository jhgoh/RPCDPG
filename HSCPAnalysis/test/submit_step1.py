#!/usr/bin/env python

cmd = "--conditions auto:run2_mc -n 1000 --era Phase2C1 --eventcontent FEVTDEBUG --mc"
cmd += " -s GEN,SIM --datatier GEN-SIM --beamspot Realistic50ns13TeVCollision"
cmd += " --customise SLHCUpgradeSimulations/Configuration/combinedCustoms.cust_2023tilted"
cmd +=               ",RPCUpgrade/HSCPAnalysis/oneHSCPAtHighEtaFilter_cff.customiseOneHSCP"
cmd +=               ",RPCUpgrade/HSCPAnalysis/Exotica_HSCP_SIM_cfi.customiseHSCPSim"
cmd += " --geometry Extended2023D1 --fileout step1.root --no_exec"

basedir = "Configuration/GenProduction/python/ThirteenTeV"
firstRun = 1

mm = [100, 200, 400, 600, 1000, 1200, 1600, 2600]
nn = [2780, 2250, 2626, 3197, 4137, 4630, 5678, 9541]
for m, n in zip(mm, nn):
    cff = "HSCPstop_M_%d_TuneCUETP8M1_13TeV_pythia8_cff" % m
    print "cmsDriver.py %s/%s " % (basedir, cff),
    print "--python_filename step1_HSCPstop_M_%d_cfg.py " % m,
    print cmd

    nJobs = 10
    nEvents = n/nJobs
    print "create-batch -n --jobName HSCPstop_M_%d --cfg step1_HSCPstop_M_%d_cfg.py --maxEvent %d --nJobs %d --firstRun %d " % (m, m, nEvents, nJobs, firstRun),
    print "--transferDest /store/user/jhgoh/RPCUpgrade/20160816_1/HSCPstop_M_%d &" % m

mm = [156, 200, 308, 494, 651, 1029, 1218, 1599]
nn = [2869, 1541, 1755, 2302, 3044, 4709, 5870, 8190]
for m, n in zip(mm, nn):
    cff = "HSCPppstau_M_%d_TuneZ2star_13TeV_pythia6_cff" % m
    print "cmsDriver.py %s/%s " % (basedir, cff),
    print "--python_filename step1_HSCPppstau_M_%d_cfg.py " % m,
    print cmd

    nJobs = 10
    nEvents = n/nJobs
    print "create-batch -n --jobName HSCPppstau_M_%d --cfg step1_HSCPppstau_M_%d_cfg.py --maxEvent %d --nJobs %d --firstRun %d " % (m, m, nEvents, nJobs, firstRun),
    print "--transferDest /store/user/jhgoh/RPCUpgrade/20160816_1/HSCPppstau_M_%d &" % m
