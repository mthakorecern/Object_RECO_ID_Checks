python3 Plotter.py \
    --inputs \
       /nfs_scratch/mithakor/ObjectReco_ID_Efficiency/AddedHEEPID/GluGlutoRadiontoHHto2B2Tau_M-1000_narrow_TuneCP5_13p6TeV_madgraph-pythia8.root  \
       /nfs_scratch/mithakor/ObjectReco_ID_Efficiency/AddedHEEPID/GluGlutoRadiontoHHto2B2Tau_M-1500_narrow_TuneCP5_13p6TeV_madgraph-pythia8.root  \
       /nfs_scratch/mithakor/ObjectReco_ID_Efficiency/AddedHEEPID/GluGlutoRadiontoHHto2B2Tau_M-2000_narrow_TuneCP5_13p6TeV_madgraph-pythia8.root  \
       /nfs_scratch/mithakor/ObjectReco_ID_Efficiency/AddedHEEPID/GluGlutoRadiontoHHto2B2Tau_M-2500_narrow_TuneCP5_13p6TeV_madgraph-pythia8.root  \
       /nfs_scratch/mithakor/ObjectReco_ID_Efficiency/AddedHEEPID/GluGlutoRadiontoHHto2B2Tau_M-3000_narrow_TuneCP5_13p6TeV_madgraph-pythia8.root  \
       /nfs_scratch/mithakor/ObjectReco_ID_Efficiency/AddedHEEPID/GluGlutoRadiontoHHto2B2Tau_M-4000_narrow_TuneCP5_13p6TeV_madgraph-pythia8.root  \
       /nfs_scratch/mithakor/ObjectReco_ID_Efficiency/AddedHEEPID/GluGlutoRadiontoHHto2B2Tau_M-4500_narrow_TuneCP5_13p6TeV_madgraph-pythia8.root  \
    --ele-skip-cuts "5,7" \
    --tau-decay-modes ""  \
    --boosted-tau-decay-modes ""\
    --require-tau-newdm-id  \
    --boosted-tau-raw-threshold 0.85    \
    --nproc 40  \
    &> plotting_eff.log &
