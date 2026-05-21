python3 Plotter.py \
    --inputs \
       /nfs_scratch/mithakor/ObjectReco_ID_Efficiency/Run2_matching_0_3/RadionTohhTohtatahbb_narrow_M-1000_TuneCP5_13TeV-madgraph-pythia8.root  \
       /nfs_scratch/mithakor/ObjectReco_ID_Efficiency/Run2_matching_0_3/RadionTohhTohtatahbb_narrow_M-1600_TuneCP5_13TeV-madgraph-pythia8.root  \
       /nfs_scratch/mithakor/ObjectReco_ID_Efficiency/Run2_matching_0_3/RadionTohhTohtatahbb_narrow_M-2000_TuneCP5_13TeV-madgraph-pythia8.root  \
       /nfs_scratch/mithakor/ObjectReco_ID_Efficiency/Run2_matching_0_3/RadionTohhTohtatahbb_narrow_M-2500_TuneCP5_13TeV-madgraph-pythia8.root  \
       /nfs_scratch/mithakor/ObjectReco_ID_Efficiency/Run2_matching_0_3/RadionTohhTohtatahbb_narrow_M-3000_TuneCP5_13TeV-madgraph-pythia8.root  \
       /nfs_scratch/mithakor/ObjectReco_ID_Efficiency/Run2_matching_0_3/RadionTohhTohtatahbb_narrow_M-4000_TuneCP5_13TeV-madgraph-pythia8.root  \
       /nfs_scratch/mithakor/ObjectReco_ID_Efficiency/Run2_matching_0_3/RadionTohhTohtatahbb_narrow_M-4500_TuneCP5_13TeV-madgraph-pythia8.root  \
    --ele-skip-cuts "7" \
    --nproc 40  \
    &> plotting_run2_eff.log &
