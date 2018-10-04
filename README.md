# TimingTriggerAnalysis

Repo for analysis of generated data. After installing production code should be placed in Delphes directory.

Simple study of background rejection can be made (from the Delphes directory) using:

```bash
root 'TimingTriggerAnalysis/simpleAnalysisCode.C("pythia_delphes.root")'
```
Followed by 

```bash
python TimingTriggerAnalysis/makeCumu.py TimingTriggerAnalysis/output.root <outDirName>
```

This produces plots showing the time distribution/background suppression 
