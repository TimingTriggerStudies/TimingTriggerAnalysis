#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif
#include "TLorentzVector.h"

double median(std::vector<double>& allTimes)
{
    size_t size = allTimes.size();

    if (size == 0)
    {
	return 0;  // Undefined, really.
    }
    else if (size == 1)
    {
	return allTimes[0];
    }
    else
    {
	sort(allTimes.begin(), allTimes.end());
	if (size % 2 == 0)
	{
	    return (allTimes[size / 2 - 1] + allTimes[size / 2]) / 2;
	}
	else 
	{
	    return allTimes[size / 2];
	}
    }
}

void getTrackVars(std::vector<Track *> trackCluster,double &maxTime,double &meanTime,double &modeTime,double&medianTime,double& rmsTime,double &clusterSize){
    maxTime = 0;
    meanTime = 0;
    modeTime = 0;
    medianTime = 0;
    rmsTime = 0;
    clusterSize = trackCluster.size();
    std::vector<double> allTimes;
    for(Int_t iTrack = 0; iTrack < trackCluster.size(); ++iTrack)
    {
	Track * track = trackCluster.at(iTrack);
	meanTime += track->T;
	if (track->T > maxTime){
	    maxTime = track->T;
	}
	allTimes.push_back(track->T);
    }
    for(Int_t iTrack = 0; iTrack < trackCluster.size(); ++iTrack)
    {
	Track * track = trackCluster.at(iTrack);
	rmsTime += 1./(clusterSize-1) * (track->T - meanTime)*(track->T - meanTime);
    }
    if (clusterSize > 1)meanTime /= clusterSize;
    if (clusterSize > 2) rmsTime = TMath::Sqrt(rmsTime);
    else rmsTime = 0;
    medianTime = median(allTimes);

}
void simpleAnalysisCode(const char *inputFile)
{
    gSystem->Load("libDelphes");

    // Create chain of root trees
    float ptThreshSeed = 10;
    float ptThresh = 1;
    float deltaRThresh = 0.4;
    float etaThresh = 3;
    TChain chain("Delphes");
    chain.Add(inputFile);

    // Create object of class ExRootTreeReader
    ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
    Long64_t numberOfEntries = treeReader->GetEntries();

    // Get pointers to branches used in this analysis
    TClonesArray *branchEvent = treeReader->UseBranch("Event");
    TClonesArray *branchTrack = treeReader->UseBranch("Track");

    // Book histograms
    TH1 *histMaxTime = new TH1F("maxTime", "time", 100, 0.0, 1E-9);
    TH1 *histMeanTime = new TH1F("meanTime", "time", 200, -1E-9, 1E-9);
    TH1 *histMeanTimeMax = new TH1F("meanTimeMax", "time", 200, -1E-9, 1E-9);
    TH1 *histMeanTimeNTracks = new TH2F("meanTimeVsNTracks", "time;ntracks", 200, -1E-9, 1E-9,100,0,100);
    TH1 *histMedianTimeNTracks = new TH2F("medianTimeVsNTracks", "time;ntracks", 200, -1E-9, 1E-9,100,0,100);
    TH1 *histMedianTime = new TH1F("medianTime", "time", 200, -1E-9, 1E-9);
    TH1 *histRMSTime = new TH1F("rmsTime", "time", 300, 0, 3E-9);
    TH2 *histMeanRMSTime = new TH2F("meanVsRMSTime",";mean;rms",200,-1E-9,1E-9,300,0,3E-9);
    TH1 *histNTrack = new TH1F("nTrack", "nTrack", 100, 0.0, 100);
    TH1 *histNJets = new TH1F("nJets", "nJets", 100, 0.0, 100);
    std::vector<std::vector<Track *>> trackClusters;

    // Loop over all events
    // for(Int_t entry = 0; entry < 10; ++entry)
    for(Int_t entry = 0; entry < numberOfEntries; ++entry)
    {
	// Load selected branches with data from specified event
	treeReader->ReadEntry(entry);

	// If event contains at least 1 track
	float maxTime = 0;
	float meanTimeMax = -999;
	std::set<Int_t> indexSet;
	if(branchTrack->GetEntries() > 0)
	{
	    for(Int_t iTrack = 0; iTrack < branchTrack->GetEntries(); ++iTrack)
	    {
		if (indexSet.find(iTrack) != indexSet.end()) continue;
		// Take first track
		Track *track = (Track*) branchTrack->At(iTrack);
		if (track->PT < ptThreshSeed || fabs(track->Eta) > etaThresh || track->T > 1) continue;
		// Plot track transverse momentum
		std::vector<Track *> trackCluster;
		// std::cout << track->VertexIndex << std::endl;
		TLorentzVector trackVector;
		trackVector.SetPtEtaPhiM(track->PT,track->Eta,track->Phi,0);
		trackCluster.push_back(track);
		indexSet.insert(iTrack);
		for(Int_t iTrack2 = 0; iTrack2 < branchTrack->GetEntries(); ++iTrack2)
		{
		    if (iTrack2 == iTrack) continue;
		    if (indexSet.find(iTrack2) != indexSet.end()) continue;
		    Track *track2 = (Track*) branchTrack->At(iTrack2);
		    if (track2->PT < ptThresh || fabs(track2->Eta) > etaThresh || track2->T > 1) continue;
		    TLorentzVector trackVector2;
		    trackVector2.SetPtEtaPhiM(track2->PT,track2->Eta,track2->Phi,0);
		    // std::cout << trackVector2.DeltaR(trackVector) << std::endl;
		    if (trackVector2.DeltaR(trackVector) < deltaRThresh)
		    {
			trackCluster.push_back(track2);
			indexSet.insert(iTrack2);
		    }
		}
		double maxTime = 0;
		double meanTime = 0;
		double rmsTime = 0;
		double clusterSize = 0;
		double medianTime = 0;
		double modeTime = 0;
		getTrackVars(trackCluster,maxTime,meanTime,modeTime,medianTime,rmsTime,clusterSize);
		if (meanTime > meanTimeMax) meanTimeMax = meanTime;
		histMaxTime->Fill(maxTime);
		histMeanTime->Fill(meanTime);
		histMedianTime->Fill(medianTime);
		histMeanRMSTime->Fill(meanTime,rmsTime);
		histMeanTimeNTracks->Fill(meanTime,clusterSize);
		histMedianTimeNTracks->Fill(medianTime,clusterSize);
		histRMSTime->Fill(rmsTime);
		histNTrack->Fill(clusterSize);
		trackClusters.push_back(trackCluster);
	    }
	    histNJets->Fill(trackClusters.size());
	    histMeanTimeMax->Fill(meanTimeMax);
	    trackClusters.clear();
	}

    }
    // Show resulting histograms
    TFile * outFile = new TFile("TimingTriggerAnalysis/output.root","RECREATE");
    outFile->cd();
    histMaxTime->Write();
    histMeanTime->Write();
    histMeanTimeMax->Write();
    histMedianTime->Write();
    histMeanRMSTime->Write();
    histMeanTimeNTracks->Write();
    histMedianTimeNTracks->Write();
    histRMSTime->Write();
    histNTrack->Write();
    histNJets->Write();
    outFile->Close();
}

