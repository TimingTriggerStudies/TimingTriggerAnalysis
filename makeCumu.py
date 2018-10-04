import ROOT as r
import sys,os

r.gROOT.SetBatch()
r.gStyle.SetOptStat(0)

if len(sys.argv) != 3:
    print "Usage: python makeCumu.py <inputName> <outputDirName>"

inputFile = r.TFile(sys.argv[1])

outputDir=sys.argv[2]

if not os.path.exists(outputDir):
    os.mkdir(outputDir)

outputFile = r.TFile("cumuOutput.root","RECREATE")

perJetHist = inputFile.Get("meanTime")
perEventHist = inputFile.Get("meanTimeMax")

outputFile.cd()
for hist in [perJetHist,perEventHist]:
    hist.SetTitle("")
    hist.GetXaxis().SetTitle("Time Delay Cut")
    tC = r.TCanvas()
    hist.Draw()
    tC.SaveAs(outputDir+"/"+hist.GetName()+".pdf")
    tC.Clear()
    totalN = hist.GetEntries()
    cumu = hist.GetCumulative(False)
    cumu.Scale(1./totalN)
    sup = cumu.Clone("backgroundSuppression_"+hist.GetName())
    sup.Divide(cumu)
    sup.Divide(cumu)
    sup.Draw()
    sup.GetXaxis().SetRangeUser(0,400E-12)
    sup.GetYaxis().SetTitle("Background Suppression")
    tC.SetLogy()
    tC.SaveAs(outputDir+"/"+sup.GetName()+".pdf")
    tC.Clear()
# outputFile.Close()
