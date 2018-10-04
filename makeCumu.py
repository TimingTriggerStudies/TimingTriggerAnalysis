import ROOT as r
r.gROOT.SetBatch()
r.gStyle.SetOptStat(0)

inputFile = r.TFile("output.root")

outputFile = r.TFile("cumuFile.root","RECREATE")

perJetHist = inputFile.Get("meanTime")
perEventHist = inputFile.Get("meanTimeMax")

outputFile.cd()
for hist in [perJetHist,perEventHist]:
    hist.SetTitle("")
    hist.GetXaxis().SetTitle("Time Delay Cut")
    tC = r.TCanvas()
    hist.Draw()
    tC.SaveAs(hist.GetName()+".pdf")
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
    tC.SaveAs(sup.GetName()+".pdf")
    tC.Clear()
# outputFile.Close()
