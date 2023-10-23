void getchargeandtime(int runNum, Double_t &total_charge_bcm1_cut,Double_t &total_charge_bcm2_cut,
Double_t &total_charge_bcm4a_cut,Double_t &total_charge_bcm4c_cut,Double_t &total_charge_bcm1_4acut,Double_t &total_charge_bcm2_4acut,Double_t &total_charge_bcm4c_4acut,Double_t &total_time_bcm1_cut,
Double_t &total_time_bcm2_cut,Double_t &total_time_bcm4a_cut,Double_t &total_time_bcm4c_cut,Double_t &total_time,Double_t &total_time_bcm1_4acut,Double_t &total_time_bcm2_4acut,Double_t &total_time_bcm4c_4acut)
{
  //Created and Read Rootfile
  
  TString filename=Form("../boiling/ROOTfiles/NPS/SCALERS/nps_replay_scalers_%d_1_-1.root", runNum); 

  TFile *file = new TFile(filename);
  if (file->IsOpen()!=1){
    return;
  }
  
  TTree *tdata = (TTree*) file->Get("T");    //data TTree
  TTree *tscal = (TTree*) file->Get("TSH");  //Scaler TTree
  if (!tdata || !tscal){
    return;
  }




  //Get Scaler Leafs
  Double_t Scal_evNum;    //event number associated with scaler reads
  tscal->SetBranchAddress("evNumber", &Scal_evNum);
  Double_t  Scal_BCM4A_charge;
  tscal->SetBranchAddress("H.BCM4A.scalerCharge",&Scal_BCM4A_charge);
  Double_t  Scal_BCM4A_current;
  tscal->SetBranchAddress("H.BCM4A.scalerCurrent",&Scal_BCM4A_current);  
  Double_t  Scal_BCM4C_charge;
  tscal->SetBranchAddress("H.BCM4C.scalerCharge",&Scal_BCM4C_charge);
  Double_t  Scal_BCM4C_current;
  tscal->SetBranchAddress("H.BCM4C.scalerCurrent",&Scal_BCM4C_current);
  Double_t  Scal_BCM1_charge;
  tscal->SetBranchAddress("H.BCM1.scalerCharge",&Scal_BCM1_charge);
  Double_t  Scal_BCM1_current;
  tscal->SetBranchAddress("H.BCM1.scalerCurrent",&Scal_BCM1_current); 
  Double_t  Scal_BCM1;
  tscal->SetBranchAddress("H.BCM1.scaler",&Scal_BCM1);
  Double_t Scal_BCM4A;
  tscal->SetBranchAddress("H.BCM4A.scaler",&Scal_BCM4A);
  Double_t  Scal_BCM2_charge;
  tscal->SetBranchAddress("H.BCM2.scalerCharge",&Scal_BCM2_charge);
  Double_t  Scal_BCM2_current;
  tscal->SetBranchAddress("H.BCM2.scalerCurrent",&Scal_BCM2_current); 

  Double_t  Scal_time;
  tscal->SetBranchAddress("H.1MHz.scalerTime",&Scal_time);
  

  //Defive Quantities To Store Previous Reads and cumulative quantities
  Double_t prev_time = 0.;
  Double_t prev_charge_bcm4a = 0.;
  Double_t prev_charge_bcm4c = 0.;
  Double_t prev_charge_bcm1 = 0.;
  Double_t prev_charge_bcm2 = 0.;
  Double_t prev_scaler_bcm1 = 0.;
  Double_t prev_scaler_bcm4a =0;

  Double_t total_charge_bcm4a = 0.;
  Double_t total_charge_bcm4c = 0.;
  Double_t total_charge_bcm1 = 0.;
  Double_t total_charge_bcm2 = 0.;


  //Loop Over Scaler Reads 
  Long64_t scal_entries = tscal->GetEntries();

  Int_t evt_flag_bcm4a[scal_entries];             //Store Flag [0 or 1], to know if scaler read passed current cut (1) or not (0)
  Int_t evt_flag_bcm4c[scal_entries];             //Store Flag [0 or 1], to know if scaler read passed current cut (1) or not (0)
  Int_t evt_flag_bcm1[scal_entries];
  Int_t evt_flag_bcm2[scal_entries]; 

  Int_t scal_evt_num[scal_entries];    //Store Event Associated with Scaler Read
  //cout << "Scaler events"<<scal_entries<<endl; 
  for (int i = 0; i < scal_entries; i++) {
    
    //**NOTE: Each scaler read is associated with as specific event number
    //        as (scaler read 1-> event 1000,  scaler read 2 -> event 2300, ...)
    //        This means events up to 1000 correspond to scaler read 1, ...
   
    tscal->GetEntry(i);
    //Save all no cut quantities.
    total_time = Scal_time;
    total_charge_bcm4a = Scal_BCM4A_charge;
    total_charge_bcm4c = Scal_BCM4C_charge;
    total_charge_bcm1 = Scal_BCM1_charge;
    total_charge_bcm2 = Scal_BCM2_charge;
    //Apply cut
    if(Scal_BCM1_current > 2){
      total_time_bcm1_cut = total_time_bcm1_cut + (Scal_time - prev_time);
	    total_charge_bcm1_cut = total_charge_bcm1_cut + (Scal_BCM1_charge - prev_charge_bcm1);
    }
    if(Scal_BCM2_current > 2){
      total_time_bcm2_cut = total_time_bcm2_cut + (Scal_time - prev_time);
	    total_charge_bcm2_cut = total_charge_bcm2_cut + (Scal_BCM2_charge - prev_charge_bcm2);
    }
    if(Scal_BCM4A_current > 2){
      total_time_bcm4a_cut = total_time_bcm4a_cut + (Scal_time - prev_time);
	    total_charge_bcm4a_cut = total_charge_bcm4a_cut + (Scal_BCM4A_charge - prev_charge_bcm4a);

      total_time_bcm1_4acut = total_time_bcm1_4acut + (Scal_time - prev_time);
      total_charge_bcm1_4acut = total_charge_bcm1_4acut + (Scal_BCM1_charge - prev_charge_bcm1);
      total_time_bcm2_4acut = total_time_bcm2_4acut + (Scal_time - prev_time);
      total_charge_bcm2_4acut = total_charge_bcm2_4acut + (Scal_BCM2_charge - prev_charge_bcm2);
      total_time_bcm4c_4acut = total_time_bcm4c_4acut + (Scal_time - prev_time);
      total_charge_bcm4c_4acut = total_charge_bcm4c_4acut + (Scal_BCM4C_charge - prev_charge_bcm4c);

    }
    if(Scal_BCM4C_current > 2){
      total_time_bcm4c_cut = total_time_bcm4c_cut + (Scal_time - prev_time);
	    total_charge_bcm4c_cut = total_charge_bcm4c_cut + (Scal_BCM4C_charge - prev_charge_bcm4c);
    }
    prev_time = Scal_time;
    prev_charge_bcm4a = Scal_BCM4A_charge;
    prev_charge_bcm4c = Scal_BCM4C_charge;
    prev_charge_bcm1 = Scal_BCM1_charge;
    prev_charge_bcm2 = Scal_BCM2_charge;


  }
}

void BCMcheck(){
  //open the runlist file
  Int_t runNUM;
  string line;
  TString filename = "BCM.dat";
  ifstream ifs;
  ifs.open(filename);

  vector <double> bcm1ratio;
  vector <double> bcm2ratio;
  vector <double> bcm4cratio;

  vector <double> bcm1ratio_4acut;
  vector <double> bcm2ratio_4acut;
  vector <double> bcm4cratio_4acut;

  vector <double> abovetime1;
  vector <double> abovetime2;
  vector <double> abovetime4c;

  vector <double> abovetime1_4acut;
  vector <double> abovetime2_4acut;
  vector <double> abovetime4c_4acut;

  vector <int> runnumlist;


  TCanvas *c_1 = new TCanvas("c_1","",800,600);
  TCanvas *c_2 = new TCanvas("c_2","",800,600);
  vector<TGraph*> graphArray;
  vector<TGraph*> graphArray_4acut;
  vector<TGraph*> graphArray_runnum;
  TGraph *g1;
  TGraph *g2;
  TGraph *g3;
  TGraph *a1;
  TGraph *a2;
  TGraph *a3;
  TGraph *b1;
  TGraph *b2;
  TGraph *b3;
  TGraph *c1;
  TGraph *c2;
  TGraph *c3;
  

  while (getline(ifs, line)){
    //convert run from string to int
    runNUM = stoi(line);
    cout << runNUM << endl;

  Double_t total_charge_bcm1_cut = 0;
  Double_t total_charge_bcm2_cut = 0;
  Double_t total_charge_bcm4a_cut = 0;
  Double_t total_charge_bcm4c_cut = 0;
  Double_t total_time_bcm1_cut = 0;
  Double_t total_time_bcm2_cut = 0;
  Double_t total_time_bcm4a_cut = 0;
  Double_t total_time_bcm4c_cut = 0;
  Double_t total_time = 0;
  Double_t total_charge_bcm1_4acut = 0;
  Double_t total_charge_bcm2_4acut = 0;
  Double_t total_charge_bcm4c_4acut = 0;
  Double_t total_time_bcm1_4acut = 0;
  Double_t total_time_bcm2_4acut = 0;
  Double_t total_time_bcm4c_4acut = 0;


    getchargeandtime(runNUM,total_charge_bcm1_cut,total_charge_bcm2_cut,total_charge_bcm4a_cut,
    total_charge_bcm4c_cut,total_charge_bcm1_4acut, total_charge_bcm2_4acut,total_charge_bcm4c_4acut,
    total_time_bcm1_cut,total_time_bcm2_cut,total_time_bcm4a_cut,total_time_bcm4c_cut,total_time,
    total_time_bcm1_4acut,total_time_bcm2_4acut,total_time_bcm4c_4acut);

    bcm1ratio.push_back(total_charge_bcm1_cut/total_charge_bcm4a_cut);
    bcm2ratio.push_back(total_charge_bcm2_cut/total_charge_bcm4a_cut);
    bcm4cratio.push_back(total_charge_bcm4c_cut/total_charge_bcm4a_cut);
    abovetime1.push_back((total_time_bcm1_cut/total_time)*100);
    abovetime2.push_back((total_time_bcm2_cut/total_time)*100);
    abovetime4c.push_back((total_time_bcm4c_cut/total_time)*100);

    bcm1ratio_4acut.push_back(total_charge_bcm1_4acut/total_charge_bcm4a_cut);
    bcm2ratio_4acut.push_back(total_charge_bcm2_4acut/total_charge_bcm4a_cut);
    bcm4cratio_4acut.push_back(total_charge_bcm4c_4acut/total_charge_bcm4a_cut);
    abovetime1_4acut.push_back((total_time_bcm1_4acut/total_time)*100);
    abovetime2_4acut.push_back((total_time_bcm2_4acut/total_time)*100);
    abovetime4c_4acut.push_back((total_time_bcm4c_4acut/total_time)*100);

    g1 = new TGraph(bcm1ratio.size(),&abovetime1[0],&bcm1ratio[0]);
    g1->SetMarkerColor(1);
    g1->SetMarkerStyle(23);
    g1->SetMarkerSize(1);
    g2 = new TGraph(bcm2ratio.size(),&abovetime2[0],&bcm2ratio[0]);
    g2->SetMarkerColor(2);
    g2->SetMarkerStyle(21);
    g2->SetMarkerSize(1);
    g3 = new TGraph(bcm4cratio.size(),&abovetime4c[0],&bcm4cratio[0]);
    g3->SetMarkerColor(3);
    g3->SetMarkerStyle(22);
    g3->SetMarkerSize(1);   

    a1 = new TGraph(bcm1ratio_4acut.size(),&abovetime1_4acut[0],&bcm1ratio_4acut[0]);
    a1->SetMarkerColor(4);
    a1->SetMarkerStyle(23);
    a1->SetMarkerSize(1);
    a2 = new TGraph(bcm2ratio_4acut.size(),&abovetime2_4acut[0],&bcm2ratio_4acut[0]);
    a2->SetMarkerColor(5);
    a2->SetMarkerStyle(21);
    a2->SetMarkerSize(1);
    a3 = new TGraph(bcm4cratio_4acut.size(),&abovetime4c_4acut[0],&bcm4cratio_4acut[0]);
    a3->SetMarkerColor(6);
    a3->SetMarkerStyle(22);
    a3->SetMarkerSize(1);

    double x = runNUM;

    b1 = new TGraph(bcm1ratio.size(),&x,&bcm1ratio[0]);
    b1->SetMarkerColor(1);
    b1->SetMarkerStyle(23);
    b1->SetMarkerSize(1);
    b2 = new TGraph(bcm2ratio.size(),&x,&bcm2ratio[0]);
    b2->SetMarkerColor(2);
    b2->SetMarkerStyle(21);
    b2->SetMarkerSize(1);
    b3 = new TGraph(bcm4cratio.size(),&x,&bcm4cratio[0]);
    b3->SetMarkerColor(3);
    b3->SetMarkerStyle(22);
    b3->SetMarkerSize(1);

    c1 = new TGraph(bcm1ratio_4acut.size(),&x,&bcm1ratio_4acut[0]);
    c1->SetMarkerColor(4);
    c1->SetMarkerStyle(23);
    c1->SetMarkerSize(1);
    c2 = new TGraph(bcm2ratio_4acut.size(),&x,&bcm2ratio_4acut[0]);
    c2->SetMarkerColor(5);
    c2->SetMarkerStyle(21);
    c2->SetMarkerSize(1);
    c3 = new TGraph(bcm4cratio_4acut.size(),&x,&bcm4cratio_4acut[0]);
    c3->SetMarkerColor(6);
    c3->SetMarkerStyle(22);
    c3->SetMarkerSize(1);

    graphArray_runnum.push_back(b1);
    graphArray_runnum.push_back(b2);
    graphArray_runnum.push_back(b3);
    graphArray_runnum.push_back(c1);
    graphArray_runnum.push_back(c2);
    graphArray_runnum.push_back(c3);
    graphArray.push_back(g1);
    graphArray.push_back(g2);
    graphArray.push_back(g3);
    graphArray_4acut.push_back(a1);
    graphArray_4acut.push_back(a2);
    graphArray_4acut.push_back(a3);
    bcm1ratio.clear();
    bcm2ratio.clear();
    bcm4cratio.clear();
    abovetime1.clear();
    abovetime2.clear();
    abovetime4c.clear();
    bcm1ratio_4acut.clear();
    bcm2ratio_4acut.clear();
    bcm4cratio_4acut.clear();
    abovetime1_4acut.clear();
    abovetime2_4acut.clear();
    abovetime4c_4acut.clear();

  }



for (int i = 0; i < graphArray.size(); i++) {
    graphArray[i]->GetYaxis()->SetRangeUser(0.9, 1.02);
    c_1->cd();
    if (i == 0) {
        graphArray[i]->GetXaxis()->SetLimits(0, 100);
        graphArray[i]->GetXaxis()->SetTitle("Percentage of good beam time");
        graphArray[i]->GetYaxis()->SetTitle("BCM*/BCM4A");
        graphArray[i]->Draw("AP");
        graphArray_4acut[i]->Draw("SP");  // Draw the first graph without "SAME"
    } else {
        graphArray_4acut[i]->Draw("SP");
        graphArray[i]->Draw("SP");  // Draw subsequent graphs with "SAME"
    }
}
for (int j = 0; j < graphArray_runnum.size(); j++){
  graphArray_runnum[j]->GetYaxis()->SetRangeUser(0.9, 1.02);
  c_2->cd();
  if (j == 0) {
    graphArray_runnum[j]->GetXaxis()->SetLimits(1000, 2000);
        graphArray_runnum[j]->GetXaxis()->SetTitle("RunNumber");
        graphArray_runnum[j]->GetYaxis()->SetTitle("BCM*/BCM4A");
        graphArray_runnum[j]->Draw("AP");
    } else {
        graphArray_runnum[j]->Draw("SP");
}
}
    TLegend *leg=new TLegend(.1,.15,.5,.4);
    leg->AddEntry(graphArray[0],"BCM1/BCM4A","p");
    leg->AddEntry(graphArray[1],"BCM2/BCM4A","p");
    leg->AddEntry(graphArray[2],"BCM4C/BCM4A","p");
    leg->AddEntry(graphArray_4acut[0],"BCM1/BCM4A_4acut","p");
    leg->AddEntry(graphArray_4acut[1],"BCM2/BCM4A_4acut","p");
    leg->AddEntry(graphArray_4acut[2],"BCM4C/BCM4A_4acut","p");
    c_1->cd();
    leg->Draw();
    TLegend *leg1=new TLegend(.1,.15,.5,.4);
    c_2->cd();
    leg1->AddEntry(graphArray_runnum[0],"BCM1/BCM4A","p");
    leg1->AddEntry(graphArray_runnum[1],"BCM2/BCM4A","p");
    leg1->AddEntry(graphArray_runnum[2],"BCM4C/BCM4A","p");
    leg1->AddEntry(graphArray_runnum[3],"BCM1/BCM4A_4acut","p");
    leg1->AddEntry(graphArray_runnum[4],"BCM2/BCM4A_4acut","p");
    leg1->AddEntry(graphArray_runnum[5],"BCM4C/BCM4A_4acut","p");
    leg1->Draw();

    string pdfname;
    c_1->Update();
    pdfname = to_string(runNUM) + ".pdf";
    c_1->SaveAs("Plot.pdf");
    c_2->Update();
    c_2->SaveAs("Plot1.pdf");
}