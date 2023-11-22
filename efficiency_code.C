#include <stdio.h> 
#include <stdlib.h> 
#include <math.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TMath.h>

void efficiency_code(){
	 	// Open the ROOT file.In this file,Monte Carlo data are stored in the form of histograms.It contain the total muons per pseudoraditity per phi and the total hits per layer per PCB
	TFile *file = new TFile("MuonTester_Eff_MC_ParticleGun.root", "READ");
	
	// Check if the file was opened successfully
    	if (!file || file->IsZombie()) {
        std::cerr << "Error opening the ROOT file." << std::endl;
        return 1;
    	}
		
	//variables
	//As you can see.i call  h_A_eta_phi6 where 6 index that all muons in the histograms correspond ti hits upper or equal to 6.Respectively for the others name(i.e. h_Efficiency_phi0_eta_A_side.) 
	TH2F *h_Efficiency_phi4_eta_A_side = new TH2F("h_Efficiency_phi4_eta_A_side", "Efficiency functions of eta,phi at A side with hits >= 4 ", 8, 1.3,2.7, 32,-3.2,3.2);
	TH2F *h_Efficiency_phi4_eta_C_side = new TH2F("h_Efficiency_phi4_eta_C_side", "Efficiency functions of eta,phi at C side with hits >= 4 ", 8, -2.7,-1.3, 32,-3.2,3.2);

	TH2F *h_A_eta_phi = (TH2F*)file->Get("h_A_eta_phi");//total muon number for A side
	TH2F *h_C_eta_phi = (TH2F*)file->Get("h_C_eta_phi");//total muon number for C side
	TH2F *h_A_eta_phi4 = (TH2F*)file->Get("h_A_eta_phi4");
	TH2F *h_C_eta_phi4 = (TH2F*)file->Get("h_C_eta_phi4");
	TMatrixD efficiency_error_Propagator_2_matrix_A_phi4(h_A_eta_phi->GetNbinsX(), h_A_eta_phi->GetNbinsY());
	TMatrixD efficiency_error_Propagator_matrix_A_phi4(h_A_eta_phi->GetNbinsX(), h_A_eta_phi->GetNbinsY());
	TMatrixD efficiency_error_Propagator_2_matrix_C_phi4(h_C_eta_phi->GetNbinsX(), h_C_eta_phi->GetNbinsY());
	TMatrixD efficiency_error_Propagator_matrix_C_phi4(h_C_eta_phi->GetNbinsX(), h_C_eta_phi->GetNbinsY());
	TMatrixD h_A_eta_phi4_matrix(h_A_eta_phi4->GetNbinsX(), h_A_eta_phi4->GetNbinsY());
	TMatrixD h_C_eta_phi4_matrix(h_C_eta_phi4->GetNbinsX(), h_C_eta_phi4->GetNbinsY());
	TMatrixD h_A_eta_phi_matrix(h_A_eta_phi->GetNbinsX(), h_A_eta_phi->GetNbinsY());
	TMatrixD h_C_eta_phi_matrix(h_C_eta_phi->GetNbinsX(), h_C_eta_phi->GetNbinsY());
    	TCanvas* canvas34 = new TCanvas("canvas34","Efficiency_A_side_phi4",1800,2000);
    	TCanvas* canvas314 = new TCanvas("canvas314","Efficiency_C_side_phi4",1800,2000);
    	
    	

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////     Efficiency calculation(Divide the histograms)       ////////////////////////////////////  	   	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	

///////////      Efficiency is define as the fraction of muons which has  hits upper or equal 4 with the total mumber	
	
	
	h_Efficiency_phi4_eta_A_side->Divide(h_A_eta_phi4, h_A_eta_phi);
	h_Efficiency_phi4_eta_C_side->Divide(h_C_eta_phi4, h_C_eta_phi);    	
    	

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
////////////////////////        Efficiency error calculation         //////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////





	
	// Fill h_A_eta_phi_matrix
	for (int binx = 1; binx <= h_A_eta_phi->GetNbinsX(); binx++) {
    		for (int biny = 1; biny <= h_A_eta_phi->GetNbinsY(); biny++) {
        		h_A_eta_phi_matrix(binx - 1, biny - 1) = h_A_eta_phi->GetBinContent(binx, biny);
    		}
	}	


	// Fill h_A_eta_phi4_matrix
	for (int binx = 1; binx <= h_A_eta_phi4->GetNbinsX(); binx++) {
    		for (int biny = 1; biny <= h_A_eta_phi4->GetNbinsY(); biny++) {
        		h_A_eta_phi4_matrix(binx - 1, biny - 1) = h_A_eta_phi4->GetBinContent(binx, biny);
    		}
	}

	// Fill h_C_eta_phi4_matrix
	for (int binx = 1; binx <= h_C_eta_phi4->GetNbinsX(); binx++) {
    		for (int biny = 1; biny <= h_C_eta_phi4->GetNbinsY(); biny++) {
        		h_C_eta_phi4_matrix(binx - 1, biny - 1) = h_C_eta_phi4->GetBinContent(binx, biny);
    		}
	}



    
//////////////////////////////////////////////////////////   efficiency error   //////////////////////////////////////////////////////////////////////////    
    

	canvas34->cd();
	h_Efficiency_phi4_eta_A_side->SetOption("colz");
	h_Efficiency_phi4_eta_A_side->GetXaxis()->SetTitle("Pseudorapidity");
	h_Efficiency_phi4_eta_A_side->GetYaxis()->SetTitle("#phi (rad)");
	h_Efficiency_phi4_eta_A_side->SetStats(0);
	h_Efficiency_phi4_eta_A_side->Draw();


	TLatex* binLabel41 = new TLatex();
	binLabel41->SetTextSize(0.02);
	binLabel41->SetTextAlign(22);
	binLabel41->SetTextColor(kBlack);

	for (int binx = 1; binx <= h_Efficiency_phi4_eta_A_side->GetNbinsX(); binx++) {
    		for (int biny = 1; biny <= h_Efficiency_phi4_eta_A_side->GetNbinsY(); biny++) {
        		double content = h_Efficiency_phi4_eta_A_side->GetBinContent(binx, biny);
        		double de_dN = 1.0 / h_A_eta_phi_matrix(binx - 1, biny - 1);
        		double delta_N = TMath::Sqrt(h_A_eta_phi4_matrix(binx - 1, biny - 1));
        		double delta_Ntotal = TMath::Sqrt(h_A_eta_phi_matrix(binx - 1, biny - 1));
        		double de_dNtotal = -h_A_eta_phi4_matrix(binx - 1, biny - 1) / (h_A_eta_phi_matrix(binx - 1, biny - 1) * h_A_eta_phi_matrix(binx - 1, biny - 1));

        		efficiency_error_Propagator_matrix_A_phi4(binx-1, biny-1) = (1/TMath::Sqrt(h_A_eta_phi_matrix(binx - 1, biny - 1)))*TMath::Sqrt(content*(1-content));

        		// Place the content directly on the bin at the center
        		double x = h_Efficiency_phi4_eta_A_side->GetXaxis()->GetBinCenter(binx);
        		double y = h_Efficiency_phi4_eta_A_side->GetYaxis()->GetBinCenter(biny);
        		binLabel41->DrawLatex(x, y, Form("%.4f +/- %.4f", content, efficiency_error_Propagator_matrix_A_phi4(binx - 1, biny - 1)));
        
        		//cout << "The efficiency at (" << binx << "," << biny << ") is: " << content << "+/-" << efficiency_error_Propagator_matrix_A_phi4(binx - 1, biny - 1) << " and the bin content of phi4A is "<<h_A_eta_phi4->GetBinContent(binx , biny) << " while the bin content of phiA is " << h_A_eta_phi->GetBinContent(binx , biny ) <<endl;
    		}
	}




	canvas314->cd();
	h_Efficiency_phi4_eta_C_side->SetOption("colz");
	h_Efficiency_phi4_eta_C_side->GetXaxis()->SetTitle("Pseudorapidity");
	h_Efficiency_phi4_eta_C_side->GetYaxis()->SetTitle("#phi (rad)");
	h_Efficiency_phi4_eta_C_side->SetStats(0);
	h_Efficiency_phi4_eta_C_side->Draw();


	TLatex* binLabel42 = new TLatex();
	binLabel42->SetTextSize(0.02);
	binLabel42->SetTextAlign(22);
	binLabel42->SetTextColor(kBlack);
		
	for (int binx = 1; binx <= h_Efficiency_phi4_eta_C_side->GetNbinsX(); binx++) {
	    for (int biny = 1; biny <= h_Efficiency_phi4_eta_C_side->GetNbinsY(); biny++) {
        	double content = h_Efficiency_phi4_eta_C_side->GetBinContent(binx, biny);

	        efficiency_error_Propagator_matrix_C_phi4(binx-1, biny-1) = (1/TMath::Sqrt(h_A_eta_phi_matrix(binx - 1, biny - 1)))*TMath::Sqrt(content*(1-content));

	        // Place the content directly on the bin at the center
	        double x = h_Efficiency_phi4_eta_C_side->GetXaxis()->GetBinCenter(binx);
	        double y = h_Efficiency_phi4_eta_C_side->GetYaxis()->GetBinCenter(biny);
	        binLabel42->DrawLatex(x, y, Form("%.4f +/- %.4f", content, efficiency_error_Propagator_matrix_C_phi4(binx - 1, biny - 1)));
	       //cout << "The efficiency at (" << binx << "," << biny << ") is: " << content << "+/-" << efficiency_error_Propagator_matrix_C_phi4(binx - 1, biny - 1) << " and the bin content of phi4C is "<<h_C_eta_phi4->GetBinContent(binx , biny ) << " while the bin content of phiC is " << h_C_eta_phi->GetBinContent(binx, biny) <<endl;

    		}
	}



/*
	//Save the Histograms as png file
	canvas1->cd();
	h_A_eta_phi->SetTitle("Total muon Number  at A Side");
	h_A_eta_phi->SetOption("colz");
	h_A_eta_phi->GetXaxis()->SetTitle("Pseudorapidiy");
	h_A_eta_phi->GetYaxis()->SetTitle("#phi (rad)");
	h_A_eta_phi->SetStats(0);
	h_A_eta_phi->Draw();
	canvas1->SaveAs("Total muon Number  at A Side.png");
	
	canvas2->cd();
	h_A_eta_phi4->SetTitle("Muon distribution for Number of Hits >=4 at A Side");
	h_A_eta_phi4->SetOption("colz");
	h_A_eta_phi4->GetXaxis()->SetTitle("Pseudorapidity");
	h_A_eta_phi4->GetYaxis()->SetTitle("#phi (rad)");
	h_A_eta_phi4->SetStats(0);
	h_A_eta_phi4->Draw();
	canvas2->SaveAs("Muon distribution for muons with hits >=4 at A Side.png");	
*/
/*
	canvas30->SaveAs("Efficiency_A_hits == 0.png");
	canvas31->SaveAs("Efficiency_A_hits >= 1.png");
	canvas32->SaveAs("Efficiency_A_hits >= 2.png");
	canvas33->SaveAs("Efficiency_A_hits >= 3.png");
	canvas34->SaveAs("Efficiency_A_hits >= 4.png");
	canvas35->SaveAs("Efficiency_A_hits >= 5.png");
	canvas36->SaveAs("Efficiency_A_hits >= 6.png");
	canvas37->SaveAs("Efficiency_A_hits >= 7.png");
	canvas38->SaveAs("Efficiency_A_hits >= 8.png");
	canvas310->SaveAs("Efficiency_C_side_with_hits == 0.png");
	canvas311->SaveAs("Efficiency_C_side_with_hits >= 1.png");
	canvas312->SaveAs("Efficiency_C_side_with_hits >= 2.png");
	canvas313->SaveAs("Efficiency_C_side_with_hits >= 3.png");
	canvas314->SaveAs("Efficiency_C_side_with_hits >= 4.png");
	canvas315->SaveAs("Efficiency_C_side_with_hits >= 5.png");
	canvas316->SaveAs("Efficiency_C_side_with_hits >= 6.png");
	canvas317->SaveAs("Efficiency_C_side_with_hits >= 7.jpg");
	canvas318->SaveAs("Efficiency_C_side_with_hits >= 8.png");
	
*/	




}
