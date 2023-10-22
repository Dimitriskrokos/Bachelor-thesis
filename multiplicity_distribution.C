#include <TGraphErrors.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <Rtypes.h>
#include <TAxis.h>
#include <TGraph2D.h>
#include <TH2Poly.h>
#include <TText.h>
#include <TMarker.h>

void multiplicity_distribution() 
{
 	// Open the ROOT file
	TFile *file = new TFile("the_differents_histograms.root", "READ");

    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open the ROOT file." << std::endl;
        return;
    }

    
    // Define the histogram names
    const char* histNames[32] = {
        "h_A_phi01_PCB_Layer",
        "h_A_phi02_PCB_Layer",
        "h_A_phi03_PCB_Layer",
        "h_A_phi04_PCB_Layer",
        "h_A_phi05_PCB_Layer",
        "h_A_phi06_PCB_Layer",
        "h_A_phi07_PCB_Layer",
        "h_A_phi08_PCB_Layer",
        "h_A_phi09_PCB_Layer",
        "h_A_phi10_PCB_Layer",
        "h_A_phi11_PCB_Layer",
        "h_A_phi12_PCB_Layer",
        "h_A_phi13_PCB_Layer",
        "h_A_phi14_PCB_Layer",
        "h_A_phi15_PCB_Layer",
        "h_A_phi16_PCB_Layer",
        "h_C_phi01_PCB_Layer",
        "h_C_phi02_PCB_Layer",
        "h_C_phi03_PCB_Layer",
        "h_C_phi04_PCB_Layer",
        "h_C_phi05_PCB_Layer",
        "h_C_phi06_PCB_Layer",
        "h_C_phi07_PCB_Layer",
        "h_C_phi08_PCB_Layer",
        "h_C_phi09_PCB_Layer",
        "h_C_phi10_PCB_Layer",
        "h_C_phi11_PCB_Layer",
        "h_C_phi12_PCB_Layer",
        "h_C_phi13_PCB_Layer",
        "h_C_phi14_PCB_Layer",
        "h_C_phi15_PCB_Layer",
        "h_C_phi16_PCB_Layer"
    };
    
    
    
    
    
    
double Multiplicity_L_A[8][8][8];//(sector,pcb,Layer)    
double Multiplicity_L_C[8][8][8];
double Multiplicity_S_A[8][8][8];//(sector,pcb,Layer)    
double Multiplicity_S_C[8][8][8];


TLatex* binLabel[32];
TH2F* hist[32];

//TCanvas *CANVAS[32];


for (int sector = 0; sector < 32; ++sector) {
//	CANVAS[sector] = new TCanvas(Form("canvas%d", sector + 1), Form("canvas%d", sector + 1), 800, 600);
	hist[sector] = dynamic_cast<TH2F*>(file->Get(histNames[sector]));
	
//	CANVAS[sector]->cd();
//	hist[sector] -> Draw();
}



///////////  Matrix initialization  ////////////////////
for(int sector = 0;sector < 8 ;sector++)
	for(int layer = 0;layer < 8 ;layer++){
		for(int PCB = 0;PCB < 8 ;PCB++){
			Multiplicity_L_A[sector][PCB][layer] = 0.;
			Multiplicity_L_C[sector][PCB][layer] = 0.;
			Multiplicity_S_A[sector][PCB][layer] = 0.;
			Multiplicity_S_C[sector][PCB][layer] = 0.;
		}
}




////////////////////////////////////////////////////////

//A_side 
for(int PCB = 1;PCB <= 8 ;PCB++){
	for(int layer = 1;layer <= 8 ;layer++){
		for(int sector = 0;sector < 16 ;sector++){
			int p = sector / 2;
        		if (sector - 2 * p == 0) {
        			double content = hist[sector]->GetBinContent(PCB, layer);
				Multiplicity_L_A[sector/2][PCB-1][layer-1] = content;	
        		}	
			else{
				double content = hist[sector]->GetBinContent(PCB, layer);
				Multiplicity_S_A[sector/2][PCB-1][layer-1] = content;
			}			
		}
	}
}
//////////////////////////////////////////////////////////////

//(sector,pcb,Layer)

////////////////////////////////////////////////////////////////////////////
/*
for(int PCB = 0;PCB < 8 ;PCB++){
	for(int sector = 0;sector < 8; sector++){
		for(int layer = 0;layer < 8 ;layer++){
			cout<<"\nA("<<"sector = "<<sector + 1<<","<<"PCB = "<<PCB + 1<<","<<"layer = "<<layer + 1<<") = "<<Multiplicity_L_A[sector][PCB][layer];
		}
	}
}

*/	
////////////////////////////////////////////////////////////////////////////

//C side
for(int PCB = 1;PCB <= 8 ;PCB++){
	for(int layer = 1;layer <= 8 ;layer++){
		for(int sector = 16;sector < 32 ;sector++){
			int p = sector / 2;
        		if (sector - 2 * p == 0) {
        			double content = hist[sector]->GetBinContent(PCB, layer);
				Multiplicity_L_C[(sector-16)/2][PCB-1][layer-1] = content;	
	        		//cout<<"\nC("<<sector-16<<","<<PCB<<","<<layer<<") = "<<Multiplicity_L_C[sector-16][PCB-1][layer-1];
        		}	
			else{
				double content = hist[sector]->GetBinContent(PCB, layer);
				Multiplicity_S_C[(sector-16)/2][PCB-1][layer-1] = content;
			}			
		}
	}
}


////////////////////////////////////////////////////////////////////////////

/*
for(int PCB = 0;PCB < 8 ;PCB++){
	for(int layer = 0;layer < 8;layer++){
		for(int sector = 0;sector < 8 ;sector++){
			cout<<"\nC("<<sector<<","<<PCB<<","<<layer<<") = "<<Multiplicity_L_C[sector][PCB][layer];
		}
	}
}
*/	
////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////// Create graph for each PCB where y = multiplicity of each layer(1 to 8) and x =layer which is group in sectors ///////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





TCanvas *canvas_L[8];
TCanvas *canvas_S[8];
TGraph *graph_L_A[8];
TGraph *graph_S_A[8];
TGraph *graph_L_C[8];
TGraph *graph_S_C[8];



Color_t colors[16] = {kBlue, kRed, kGreen, kOrange, kMagenta, kCyan, kYellow, kViolet, kPink, kAzure, kTeal, kSpring, kGray, kBlack, kRed + 2, kBlue + 2};

TLegend *legend_L[8];
TLegend *legend_S[8];


for(int PCB = 0;PCB<8;PCB++){
	legend_S[PCB] = new TLegend(0.8, 0.85, 0.9, 0.9);
	//legend_S[PCB]->SetHeader("Legend");
	
	legend_L[PCB] = new TLegend(0.8, 0.85, 0.9, 0.9);
	//legend_L[PCB]->SetHeader("Legend");

	
	TString canvasName_L = Form("Multiplicity_per_layer_for_PCB%d_at_L.png", PCB + 1);
	TString canvasName_S = Form("Multiplicity_per_layer_for_PCB%d_at_S.png", PCB+1);
	
	canvas_L[PCB] = new TCanvas(canvasName_L, Form("Multiplicity_per_layer_for_PCB%d_at_L.png", PCB + 1), 800, 600);
	canvas_S[PCB] = new TCanvas(canvasName_S, Form("Multiplicity_per_layer_for_PCB%d_at_S.png", PCB + 1), 800, 600);
	graph_L_A[PCB] = new TGraph();
	graph_L_C[PCB] = new TGraph();
	graph_S_A[PCB] = new TGraph();
	graph_S_C[PCB] = new TGraph();
	
	
	
	
		//cout<<"           For PCB = "<<PCB<<"                 \n\n\n\n";
		Int_t pointIndex = 0;
		int layer = 1;
		for(int sector = 0;sector<8;sector++){
			
			for(int LAYER = 0;LAYER<8;LAYER++){
				graph_L_A[PCB]->SetPoint(graph_L_A[PCB]->GetN(),layer,Multiplicity_L_A[sector][PCB][LAYER]);
				graph_L_C[PCB]->SetPoint(graph_L_C[PCB]->GetN(),layer,Multiplicity_L_C[sector][PCB][LAYER]);
				graph_S_A[PCB]->SetPoint(graph_S_A[PCB]->GetN(),layer,Multiplicity_S_A[sector][PCB][LAYER]);
				graph_S_C[PCB]->SetPoint(graph_S_C[PCB]->GetN(),layer,Multiplicity_S_C[sector][PCB][LAYER]);
				//graph_L[PCB]->SetPoint(graph_L[PCB]->GetN(),layer,Multiplicity_L_C[sector][PCB][LAYER]);
				//graph_L[PCB]->SetPointColor(pointIndex, kBlue + sector);
				layer++;
				pointIndex++;
			}
			layer+=20;
		}

	
	
	
	
	
	
	graph_L_A[PCB]->SetMarkerStyle(20);
    	graph_L_A[PCB]->SetTitle(Form("Multiplicity_per_layer_for_PCB%d_at_L", PCB + 1));
    	graph_L_A[PCB]->SetMarkerColor(colors[0]);
    	graph_L_A[PCB]->GetXaxis()->SetTitle("layer");
	graph_L_A[PCB]->GetYaxis()->SetTitle("Multiplicity");
	graph_S_A[PCB]->SetMarkerStyle(20);
    	graph_S_A[PCB]->SetTitle(Form("Multiplicity_per_layer_for_PCB%d_at_S", PCB + 1));
    	graph_S_A[PCB]->SetMarkerColor(colors[1]);
    	graph_S_A[PCB]->GetXaxis()->SetTitle("layer");
	graph_S_A[PCB]->GetYaxis()->SetTitle("Multiplicity");
	graph_L_C[PCB]->SetMarkerStyle(20);
    	graph_L_C[PCB]->SetTitle(Form("Multiplicity_per_layer_for_PCB%d_at_L", PCB + 1));
    	graph_L_C[PCB]->SetMarkerColor(colors[2]);
    	graph_L_C[PCB]->GetXaxis()->SetTitle("layer");
	graph_L_C[PCB]->GetYaxis()->SetTitle("Multiplicity");
	graph_S_C[PCB]->SetMarkerStyle(20);
    	graph_S_C[PCB]->SetTitle(Form("Multiplicity_per_layer_for_PCB%d_at_S", PCB + 1));
    	graph_S_C[PCB]->SetMarkerColor(colors[13]);
    	graph_S_C[PCB]->GetXaxis()->SetTitle("layer");
	graph_S_C[PCB]->GetYaxis()->SetTitle("Multiplicity");

	
	
	canvas_L[PCB] ->cd();
	graph_L_A[PCB]->Draw("AP");
	graph_L_C[PCB]->Draw("Psame");
	
	legend_L[PCB]->AddEntry(graph_L_A[PCB], "sectors of A side", "p");
	legend_L[PCB]->AddEntry(graph_L_C[PCB], "sectors of C side", "p");
	legend_L[PCB]->Draw();
	
	canvas_S[PCB] ->cd();
	graph_S_A[PCB]->Draw("AP");
	graph_S_C[PCB]->Draw("Psame");

	legend_S[PCB]->AddEntry(graph_S_A[PCB], "sectors of A side", "p");
	legend_S[PCB]->AddEntry(graph_S_C[PCB], "sectors of C side", "p");
	
	legend_S[PCB]->Draw();

	
	canvas_L[PCB] ->SaveAs(canvasName_L);
	canvas_S[PCB] ->SaveAs(canvasName_S);

	canvas_L[PCB] ->Close();
	canvas_S[PCB] ->Close();


}


}








