#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TVectorT.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <string>

using namespace std;

// ================================================================================================
int main(int argc, char ** argv)
{
	if ( argc<2 )
	{
		cerr << "Wrong number of arguments. Instead use\n"
			<< "   combiner [list of runs...]\n\n";
		return -1;
	}

	TVectorT<double> kinematics(5);
	TVectorT<double> totalCharge(1);
	totalCharge[0] = 0.;

	// Establish the output file
	TFile * outfile = new TFile("combine_out.root","RECREATE");
	TTree * outtree = new TTree("sk","skimmed tree");
	
	ofstream outtxt ("eventCounts.txt");

	double 	sk_e_cer		,
		sk_e_prl1		,	sk_e_prl2		,
		sk_e_ntrk		,	
		sk_e_ytar		,	
		sk_e_yptar		,	
		sk_gold_e_yptar		,	
		sk_e_xptar		,	
		sk_gold_e_xptar		,	
		sk_e_delta		,	
		sk_gold_e_delta		,	
		sk_e_mom		,	
		sk_gold_e_mom		,	
		sk_e_beta		,	
		sk_gold_e_beta		,	
		sk_e_x			,	
		sk_e_y			,	
		sk_e_z			,	
		sk_e_px			,	sk_e_py			,
		sk_e_pz			,	
		sk_gold_e_px		,	
		sk_gold_e_py		,	
		sk_gold_e_pz		,	
		sk_t1			,	sk_t4			,
		sk_tcoinc		,
		sk_Q2			,	sk_W2			,
		sk_Nu			,	sk_ph_q			,
		sk_th_q			,	sk_xB			,
		sk_q3m			,	sk_q_x			,
		sk_q_y			,	sk_q_z			,
		sk_e_th			,
		sk_rast_Q2		,	sk_rast_W2		,
		sk_rast_Nu		,	sk_rast_ph_q		,
		sk_rast_th_q		,	sk_rast_xB		,
		sk_rast_q3m		,	sk_rast_q_x		,
		sk_rast_q_y		,	sk_rast_q_z		,		
		sk_e_rast_th		,
		sk_exRa_Q2		,	sk_exRa_W2		,
		sk_exRa_Nu		,	sk_exRa_ph_q		,
		sk_exRa_th_q		,	sk_exRa_xB		,
		sk_exRa_q3m		,	sk_exRa_q_x		,
		sk_exRa_q_y		,	sk_exRa_q_z		,
		sk_e_exRa_th		,
		sk_e_ext_deltaDp	,	sk_e_ext_deltaP		,
		sk_e_ext_deltaTh	,	sk_e_ext_delta		,
		sk_e_ext_mom		,	sk_e_ext_yptar		,
		sk_e_ext_xptar		,	sk_e_ext_px		,
		sk_e_ext_py		,	sk_e_ext_pz		,		
		sk_kin_eThetaC		,	sk_kin_pThetaC		,
		sk_kin_eMomC		,	sk_kin_pMomC		,
		sk_kin_Ebeam		,	sk_kin_Q		,
		sk_BCMcharge		,	sk_BCMcurr		,
		sk_BCMrenew 		,	sk_evtypebits		;


	outtree -> Branch("L_prl1"		,&sk_e_prl1		,"L_prl1/D"		);
	outtree -> Branch("L_prl2"		,&sk_e_prl2		,"L_prl2/D"		);

	outtree -> Branch("L_cer"		,&sk_e_cer		,"L_cer/D"		);

	outtree -> Branch("L_ntrk"		,&sk_e_ntrk		,"L_ntrk/D"		);

	outtree -> Branch("L_ytar"		,&sk_e_ytar		,"L_ytar/D"		);

	outtree -> Branch("L_yptar"		,&sk_e_yptar		,"L_yptar/D"		);
	outtree -> Branch("L_gold_yptar"	,&sk_gold_e_yptar	,"L_gold_yptar/D"	);

	outtree -> Branch("L_xptar"		,&sk_e_xptar		,"L_xptar/D"		);
	outtree -> Branch("L_gold_xptar"	,&sk_gold_e_xptar	,"L_gold_xptar/D"	);
	
	outtree -> Branch("L_dp"		,&sk_e_delta		,"L_dp/D"		);
	outtree -> Branch("L_gold_dp"		,&sk_gold_e_delta	,"L_gold_dp/D"		);
	
	outtree -> Branch("L_mom"		,&sk_e_mom		,"L_mom/D"		);
	outtree -> Branch("L_gold_mom"		,&sk_gold_e_mom		,"L_gold_mom/D"		);

	outtree -> Branch("L_beta"		,&sk_e_beta		,"L_beta/D"		);
	outtree -> Branch("L_gold_beta"		,&sk_gold_e_beta	,"L_gold_beta/D"	);

	outtree -> Branch("L_vx"		,&sk_e_x		,"L_vx/D"		);
	outtree -> Branch("L_vy"		,&sk_e_y		,"L_vy/D"		);
	outtree -> Branch("L_vz"		,&sk_e_z		,"L_vz/D"		);
	
	outtree -> Branch("L_px"		,&sk_e_px		,"L_px/D"		);
	outtree -> Branch("L_py"		,&sk_e_py		,"L_py/D"		);
	outtree -> Branch("L_pz"		,&sk_e_pz		,"L_pz/D"		);
	
	outtree -> Branch("L_gold_px"		,&sk_gold_e_px		,"L_gold_px/D"		);
	outtree -> Branch("L_gold_py"		,&sk_gold_e_py		,"L_gold_py/D"		);
	outtree -> Branch("L_gold_pz"		,&sk_gold_e_pz		,"L_gold_pz/D"		);

	outtree -> Branch("t1"			,&sk_t1			,"t1/D"			);
	outtree -> Branch("t4"			,&sk_t4			,"t4/D"			);
	outtree -> Branch("tcoinc"		,&sk_tcoinc		,"tcoinc/D"		);
	
	
	outtree -> Branch("Q2"			,&sk_Q2			,"Q2/D"			);
	outtree -> Branch("W2"			,&sk_W2			,"W2/D"			);
	outtree -> Branch("Nu"			,&sk_Nu			,"Nu/D"			);
	outtree -> Branch("ph_q"		,&sk_ph_q		,"ph_q/D"		);
	outtree -> Branch("th_q"		,&sk_th_q		,"th_q/D"		);
	outtree -> Branch("xB"			,&sk_xB			,"xB/D"			);
	outtree -> Branch("q3m"			,&sk_q3m		,"q3m/D"		);
	outtree -> Branch("q_x"			,&sk_q_x		,"q_x/D"		);
	outtree -> Branch("q_y"			,&sk_q_y		,"q_y/D"		);
	outtree -> Branch("q_z"			,&sk_q_z		,"q_z/D"		);
	outtree -> Branch("L_theta"		,&sk_e_th		,"L_theta/D"		);

	outtree -> Branch("rast_Q2"		,&sk_rast_Q2		,"rast_Q2/D"		);
	outtree -> Branch("rast_W2"		,&sk_rast_W2		,"rast_W2/D"		);
	outtree -> Branch("rast_Nu"		,&sk_rast_Nu		,"rast_Nu/D"		);
	outtree -> Branch("rast_ph_q"		,&sk_rast_ph_q		,"rast_ph_q/D"		);
	outtree -> Branch("rast_th_q"		,&sk_rast_th_q		,"rast_th_q/D"		);
	outtree -> Branch("rast_xB"		,&sk_rast_xB		,"rast_xB/D"		);
	outtree -> Branch("rast_q3m"		,&sk_rast_q3m		,"rast_q3m/D"		);
	outtree -> Branch("rast_q_x"		,&sk_rast_q_x		,"rast_q_x/D"		);
	outtree -> Branch("rast_q_y"		,&sk_rast_q_y		,"rast_q_y/D"		);
	outtree -> Branch("rast_q_z"		,&sk_rast_q_z		,"rast_q_z/D"		);
	outtree -> Branch("rast_L_theta"	,&sk_e_rast_th		,"rast_L_theta/D"	);

	outtree -> Branch("exRa_Q2"		,&sk_exRa_Q2		,"exRa_Q2/D"		);
	outtree -> Branch("exRa_W2"		,&sk_exRa_W2		,"exRa_W2/D"		);
	outtree -> Branch("exRa_Nu"		,&sk_exRa_Nu		,"exRa_Nu/D"		);
	outtree -> Branch("exRa_ph_q"		,&sk_exRa_ph_q		,"exRa_ph_q/D"		);
	outtree -> Branch("exRa_th_q"		,&sk_exRa_th_q		,"exRa_th_q/D"		);
	outtree -> Branch("exRa_xB"		,&sk_exRa_xB		,"exRa_xB/D"		);
	outtree -> Branch("exRa_q3m"		,&sk_exRa_q3m		,"exRa_q3m/D"		);
	outtree -> Branch("exRa_q_x"		,&sk_exRa_q_x		,"exRa_q_x/D"		);
	outtree -> Branch("exRa_q_y"		,&sk_exRa_q_y		,"exRa_q_y/D"		);
	outtree -> Branch("exRa_q_z"		,&sk_exRa_q_z		,"exRa_q_z/D"		);
	outtree -> Branch("exRa_L_theta"	,&sk_e_exRa_th		,"exRa_L_theta/D"	);

	outtree -> Branch("L_ext_delta_dp"	,&sk_e_ext_deltaDp	,"L_ext_delta_dp/D"	);
	outtree -> Branch("L_ext_delta_p"	,&sk_e_ext_deltaP	,"L_ext_delta_p/D"	);
	outtree -> Branch("L_ext_delta_yptar"	,&sk_e_ext_deltaTh	,"L_ext_delta_yptar/D"	);
	outtree -> Branch("L_ext_dp"		,&sk_e_ext_delta	,"L_ext_dp/D"		);
	outtree -> Branch("L_ext_mom"		,&sk_e_ext_mom		,"L_ext_mom/D"		);
	outtree -> Branch("L_ext_yptar"		,&sk_e_ext_yptar	,"L_ext_yptar/D"	);
	outtree -> Branch("L_ext_xptar"		,&sk_e_ext_xptar	,"L_ext_xptar/D"	);
	outtree -> Branch("L_ext_px"		,&sk_e_ext_px		,"L_ext_px/D"		);
	outtree -> Branch("L_ext_py"		,&sk_e_ext_py		,"L_ext_py/D"		);
	outtree -> Branch("L_ext_pz"		,&sk_e_ext_pz		,"L_ext_pz/D"		);	


	outtree -> Branch("Kin_L_thetaC"	,&sk_kin_eThetaC	,"Kin_L_thetaC/D"	);
	outtree -> Branch("Kin_R_thetaC"	,&sk_kin_pThetaC	,"Kin_R_thetaC/D"	);
	outtree -> Branch("Kin_L_momC"		,&sk_kin_eMomC		,"Kin_L_momC/D"		);
	outtree -> Branch("Kin_R_momC"		,&sk_kin_pMomC		,"Kin_R_momC/D"		);
	outtree -> Branch("Kin_Ebeam"		,&sk_kin_Ebeam		,"Kin_Ebeam/D"		);
	outtree -> Branch("Kin_Q"		,&sk_kin_Q		,"Kin_Q/D"		);
	
	outtree -> Branch("BCM_curr"		,&sk_BCMcurr		,"BCM_curr/D"		);
	outtree -> Branch("BCM_charge"		,&sk_BCMcharge		,"BCM_charge/D"		);
	outtree -> Branch("BCM_isrenew"		,&sk_BCMrenew		,"BCM_isrenew/D"	);
	outtree -> Branch("evtypebits"          ,&sk_evtypebits         ,"evtypebits/D"         );

	int totEv 	   = 0; // number of events for individual file
	int totEv_050to060 = 0;
	int totEv_060to070 = 0;
	int totEv_070to080 = 0;
	int totEv_080to090 = 0;
	int totEv_090to100 = 0;
	int totEv_100to110 = 0;
	int totEv_110to120 = 0;
	int totEv_120to130 = 0;
	int totEv_130to140 = 0;
	int totEv_140to150 = 0;
	int totEv_150to160 = 0;
	int totEv_160to170 = 0;
	int totEv_170to180 = 0;
	int totEv_180to190 = 0;
	int totEv_190to200 = 0;
	int totEv_200to210 = 0;
	int totEv_210to220 = 0;
	int totEv_220to230 = 0;
	int totEv_230to240 = 0;
	int totEv_240to250 = 0;
	int totEv_above250 = 0;
	// Loop over input files
	for (int i=1 ; i<argc ; i++)
	{
		//cout << "Working on file " << argv[i] << "\n";
		TFile * thisFile = new TFile(argv[i]);

		if (i==1) // Copy over the kinematics
		{
			TVectorT<double> * kin = (TVectorT<double>*)thisFile->Get("kinematics");
			kinematics[0] = (*kin)[0];
			kinematics[1] = (*kin)[1];
			kinematics[2] = (*kin)[2];
			kinematics[3] = (*kin)[3];
			kinematics[4] = (*kin)[4];
		}

		// Copy over the charge
		TVectorT<double> * charge = (TVectorT<double>*)thisFile->Get("totalCharge");
		totalCharge[0] += (*charge)[0];
		// Set up the tree and branches
		TTree * thisTree = (TTree*) thisFile->Get("sk");


		thisTree -> SetBranchAddress("L_prl1"		,&sk_e_prl1		);
		thisTree -> SetBranchAddress("L_prl2"		,&sk_e_prl2		);

		thisTree -> SetBranchAddress("L_cer"		,&sk_e_cer		);

		thisTree -> SetBranchAddress("L_ntrk"		,&sk_e_ntrk		);

		thisTree -> SetBranchAddress("L_ytar"		,&sk_e_ytar		);

		thisTree -> SetBranchAddress("L_yptar"		,&sk_e_yptar		);
		thisTree -> SetBranchAddress("L_gold_yptar"	,&sk_gold_e_yptar	);

		thisTree -> SetBranchAddress("L_xptar"		,&sk_e_xptar		);
		thisTree -> SetBranchAddress("L_gold_xptar"	,&sk_gold_e_xptar	);
		
		thisTree -> SetBranchAddress("L_dp"		,&sk_e_delta		);
		thisTree -> SetBranchAddress("L_gold_dp"	,&sk_gold_e_delta	);
		
		thisTree -> SetBranchAddress("L_mom"		,&sk_e_mom		);
		thisTree -> SetBranchAddress("L_gold_mom"	,&sk_gold_e_mom		);

		thisTree -> SetBranchAddress("L_beta"		,&sk_e_beta		);
		thisTree -> SetBranchAddress("L_gold_beta"	,&sk_gold_e_beta	);

		thisTree -> SetBranchAddress("L_vx"		,&sk_e_x		);
		thisTree -> SetBranchAddress("L_vy"		,&sk_e_y		);
		thisTree -> SetBranchAddress("L_vz"		,&sk_e_z		);
		
		thisTree -> SetBranchAddress("L_px"		,&sk_e_px		);
		thisTree -> SetBranchAddress("L_py"		,&sk_e_py		);
		thisTree -> SetBranchAddress("L_pz"		,&sk_e_pz		);
		
		thisTree -> SetBranchAddress("L_gold_px"	,&sk_gold_e_px		);
		thisTree -> SetBranchAddress("L_gold_py"	,&sk_gold_e_py		);
		thisTree -> SetBranchAddress("L_gold_pz"	,&sk_gold_e_pz		);

		thisTree -> SetBranchAddress("t1"		,&sk_t1			);
		thisTree -> SetBranchAddress("t4"		,&sk_t4			);
		thisTree -> SetBranchAddress("tcoinc"		,&sk_tcoinc		);
		
		
		thisTree -> SetBranchAddress("Q2"		,&sk_Q2			);
		thisTree -> SetBranchAddress("W2"		,&sk_W2			);
		thisTree -> SetBranchAddress("Nu"		,&sk_Nu			);
		thisTree -> SetBranchAddress("ph_q"		,&sk_ph_q		);
		thisTree -> SetBranchAddress("th_q"		,&sk_th_q		);
		thisTree -> SetBranchAddress("xB"		,&sk_xB			);
		thisTree -> SetBranchAddress("q3m"		,&sk_q3m		);
		thisTree -> SetBranchAddress("q_x"		,&sk_q_x		);
		thisTree -> SetBranchAddress("q_y"		,&sk_q_y		);
		thisTree -> SetBranchAddress("q_z"		,&sk_q_z		);
		thisTree -> SetBranchAddress("L_theta"		,&sk_e_th		);

		thisTree -> SetBranchAddress("rast_Q2"		,&sk_rast_Q2		);
		thisTree -> SetBranchAddress("rast_W2"		,&sk_rast_W2		);
		thisTree -> SetBranchAddress("rast_Nu"		,&sk_rast_Nu		);
		thisTree -> SetBranchAddress("rast_ph_q"	,&sk_rast_ph_q		);
		thisTree -> SetBranchAddress("rast_th_q"	,&sk_rast_th_q		);
		thisTree -> SetBranchAddress("rast_xB"		,&sk_rast_xB		);
		thisTree -> SetBranchAddress("rast_q3m"		,&sk_rast_q3m		);
		thisTree -> SetBranchAddress("rast_q_x"		,&sk_rast_q_x		);
		thisTree -> SetBranchAddress("rast_q_y"		,&sk_rast_q_y		);
		thisTree -> SetBranchAddress("rast_q_z"		,&sk_rast_q_z		);
		thisTree -> SetBranchAddress("rast_L_theta"	,&sk_e_rast_th		);

		thisTree -> SetBranchAddress("exRa_Q2"		,&sk_exRa_Q2		);
		thisTree -> SetBranchAddress("exRa_W2"		,&sk_exRa_W2		);
		thisTree -> SetBranchAddress("exRa_Nu"		,&sk_exRa_Nu		);
		thisTree -> SetBranchAddress("exRa_ph_q"	,&sk_exRa_ph_q		);
		thisTree -> SetBranchAddress("exRa_th_q"	,&sk_exRa_th_q		);
		thisTree -> SetBranchAddress("exRa_xB"		,&sk_exRa_xB		);
		thisTree -> SetBranchAddress("exRa_q3m"		,&sk_exRa_q3m		);
		thisTree -> SetBranchAddress("exRa_q_x"		,&sk_exRa_q_x		);
		thisTree -> SetBranchAddress("exRa_q_y"		,&sk_exRa_q_y		);
		thisTree -> SetBranchAddress("exRa_q_z"		,&sk_exRa_q_z		);
		thisTree -> SetBranchAddress("exRa_L_theta"	,&sk_e_exRa_th		);

		thisTree -> SetBranchAddress("L_ext_delta_dp"	,&sk_e_ext_deltaDp	);
		thisTree -> SetBranchAddress("L_ext_delta_p"	,&sk_e_ext_deltaP	);
		thisTree -> SetBranchAddress("L_ext_delta_yptar",&sk_e_ext_deltaTh	);
		thisTree -> SetBranchAddress("L_ext_dp"		,&sk_e_ext_delta	);
		thisTree -> SetBranchAddress("L_ext_mom"	,&sk_e_ext_mom		);
		thisTree -> SetBranchAddress("L_ext_yptar"	,&sk_e_ext_yptar	);
		thisTree -> SetBranchAddress("L_ext_xptar"	,&sk_e_ext_xptar	);
		thisTree -> SetBranchAddress("L_ext_px"		,&sk_e_ext_px		);
		thisTree -> SetBranchAddress("L_ext_py"		,&sk_e_ext_py		);
		thisTree -> SetBranchAddress("L_ext_pz"		,&sk_e_ext_pz		);	


		thisTree -> SetBranchAddress("Kin_L_thetaC"	,&sk_kin_eThetaC	);
		thisTree -> SetBranchAddress("Kin_R_thetaC"	,&sk_kin_pThetaC	);
		thisTree -> SetBranchAddress("Kin_L_momC"	,&sk_kin_eMomC		);
		thisTree -> SetBranchAddress("Kin_R_momC"	,&sk_kin_pMomC		);
		thisTree -> SetBranchAddress("Kin_Ebeam"	,&sk_kin_Ebeam		);
		thisTree -> SetBranchAddress("Kin_Q"		,&sk_kin_Q		);
		
		thisTree -> SetBranchAddress("BCM_curr"		,&sk_BCMcurr		);
		thisTree -> SetBranchAddress("BCM_charge"	,&sk_BCMcharge		);
		thisTree -> SetBranchAddress("BCM_isrenew"	,&sk_BCMrenew		);
		thisTree -> SetBranchAddress("evtypebits"	,&sk_evtypebits		);

		// Loop over events
		int indEv = 0; // number of events for individual file
		int indEv_050to060 = 0;
		int indEv_060to070 = 0;
		int indEv_070to080 = 0;
		int indEv_080to090 = 0;
		int indEv_090to100 = 0;
		int indEv_100to110 = 0;
		int indEv_110to120 = 0;
		int indEv_120to130 = 0;
		int indEv_130to140 = 0;
		int indEv_140to150 = 0;
		int indEv_150to160 = 0;
		int indEv_160to170 = 0;
		int indEv_170to180 = 0;
		int indEv_180to190 = 0;
		int indEv_190to200 = 0;
		int indEv_200to210 = 0;
		int indEv_210to220 = 0;
		int indEv_220to230 = 0;
		int indEv_230to240 = 0;
		int indEv_240to250 = 0;
		int indEv_above250 = 0;

		double nEvents = thisTree->GetEntries();
		//cout << "Beginning to loop over " << nEvents << " events.\n";
		for (int event=0 ; event < nEvents ; event++)
		{
			thisTree->GetEvent(event);
			outtree->Fill();
			
			// Now let's calculate our current status on numbers
			if( (( int(sk_evtypebits) >> 1 ) & 1)				){
			if ( ((sk_e_prl1 + sk_e_prl2)/1000.)/(sk_e_ext_mom) > 0.8 && 
			     ((sk_e_prl1 + sk_e_prl2)/1000.)/(sk_e_ext_mom) < 1.3 	){
			if( pow(((sk_e_ytar+1.5*sk_e_ext_yptar)/0.08),2) + pow(((1.5*sk_e_ext_xptar)/0.08),2) <= 1 ){
			if( abs(sk_e_ext_delta)<0.045 ) {
			if( sk_e_z > -0.1 && sk_e_z < 0.11 ){
			if( sk_exRa_Q2 > 1.6){
			
				indEv ++;
				if( ( sk_exRa_xB >= 0.50) && (sk_exRa_xB < 0.60) ) indEv_050to060++;
				if( ( sk_exRa_xB >= 0.60) && (sk_exRa_xB < 0.70) ) indEv_060to070++;
				if( ( sk_exRa_xB >= 0.70) && (sk_exRa_xB < 0.80) ) indEv_070to080++;
				if( ( sk_exRa_xB >= 0.80) && (sk_exRa_xB < 0.90) ) indEv_080to090++;
				if( ( sk_exRa_xB >= 0.90) && (sk_exRa_xB < 1.00) ) indEv_090to100++;
				if( ( sk_exRa_xB >= 1.00) && (sk_exRa_xB < 1.10) ) indEv_100to110++;
				if( ( sk_exRa_xB >= 1.10) && (sk_exRa_xB < 1.20) ) indEv_110to120++;
				if( ( sk_exRa_xB >= 1.20) && (sk_exRa_xB < 1.30) ) indEv_120to130++;
				if( ( sk_exRa_xB >= 1.30) && (sk_exRa_xB < 1.40) ) indEv_130to140++;
				if( ( sk_exRa_xB >= 1.40) && (sk_exRa_xB < 1.50) ) indEv_140to150++;
				if( ( sk_exRa_xB >= 1.50) && (sk_exRa_xB < 1.60) ) indEv_150to160++;
				if( ( sk_exRa_xB >= 1.60) && (sk_exRa_xB < 1.70) ) indEv_160to170++;
				if( ( sk_exRa_xB >= 1.70) && (sk_exRa_xB < 1.80) ) indEv_170to180++;
				if( ( sk_exRa_xB >= 1.80) && (sk_exRa_xB < 1.90) ) indEv_180to190++;
				if( ( sk_exRa_xB >= 1.90) && (sk_exRa_xB < 2.00) ) indEv_190to200++;
				if( ( sk_exRa_xB >= 2.00) && (sk_exRa_xB < 2.10) ) indEv_200to210++;
				if( ( sk_exRa_xB >= 2.10) && (sk_exRa_xB < 2.20) ) indEv_210to220++;
				if( ( sk_exRa_xB >= 2.20) && (sk_exRa_xB < 2.30) ) indEv_220to230++;
				if( ( sk_exRa_xB >= 2.30) && (sk_exRa_xB < 2.40) ) indEv_230to240++;
				if( ( sk_exRa_xB >= 2.40) && (sk_exRa_xB < 2.50) ) indEv_240to250++;
				if( ( sk_exRa_xB >= 2.50) ) 	    		   indEv_above250++;
			
			}
			}
			}
			}
			}
			}
			

		}
	
		outtxt        << "\n\n\tFile: " << argv[i] << "\n\tCharge       : " << (*charge)[0] 
					      << "\n\tindEv         : " << indEv  
					      << "\n\tindEv_050to060: " << indEv_050to060
					      << "\n\tindEv_060to070: " << indEv_060to070
					      << "\n\tindEv_070to080: " << indEv_070to080
					      << "\n\tindEv_080to090: " << indEv_080to090
					      << "\n\tindEv_090to100: " << indEv_090to100
					      << "\n\tindEv_100to110: " << indEv_100to110
					      << "\n\tindEv_110to120: " << indEv_110to120
					      << "\n\tindEv_120to130: " << indEv_120to130
					      << "\n\tindEv_130to140: " << indEv_130to140
					      << "\n\tindEv_140to150: " << indEv_140to150
					      << "\n\tindEv_150to160: " << indEv_150to160
					      << "\n\tindEv_160to170: " << indEv_160to170
					      << "\n\tindEv_170to180: " << indEv_170to180
					      << "\n\tindEv_180to190: " << indEv_180to190
					      << "\n\tindEv_190to200: " << indEv_190to200
					      << "\n\tindEv_200to210: " << indEv_200to210
					      << "\n\tindEv_210to220: " << indEv_210to220
					      << "\n\tindEv_220to230: " << indEv_220to230
					      << "\n\tindEv_230to240: " << indEv_230to240
					      << "\n\tindEv_240to250: " << indEv_240to250
					      << "\n\tindEv_above250: " << indEv_above250;

		totEv 		+= indEv 	 ; // number of events for individual file
		totEv_050to060  += indEv_050to060 ;
		totEv_060to070  += indEv_060to070 ;
		totEv_070to080  += indEv_070to080 ;
		totEv_080to090  += indEv_080to090 ;
		totEv_090to100  += indEv_090to100 ;
		totEv_100to110  += indEv_100to110 ;
		totEv_110to120  += indEv_110to120 ;
		totEv_120to130  += indEv_120to130 ;
		totEv_130to140  += indEv_130to140 ;
		totEv_140to150  += indEv_140to150 ;
		totEv_150to160  += indEv_150to160 ;
		totEv_160to170  += indEv_160to170 ;
		totEv_170to180  += indEv_170to180 ;
		totEv_180to190  += indEv_180to190 ;
		totEv_190to200  += indEv_190to200 ;
		totEv_200to210  += indEv_200to210 ;
		totEv_210to220  += indEv_210to220 ;
		totEv_220to230  += indEv_220to230 ;
		totEv_230to240  += indEv_230to240 ;
		totEv_240to250  += indEv_240to250 ;
		totEv_above250  += indEv_above250 ;

		thisFile->Close();
	}


	cout << "\n=================================================== \n";
	cout 			      << "\n\tTotal Charge    : " << totalCharge[0] 
				      << "\n\tTotal Events         : " << totEv 		
				      << "\n\tTotal Events 050to060: " << totEv_050to060  
				      << "\n\tTotal Events 060to070: " << totEv_060to070  
				      << "\n\tTotal Events 070to080: " << totEv_070to080  
				      << "\n\tTotal Events 080to090: " << totEv_080to090  
				      << "\n\tTotal Events 090to100: " << totEv_090to100  
				      << "\n\tTotal Events 100to110: " << totEv_100to110  
				      << "\n\tTotal Events 110to120: " << totEv_110to120  
				      << "\n\tTotal Events 120to130: " << totEv_120to130  
				      << "\n\tTotal Events 130to140: " << totEv_130to140  
				      << "\n\tTotal Events 140to150: " << totEv_140to150  
				      << "\n\tTotal Events 150to160: " << totEv_150to160  
				      << "\n\tTotal Events 160to170: " << totEv_160to170  
				      << "\n\tTotal Events 170to180: " << totEv_170to180  
				      << "\n\tTotal Events 180to190: " << totEv_180to190  
				      << "\n\tTotal Events 190to200: " << totEv_190to200  
				      << "\n\tTotal Events 200to210: " << totEv_200to210  
				      << "\n\tTotal Events 210to220: " << totEv_210to220  
				      << "\n\tTotal Events 220to230: " << totEv_220to230  
				      << "\n\tTotal Events 230to240: " << totEv_230to240  
				      << "\n\tTotal Events 240to250: " << totEv_240to250  
				      << "\n\tTotal Events above250: " << totEv_above250  ;
	cout << "\n=================================================== \n";


	outfile->cd   ();
	outtree->Write();
	totalCharge.Write("totalCharge"); 
	kinematics.Write("kinematics");
	outfile->Close();
	outtxt.close();
	return 0;
}
