#include "TMath.h"
#include "TFile.h"
#include "TTree.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <string>

//string inputname = "../../Simulation_Files/double_beta_decay_116Cd.root"
//string savefilename = "../../Calculations/entries_116Cd.txt"


void SaveEnergyEntries(){
    std::vector<std::string> inputname;
    std::vector<string> savefile;

    inputname = {
                "../rootfiles/double_beta_decay_114Cd.root",
                "../rootfiles/double_beta_decay_116Cd.root",
                "../rootfiles/double_beta_decay_128Te.root",
                "../rootfiles/double_beta_decay_130Te.root",
                "../rootfiles/double_beta_decay_70Zn.root", 
                "../rootfiles/chain_decay_40K_epoxy_coating.root", 
                "../rootfiles/chain_decay_40K_glyptal_coating.root", 
                "../rootfiles/chain_decay_232Th_epoxy_coating.root",
                "../rootfiles/chain_decay_232Th_glyptal_coating.root",
                "../rootfiles/chain_decay_238U_epoxy_coating.root",
                "../rootfiles/chain_decay_238U_glyptal_coating.root",
                "../rootfiles/shielding_stdPb_210Pb_chain.root",
                "../rootfiles/shielding_ULAPb_210Pb_chain.root"
                };
    savefile = {
                "../rootfiles/edep_entries/entries_114Cd.txt", 
                "../rootfiles/edep_entries/entries_116Cd.txt", 
                "../rootfiles/edep_entries/entries_128Te.txt", 
                "../rootfiles/edep_entries/entries_130Te.txt", 
                "../rootfiles/edep_entries/entries_70Zn.txt", 
                "../rootfiles/edep_entries/entries_40K_epoxy.txt", 
                "../rootfiles/edep_entries/entries_40K_glyptal.txt", 
                "../rootfiles/edep_entries/entries_232Th_epoxy.txt",
                "../rootfiles/edep_entries/entries_232Th_glyptal.txt",
                "../rootfiles/edep_entries/entries_238U_epoxy.txt",
                "../rootfiles/edep_entries/entries_238U_glyptal.txt",
                "../rootfiles/edep_entries/entries_shielding_stdPb_210Pb_chain.txt",
                "../rootfiles/edep_entries/entries_shielding_ULAPb_210Pb_chain.txt"
                };

    ofstream myfile;

     for (unsigned int i=0; i<inputname.size(); i++){         
         TFile *input_file = new TFile(inputname[i].c_str(), "READ");
         TTree *input_tree = 
           dynamic_cast<TTree*>(input_file->Get("array64SensitiveCrystal_sd_data/array64SensitiveCrystalDataTree"));
         // Check if Tree exists
         if(input_tree == NULL){
             std::cout << "Error, Tree not valid!" << std::endl;
         }
         else{
             // Further I/O
             std::vector<Double_t> *edep = NULL;
             std::vector<Int_t> *crystal_id = NULL;
             std::vector<Int_t> *subcrystal_id = NULL;
             std::vector<Int_t> *particleID = NULL;

             input_tree->SetBranchStatus("*", false);
             input_tree->SetBranchStatus("crystal_id", true);
             input_tree->SetBranchStatus("subcrystal_id", true);
             input_tree->SetBranchStatus("edep", true);
             input_tree->SetBranchStatus("particleID", true);

             input_tree->SetBranchAddress("crystal_id", &crystal_id);
             input_tree->SetBranchAddress("subcrystal_id", &subcrystal_id);
             input_tree->SetBranchAddress("edep", &edep);
             input_tree->SetBranchAddress("particleID", &particleID);

             // Length of loop over all events
             unsigned int signal_counts = input_tree->GetEntries();

             
             myfile.open(savefile[i].c_str(), ios::out);
             if(myfile.is_open()){
                 unsigned int i_entry = 0;
                 std::cout << "Writing edep-Entries to file '" << savefile[i] << "'" << std::endl << std::endl;
                 time_t now = time(0);
                 char* date_time = ctime(&now);
                 myfile << "# Date: " << date_time << std::endl;
                 myfile << "# Energy in MeV" << "\t" << "crystal_id" << "\t" << "subcrystal_id" << "\t"  << "instance" << "\t" << "particleID" << std::endl;
                 for(; i_entry < signal_counts; i_entry++){
                     input_tree->GetEntry(i_entry);
                     myfile << edep->front() << "\t" << crystal_id->front() << "\t" << subcrystal_id->front() << "\t"  << 0 << "\t" << particleID->front() << std::endl;
                     unsigned int edep_len = edep->size();
                     if(edep_len>1){
                        for(unsigned int j_entries=1; j_entries<edep_len; j_entries++){  
                            myfile << edep->at(j_entries) << "\t" << crystal_id->at(j_entries) << "\t" << subcrystal_id->at(j_entries) << "\t" << j_entries << "\t" << particleID->at(j_entries) << std::endl; 
                        }
                     }
                 }
                 //std::cout << i_entry << endl;
                 myfile.close();
             }
             else{
                 std::cout << "Unable to open file " << savefile[i] << std::endl;
             }
         }
     }
}

void ausfuehren(){
    SaveEnergyEntries();
}
