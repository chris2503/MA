#include "TMath.h"
#include "TFile.h"
#include "TTree.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <string>
#include <dirent.h>


void SaveEnergyEntries(std::string dir){
    std::vector<std::string> inputname;
    std::vector<string> savefile;

    inputname = {
                "../" + dir + "/double_beta_decay_114Cd.root",
                "../" + dir + "/double_beta_decay_116Cd.root",
                "../" + dir + "/double_beta_decay_128Te.root",
                "../" + dir + "/double_beta_decay_130Te.root",
                "../" + dir + "/double_beta_decay_70Zn.root", 
                "../" + dir + "/chain_decay_40K_epoxy_coating.root", 
                "../" + dir + "/chain_decay_40K_glyptal_coating.root", 
                "../" + dir + "/chain_decay_232Th_epoxy_coating.root",
                "../" + dir + "/chain_decay_232Th_glyptal_coating.root",
                "../" + dir + "/chain_decay_238U_epoxy_coating.root",
                "../" + dir + "/chain_decay_238U_glyptal_coating.root",
                "../" + dir + "/chain_decay_232Th_delrin_layer.root",
                "../" + dir + "/chain_decay_238U_delrin_layer.root",
                "../" + dir + "/chain_decay_232Th_delrin_screws.root",
                "../" + dir + "/chain_decay_238U_delrin_screws.root"
                };

    savefile = {
                "../" + dir + "/edep_entries/entries_114Cd.txt", 
                "../" + dir + "/edep_entries/entries_116Cd.txt", 
                "../" + dir + "/edep_entries/entries_128Te.txt", 
                "../" + dir + "/edep_entries/entries_130Te.txt", 
                "../" + dir + "/edep_entries/entries_70Zn.txt", 
                "../" + dir + "/edep_entries/entries_40K_epoxy.txt", 
                "../" + dir + "/edep_entries/entries_40K_glyptal.txt", 
                "../" + dir + "/edep_entries/entries_232Th_epoxy.txt",
                "../" + dir + "/edep_entries/entries_232Th_glyptal.txt",
                "../" + dir + "/edep_entries/entries_238U_epoxy.txt",
                "../" + dir + "/edep_entries/entries_238U_glyptal.txt",
                "../" + dir + "/edep_entries/entries_232Th_dplate.txt",
                "../" + dir + "/edep_entries/entries_238U_dplate.txt",
                "../" + dir + "/edep_entries/entries_232Th_dscrew.txt",
                "../" + dir + "/edep_entries/entries_238U_dscrew.txt"
                };
    ofstream myfile;

     for (unsigned int i=0; i<inputname.size(); i++){         
         TFile *input_file = new TFile(inputname[i].c_str(), "READ");
         TTree *input_tree = 
           dynamic_cast<TTree*>(input_file->Get("XDEM_sd_data/XDEMDataTree"));
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
                    if(edep->front()==0){
                        ;
                    }
                    else{
                        myfile << edep->front() << "\t" << crystal_id->front() << "\t" << subcrystal_id->front() << "\t"  << 0 << "\t" << particleID->front() << std::endl;
                    }
                    unsigned int edep_len = edep->size();
                    if(edep_len>1){
                        for(unsigned int j_entries=1; j_entries<edep_len; j_entries++){
                            if(edep->at(j_entries)==0){
                                continue;
                            }
                            else{
                            myfile << edep->at(j_entries) << "\t" << crystal_id->at(j_entries) << "\t" << subcrystal_id->at(j_entries) << "\t" << j_entries << "\t" << particleID->at(j_entries) << std::endl; 
                            }
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
    std::vector<std::string> dir;
    dir = {
            "rootfiles_noguard_06_09_2018",
            "rootfiles_zcut_only_06_09_2018",
            "rootfiles_guard_only_06_09_2018",
            "rootfiles_guard+zcut_06_09_2018"
            };

    for(unsigned int i=0; i<dir.size(); i++){
        SaveEnergyEntries(dir[i]);
    }
}
