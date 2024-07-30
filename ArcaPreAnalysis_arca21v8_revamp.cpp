// Script for ARCA analysis writen by G.Zarpapis 2023, based on the scripts of A.Sinopoulou, K.Pikounis and K.Tzamariudaki///
// Small changes for snakemake framework by V.Tsourapis 2024
// Revert to v8 and back header reading (July 2024)
// Revised script for all ARCA analysis, June 2024

// C++ Libraries
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// My Libraries
#include "Tools.hh"

// ROOT Libraries
#include "TCanvas.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLeaf.h"
#include "TMath.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

// AANET Libraries
#include <fluxfunctions.h>

//#include "AAObject.hh"
#include "Det.hh"
#include "Evt.hh"
#include "Hit.hh"
#include "Trk.hh"
#include "Trk_util.hh"
#include "Vec.hh"
#include "externals/km3net-dataformat/offline/io_ascii.hh"
#include "externals/km3net-dataformat/tools/reconstruction.hh"

#include "Jeep/JParser.hh"

class ZZPrint {
private:
    int verbosity = 0;
    int persistent_verbosity = 0;
    int cutoff_verbosity_level = 1;

    int* ptr; // This is for differentiating uses of print(num) and print()

public:


    void endl() {
        if (*(this->ptr) <= this->cutoff_verbosity_level) std::cout << "\n";
        std::cout.flush();
    }
    void flush() {
        std::cout.flush();
    }
    void setCutoff(const int& cutoff) { this->cutoff_verbosity_level = cutoff; }
    void setVerbosity(const int& verbosity) {
        this->verbosity = verbosity;
        this->persistent_verbosity = verbosity;
    }


    template <typename T>
    ZZPrint& operator<<(const T& message) {
        if (*(this->ptr) < this->cutoff_verbosity_level) {
            std::cout << message;
        }
        return *this;
    }

    ZZPrint& operator()(const int& verbosity) {
        setVerbosity(verbosity);
        this->ptr = &(this->verbosity);
        return *this;
    }

    ZZPrint& operator()() {
        this->ptr = &(this->persistent_verbosity);
        return *this;
    }

};

std::map<int, std::map<int, std::map<int, vector<Hit>>>> makeHitsMap(const vector<Hit>& hits, Det& det);
//std::vector<float> find_lambda(const ncsr::Vector& muPos, const ncsr::Vector& muDir, const ncsr::Vector& omPos, const float& theta_cherenkov);
//std::vector<float> find_lambda(const Vec& muPos, const Vec& muDir, const Vec& omPos, const float& theta_cherenkov);

bool cherenkov_condition(const double& cos_angle, const double& d_closest, const double& tres);
bool cherenkov_condition(const Trk& track, const Hit& hit);
bool cherenkov_condition_cascade(const Trk& track, const Hit& hit);

float icecube_HESE_flux(float E);
float random_flux(float E);
float icecube_flux_diff2019(float E);
float LoI_flux(float E);
float Antares_flux(float E);
float new_flux(float E);

int main(int argc, char** argv) {
    // Benchmarking
    zarp::Benchmarking clock;
    clock.start();

    //  Cherenkov Constants
    const float dndl = 0.0298;
    const float water_index = 1.3499;
    //const float theta_cher = 42.153902;    // Unused?
    //const float c_water = 0.2222495796574987;    // Unused?
    const float c_light = 299792458 * 1e-9;
    const float v_light = c_light / (water_index + dndl);
    //const float cos_cherenkov_angle = 1.0 / water_index;    // Unused?
    const double reco_Vert_dist_for_coinc = 120;

    // Program input
    std::string inputFileName;
    std::string detectorFile;
    std::string outputFile;
    std::string special_text;
    int run_number;
    int run_subNumber;
    int n_files;
    int verbosity_cutoff;

    // Assign program input to variables
    try {
        JParser<> zap("Arca Anything Software of NCSR Demokritos");

        zap['f'] = make_field(inputFileName) = "";
        zap['a'] = make_field(detectorFile) = "";
        zap['o'] = make_field(outputFile) = (std::string(argv[0]) + ".root").c_str();
        zap['r'] = make_field(run_number) = -1;
        zap['s'] = make_field(run_subNumber) = -1;
        zap['n'] = make_field(n_files) = -1;
        zap['t'] = make_field(verbosity_cutoff) = 1;
        zap['x'] = make_field(special_text) = "arca";

        zap(argc, argv);
    }
    catch (const exception& error) {
        FATAL(error.what() << std::endl);
    }

    // std::cout but you can select what shows and what not
    ZZPrint print_mode;
    print_mode.setVerbosity(0);
    print_mode.setCutoff(verbosity_cutoff);

    // Load Detector
    Det det = Det(detectorFile);

    /*
        1) Extract lines and DOMs from detector file
        2) Calculate detector's center
    */
    int number_of_doms_in_det = 0;
    int number_of_lines_in_det = 0;
    std::vector<int> doms_in_det;
    std::vector<int> lines_in_det;
    std::map<int, std::map<int, Dom*>> dom_map = det.lines_floor_map();  // This has the full information
    Vec mean_detector_pos(0, 0, 0);

    std::vector<int> du_vec;  // Not mine, same as lines_in_det, host it for now ~ G.Z

    for (std::map<int, std::map<int, Dom*>>::iterator i = dom_map.begin(); i != dom_map.end(); ++i) {
        ++number_of_lines_in_det;
        lines_in_det.push_back(i->first);
        du_vec.push_back(i->first);  // Hosting ~ G.Z
        for (std::map<int, Dom*>::iterator k = i->second.begin(); k != i->second.end(); ++k) {
            if (k->first == 0) continue;  // Skip the base module (index = 0, DOMs are 1->18)

            ++number_of_doms_in_det;
            doms_in_det.push_back(k->second->id);

            // Calculate the sum of dom positions
            mean_detector_pos += {k->second->pos.x, k->second->pos.y, k->second->pos.z};
        }
    }

    // Calculate the mean of DOM positions, basically the detector's center
    mean_detector_pos = mean_detector_pos / number_of_doms_in_det;


    Vec det_dom_pos(0, 0, 0);
    double range_detector_max;
    double range_detector_min;
    bool first_search_of_detector_ranges = true;

    for (std::map<int, std::map<int, Dom*>>::iterator i = dom_map.begin(); i != dom_map.end(); ++i) {
        for (std::map<int, Dom*>::iterator k = i->second.begin(); k != i->second.end(); ++k) {
            if (k->first == 0) continue;  // Skip the base module (index = 0, DOMs are 1->18)

            // Initialize max and min ranges of doms from detector
            if (first_search_of_detector_ranges){
                det_dom_pos = Vec(k->second->pos.x, k->second->pos.y, k->second->pos.z);
                range_detector_max = (det_dom_pos - mean_detector_pos).len();
                range_detector_min = (det_dom_pos - mean_detector_pos).len();
                first_search_of_detector_ranges = false;
                continue;
            }
            else{
                det_dom_pos = Vec(k->second->pos.x, k->second->pos.y, k->second->pos.z);
                if ((det_dom_pos - mean_detector_pos).len() > range_detector_max) range_detector_max = (det_dom_pos - mean_detector_pos).len();
                if ((det_dom_pos - mean_detector_pos).len() < range_detector_min) range_detector_min = (det_dom_pos - mean_detector_pos).len();
            }
        }
    }

    // Welcome message
    if (true) {  // Remove detector's inital message
        print_mode(0) << "\033[1A\033[2K\033[1A\033[2K\r";
    }

    print_mode(0) << "\n\t\t\033[1mARCA-" << number_of_lines_in_det << " Pre-Processing\033[22m\n\n";

    // Load input file and check its availability
    TFile* inputFile = new TFile();
    if (!gSystem->AccessPathName(inputFileName.c_str())) {
        std::cout << "Hello1" << std::endl;
        inputFile = TFile::Open(inputFileName.c_str());
        print_mode(0) << "Opening file:\n   \033[1m" << inputFile->GetName() << "\033[22m\n";
    }
    else {
        std::cout << "Hello2" << std::endl;
        print_mode(0) << " - \033[31mERROR\033[0m: there is no file " << inputFile->GetName() << "\n";
        return 0;
    }

    if (!inputFile) {
        std::cout << "Hello3" << std::endl;
        print_mode(0) << " - \033[31mERROR\033[0m: could not open file.\n";
        return 0;
    }
    else {
        std::cout << "Hello4" << std::endl;
        print_mode(0) << "   File is \033[32mopen\033[0m!\n";
    }
    print_mode(0).endl();
    std::cout << "Hello5" << std::endl;

    // Get Event branch from Tree "E"
    Evt* event_pointer = new Evt();
    TTree* inputTree = static_cast<TTree*>(inputFile->Get("E"));
    inputTree->SetBranchAddress("Evt", &event_pointer);
    int number_of_events = inputTree->GetEntries();
    cout << " number_of_events " << number_of_events << endl;

    // ** File from Header ** //

    struct JCan {
        float Zmin;
        float Zmax;
        float Radius;
    };
    JCan Can;

    TString fname = inputFileName;
    std::string can;       // info in the tag  "can" as described in ANTRS_SOFT/2002-004
    std::string DAQ;       // in MUPAGE info in the tag "livetime"
    std::string genvol;    // info in the tag "genvol"
    std::string livetime;  // in MUPAGE info in the tag "livetime"      Same Comment?? ~G.Z.
    std::string spectrum;  // info in the tag "spectrum" ->spectral index of generated neutrino spectrum

    float index;                        // index ->spectral index of generated neutrino spectrum
    float zmin, zmax, rad, gvol, Ngen;  // info in the tag "genvol" as described in ANTRS_SOFT/2002-004
    float LiveT, ELiveT;                // info in the tag  "livetime"-> mupage live time [s] and error [s]

    bool mupage_bool = false;
    bool data_bool = false;
    bool nu_bool = false;
    bool k40_bool = false;

    Head* h = (Head*)inputFile->Get("Head");
    print_mode(0) << "Opening file's header:\n\033[1mkey : value\033[22m\n";
    foreach_map(key, value, *h) {
        print_mode(0) << key << " : " << value << "\n";
	if (key == "genvol") {
            genvol = value;
            nu_bool = true;
        }
	if (key == "K40") {
            k40_bool = true;
        }
        if (key == "spectrum") {
            spectrum = value;
	    //nu_bool = true;
	}
        if (key == "can")
            can = value;
        if (key == "DAQ") {
            DAQ = value;
            data_bool = true;
        }
        if (key == "livetime") {
            livetime = value;
            mupage_bool = true;
        }
    }
    print_mode(0).flush();

    // If file is MUPAGE, can't be data
    if (mupage_bool) data_bool = false;

     // Extract values from head
    std::stringstream ss;
    ss.str(genvol);
    ss >> zmin >> zmax >> rad >> gvol >> Ngen;
    ss.clear(); // Clears the stringstream flags, if fail flags are positive no more input is accepted
    ss.flush(); // Clears the output content

    ss.str(spectrum);
    ss >> index;
    ss.clear();
    ss.flush();

    ss.str(can);
    ss >> Can.Zmin >> Can.Zmax >> Can.Radius;
    ss.clear();
    ss.flush();

    // Extract livetime in cases of mupage or data
    if (mupage_bool) {
        ss.str(livetime);
        ss >> LiveT >> ELiveT;  // live_time + error in seconds
        ss.clear();
        ss.flush();

        for (unsigned i = 0; i < livetime.size(); i++) {
            print_mode(0) << livetime[i];
            print_mode(0).endl();
        }
    }
    else {
        if (data_bool) {
            ss.str(DAQ);
            ss >> LiveT;  // live_time in seconds
            ss.clear();
            ss.flush();
        }
        else {
            LiveT = -1;
            ELiveT = -1;
        }
    }

    //==================== Define Variables of the output==========================================================

    // Generic
    bool initialized_flag; // Helps with first initialization

    // Event
    int pseudo_runid;
    int pseudo_subRunid;
    float pseudo_livetime;
    int evt_id, evt_num;
    int number_of_L0_events = 0;
    int number_of_L1_events = 0;
    int number_of_anti_noise = 0;
    int number_of_reco_tracks = 0;
    float evt_overlays;
    int events_showers = 0;
    int events_showers_MX = 0;
    int events_muons = 0;
    int events_factory_limit = 0;
    bool flag_muon_3D;
    bool flag_shower_MX;
    bool flag_shower_3D;
    bool flag_factory_limit;


    // True Tracks
    float mu_pos[3];
    float mu_dir[3];
    float nu_pos[3];
    float nu_dir[3];
    float E_mu_max, E_mu_depos, E_bundle, E_visible;
    float cos_zen_mu, zenith_mu, logE_mu, logEbundle, logEdepos;
    float logE_nu;
    float E_nu, cos_zen_nu, zenith_nu;
    bool mu_vertex_in, neut_vertex_in;
    float R_true_muon, R_true_nu;


    // Neutrino Weights
    float w1;
    float w2;
    float w3;
    float w3atm;
    float w3iceHESE;
    float w3loi;
    float w3random;
    float w3ice;
    float w3antares;
    float w3new;
    float neutrino_ngen;
    int neutrino_type;
    int neutrino_interaction; //CC or NC


    // Reco Muon Tracks
    float best_trk_pos[3];
    float best_trk_dir[3];
    float best_trk_pos_shower[3];
    float best_trk_dir_shower[3];
    float Ereco_shower, logEreco_shower, Ereco_shower_corrected, Nhit_shower;
    float Ereco, Ereco2;
    float cos_zen, zenith;
    float jlik, jbeta0, jbeta0_deg, logbeta0, GChi, GNhit, EChi, GLik, GNit, Snpe, SnpeT, Slen, ENDF, ENhit, Elen, fi, fi_rad;
    float Q1, Q1new, R_shower;
    float logEreco, logEreco2;
    float multiplicity_muon;      // How many muons per event
    float multiplicity_neutrino;  // How many neutrinos per event
    float R;
    float reco_type;
    float Q, costheta_max, theta_max, costheta_min, theta_min, diff_theta;
    float num_of_good_sol, downsol, upsol;
    float max_lik_down, max_lik_up, lik_sol_down, lik_sol_up;
    float dlik, normdlik;
    float bjorken_y;

    // Reco Shower Tracks
    float zenith_shower, delta_zenith_shower_mu, delta_zenith_shower_nu, lik_shower, beta0_shower, beta0_shower_deg;

    // Reco Track Conditions
    bool found_reco_track;
    bool is_muon_event;
    bool is_shower_event;
    //bool jgandalf_reco_evt; // Unused
    bool LooseShower;
    bool all_in;
    bool reco_vertex_in;


    // Mixed Tracks
    float diffangle, delta_zenith;
    float diffangle_track_shower, delta_zenith_track_shower;
    float diffangle_shower_nu, diffangle_shower_mu;


    // Hits
    float max_ToT_trig, ToT_trig;
    float D, D3d;

    int num_triggered_hits;
    int num_cherenkov_hits;
    int num_cascade_hits;
    std::vector<const Hit*> triggered_hits;
    std::vector<const Hit*> cherenkov_hits;
    std::vector<const Hit*> cascade_hits;
    std::vector<int> triggered_hits_per_event;
    std::vector<int> cherenkov_hits_per_event;
    std::vector<int> cascade_hits_per_event;

    // Maps' first element is the key for DU (a.k. line), (detector) DOM_ID or (global) PMT_ID.
    // Maps' second element is a counter -> how many triggered / cherenkov / cascade hits were recorded.
    int num_triggered_lines;
    int num_triggered_doms;
    int num_triggered_pmts;
    std::map<int, int> triggered_lines;
    std::map<int, int> triggered_doms;
    std::map<int, int> triggered_pmts;
    std::vector<int> triggered_lines_per_event;
    std::vector<int> triggered_doms_per_event;
    std::vector<int> triggered_pmts_per_event;

    int num_cherenkov_lines;
    int num_cherenkov_doms;
    int num_cherenkov_pmts;
    std::map<int, int> cherenkov_lines;
    std::map<int, int> cherenkov_doms;
    std::map<int, int> cherenkov_pmts;
    std::vector<int> cherenkov_lines_per_event;
    std::vector<int> cherenkov_doms_per_event;
    std::vector<int> cherenkov_pmts_per_event;

    int num_cascade_lines;
    int num_cascade_doms;
    int num_cascade_pmts;
    std::map<int, int> cascade_lines;
    std::map<int, int> cascade_doms;
    std::map<int, int> cascade_pmts;
    std::vector<int> cascade_lines_per_event;
    std::vector<int> cascade_doms_per_event;
    std::vector<int> cascade_pmts_per_event;

    // Ratios 
    float ratio_cherenkov_hits; // Cherenkov / Triggered
    float ratio_cherenkov_lines;
    float ratio_cherenkov_doms;
    float ratio_cherenkov_pmts;
    float ratio_cascade_hits; // Cascade / Triggered
    float ratio_cascade_lines;
    float ratio_cascade_doms;
    float ratio_cascade_pmts;

    // Border Hits
    int int_Nborder_hits;
    int int_Nborder_cherenkov_hits;
    int int_Nborder_dtres_hits;
    int int_Nborder_DOMs;
    int int_Nhits_upper;
    int int_Nhits_lower;
    int int_Nhits_border_upper;
    int int_Nhits_border_lower;
    int int_Nhits_cherenkov_upper;
    int int_Nhits_cherenkov_lower;
    int int_Nhits_border_cherenkov_upper;
    int int_Nhits_border_cherenkov_lower;

    float NtrackIT;
    float NtrackEarly;
    float Nallbehind;
    float NtrackLate;
    float Nallbehind2, Nallbehind3, NtrackIT30, NtrackIT10;
    float TrLengthIT;
    float lambdaIT_max;
    float ToT_IT;
    float ToT_allIT;
    float strict_lambda_min;
    float lambdaCH_max;
    float maxML1lamda;
    float minML1lamda;
    float TrLengthIT_2;
    float TrLengthIT_3;
    float Ntrack_all;
    float Ntrack_all2;
    float itoverlen;

    float sum_ToT_mu;   // Sum of the ToT for hits that fulfill the Cerenkov(muon) hypothesis using JGandalf’s best solution.
    float sum_ToT_casc_mu;  // Sum of the ToT for hits that fulfill the Cerenkov(cascade) hypothesis using JGandalf’s best solution
    float sum_ToT_casc;       // Sum of the ToT for hits that fulfill the Cerenkov(cascademu) hypothesis using JGandalf’s best solution.

    float ratio_closehits_muon;
    float ratio_closehits_cascmuon;
    float ratio_closehits_casc;
    float redToT_muon;
    float redToT_cascmuon;
    float redToT_casc;
    float diff_dist_mu;
    float diff_dist_casc_mu;
    float diff_dist_casc;

    float min_diff_sol;
    float max_diff_sol;
    float min_diff_sollik;
    float max_diff_sollik;
    float min_zen_sol;
    float max_zen_sol;
    float nhits_casc_mu, nhits_casc, nhits_mu;
    float nhits_casc_mu_100, nhits_casc_100, nhits_mu_100;
    float min_dist_casc_mu;
    float min_dist_casc;
    float min_dist_mu;
    float max_dist_casc_mu, max_dist_mu, max_dist_casc;
    float ToT_border_casc;
    float ToT_border_mu;
    float ToT_border_cascmu;

    float ratio110;
    float ratio410;
    float ratio130;
    float ratio430;
    float ratio1;
    float ratio2;
    float ratio4;
    float ratio5;
    float ratio6;
    float ratiol;
    float ratiol_trig;
    float ratio4_trig;
    float ratio3;
    float ratio310;
    float ratio330;
    float mean_tres_it;
    float nhits_casc_mu_50;
    float nhits_casc_50;
    float nhits_mu_50;
    float nhits_casc_mu_30;
    float nhits_casc_30;
    float nhits_mu_30;
    float myratio50_muon;
    float myratio50_casc;
    float myratio50_cascmuon;
    float myratio30_muon;
    float myratio30_casc;
    float myratio30_cascmuon;

    float myratio50_casc_over_mu;
    float myratio30_casc_over_mu;
    float myratio50_cascmuon_over_mu;
    float myratio30_cascmuon_over_mu;
    float ratio_closehits_casc_over_mu;
    float ratio_closehits_cascmuon_over_mu;
    float redToT_casc_over_mu;
    float redToT_cascmuon_over_mu;  

    float Nborder_hits;
    float Nborder_cherenkov_hits;
    float Nborder_cherenkov_2lines_hits;
    float Nborder_dtres_hits;
    float Nborder_DOMs;
    float Nhits_upper;
    float Nhits_lower;
    float Nhits_border_upper;
    float Nhits_border_lower;
    float Nhits_cherenkov_upper;
    float Nhits_cherenkov_lower;
    float Nhits_border_cherenkov_upper;
    float Nhits_border_cherenkov_lower;

    float NtrackIT50;

    float NtrackIT50_2;
    float NtrackIT30_2;
    float NtrackIT10_2;

    float NtrackIT50_3;
    float NtrackIT30_3;
    float NtrackIT10_3;


    // Cherenkov Pars for triggered hits
    double reco_time;
    double reco_tres;
    double reco_d_track;
    double reco_d_photon;
    double reco_d_closest;
    double reco_cos_angle;

    //Lilia's variables
    float Log_No_PMTs_with_dist_weights_MAXhits; //var1
    float No_OMs_length_potential_length_chosen; //var2
    float PMTs_hit_to_nohit; //var4
    float MaxLenPos_OvrMaxDist; //var5
    float TrLenIT_3_OvrMaxDist; //var6


    // Cuts: true -> passed, false -> failed
    bool L0_test;
    bool anti_noise;
    bool L1_test;


    // Association of TTrees in Post-Processing
    int associated_n_events = 0;
    int associated_n_hits;


    // Output
    TFile* file_out = new TFile(outputFile.c_str(), "RECREATE");

    /*
        ====================================
        |                                  |
        |           Histograms             |
        |                                  |
        ====================================
    */


    /*
        ====================================
        |                                  |
        |            Branches              |
        |                                  |
        ====================================
    */

    // Quality Parameters Tree
    TTree* qualityParametersTree = new TTree("QualityParameters", "Quality Parameters Tree");

    // Quality Parameters Branches
    qualityParametersTree->Branch("run_number", &run_number, "run_number/I");
    qualityParametersTree->Branch("run_subNumber", &run_subNumber, "run_subNumber/I");
    qualityParametersTree->Branch("livetime", &LiveT, "livetime/F");
    qualityParametersTree->Branch("associated_n_events", &associated_n_events, "associated_n_events/I");

    // Processed Events Tree
    TTree* processedEventsTree = new TTree("ProcessedEvents", "Processed Events Tree");

    // Processed Events Branches
    processedEventsTree->Branch("pseudo_runid", &pseudo_runid, "run_number/I");
    processedEventsTree->Branch("pseudo_subRunid", &pseudo_subRunid, "run_subNumber/I");
    processedEventsTree->Branch("pseudo_livetime", &pseudo_livetime, "pseudo_livetime/F"); 
    processedEventsTree->Branch("evt_id", &evt_id, "evt_id/I"); 
    processedEventsTree->Branch("evt_num", &evt_num, "evt_num/I"); 
    processedEventsTree->Branch("evt_overlays", &evt_overlays, "evt_overlays/F"); 
    processedEventsTree->Branch("reco_type", &reco_type, "reco_type/F"); 
    processedEventsTree->Branch("jbeta0", &jbeta0, "jbeta0/F"); 
    processedEventsTree->Branch("jbeta0_deg", &jbeta0_deg, "jbeta0_deg/F"); 
    processedEventsTree->Branch("D", &D, "D/F"); 
    processedEventsTree->Branch("D3d", &D3d, "D3d/F"); 
    processedEventsTree->Branch("cos_zen", &cos_zen, "cos_zen/F"); 
    processedEventsTree->Branch("fi", &fi, "fi/F"); 
    processedEventsTree->Branch("R", &R, "R/F"); 
    processedEventsTree->Branch("mean_tres_it", &mean_tres_it, "mean_tres_it/F"); 
    processedEventsTree->Branch("zenith", &zenith, "zenith/F"); 
    processedEventsTree->Branch("logEreco", &logEreco, "logEreco/F"); 
    processedEventsTree->Branch("Ereco", &Ereco, "Ereco/F"); 
    processedEventsTree->Branch("E_mu_max", &E_mu_max, "E_mu_max/F");
    processedEventsTree->Branch("E_mu_depos", &E_mu_depos, "E_mu_depos/F");
    processedEventsTree->Branch("E_bundle", &E_bundle, "E_bundle/F");
    processedEventsTree->Branch("E_visible", &E_visible, "E_visible/F");
    processedEventsTree->Branch("E_nu", &E_nu, "E_nu/F");
    processedEventsTree->Branch("multiplicity_muon", &multiplicity_muon, "multiplicity_muon/F");
    processedEventsTree->Branch("multiplicity_neutrino",&multiplicity_neutrino,"multiplicity_neutrino/F");    
    processedEventsTree->Branch("logEreco2", &logEreco2, "logEreco2/F"); 
    processedEventsTree->Branch("logbeta0", &logbeta0, "logbeta0/F"); 
    processedEventsTree->Branch("jlik", &jlik, "jlik/F"); 
    processedEventsTree->Branch("GNhit", &GNhit, "GNhit/F"); 
    processedEventsTree->Branch("Snpe", &Snpe, "Snpe/F"); 
    processedEventsTree->Branch("SnpeT", &SnpeT, "SnpeT/F"); 
    processedEventsTree->Branch("Slen", &Slen, "Slen/F");  // jgandalf 
    processedEventsTree->Branch("Elen", &Elen, "Elen/F");  // jenergy 
    processedEventsTree->Branch("best_trk_pos_x", &best_trk_pos[0], "best_trk_pos_x/F"); 
    processedEventsTree->Branch("best_trk_pos_y", &best_trk_pos[1], "best_trk_pos_y/F"); 
    processedEventsTree->Branch("best_trk_pos_z", &best_trk_pos[2], "best_trk_pos_z/F"); 
    processedEventsTree->Branch("best_trk_dir_x", &best_trk_dir[0], "best_trk_dir_x/F");
    processedEventsTree->Branch("best_trk_dir_y", &best_trk_dir[1], "best_trk_dir_y/F");
    processedEventsTree->Branch("best_trk_dir_z", &best_trk_dir[2], "best_trk_dir_z/F");
    processedEventsTree->Branch("Q1", &Q1, "Q1/F"); 
    processedEventsTree->Branch("Q1new", &Q1new, "Q1new/F");  // lik/ndof 
    processedEventsTree->Branch("diffangle", &diffangle, "diffangle/F"); 
    processedEventsTree->Branch("delta_zenith", &delta_zenith, "delta_zenith/F"); 
    processedEventsTree->Branch("mu_pos_x", &mu_pos[0], "mu_pos_x/F"); 
    processedEventsTree->Branch("mu_pos_y", &mu_pos[1], "mu_pos_y/F"); 
    processedEventsTree->Branch("mu_pos_z", &mu_pos[2], "mu_pos_z/F"); 
    processedEventsTree->Branch("mu_dir_x", &mu_dir[0], "mu_dir_x/F"); 
    processedEventsTree->Branch("mu_dir_y", &mu_dir[1], "mu_dir_y/F"); 
    processedEventsTree->Branch("mu_dir_z", &mu_dir[2], "mu_dir_z/F"); 
    processedEventsTree->Branch("nu_pos_x", &nu_pos[0], "nu_pos_x/F"); 
    processedEventsTree->Branch("nu_pos_y", &nu_pos[1], "nu_pos_y/F"); 
    processedEventsTree->Branch("nu_pos_z", &nu_pos[2], "nu_pos_z/F");
    processedEventsTree->Branch("nu_dir_x", &nu_dir[0], "nu_dir_x/F"); 
    processedEventsTree->Branch("nu_dir_y", &nu_dir[1], "nu_dir_y/F"); 
    processedEventsTree->Branch("nu_dir_z", &nu_dir[2], "nu_dir_z/F"); 
    processedEventsTree->Branch("cos_zen_mu", &cos_zen_mu, "cos_zen_mu/F"); 
    processedEventsTree->Branch("zenith_mu", &zenith_mu, "zenith_mu/F"); 
    processedEventsTree->Branch("logE_mu", &logE_mu, "logE_mu/F"); 
    processedEventsTree->Branch("logE_nu", &logE_nu, "logE_nu/F"); 
    processedEventsTree->Branch("logEbundle", &logEbundle, "logEbundle/F"); 
    processedEventsTree->Branch("logEdepos", &logEdepos, "logEdepos/F"); 
    processedEventsTree->Branch("bjorken_y", &bjorken_y, "bjorken_y/F"); 
    processedEventsTree->Branch("cos_zen_nu", &cos_zen_nu, "cos_zen_nu/F"); 
    processedEventsTree->Branch("zenith_nu", &zenith_nu, "zenith_nu/F"); 
    processedEventsTree->Branch("w1", &w1, "w1/F");
    processedEventsTree->Branch("w2", &w2, "w2/F");
    processedEventsTree->Branch("w3atm", &w3atm, "w3atm/F"); 
    processedEventsTree->Branch("w3iceHESE", &w3iceHESE, "w3iceHESE/F"); 
    processedEventsTree->Branch("w3loi", &w3loi, "w3loi/F"); 
    processedEventsTree->Branch("w3random", &w3random, "w3random/F"); 
    processedEventsTree->Branch("w3ice", &w3ice, "w3ice/F");
    processedEventsTree->Branch("w3antares", &w3antares, "w3antares/F");
    processedEventsTree->Branch("w3new",&w3new,"w3new/F");
    processedEventsTree->Branch("neutrino_ngen",&neutrino_ngen,"neutrino_ngen/F");
    processedEventsTree->Branch("neutrino_type", &neutrino_type, "neutrino_type/I");
    processedEventsTree->Branch("neutrino_interaction",&neutrino_interaction,"neutrino_interaction/I");
    processedEventsTree->Branch("Ntrack_all", &Ntrack_all, "Ntrack_all/F"); 
    processedEventsTree->Branch("Nallbehind", &Nallbehind, "Nallbehind/F"); 
    processedEventsTree->Branch("Nallbehind2", &Nallbehind2, "Nallbehind2/F"); 
    processedEventsTree->Branch("Nallbehind3", &Nallbehind3, "Nallbehind3/F"); 
    processedEventsTree->Branch("NtrackIT10", &NtrackIT10, "NtrackIT10/F");  // in time with the track 
    processedEventsTree->Branch("NtrackIT30", &NtrackIT30, "NtrackIT30/F"); 
    processedEventsTree->Branch("NtrackIT", &NtrackIT, "NtrackIT/F"); 
    processedEventsTree->Branch("NtrackEarly", &NtrackEarly, "NtrackEarly/F"); 
    processedEventsTree->Branch("NtrackLate", &NtrackLate, "NtrackLate/F"); 
    processedEventsTree->Branch("TrLengthIT", &TrLengthIT, "TrLengthIT/F"); 
    processedEventsTree->Branch("TrLengthIT_2", &TrLengthIT_2, "TrLengthIT_2/F"); 
    processedEventsTree->Branch("TrLengthIT_3", &TrLengthIT_3, "TrLengthIT_3/F"); 
    processedEventsTree->Branch("ratio130", &ratio130, "ratio130/F"); 
    processedEventsTree->Branch("ratio430", &ratio430, "ratio430/F"); 
    processedEventsTree->Branch("ratio330", &ratio330, "ratio330/F"); 
    processedEventsTree->Branch("ratio110", &ratio110, "ratio110/F"); 
    processedEventsTree->Branch("ratio410", &ratio410, "ratio410/F"); 
    processedEventsTree->Branch("ratio310", &ratio310, "ratio310/F");
    processedEventsTree->Branch("ratio1", &ratio1, "ratio1/F"); 
    processedEventsTree->Branch("ratio2", &ratio2, "ratio2/F"); 
    processedEventsTree->Branch("ratio3", &ratio3, "ratio3/F"); 
    processedEventsTree->Branch("ratio4", &ratio4, "ratio4/F"); 
    processedEventsTree->Branch("ratio5", &ratio5, "ratio5/F"); 
    processedEventsTree->Branch("ratio6", &ratio6, "ratio6/F"); 
    processedEventsTree->Branch("ratiol", &ratiol, "ratiol/F"); 
    processedEventsTree->Branch("ratiol_trig", &ratiol_trig, "ratiol_trig/F"); 
    processedEventsTree->Branch("ratio4_trig", &ratio4_trig, "ratio4_trig/F"); 
    processedEventsTree->Branch("myratio50_muon", &myratio50_muon, "myratio50_muon/F");  // first hit of track 
    processedEventsTree->Branch("myratio30_muon", &myratio30_muon, "myratio30_muon/F"); 
    processedEventsTree->Branch("myratio50_cascmuon", &myratio50_cascmuon, "myratio50_cascmuon/F"); 
    processedEventsTree->Branch("myratio30_cascmuon", &myratio30_cascmuon, "myratio30_cascmuon/F"); 
    processedEventsTree->Branch("myratio50_casc", &myratio50_casc, "myratio50_casc/F"); 
    processedEventsTree->Branch("myratio30_casc", &myratio30_casc, "myratio30_casc/F");
    processedEventsTree->Branch("myratio50_cascmuon_over_mu",&myratio50_cascmuon_over_mu, "myratio50_cascmuon_over_mu/F");
    processedEventsTree->Branch("myratio30_cascmuon_over_mu",&myratio30_cascmuon_over_mu, "myratio30_cascmuon_over_mu/F");
    processedEventsTree->Branch("myratio50_casc_over_mu",&myratio50_casc_over_mu, "myratio50_casc_over_mu/F");
    processedEventsTree->Branch("myratio30_casc_over_mu",&myratio30_casc_over_mu, "myratio30_casc_over_mu/F");
    processedEventsTree->Branch("redToT_cascmuon_over_mu", &redToT_cascmuon_over_mu, "redToT_cascmuon_over_mu/F");
    processedEventsTree->Branch("redToT_casc_over_mu", &redToT_casc_over_mu, "redToT_casc_over_mu/F");
    processedEventsTree->Branch("ratio_closehits_cascmuon_over_mu", &ratio_closehits_cascmuon_over_mu, "ratio_closehits_cascmuon_over_mu/F");
    processedEventsTree->Branch("ratio_closehits_casc_over_mu", &ratio_closehits_casc_over_mu, "ratio_closehits_casc_over_mu/F");
    processedEventsTree->Branch("diff_theta", &diff_theta, "diff_theta/F"); 
    processedEventsTree->Branch("ratio_closehits_muon", &ratio_closehits_muon, "ratio_closehits_muon/F"); 
    processedEventsTree->Branch("ratio_closehits_cascmuon", &ratio_closehits_cascmuon, "ratio_closehits_cascmuon/F"); 
    processedEventsTree->Branch("ratio_closehits_casc", &ratio_closehits_casc, "ratio_closehits_casc/F"); 
    processedEventsTree->Branch("redToT_muon", &redToT_muon, "redToT_muon/F"); 
    processedEventsTree->Branch("redToT_cascmuon", &redToT_cascmuon, "redToT_cascmuon/F"); 
    processedEventsTree->Branch("redToT_casc", &redToT_casc, "redToT_casc/F"); 
    processedEventsTree->Branch("diff_dist_mu", &diff_dist_mu, "diff_dist_mu/F"); 
    processedEventsTree->Branch("diff_dist_casc_mu", &diff_dist_casc_mu, "diff_dist_casc_mu/F"); 
    processedEventsTree->Branch("diff_dist_casc", &diff_dist_casc, "diff_dist_casc/F"); 
    processedEventsTree->Branch("max_lik_down", &max_lik_down, "max_lik_down/F"); 
    processedEventsTree->Branch("max_lik_up", &max_lik_up, "max_lik_up/F"); 
    processedEventsTree->Branch("costheta_min", &costheta_min, "costheta_min/F"); 
    processedEventsTree->Branch("costheta_max", &costheta_max, "costheta_max/F"); 
    processedEventsTree->Branch("num_of_good_sol", &num_of_good_sol, "num_of_good_sol/F"); 
    processedEventsTree->Branch("downsol", &downsol, "downsol/F"); 
    processedEventsTree->Branch("upsol", &upsol, "upsol/F"); 
    processedEventsTree->Branch("nhits_casc", &nhits_casc, "nhits_casc/F"); 
    processedEventsTree->Branch("nhits_casc_100", &nhits_casc_100, "nhits_casc_100/F"); 
    processedEventsTree->Branch("sum_ToT_casc", &sum_ToT_casc, "sum_ToT_casc/F"); 
    processedEventsTree->Branch("min_dist_casc", &min_dist_casc, "min_dist_casc/F"); 
    processedEventsTree->Branch("max_dist_casc", &max_dist_casc, "max_dist_casc/F"); 
    processedEventsTree->Branch("nhits_mu", &nhits_mu, "nhits_mu/F");
    processedEventsTree->Branch("nhits_mu_100", &nhits_mu_100, "nhits_mu_100/F");
    processedEventsTree->Branch("sum_ToT_mu", &sum_ToT_mu, "sum_ToT_mu/F"); 
    processedEventsTree->Branch("min_dist_mu", &min_dist_mu, "min_dist_mu/F"); 
    processedEventsTree->Branch("max_dist_mu", &max_dist_mu, "max_dist_mu/F"); 
    processedEventsTree->Branch("nhits_casc_mu", &nhits_casc_mu, "nhits_casc_mu/F");
    processedEventsTree->Branch("nhits_casc_mu_100", &nhits_casc_mu_100, "nhits_casc_mu_100/F");
    processedEventsTree->Branch("sum_ToT_casc_mu", &sum_ToT_casc_mu, "sum_ToT_casc_mu/F");
    processedEventsTree->Branch("min_dist_casc_mu", &min_dist_casc_mu, "min_dist_casc_mu/F");
    processedEventsTree->Branch("max_dist_casc_mu", &max_dist_casc_mu, "max_dist_casc_mu/F");
    processedEventsTree->Branch("min_diff_sollik", &min_diff_sollik, "min_diff_sollik/F"); 
    processedEventsTree->Branch("max_diff_sollik", &max_diff_sollik, "max_diff_sollik/F"); 
    processedEventsTree->Branch("min_diff_sol", &min_diff_sol, "min_diff_sol/F"); 
    processedEventsTree->Branch("max_diff_sol", &max_diff_sol, "max_diff_sol/F"); 
    processedEventsTree->Branch("min_zen_sol", &min_zen_sol, "min_zen_sol/F"); 
    processedEventsTree->Branch("max_zen_sol", &max_zen_sol, "max_zen_sol/F"); 
    processedEventsTree->Branch("diffangle_shower_nu", &diffangle_shower_nu, "diffangle_shower_nu/F"); 
    processedEventsTree->Branch("diffangle_shower_mu", &diffangle_shower_mu, "diffangle_shower_mu/F"); 
    processedEventsTree->Branch("delta_zenith_shower_mu", &delta_zenith_shower_mu, "delta_zenith_shower_mu/F"); 
    processedEventsTree->Branch("delta_zenith_shower_nu", &delta_zenith_shower_nu, "delta_zenith_shower_nu/F"); 
    processedEventsTree->Branch("zenith_shower", &zenith_shower, "zenith_shower/F"); 
    processedEventsTree->Branch("beta0_shower", &beta0_shower, "beta0_shower/F"); 
    processedEventsTree->Branch("beta0_shower_deg", &beta0_shower_deg, "beta0_shower_deg/F"); 
    processedEventsTree->Branch("lik_shower", &lik_shower, "lik_shower/F"); 
    processedEventsTree->Branch("Nhit_shower", &Nhit_shower, "Nhit_shower/F"); 
    processedEventsTree->Branch("best_trk_pos_shower_x", &best_trk_pos_shower[0], "best_trk_pos_shower_x/F"); 
    processedEventsTree->Branch("best_trk_pos_shower_y", &best_trk_pos_shower[1], "best_trk_pos_shower_y/F"); 
    processedEventsTree->Branch("best_trk_pos_shower_z", &best_trk_pos_shower[2], "best_trk_pos_shower_z/F"); 
    processedEventsTree->Branch("best_trk_dir_shower_x", &best_trk_dir_shower[0], "best_trk_dir_shower_x/F");
    processedEventsTree->Branch("best_trk_dir_shower_y", &best_trk_dir_shower[1], "best_trk_dir_shower_y/F");
    processedEventsTree->Branch("best_trk_dir_shower_z", &best_trk_dir_shower[2], "best_trk_dir_shower_z/F");
    processedEventsTree->Branch("Ereco_shower", &Ereco_shower, "Ereco_shower/F"); 
    processedEventsTree->Branch("Ereco_shower_corrected", &Ereco_shower_corrected, "Ereco_shower_corrected/F"); 
    processedEventsTree->Branch("logEreco_shower", &logEreco_shower, "logEreco_shower/F"); 
    processedEventsTree->Branch("dlik", &dlik, "dlik/F"); 
    processedEventsTree->Branch("normdlik", &normdlik, "normdlik/F"); 
    processedEventsTree->Branch("itoverlen", &itoverlen, "itoverlen/F"); 
    processedEventsTree->Branch("delta_zenith_track_shower", &delta_zenith_track_shower, "delta_zenith_track_shower/F"); 
    processedEventsTree->Branch("diffangle_track_shower", &diffangle_track_shower, "diffangle_track_shower/F"); 
    processedEventsTree->Branch("flag_muon_3D", &flag_muon_3D, "flag_muon_3D/O"); 
    processedEventsTree->Branch("flag_shower_3D", &flag_shower_3D, "flag_shower_3D/O"); 
    processedEventsTree->Branch("flag_shower_MX", &flag_shower_MX, "flag_shower_MX/O"); 
    processedEventsTree->Branch("ToT_border_mu", &ToT_border_mu, "ToT_border_mu/F"); 
    processedEventsTree->Branch("ToT_border_casc", &ToT_border_casc, "ToT_border_casc/F"); 
    processedEventsTree->Branch("ToT_border_cascmu", &ToT_border_cascmu, "ToT_border_cascmu/F"); 
    processedEventsTree->Branch("ToT_trig", &ToT_trig, "ToT_trig/F"); 
    processedEventsTree->Branch("max_ToT_trig", &max_ToT_trig, "max_ToT_trig/F"); 
    processedEventsTree->Branch("ToT_IT", &ToT_IT, "ToT_IT/F"); 
    processedEventsTree->Branch("ToT_allIT", &ToT_allIT, "ToT_allIT/F"); 
    processedEventsTree->Branch("Nborder_hits", &Nborder_hits, "Nborder_hits/F"); 
    processedEventsTree->Branch("Nborder_cherenkov_hits", &Nborder_cherenkov_hits, "Nborder_cherenkov_hits/F"); 
    processedEventsTree->Branch("Nborder_dtres_hits", &Nborder_dtres_hits, "Nborder_dtres_hits/F"); 
    processedEventsTree->Branch("Nborder_DOMs", &Nborder_DOMs, "Nborder_DOMs/F"); 
    processedEventsTree->Branch("Nhits_upper", &Nhits_upper, "Nhits_upper/F"); 
    processedEventsTree->Branch("Nhits_lower", &Nhits_lower, "Nhits_lower/F"); 
    processedEventsTree->Branch("Nhits_border_upper", &Nhits_border_upper, "Nhits_border_upper/F"); 
    processedEventsTree->Branch("Nhits_border_lower", &Nhits_border_lower, "Nhits_border_lower/F"); 
    processedEventsTree->Branch("Nhits_cherenkov_upper", &Nhits_cherenkov_upper, "Nhits_cherenkov_upper/F"); 
    processedEventsTree->Branch("Nhits_cherenkov_lower", &Nhits_cherenkov_lower, "Nhits_cherenkov_lower/F"); 
    processedEventsTree->Branch("Nhits_border_cherenkov_upper", &Nhits_border_cherenkov_upper, "Nhits_border_cherenkov_upper/F"); 
    processedEventsTree->Branch("Nhits_border_cherenkov_lower", &Nhits_border_cherenkov_lower, "Nhits_border_cherenkov_lower/F"); 
    processedEventsTree->Branch("NtrackIT50", &NtrackIT50, "NtrackIT50/F");
    processedEventsTree->Branch("NtrackIT50_2", &NtrackIT50_2, "NtrackIT50_2/F");
    processedEventsTree->Branch("NtrackIT30_2", &NtrackIT30_2, "NtrackIT30_2/F");
    processedEventsTree->Branch("NtrackIT10_2", &NtrackIT10_2, "NtrackIT10_2/F");
    processedEventsTree->Branch("NtrackIT50_3", &NtrackIT50_3, "NtrackIT50_3/F");
    processedEventsTree->Branch("NtrackIT30_3", &NtrackIT30_3, "NtrackIT30_3/F");
    processedEventsTree->Branch("NtrackIT10_3", &NtrackIT10_3, "NtrackIT10_3/F");
    processedEventsTree->Branch("num_triggered_hits", &num_triggered_hits, "num_triggered_hits/I"); 
    processedEventsTree->Branch("num_triggered_lines", &num_triggered_lines, "num_triggered_lines/I");  
    processedEventsTree->Branch("num_triggered_doms", &num_triggered_doms, "num_triggered_doms/I");  
    processedEventsTree->Branch("num_triggered_pmts", &num_triggered_pmts, "num_triggered_pmts/I");  
    processedEventsTree->Branch("num_cherenkov_hits", &num_cherenkov_hits, "num_cherenkov_hits/I"); 
    processedEventsTree->Branch("num_cherenkov_lines", &num_cherenkov_lines, "num_cherenkov_lines/I");  
    processedEventsTree->Branch("num_cherenkov_doms", &num_cherenkov_doms, "num_cherenkov_doms/I");  
    processedEventsTree->Branch("num_cherenkov_pmts", &num_cherenkov_pmts, "num_cherenkov_pmts/I");  
    processedEventsTree->Branch("num_cascade_hits", &num_cascade_hits, "num_cascade_hits/I"); 
    processedEventsTree->Branch("num_cascade_lines", &num_cascade_lines, "num_cascade_lines/I");  
    processedEventsTree->Branch("num_cascade_doms", &num_cascade_doms, "num_cascade_doms/I");  
    processedEventsTree->Branch("num_cascade_pmts", &num_cascade_pmts, "num_cascade_pmts/I");  
    processedEventsTree->Branch("ratio_cherenkov_hits", &ratio_cherenkov_hits, "ratio_cherenkov_hits/F"); 
    processedEventsTree->Branch("ratio_cherenkov_lines", &ratio_cherenkov_lines, "ratio_cherenkov_lines/F");  
    processedEventsTree->Branch("ratio_cherenkov_doms", &ratio_cherenkov_doms, "ratio_cherenkov_doms/F");  
    processedEventsTree->Branch("ratio_cherenkov_pmts", &ratio_cherenkov_pmts, "ratio_cherenkov_pmts/F");  
    processedEventsTree->Branch("ratio_cascade_hits", &ratio_cascade_hits, "ratio_cascade_hits/F"); 
    processedEventsTree->Branch("ratio_cascade_lines", &ratio_cascade_lines, "ratio_cascade_lines/F"); 
    processedEventsTree->Branch("ratio_cascade_doms", &ratio_cascade_doms, "ratio_cascade_doms/F");  
    processedEventsTree->Branch("ratio_cascade_pmts", &ratio_cascade_pmts, "ratio_cascade_pmts/F");
    processedEventsTree->Branch("associated_n_hits", &associated_n_hits, "associated_n_hits/I");
    processedEventsTree->Branch("Log_No_PMTs_with_dist_weights_MAXhits", &Log_No_PMTs_with_dist_weights_MAXhits, "Log_No_PMTs_with_dist_weights_MAXhits/F");
    processedEventsTree->Branch("No_OMs_length_potential_length_chosen", &No_OMs_length_potential_length_chosen, "No_OMs_length_potential_length_chosen/F");
    processedEventsTree->Branch("PMTs_hit_to_nohit", &PMTs_hit_to_nohit, "PMTs_hit_to_nohit/F");
    processedEventsTree->Branch("MaxLenPos_OvrMaxDist", &MaxLenPos_OvrMaxDist, "MaxLenPos_OvrMaxDist/F");
    processedEventsTree->Branch("TrLenIT_3_OvrMaxDist", &TrLenIT_3_OvrMaxDist, "TrLenIT_3_OvrMaxDist/F");
     


    //=============================================================================================================================

    /*
        ====================================
        |                                  |
        |           Event Loop             |
        |                                  |
        ====================================
    */

    // Event Loop
    for (int event_num = 0; event_num < number_of_events; ++event_num) {

		
    	if (event_num >= 30) print_mode.setVerbosity(1);
        

        // Testing Guard Clause: Check first 10 events
        //if (true && event_num > 20) break;

        //
        inputTree->GetEntry(event_num);
        const Evt event = *event_pointer;

	/*These can't be run in case of data
	  nu_bool = event.mc_trks[0].is_neutrino();
	  mupage_bool = event.mc_trks[0].is_muon();*/

        // Print event information
        print_mode() << "\n\n\033[1mStarting Event " << event_num << "\033[22m : ";
        //event.print();

        /*
            ~ Explanation of bitwise operators in C++ ~

            Bitwise AND: &
                compares variables bit to bit, for example:
                1001 & 0110 ==> 0000 false
                0100 & 1111 ==> 0100 true
                0011 & 0111 ==> 0011 true

            Bitwise Left Shift: <<
                moves bits to the left, destroying far left bits and adding 0 bits at the far right.
                0110 << 1 ==> 1100 (moved left by 1).
                10010110 << 2 ==> 01011000 (moved left by 2, destroyed the 2 far left bits that had value 10, added 00 at far right).

            In binary, integers are represented by 4 bytes or 32 bits

            1       == 00000000 00000000 00000000 00000001
            1 << 1  == 00000000 00000000 00000000 00000010
            1 << 2  == 00000000 00000000 00000000 00000100
            1 << 4  == 00000000 00000000 00000000 00010000
            1 << 31 == 10000000 00000000 00000000 00000000

            ~ G.Z.
        */

        flag_shower_3D = event.trigger_mask & (1 << 1);
        flag_shower_MX = event.trigger_mask & (1 << 2);
        flag_muon_3D = event.trigger_mask & (1 << 4);
        flag_factory_limit = event.trigger_mask & (1 << 31);

        // If true adds 1, if false adds 0.
        events_showers += flag_shower_3D;
        events_showers_MX += flag_shower_MX;
        events_muons += flag_muon_3D;
        events_factory_limit += flag_factory_limit;

        evt_num = event.frame_index;
        evt_id = event.id;
        evt_overlays = event.overlays;

        // Transfer out of Loop
        // print_mode() << " flag_shower_3D : " << flag_shower_3D << " flag_shower_MX : " << flag_shower_MX << " flag_muon_3D : " << flag_muon_3D << " flag_factory_limit : " << flag_factory_limit << std::endl;

        pseudo_livetime = LiveT;
        pseudo_runid = run_number;
	pseudo_subRunid = run_subNumber;

        // Per-Event Initialization to 0
        max_ToT_trig = logE_nu = E_nu = cos_zen_nu = zenith_nu = 0;

        mu_vertex_in = neut_vertex_in = reco_vertex_in = false;

        // Neutrino Weights;
        w3 = w3atm = w3iceHESE = w3loi = w3antares = w3random = w3ice = w3new = w1 = w2 = 0;
	neutrino_ngen = 0;
	neutrino_type = 0;
	neutrino_interaction=1;

        cos_zen = zenith = bjorken_y = logEreco = logEreco2 = multiplicity_neutrino = multiplicity_muon = evt_overlays = NtrackIT = NtrackEarly = Nallbehind = NtrackLate = Nallbehind2 = Nallbehind3 = NtrackIT30 = NtrackIT10 = TrLengthIT = lambdaIT_max = strict_lambda_min = lambdaCH_max = maxML1lamda = TrLengthIT_2 = TrLengthIT_3 = Ntrack_all = Ntrack_all2 = jlik = jbeta0 = jbeta0_deg = logbeta0 = GChi = GNhit = EChi = GLik = GNit = Snpe = SnpeT = Slen = ENDF = ENhit = Elen = fi = fi_rad = Ereco_shower = Ereco_shower_corrected = Nhit_shower = Ereco = Ereco2 = Q1 = Q1new = ratio_closehits_muon = ratio_closehits_cascmuon = ratio_closehits_casc = redToT_muon = redToT_cascmuon = redToT_casc = diff_dist_mu = diff_dist_casc_mu = diff_dist_casc = min_diff_sol = max_diff_sol = min_diff_sollik = max_diff_sollik = min_zen_sol = max_zen_sol = nhits_casc_mu = nhits_casc = nhits_mu = nhits_casc_mu_100 = nhits_casc_100 = nhits_mu_100 = max_dist_casc_mu = max_dist_mu = max_dist_casc = sum_ToT_casc_mu = sum_ToT_mu = sum_ToT_casc = Q = costheta_max = theta_max = costheta_min = theta_min = diff_theta = num_of_good_sol = downsol = upsol = max_lik_down = max_lik_up = lik_sol_down = lik_sol_up = ratio110 = ratio410 = ratio130 = ratio430 = ratio1 = ratio2 = ratio4 = ratio5 = ratio6 = ratiol = ratiol_trig = ratio4_trig = ratio3 = ratio310 = ratio330 = mean_tres_it = R = R_true_muon = R_true_nu = nhits_casc_mu_50 = nhits_casc_50 = nhits_mu_50 = nhits_casc_mu_30 = nhits_casc_30 = nhits_mu_30 = myratio50_muon = myratio50_casc = myratio50_cascmuon = myratio30_muon = myratio30_casc = myratio30_cascmuon = D = D3d = itoverlen = 0;
	Log_No_PMTs_with_dist_weights_MAXhits = No_OMs_length_potential_length_chosen = PMTs_hit_to_nohit = MaxLenPos_OvrMaxDist = TrLenIT_3_OvrMaxDist = 0.0;

        mu_pos[0] = mu_pos[1] = mu_pos[2] = mu_dir[0] = mu_dir[1] = mu_dir[2] = nu_pos[0] = nu_pos[1] = nu_pos[2] = nu_dir[0] = nu_dir[1] = nu_dir[2] = 0;

        /* BDT unused
        downLikMax = upLikMax = cosZenith_min = cosZenith_max = prefitN_1degree = prefitN_down = prefitN_up = cascNHits = cascNHits_100m = cascDistanceMin = cascDistanceMax = nHits = nHits_100m = distanceMin = distanceMax = cascmuNHits = cascmuNHits_100m = cascmuDistanceMin = cascmuDistanceMax = thetaDiff = 0;
        */

        Nborder_cherenkov_hits = Nborder_hits = Nborder_cherenkov_2lines_hits = Nborder_dtres_hits = Nborder_DOMs = Nhits_upper = Nhits_lower = Nhits_border_upper = Nhits_border_lower = Nhits_cherenkov_upper = Nhits_cherenkov_lower = Nhits_border_cherenkov_upper = Nhits_border_cherenkov_lower = 0.0;

        best_trk_pos[0] = best_trk_pos[1] = best_trk_pos[2] = best_trk_dir[0] = best_trk_dir[1] = best_trk_dir[2] = best_trk_pos_shower[0] = best_trk_pos_shower[1] = best_trk_pos_shower[2] = best_trk_dir_shower[0] = best_trk_dir_shower[1] = best_trk_dir_shower[2] = 0;

        int_Nborder_cherenkov_hits = int_Nborder_hits = int_Nborder_dtres_hits = int_Nborder_DOMs = int_Nhits_upper = int_Nhits_lower = int_Nhits_border_upper = int_Nhits_border_lower = int_Nhits_cherenkov_upper = int_Nhits_cherenkov_lower = int_Nhits_border_cherenkov_upper = int_Nhits_border_cherenkov_lower = NtrackIT50 = NtrackIT50_2 = NtrackIT30_2 = NtrackIT10_2 = NtrackIT50_3 = NtrackIT30_3 = NtrackIT10_3 = 0;

        ToT_border_mu = 0;
        ToT_border_casc = 0;
        ToT_border_cascmu = 0;
        ToT_trig = 0;
        ToT_IT = 0;
        ToT_allIT = 0;

        associated_n_hits = 0;
        // Initialization End //


        // Reset Conditions



        // Per-Event Initialization to numbers
        min_dist_casc_mu = min_dist_casc = min_dist_mu = 99999.9;
        minML1lamda = 9999;




        // Per-Event Initialization Complete

        /*
            ====================================
            |                                  |
            |       True Tracks Analysis       |
            |                                  |
            ====================================
        */
        print_mode() << "\033[1;35m[True Tracks Analysis]\033[22;39m\n";

        // Block: MC Track Loop, get MC muon / neutrino tracks // Platinum
        Trk mu;
        Trk nu;
        E_mu_max = 0;
        E_mu_depos = 0;
        E_bundle = 0;
	E_visible = 0;
        cos_zen_mu = 0;
        zenith_mu = 0;
        logE_mu = 0;
        logEbundle = 0; logEdepos = 0;

        for (unsigned i_track = 0; i_track < event.mc_trks.size(); ++i_track) {

            // Check if track belongs to neutrino
            if (event.mc_trks[i_track].is_neutrino()) {
                ++multiplicity_neutrino;
                nu = event.mc_trks[0]; //wrong in case of mupage file, but we don't care in that case
            }

            // Check if track belongs to muon
            if (event.mc_trks[i_track].is_muon()) {
                ++multiplicity_muon;
                E_bundle = E_bundle + event.mc_trks[i_track].E;

                if (event.mc_trks[i_track].haveusr("energy_lost_in_can")) {
                    E_mu_depos = E_mu_depos + event.mc_trks[i_track].getusr(mc_usr_keys::energy_lost_in_can);
                }

                // Final muon track is of the one with the highest energy
                if (event.mc_trks[i_track].E >= E_mu_max) {
                    mu = event.mc_trks[i_track];
                    E_mu_max = event.mc_trks[i_track].E;  // MC muon energy is in GeV
                }
            }
        }

        // Check if tracks have been filled
        bool mu_mc_track_found = mu.dir.len() != 0 && mu.pos.len() != 0;
        bool nu_mc_track_found = nu.dir.len() != 0 && nu.pos.len() != 0;

        // Extract muon MC variables
        if (mu_mc_track_found) {
            cos_zen_mu = -mu.dir.z / mu.dir.len();
            zenith_mu = acos(cos_zen_mu) * TMath::RadToDeg();
            logE_mu = TMath::Log10(mu.E);
            logEbundle = TMath::Log10(E_bundle);
            logEdepos = TMath::Log10(E_mu_depos);

            mu_pos[0] = mu.pos.x;
            mu_pos[1] = mu.pos.y;
            mu_pos[2] = mu.pos.z;

            mu_dir[0] = mu.dir.x;
            mu_dir[1] = mu.dir.y;
            mu_dir[2] = mu.dir.z;

            R_true_muon = (mu.pos - mean_detector_pos).lenxy();
            mu_vertex_in = R_true_muon < 200 && mu.pos.z < 650 && mu.pos.z > 120; // Maybe outdated? ~ G.Z.
        }

        // Extract neutrino MC variables
        if (nu_mc_track_found) {
            E_nu = nu.E;
            cos_zen_nu = -nu.dir.z / nu.dir.len();
            zenith_nu = std::acos(cos_zen_nu) * TMath::RadToDeg();
            logE_nu = TMath::Log10(nu.E);

            nu_pos[0] = nu.pos.x;
            nu_pos[1] = nu.pos.y;
            nu_pos[2] = nu.pos.z;

            nu_dir[0] = nu.dir.x;
            nu_dir[1] = nu.dir.y;
            nu_dir[2] = nu.dir.z;

	    R_true_nu = (nu.pos - mean_detector_pos).lenxy();
	    neut_vertex_in = R_true_nu < 200 && nu.pos.z < 650 && nu.pos.z > 120; // Maybe outdated?
        }

        // Extract bjorken_y
        if (mu_mc_track_found && nu_mc_track_found) bjorken_y = (1.0 - (mu.E / nu.E));

	//=========================================================================
	//VISIBLE ENERGY (6/3/24 -> update 20/6/24)
	for (unsigned i_track = 0; i_track < event.mc_trks.size(); ++i_track){

	  Trk temp = event.mc_trks[i_track];
	  //instrumented volume
	  bool inVolume = (temp.pos - mean_detector_pos).lenxy() < 200 && temp.pos.z < 690 && temp.pos.z > 64;
	  bool tauFromNuTau = TMath::Abs(event.mc_trks[i_track].type) == 15 && TMath::Abs(event.mc_trks[0].type) == 16;

	  if(temp.status == 1){//compute ONLY for final particles
	    
	    if( inVolume ){

	      if (event.mc_trks[i_track].haveusr("energy_lost_in_can")) {
		E_visible += event.mc_trks[i_track].getusr(mc_usr_keys::energy_lost_in_can);
	      }//if muon keep lost in can
	      else{
		if( !event.mc_trks[i_track].is_neutrino() && !tauFromNuTau ){//if tau from nutau decay get daughters
		  E_visible += event.mc_trks[i_track].E;
		}//if not nu (nu don't leave Energy) keep whole E
	      }//end else

	    }//end if mc_track in

	    else{

	      if (event.mc_trks[i_track].haveusr("energy_lost_in_can")) {
		E_visible += event.mc_trks[i_track].getusr(mc_usr_keys::energy_lost_in_can);
	      }//if muon keep lost in can
	    
	    }//end if mc_track out

	  }//end if final state

	}//end for
	//===============================================================================

        // Print muon track results
        print_mode() << "\033[1m[True Muon Track]\033[22m : ";
        if (mu.dir.len() != 0 && mu.pos.len() != 0) {
            print_mode() << "\033[32mFound\033[39m\n";
            print_mode() << " Pos: (" << mu.pos.x << ", " << mu.pos.y << ", " << mu.pos.z << ")\n";
            print_mode() << "  Vertex is " << (mu_vertex_in ? "\033[32minside" : "\033[31moutside") << "\033[39m detector volume.\n";
            print_mode() << " Dir: (" << mu.dir.x << ", " << mu.dir.y << ", " << mu.dir.z << ")\n";
            print_mode() << " Energy: " << mu.E << " GeV\n";
            print_mode() << " Bundle Energy (" << multiplicity_muon << " muons): " << E_bundle << " GeV\n";
            print_mode() << " Deposited Energy: " << E_mu_depos << " GeV\n";
        }
        else print_mode() << "\033[31mNot Found\033[39m\n";

        // Print neutrino track results
        print_mode() << "\033[1m[True Neutrino Track]\033[22m : ";
        if (nu.dir.len() != 0 && nu.pos.len() != 0) {
            print_mode() << "\033[32mFound\033[39m\n";
            print_mode() << " Pos: (" << nu.pos.x << ", " << nu.pos.y << ", " << nu.pos.z << ")\n";
	    print_mode() << "  Vertex is " << (neut_vertex_in ? "\033[32minside" : "\033[31moutside") << "\033[39m detector volume.\n";
            print_mode() << " Dir: (" << nu.dir.x << ", " << nu.dir.y << ", " << nu.dir.z << ")\n";
            print_mode() << " Energy: " << nu.E << "\n";
	    print_mode() << " E visible: " << E_visible << "\n";
            if (multiplicity_neutrino > 1) print_mode() << " Neutrino Multiplicity: " << multiplicity_neutrino << "\n";
        }
        else print_mode() << "\033[31mNot Found\033[39m\n";
	// Block End //


        // Block: Neutrino Weights // Gold
        if (event.mc_trks.size() > 0 && nu_bool) {

            // Fatal Exit
            if (n_files == -1) {
                print_mode() << "Fatal Error: n_files variable is -1. We need to know how many neutrino files are used.\n";
                print_mode() << "Add this in the arguments: -n <n_files>";
                print_mode().endl();
                exit(1);
            }

            if (multiplicity_neutrino < 1) {
                print_mode() << "Z-Warning: nu_bool is true but no neutrino tracks";
                print_mode().endl();
            }

	    
            neutrino_type = nu.type; //neutrino type to keep for branches
	    neutrino_interaction = event.w2list[10]; // CC -> 2 , NC -> 3 //v9

	    /*CASES::
	      nue CC: 12 2
	      nue NC: 12 3
	      
	      anue CC: -12 2
	      anue CC GLRES: -12 2 ?
	      anue NC: -12 3

	      numu CC: 14 2
	      numu NC: 14 3

	      anumu CC: -14 2
	      anumu NC: -14 3

	      nutau CC shower: 16 2 ?
	      nutau CC muon: 16 2 ?
	      nutau NC: 16 3

	      anutau CC shower: -16 2 ?
	      anutau CC muon: -16 2 ?
	      anutau NC: -16 3
	     */
	    
	    //neutrino_ngen = 1./event.w[3]; //old: events_per_file //v9
	    neutrino_ngen = Ngen; //revamped v8
            float number_of_gen_events = neutrino_ngen * n_files;
            float number_of_years = 1.0;
            float number_ofARCA_blocks = 2.0;
            float rate_atm = number_of_years * number_ofARCA_blocks * w3;  // Unused
            float coszen_nu_for_w3 = -nu.dir.z / sqrt(nu.dir.x * nu.dir.x + nu.dir.y * nu.dir.y + nu.dir.z * nu.dir.z);
            float zen_nu_rad_for_w3 = acos(coszen_nu_for_w3);

            // Atm.nu flux checked with Rosa. Should have the zen in rad!!!
	    w1 = event.w[1];
            w2 = event.w[1] / (number_of_gen_events);
            float weight_from_w2 = w2 * (getFluxConventional(nu.type, nu.E, zen_nu_rad_for_w3) + getFluxPrompt(nu.type, nu.E, zen_nu_rad_for_w3));
            w3atm = weight_from_w2 * correctKnee(nu.E);

            // Astrophysical fluxes. One used from lilia, LoI's, up to date IC!!
            // float rate_LoI_flux = number_of_years*number_ofARCA_blocks*w2*LoI_flux(nu.E);

            w3iceHESE = w2 * icecube_HESE_flux(nu.E);
            w3loi = w2 * LoI_flux(nu.E);
            w3ice = w2 * icecube_flux_diff2019(nu.E);
            w3random = w2 * random_flux(nu.E);
            w3antares = w2 * Antares_flux(nu.E);
	    w3new = w2 * new_flux(nu.E);

        }
        // Block End //

        /*
            ====================================
            |                                  |
            |       Reco Tracks Analysis       |
            |                                  |
            ====================================
        */
        print_mode() << "\033[1;35m[Reconstructed Tracks Analysis]\033[22;39m\n";



        // Block: Get Best Tracks //
        Trk best_trk;
        Trk best_trk_shower;

        // Conditions
        found_reco_track = event.trks.size() > 0 ? true : false;
        is_muon_event = has_reconstructed_jppmuon(event);
        is_shower_event = has_reconstructed_aashower(event);
        LooseShower = false;

        // Assign best track from JPP function
        if (found_reco_track && is_muon_event) best_trk = get_best_reconstructed_track<JPP_RECONSTRUCTION_TYPE>(event);
        if (found_reco_track && is_shower_event) best_trk_shower = get_best_reconstructed_aashower(event);

        // Check conditions
        //jgandalf_reco_evt = has_jppmuon_gandalf(best_trk); // Unused
        all_in = found_reco_track && is_muon_event && best_trk.rec_stages.size() > 4;

        // Print conditions for reconstructed track
        print_mode() << "\033[1m[Found Reconstructed Track]\033[22m : " << (found_reco_track ? "\033[32mTrue\033[39m\n" : "\033[31mFalse\033[39m\n");
        print_mode() << "\033[1m[Has Reconstructed JPPMuon]\033[22m : " << (is_muon_event ? "\033[32mTrue\033[39m\n" : "\033[31mFalse\033[39m\n");
        print_mode() << "\033[1m[Has Reconstructed AAShower]\033[22m : " << (is_shower_event ? "\033[32mTrue\033[39m\n" : "\033[31mFalse\033[39m\n");

        // Print how many muon reconstructed track stages passed
        print_mode() << "\033[1m[Best Muon Track Reco Stages]\033[22m : ";
        print_mode() << (best_trk.rec_stages.size() > 4 ? "\033[32m" : "\033[31m");
        print_mode() << best_trk.rec_stages.size() << "\033[39m\n";

        // Print how many shower reconstructed track stages passed
        print_mode() << "\033[1m[Best Shower Track Reco Stages]\033[22m : " << best_trk_shower.rec_stages.size() << "\n";

        if (found_reco_track) number_of_reco_tracks++;


        // Guard Clause: Skip event if not all_in
        if (!all_in) {
            print_mode() << "\033[33mSkipping this event.\033[39m\n";
            continue;
        }

        // Extract best track information
        jbeta0 = best_trk.fitinf[0];  // error on the estimated direction [rad] JGand
        GChi = best_trk.fitinf[2];  // chi2 from JGandalf.cc      reco_fitinf_lik = GChi
        GNhit = best_trk.fitinf[3];  // number of hits from JGandalf.cc     Nhits = GNhit
        Ereco = best_trk.fitinf[4];  // uncorrected energy
        EChi = best_trk.fitinf[5];  // chi2 from JEnergy.cc
        GLik = best_trk.fitinf[6];  // control parameter from JGandalf.cc
        GNit = best_trk.fitinf[7];  // number of iteration from JGandalf.cc
        Snpe = best_trk.fitinf[8];  // number of photo-electrons up to the barycentre from JStart.cc
        SnpeT = best_trk.fitinf[9];  // number of photo-electrons along the whole track from JStart.cc
        Slen = best_trk.fitinf[10];  // distance between first and last hits in metres from JStart.cc
        Elen = best_trk.fitinf[13];  // range of a muon with the reconstructed energy [m] from JEnergy.cc
        ENDF = best_trk.fitinf[15];  // number of degrees of freedom from JEnergy.cc
        ENhit = best_trk.fitinf[16];  // number of hits fron JEnergy

        best_trk_pos[0] = best_trk.pos.x;
        best_trk_pos[1] = best_trk.pos.y;
        best_trk_pos[2] = best_trk.pos.z;

        best_trk_dir[0] = best_trk.dir.x;
        best_trk_dir[1] = best_trk.dir.y;
        best_trk_dir[2] = best_trk.dir.z;

        jlik = best_trk.lik;
        Ereco2 = best_trk.E;
        fi_rad = best_trk.dir.phi();
        fi = fi_rad * TMath::RadToDeg();

        cos_zen = -best_trk.dir.z / best_trk.dir.len();
        zenith = acos(-best_trk.dir.z) * TMath::RadToDeg();
        Q = -jlik / ENDF;

        if (jbeta0 > 0) {
            jbeta0_deg = jbeta0 * 180.0 / TMath::Pi();
            logbeta0 = log10(jbeta0);
        }
        else {
            jbeta0_deg = 1000.0; // Impossible value for debugging
        }

        if (Ereco > 0 && Ereco2 > 0) {
            logEreco = log10(Ereco);
            logEreco2 = log10(Ereco2);
        }

        if (GNhit > 0) {
            Q1 = -jlik / GNhit;
            Q1new = -GChi / GNhit;
        }

        // Check if best_trk.pos (vertex) is inside the detector volume
        R = (best_trk.pos - mean_detector_pos).lenxy();
        reco_vertex_in = R < range_detector_max  && best_trk.pos.z < 650 && best_trk.pos.z > 120;

        // Extract best track information for shower
        zenith_shower = lik_shower = beta0_shower = beta0_shower_deg = -1000;
        Nhit_shower = Ereco_shower = logEreco_shower = Ereco_shower_corrected = R_shower = -1000;
        if (is_shower_event) {

            best_trk_pos_shower[0] = best_trk_shower.pos.x;
            best_trk_pos_shower[1] = best_trk_shower.pos.y;
            best_trk_pos_shower[2] = best_trk_shower.pos.z;

            best_trk_dir_shower[0] = best_trk_shower.dir.x;
            best_trk_dir_shower[1] = best_trk_shower.dir.y;
            best_trk_dir_shower[2] = best_trk_shower.dir.z;

            zenith_shower = acos(-best_trk_shower.dir.z) * TMath::RadToDeg();
            lik_shower = best_trk_shower.lik;

            // from Javier at http://git.km3net.de/jbarriosmarti/ConvertToAanet/blob/JStart/inc/AAshowerfit_file.cc
            beta0_shower = best_trk_shower.error_matrix[24] + best_trk_shower.dir.z * best_trk_shower.dir.z * best_trk_shower.error_matrix[32];
            beta0_shower_deg = (beta0_shower > 0 ? beta0_shower * TMath::RadToDeg() : 1000); // 1000 -> impossible value for debugging

            Nhit_shower = best_trk_shower.fitinf[1];
            Ereco_shower = best_trk_shower.E;
            logEreco_shower = TMath::Log10(Ereco_shower);
            Ereco_shower_corrected = Ereco_shower * 1.05;  // should it be used for v5.1 ?

            // Check if shower is loose
            R_shower = best_trk_shower.pos.lenxy();
            LooseShower = R_shower < 250000.0 && best_trk_shower.pos.z < 600.0 && lik_shower > 2200.0;
        }

        // Print muon results
        print_mode() << "\033[1m[Loose Shower]\033[22m : " << (LooseShower ? "True\n" : "False\n");
        print_mode() << "\033[1m[Best Muon Track]\033[22m\n";
        print_mode() << " Pos: (" << best_trk.pos.x << ", " << best_trk.pos.y << ", " << best_trk.pos.z << ")\n";
        print_mode() << "  Vertex is " << (reco_vertex_in ? "\033[32minside" : "\033[31moutside") << "\033[39m detector volume.\n";
        print_mode() << " Dir: (" << best_trk.dir.x << ", " << best_trk.dir.y << ", " << best_trk.dir.z << ")\n";
        print_mode() << " Ereco: " << Ereco << "\n logEreco: " << logEreco << "\n";
        print_mode() << " Ereco2: " << Ereco2 << "\n logEreco2: " << logEreco2 << "\n";
        print_mode() << " zenith: " << zenith << "\n cos(zenith) : " << cos_zen << "\n";
        print_mode() << " beta0: " << jbeta0_deg << "\n";

        // Print shower results
        print_mode() << "\033[1m[Best Shower Track]\033[22m\n";
        print_mode() << " Pos: (" << best_trk_shower.pos.x << ", " << best_trk_shower.pos.y << ", " << best_trk_shower.pos.z << ")\n";
        print_mode() << " Dir: (" << best_trk_shower.dir.x << ", " << best_trk_shower.dir.y << ", " << best_trk_shower.dir.z << ")\n";
        print_mode() << " zenith: " << zenith_shower << "\n";
        print_mode() << " beta0: " << beta0_shower_deg << "\n";
        // Block End //


        // Block: Find Solutions // Gold
        min_diff_sollik = min_diff_sol = min_zen_sol = 1000;
        max_diff_sollik = max_diff_sol = max_zen_sol = -1000;
        costheta_max = -1000;
        costheta_min = 1000;
        theta_max = -1000;
        theta_min = 1000;
        downsol = upsol = 0;
        lik_sol_down = lik_sol_up = 0;
        max_lik_down = max_lik_up = 0;
        num_of_good_sol = 0;

        for (unsigned int i_tracks = 0; i_tracks < event.trks.size(); i_tracks++) {
            const Trk* current_track = &event.trks[i_tracks];

            // [ ! ] So it's fine on any other case?
            if (current_track->rec_type == 4000 && current_track->rec_stages.size() < 3) continue;
            if (current_track->rec_type == 101) continue;

            if ((event.trks.size() > 1) && (current_track->lik != jlik) && (current_track->lik > 0) && (jlik > 0)) {

                float diffangle_sol = acos(current_track->dir.dot(best_trk.dir)) * TMath::RadToDeg();
                float zenith_sol = acos(-current_track->dir.z) * TMath::RadToDeg();
                //float cos_zen_sol = -current_track->dir.z / current_track->dir.len();  // Unused
                float lik_sol = current_track->lik;
                //float lik_per = fabs((lik_sol - jlik) / jlik);  // Unused

                // Number of good solutions (very small diffangle)
                if (diffangle_sol < 1.0) num_of_good_sol++;


                // Min - Max comparisons
                if (costheta_max < -current_track->dir.z) {
                    costheta_max = -current_track->dir.z;
                    theta_max = acos(-current_track->dir.z) * TMath::RadToDeg();
                }
                if (costheta_min > -current_track->dir.z) {
                    costheta_min = -current_track->dir.z;
                    theta_min = acos(-current_track->dir.z) * TMath::RadToDeg();
                }


                // Downgoing -> direction is negative
                if (current_track->dir.z < 0) {
                    downsol++;
                    if (max_lik_down < current_track->lik) {
                        max_lik_down = current_track->lik;
                    }
                }
                else {
                    upsol++;
                    if (max_lik_up < current_track->lik) {
                        max_lik_up = current_track->lik;
                    }
                }

                if (diffangle_sol <= min_diff_sollik) {
                    if (fabs((lik_sol - jlik) / jlik) <= 0.1) {
                        min_diff_sollik = diffangle_sol;
                    }
                }

                if (diffangle_sol > max_diff_sollik) {
                    if (fabs((lik_sol - jlik) / jlik) <= 0.1) {
                        max_diff_sollik = diffangle_sol;
                    }
                }

                if (diffangle_sol <= min_diff_sol) min_diff_sol = diffangle_sol;
                if (diffangle_sol > max_diff_sol) max_diff_sol = diffangle_sol;
                if (zenith_sol <= min_zen_sol) min_zen_sol = zenith_sol;
                if (zenith_sol > max_zen_sol) max_zen_sol = zenith_sol;
            }
        }

        dlik = max_lik_up - max_lik_down;
        if (jlik > 0) normdlik = dlik / jlik;
        // Block End //


        // Block: Compare best tracks of muon and shower // Platinum
        // Initialize with impossible values
        diffangle_track_shower = 1000;
        delta_zenith_track_shower = 1000;

        if (is_shower_event) {

            // Tracks' angle difference between best muon track and best shower track
            diffangle_track_shower = acos(best_trk_shower.dir.dot(best_trk.dir)) * TMath::RadToDeg();

            // Zenith angle difference between best muon track and best shower track
            delta_zenith_track_shower = zenith - zenith_shower;
        }
        print_mode() << "\033[1m[Best Muon, Best Shower]\033[22m\n # Diffangle: " << diffangle_track_shower << "\n # Delta Zenith: " << delta_zenith_track_shower << "\n";
        // Block End //


        /*
            ====================================
            |                                  |
            |      Mixed Tracks Analysis       |
            |                                  |
            ====================================
        */
        print_mode() << "\033[1;35m[True-Reco Tracks Comparison]\033[22;39m\n";


        // Block: Compare Reconstructed tracks to True tracks and get angles // Platinum
        // Initialize with impossible values
        diffangle = 1000;
        delta_zenith = 1000;
        diffangle_shower_mu = 1000;
        diffangle_shower_nu = 1000;
        delta_zenith_shower_mu = 1000;
        delta_zenith_shower_nu = 1000;

        if (event.mc_trks.size() > 0) {

            // Tracks' angle difference between best track and mu track
            diffangle = acos(mu.dir.dot(best_trk.dir)) * TMath::RadToDeg();

            // Zenith angle difference between best track and mu track
            delta_zenith = zenith_mu - zenith;

            // If shower exists
            if (is_shower_event) {

                // Tracks' angle difference between best shower track and (mu/nu) track
                diffangle_shower_mu = acos(mu.dir.dot(best_trk_shower.dir)) * TMath::RadToDeg();
                diffangle_shower_nu = acos(nu.dir.dot(best_trk_shower.dir)) * TMath::RadToDeg();

                // Zenith angle difference between best shower track and (mu/nu) track
                delta_zenith_shower_mu = zenith_mu - zenith_shower;
                delta_zenith_shower_nu = zenith_nu - zenith_shower;
            }
        }

        // Print Results
        print_mode() << "\033[1m[Best Muon, True Muon]\033[22m\n # Diffangle: " << diffangle << "\n # Delta Zenith: " << delta_zenith << "\n";
        print_mode() << "\033[1m[Best Shower, True Muon]\033[22m\n # Diffangle: " << diffangle_shower_mu << "\n # Delta Zenith: " << delta_zenith_shower_mu << "\n";
        print_mode() << "\033[1m[Best Shower, True Neutrino]\033[22m\n # Diffangle: " << diffangle_shower_nu << "\n # Delta Zenith: " << delta_zenith_shower_nu << "\n";
        // Block End //


        /*
            ====================================
            |                                  |
            |          Hits Analysis           |
            |                                  |
            ====================================
        */
        print_mode() << "\033[1;35m[Reconstructed Hits Analysis]\033[22;39m\n";


        // Block: Fill vectors of pointers to triggered, cherenkov and cascade hits // Platinum
        // Initialize
        num_triggered_hits = 0;
        num_cherenkov_hits = 0;
        num_cascade_hits = 0;
        triggered_hits.clear();
        cherenkov_hits.clear();
        cascade_hits.clear();

        for (unsigned i_hit = 0; i_hit < event.hits.size(); ++i_hit) {

            // Guard Clause: Trigger
            if (!event.hits[i_hit].trig > 0) continue;

            // Increment triggered hits
            ++num_triggered_hits;

            // Get Cherenkov output
            cherenkov_pars(best_trk, event.hits[i_hit].pos, event.hits[i_hit].dir, reco_time, reco_d_track, reco_d_photon, reco_d_closest, reco_cos_angle);

            reco_tres = event.hits[i_hit].t - reco_time;

            if (cherenkov_condition(reco_cos_angle, reco_d_closest, reco_tres)) {

                // Increment Cherenkov hits
                ++num_cherenkov_hits;
            }

            if (cherenkov_condition_cascade(best_trk, event.hits[i_hit])) {

                // Increment Cascade hits
                ++num_cascade_hits;
            }
        }

        // reserve size for vector of triggered hits. It's more efficient than push_back in previous loop
        triggered_hits.reserve(num_triggered_hits);
        cherenkov_hits.reserve(num_cherenkov_hits);
        cascade_hits.reserve(num_cascade_hits);

        for (unsigned i_hit = 0; i_hit < event.hits.size(); ++i_hit) {

            // Guard Clause: Trigger
            if (!event.hits[i_hit].trig > 0) continue;

            // Save pointer to triggered hit
            triggered_hits.push_back(&event.hits[i_hit]);

            // Get Cherenkov output
            cherenkov_pars(best_trk, event.hits[i_hit].pos, event.hits[i_hit].dir, reco_time, reco_d_track, reco_d_photon, reco_d_closest, reco_cos_angle);

            reco_tres = event.hits[i_hit].t - reco_time;

            if (cherenkov_condition(reco_cos_angle, reco_d_closest, reco_tres)) {

                // Save pointer to Cherenkov hit
                cherenkov_hits.push_back(&event.hits[i_hit]);
            }

            if (cherenkov_condition_cascade(best_trk, event.hits[i_hit])) {

                // Save pointer to cascade hit
                cascade_hits.push_back(&event.hits[i_hit]);
            }

        }

        // Extract info
        triggered_hits_per_event.push_back(num_triggered_hits);
        cherenkov_hits_per_event.push_back(num_cherenkov_hits);
        cascade_hits_per_event.push_back(num_cascade_hits);
        // Block End //



        // Block: Get triggered / Cherenkov, Doms / Pmts // Gold+
        // Generic iterator
        std::map<int, int>::iterator it;
        Pmt current_pmt;
        Dom current_dom;
        int line;

        // Clear maps
        triggered_lines.clear();
        triggered_doms.clear();
        triggered_pmts.clear();
        cherenkov_lines.clear();
        cherenkov_doms.clear();
        cherenkov_pmts.clear();
        cascade_lines.clear();
        cascade_doms.clear();
        cascade_pmts.clear();

        // Triggered Hits Loop
        for (unsigned i_hit = 0; i_hit < triggered_hits.size(); ++i_hit) {

            // Extract line
            current_pmt = det.get_pmt(triggered_hits[i_hit]->dom_id, triggered_hits[i_hit]->channel_id);
            current_dom = det.get_dom(current_pmt);
            line = current_dom.line_id;

            // Increment hit count in maps
            triggered_lines[line]++;
            triggered_doms[triggered_hits[i_hit]->dom_id]++;
            triggered_pmts[triggered_hits[i_hit]->pmt_id]++;
        }

        // Cherenkov Hits Loop
        for (unsigned i_hit = 0; i_hit < cherenkov_hits.size(); ++i_hit) {

            // Extract line
            current_pmt = det.get_pmt(cherenkov_hits[i_hit]->dom_id, cherenkov_hits[i_hit]->channel_id);
            current_dom = det.get_dom(current_pmt);
            line = current_dom.line_id;

            // Increment hit count in maps
            cherenkov_lines[line]++;
            cherenkov_doms[cherenkov_hits[i_hit]->dom_id]++;
            cherenkov_pmts[cherenkov_hits[i_hit]->pmt_id]++;
        }

        // Cascade Hits Loop
        for (unsigned i_hit = 0; i_hit < cascade_hits.size(); ++i_hit) {

            // Extract line
            current_pmt = det.get_pmt(cascade_hits[i_hit]->dom_id, cascade_hits[i_hit]->channel_id);
            current_dom = det.get_dom(current_pmt);
            line = current_dom.line_id;

            // Increment hit count in maps
            cascade_lines[line]++;
            cascade_doms[cascade_hits[i_hit]->dom_id]++;
            cascade_pmts[cascade_hits[i_hit]->pmt_id]++;
        }

        // Extract info
        num_triggered_lines = triggered_lines.size();
        num_triggered_doms = triggered_doms.size();
        num_triggered_pmts = triggered_pmts.size();
        num_cherenkov_lines = cherenkov_lines.size();
        num_cherenkov_doms = cherenkov_doms.size();
        num_cherenkov_pmts = cherenkov_pmts.size();
        num_cascade_lines = cascade_lines.size();
        num_cascade_doms = cascade_doms.size();
        num_cascade_pmts = cascade_pmts.size();
        triggered_lines_per_event.push_back(num_triggered_lines);
        triggered_doms_per_event.push_back(num_triggered_doms);
        triggered_pmts_per_event.push_back(num_triggered_pmts);
        cherenkov_lines_per_event.push_back(num_cherenkov_lines);
        cherenkov_doms_per_event.push_back(num_cherenkov_doms);
        cherenkov_pmts_per_event.push_back(num_cherenkov_pmts);
        cascade_lines_per_event.push_back(num_cascade_lines);
        cascade_doms_per_event.push_back(num_cascade_doms);
        cascade_pmts_per_event.push_back(num_cascade_pmts);

        // If denominator equals 0 (bad division), assign -1 to let user know.
        ratio_cherenkov_hits = (num_triggered_hits > 0) ? num_cherenkov_hits * 1.0 / num_triggered_hits : -1;
        ratio_cherenkov_lines = (num_triggered_lines > 0) ? num_cherenkov_lines * 1.0 / num_triggered_lines : -1;
        ratio_cherenkov_doms = (num_triggered_doms > 0) ? num_cherenkov_doms * 1.0 / num_triggered_doms : -1;
        ratio_cherenkov_pmts = (num_triggered_pmts > 0) ? num_cherenkov_pmts * 1.0 / num_triggered_pmts : -1;

        ratio_cascade_hits = (num_triggered_hits > 0) ? num_cascade_hits * 1.0 / num_triggered_hits : -1;
        ratio_cascade_lines = (num_triggered_lines > 0) ? num_cascade_lines * 1.0 / num_triggered_lines : -1;
        ratio_cascade_doms = (num_triggered_doms > 0) ? num_cascade_doms * 1.0 / num_triggered_doms : -1;
        ratio_cascade_pmts = (num_triggered_pmts > 0) ? num_cascade_pmts * 1.0 / num_triggered_pmts : -1;

        // Print output for triggered counters
        print_mode() << "\033[1m[Triggered Counters]\033[22m\n";
        print_mode() << " # Hits : " << num_triggered_hits << "\n # PMTs : " << num_triggered_pmts;
        print_mode() << "\n # DOMs : " << num_triggered_doms << "\n # DUs  : " << num_triggered_lines << "\n";

        // Print output for cherenkov counters
        print_mode() << "\033[1m[Cherenkov Counters]\033[22m\n";
        print_mode() << " # Hits : " << num_cherenkov_hits << "\n # PMTs : " << num_cherenkov_pmts;
        print_mode() << "\n # DOMs : " << num_cherenkov_doms << "\n # DUs  : " << num_cherenkov_lines << "\n";

        // Print output for cascade counters
        print_mode() << "\033[1m[Cascade Counters]\033[22m\n";
        print_mode() << " # Hits : " << num_cascade_hits << "\n # PMTs : " << num_cascade_pmts;
        print_mode() << "\n # DOMs : " << num_cascade_doms << "\n # DUs  : " << num_cascade_lines << "\n";

        // Print output for cherenkov / trigger ratios
        print_mode() << "\033[1m[Cherenkov Ratios]\033[22m\n";
        print_mode() << " # Hits : " << ratio_cherenkov_hits << "\n # PMTs : " << ratio_cherenkov_pmts;
        print_mode() << "\n # DOMs : " << ratio_cherenkov_doms << "\n # DUs  : " << ratio_cherenkov_lines << "\n";

        // Print output for cascade / trigger ratios
        print_mode() << "\033[1m[Cascade Ratios]\033[22m\n";
        print_mode() << " # Hits : " << ratio_cascade_hits << "\n # PMTs : " << ratio_cascade_pmts;
        print_mode() << "\n # DOMs : " << ratio_cascade_doms << "\n # DUs  : " << ratio_cascade_lines << "\n";

        // Block End //



        // Block: Pikounis LOM Finder revised by Zarpapis // Platinum
        // Initialize
        const Hit* largest_hit = 0;
        int largest_dom;
        //int largest_pmt;  // Unused
        std::map<int, unsigned> tot_per_dom;
        std::map<int, unsigned> tot_per_pmt;
        std::map<int, unsigned>::iterator largest_tot;
        std::vector<const Hit*> hit_candidates;

        // Get sum of tot for each DOM, PMT (referrenced to by their IDs)
        for (int i_hit = 0; i_hit < num_triggered_hits; ++i_hit) {

            tot_per_dom[triggered_hits[i_hit]->dom_id] += triggered_hits[i_hit]->tot;
            tot_per_pmt[triggered_hits[i_hit]->pmt_id] += triggered_hits[i_hit]->tot;
        }

        // Find largest DOM
        for (std::map<int, unsigned>::iterator i = tot_per_dom.begin(); i != tot_per_dom.end(); ++i) {

            // Initialize
            if (i == tot_per_dom.begin()) largest_tot = i;

            // Compare
            if (largest_tot->second < i->second) largest_tot = i;
        }

        // Get ID of DOM with largest total in .tot
        largest_dom = largest_tot->first;

        // Reserve size for vector of Hit pointers, by using the counted hits in triggered_doms map
        hit_candidates.reserve(triggered_doms[largest_dom]);

        // Find the hits that belong to the largest DOM
        for (int i_hit = 0; i_hit < num_triggered_hits; ++i_hit) {
            if (triggered_hits[i_hit]->dom_id == largest_dom) hit_candidates.push_back(triggered_hits[i_hit]);
        }

        // Scan through found hits to get the largest.
        for (int i_hit = 0; i_hit < triggered_doms[largest_dom]; ++i_hit) {

            // Initialize
            if (!largest_hit) largest_hit = hit_candidates[i_hit];

            // Compare
            if (largest_hit->tot < hit_candidates[i_hit]->tot) largest_hit = hit_candidates[i_hit];

            // If equal tot, get hit with earlier time
            else if (largest_hit->tot == hit_candidates[i_hit]->tot && largest_hit->t > hit_candidates[i_hit]->t) largest_hit = hit_candidates[i_hit];
        }

        // Extract info
        float time_of_largest_pulse = largest_hit->t;
        float LOMx = largest_hit->pos.x;
        float LOMy = largest_hit->pos.y;
        float LOMz = largest_hit->pos.z;

        // Print Results
        // Print results
        if (true) {
            print_mode() << "\033[1m[Largest OM, Largest Hit]\033[22m : ";
            if (!largest_hit) {
                print_mode() << "\033[31mNot Found.\n  \033[39m There were no hits in this event.\n";
            }
            else {
                print_mode() << "\033[32mFound\033[39m\n ";
                //largest_hit->print();
                print_mode() << "\n";
            }
        }
        // Block End //



        // Block: d_track && d_closest to find first photon position // Platinum
        // Initialization
        double* reco_d_track_min_ptr = 0; // NULL pointer
        double* reco_d_track_max_ptr = 0; // NULL pointer
        double reco_d_track_min;
        double reco_d_track_max;
        double reco_d_track_mean = 0;
        double reco_d_track_rms = 0;

        double* reco_d_closest_min_ptr = 0; // NULL pointer
        double* reco_d_closest_max_ptr = 0; // NULL pointer
        double reco_d_closest_min;  // Unused
        double reco_d_closest_max;  // Unused
        double reco_d_closest_mean = 0;
        double reco_d_closest_rms = 0;

        std::vector<double> d_track_vec;
        std::vector<double> d_closest_vec;
        std::vector<double> d_track_diff; // = d_track_vec[i] - reco_d_track_mean
        std::vector<double> d_closest_diff; // = d_closest_vec[i] - reco_d_closest_mean

        for (unsigned i_hit = 0; i_hit < event.hits.size(); ++i_hit) {
            const Hit* current_hit = &event.hits[i_hit];

            // Guard Clause
            if (!current_hit->trig > 0) continue;

            // Get Cherenkov output
            cherenkov_pars(best_trk, current_hit->pos, current_hit->dir, reco_time, reco_d_track, reco_d_photon, reco_d_closest, reco_cos_angle);

            // Time resolution
            reco_tres = current_hit->t - reco_time;

            // Guard Clause: Cherenkov Condition
            if (!cherenkov_condition(reco_cos_angle, reco_d_closest, reco_tres)) continue;

            // Include this for next block of code
            d_track_vec.push_back(reco_d_track);
            d_closest_vec.push_back(reco_d_closest);
        }

        // Initialize all NULL pointers with first value, if non-empty vector
        if (d_track_vec.size() > 0) {
            reco_d_track_min_ptr = &d_track_vec[0];
            reco_d_track_max_ptr = &d_track_vec[0];
            reco_d_closest_min_ptr = &d_closest_vec[0];
            reco_d_closest_max_ptr = &d_closest_vec[0];
        }

        for (unsigned i = 0; i < d_track_vec.size(); ++i) {

            // Make comparisons for min / max
            if (*reco_d_track_min_ptr > d_track_vec[i]) reco_d_track_min_ptr = &d_track_vec[i];
            if (*reco_d_track_max_ptr < d_track_vec[i]) reco_d_track_max_ptr = &d_track_vec[i];
            if (*reco_d_closest_min_ptr > d_closest_vec[i]) reco_d_closest_min_ptr = &d_closest_vec[i];
            if (*reco_d_closest_max_ptr < d_closest_vec[i]) reco_d_closest_max_ptr = &d_closest_vec[i];

            // Calculate sums for means
            reco_d_track_mean += d_track_vec[i];
            reco_d_closest_mean += d_closest_vec[i];

            // Calculate sums of squares for rms
            reco_d_track_rms += d_track_vec[i] * d_track_vec[i];
            reco_d_closest_rms += d_closest_vec[i] * d_closest_vec[i];
        }

        // Calculate means
        reco_d_track_mean /= d_track_vec.size();
        reco_d_closest_mean /= d_closest_vec.size();

        // Calculate rms
        reco_d_track_rms /= d_track_vec.size();
        reco_d_track_rms = std::sqrt(reco_d_track_rms);
        reco_d_closest_rms /= d_closest_vec.size();
        reco_d_closest_rms = std::sqrt(reco_d_closest_rms);

        // Prepare differences
        d_track_diff = d_track_vec;
        d_closest_diff = d_closest_vec;

        // Calculate differences
        for (unsigned i = 0; i < d_track_vec.size(); ++i) {
            d_track_diff[i] -= reco_d_track_mean;
            d_closest_diff[i] -= reco_d_closest_mean;
        }

        // De-reference
        if (d_track_vec.size() > 0) {
            reco_d_track_min = *reco_d_track_min_ptr;
            reco_d_track_max = *reco_d_track_max_ptr;
        }
        else {
            reco_d_track_min = 0;
            reco_d_track_max = 0;
        }

        reco_d_closest_min = *reco_d_closest_min_ptr; // Unused
        reco_d_closest_max = *reco_d_closest_max_ptr; // Unused

        // Position of first photon
        Vec d_track_min_vec = best_trk.pos + best_trk.dir * reco_d_track_min;
        // Block End //


        // Block: Find first hit (by time) // Platinum
        // Initialize
        const Hit* first_hit_candidate = 0;  // NULL pointer
        const Hit* first_hit = 0;            // NULL pointer
        double first_d_track;

        for (unsigned i_hit = 0; i_hit < event.hits.size(); ++i_hit) {
            first_hit_candidate = &event.hits[i_hit];
	    if( (first_hit_candidate->t > 0.0 ) && (first_hit_candidate->trig > 0) ){
                // if first_hit pointer is uninitialized (NULL), initialize with candidate pointer
                if (first_hit == 0)
                    first_hit = first_hit_candidate;

                // first_hit pointer is surely initialized, so its time is compared with candidate's
                if (first_hit_candidate->t < first_hit->t)
                    first_hit = first_hit_candidate;
            }
        }

        // Get d_track of first hit
        cherenkov_pars(best_trk, first_hit->pos, first_hit->dir, reco_time, first_d_track, reco_d_photon, reco_d_closest, reco_cos_angle);

        // Print results
        if (true) {
            print_mode() << "\033[1m[First Hit]\033[22m : ";
            if (first_hit_candidate == 0) {
                print_mode() << "\033[31mNot Found.\n  \033[39m There were no hits in this event.\n";
            }
            else {
                if (first_hit == 0) {
                    print_mode() << "\033[31mNot Found.\n  \033[39m Found hits but not \033[3mtriggered\033[23m hits in this event.\n";
                }
                else {
                    print_mode() << "\033[32mFound\033[39m\n ";
                    //first_hit->print();
                    print_mode() << "\n";
                }
            }
        }
        // Block End //


        // Block: Minimum distance from track to DOM ???? Why???
        // Initialize
        double d_closest_min;
        double ag = 0;
        initialized_flag = false; // Helps with first initialization

        for (unsigned i_hit = 0; i_hit < triggered_hits.size(); ++i_hit) {
            const Hit* current_hit = triggered_hits[i_hit];

            // Guard Clause
            if (!current_hit->trig > 0) continue;

            // Get Cherenkov output
            cherenkov_pars(best_trk, current_hit->pos, current_hit->dir, reco_time, reco_d_track, reco_d_photon, reco_d_closest, reco_cos_angle);

            // Initialize d_closest_min
            if (!initialized_flag) {
                d_closest_min = reco_d_closest;
                initialized_flag = true;
            }

            // Calculate
            if (d_closest_min > reco_d_closest) {
                d_closest_min = reco_d_closest;
                ag = (current_hit->pos - best_trk.pos).dot(best_trk.dir);
            }
        }

        Vec closest_point = best_trk.pos + best_trk.dir * ag;
        D3d = closest_point.len();
        D = closest_point.lenxy();
        // Block End //


        // Block: Jannik's proposal, distance of first / last photon to the track position // Platinum
        Vec trkPos_FirstPhoton = best_trk.pos + best_trk.dir * reco_d_track_min - mean_detector_pos;
        Vec trkPos_LastPhoton = best_trk.pos + best_trk.dir * reco_d_track_max - mean_detector_pos;
        Vec trkPos_MeanPhoton = best_trk.pos + best_trk.dir * reco_d_track_mean - mean_detector_pos;

        // Print Results
        print_mode() << "\033[1m[Yannik's Proposal]\033[22m\n";
        print_mode() << " First Photon's pos: (" << trkPos_FirstPhoton.x << ", " << trkPos_FirstPhoton.y << ", " << trkPos_FirstPhoton.z << ")\n";
        print_mode() << " Last Photon's pos : (" << trkPos_LastPhoton.x << ", " << trkPos_LastPhoton.y << ", " << trkPos_LastPhoton.z << ")\n";
        print_mode() << " Mean Photon's pos : (" << trkPos_MeanPhoton.x << ", " << trkPos_MeanPhoton.y << ", " << trkPos_MeanPhoton.z << ")\n";
        // Block End //


        // Block: Under construction //
        print_mode() << "\033[1m[NTrack Analysis]\033[22m\n";
        bool ismuon = false;
        bool iscasc = false;

        float sum_pos_z_it = 0.0;
        float sum_pos_z_trig = 0.0;
        float sum_tres = 0.0;

        // parametrization of the shower maximum //number extracted from $JPP/software/JPhysics/JGeanz.hh
        const float a0 = 1.85;
        const float a1 = 0.62;
        const float b = 0.54;

        // MC Maps
        double casc_time = -999.9;
        double d_shower_photon = 0.0;
        double d_shower = 0;

        /*
        for (unsigned i_hit = 0; i_hit < event.hits.size(); ++i_hit){
            const Hit* current_hit = &event.hits[i_hit];

            // Guard Clause
            if (!current_hit->trig > 0) continue;

            // Get Cherenkov output
            cherenkov_pars(best_trk, current_hit->pos, current_hit->dir, reco_time, reco_d_track, reco_d_photon, reco_d_closest, reco_cos_angle);

            current_pmt = det.get_pmt(current_hit->dom_id, current_hit->channel_id);
            current_dom = det.get_dom(current_pmt);
            line = current_dom.line_id;

            // Cherenkov Condition;
            if (reco_cos_angle < 0 && reco_d_closest < reco_Vert_dist_for_coinc && fabs(reco_tres) < 10.0) {

            }

            if (best_trk.rec_type == JPP_RECONSTRUCTION_TYPE && best_trk.E > TMath::Exp(-a0 / a1)) {
                d_shower = (a0 + a1 * TMath::Log(best_trk.E) - 1) * b;
            }

        }
        */

        for (unsigned i_hit = 0; i_hit < event.hits.size(); ++i_hit) {
            Hit reco_current_hit = event.hits[i_hit];
            int trig = reco_current_hit.trig;
            int current_dom_id = reco_current_hit.dom_id;
            int current_channel_id = reco_current_hit.channel_id;
            Pmt current_pmt = det.get_pmt(current_dom_id, current_channel_id);


            float distance_from_LOM = sqrt((LOMx - reco_current_hit.pos.x) * (LOMx - reco_current_hit.pos.x) + (LOMy - reco_current_hit.pos.y) * (LOMy - reco_current_hit.pos.y) + (LOMz - reco_current_hit.pos.z) * (LOMz - reco_current_hit.pos.z));

            cherenkov_pars(best_trk, reco_current_hit.pos, reco_current_hit.dir, reco_time, reco_d_track, reco_d_photon, reco_d_closest, reco_cos_angle);
            reco_tres = reco_current_hit.t - reco_time;  // true time - expected by Chere time

            if (trig > 0) {


                int reco_current_dom_id = reco_current_hit.dom_id;
                int reco_current_channel_id = reco_current_hit.channel_id;
                Pmt reco_current_pmt = det.get_pmt(reco_current_dom_id, reco_current_channel_id);
                Dom reco_current_dom = det.get_dom(reco_current_pmt);
                int reco_floor = reco_current_dom.floor_id;
                double reco_hit_tot = reco_current_hit.tot;

                if (best_trk.rec_type == JPP_RECONSTRUCTION_TYPE) {
                    if (best_trk.E > TMath::Exp(-a0 / a1)) {
                        d_shower = (a0 + a1 * TMath::Log(best_trk.E) - 1) * b;
                    }
                }  // xrisimopoiei tin corrected energy


                sum_pos_z_trig += reco_current_hit.pos.z;


                // Casc Cherenkov//
                Vec v = reco_current_hit.pos - (best_trk.pos + best_trk.dir * d_shower);
                d_shower_photon = TMath::Sqrt(v.dot(v));
                casc_time = best_trk.t + d_shower / c_light + d_shower_photon / v_light;
                /////////////////

                double reco_casc_tres = reco_current_hit.t - casc_time;
                sum_tres = sum_tres + reco_tres;
                float ML1lambda = reco_d_track;
                float ML1lambda_2 = reco_d_track - reco_d_track_min;

                // temp solution
                double temp_time, temp_d_track, temp_d_photon, temp_d_closest, temp_cos_angle;
                cherenkov_pars(best_trk, first_hit->pos, first_hit->dir, temp_time, temp_d_track, temp_d_photon, temp_d_closest, temp_cos_angle);

                float ML1lambda_3 = reco_d_track - temp_d_track;
                //float ML1lambda_2 = reco_d_track - minML1lamda;

                // ----- Dtres ---- //
                if (reco_d_closest < reco_Vert_dist_for_coinc && fabs(reco_tres) < 10.0) {

                    if (reco_floor == 18 || reco_floor == 17) {  // || reco_floor == 16){
                        int_Nborder_dtres_hits++;
                        Nborder_dtres_hits = float(int_Nborder_dtres_hits);
                    }
                }

                //--Cherenkov hypothesis--//
                if (reco_cos_angle < 0 && reco_d_closest < reco_Vert_dist_for_coinc && fabs(reco_tres) < 10.0) {
                    ismuon = 1;
                    if (reco_current_hit.channel_id < 12) {  // pmt hits in the upper hemishpere of DOM
                        int_Nhits_cherenkov_upper++;
                        Nhits_cherenkov_upper = float(int_Nhits_cherenkov_upper);
                    }
                    if (reco_current_hit.channel_id > 11) {
                        int_Nhits_cherenkov_lower++;
                        Nhits_cherenkov_lower = float(int_Nhits_cherenkov_lower);
                    }

                    if (reco_floor == 18 || reco_floor == 17) {  // || reco_floor == 16){
                        int_Nborder_cherenkov_hits++;
                        Nborder_cherenkov_hits = float(int_Nborder_cherenkov_hits);
                        if (reco_current_hit.channel_id < 12) {
                            int_Nhits_border_cherenkov_upper++;
                            Nhits_border_cherenkov_upper = float(int_Nhits_border_cherenkov_upper);
                        }
                        if (reco_current_hit.channel_id > 11) {
                            int_Nhits_border_cherenkov_lower++;
                            Nhits_border_cherenkov_lower = float(int_Nhits_border_cherenkov_lower);
                        }
                    }  // end of border doms
                }      // end of Cherenkov hypothesis

                // all (not only cherenkov)
                if (reco_current_hit.channel_id < 12) {
                    int_Nhits_upper++;
                    Nhits_upper = float(int_Nhits_upper);
                }
                if (reco_current_hit.channel_id > 11) {
                    int_Nhits_lower++;
                    Nhits_lower = float(int_Nhits_lower);
                }
                if (reco_floor == 18 || reco_floor == 17) {  // || reco_floor == 16){
                    int_Nborder_hits++;
                    Nborder_hits = float(int_Nborder_hits);
                    if (reco_current_hit.channel_id < 12) {
                        int_Nhits_border_upper++;
                        Nhits_border_upper = float(int_Nhits_border_upper);
                    }
                    if (reco_current_hit.channel_id > 11) {
                        int_Nhits_border_lower++;
                        Nhits_border_lower = float(int_Nhits_border_lower);
                    }
                }

                if (d_shower_photon < reco_Vert_dist_for_coinc && fabs(reco_casc_tres) < 10.0) {
                    iscasc = 1;
                }

                // Alfonsos vars
                Vec vdist = reco_current_hit.pos - best_trk.pos;
                float dist = vdist.len(); // TMath::Sqrt(vdist.dot(vdist));
                double true_dist = (reco_current_hit.pos - mu.pos).len();


                Vec l_min_lambda = best_trk.pos + (best_trk.dir * reco_d_track_min);
                Vec lmin_pos1g = best_trk.pos + (best_trk.dir * temp_d_track);


                double mydist = (reco_current_hit.pos - l_min_lambda).len();
                double mydist1g = (reco_current_hit.pos - lmin_pos1g).len();



                // print_mode()<<" dist : "<<dist<<std::endl;

                if (ismuon && iscasc) {
                    nhits_casc_mu++;
                    if (std::isnan(dist)) {

                    }
                    else {
                        if (dist < 100)
                            nhits_casc_mu_100++;
                        if (dist < min_dist_casc_mu)
                            min_dist_casc_mu = dist;
                        if (dist > max_dist_casc_mu)
                            max_dist_casc_mu = dist;
                    }
                    sum_ToT_casc_mu += reco_hit_tot;
                    if (float(reco_floor) > 16.5) {
                        ToT_border_cascmu += reco_hit_tot;
                    }
                }

                if (ismuon && !iscasc) {
                    nhits_mu++;
                    if (std::isnan(dist)) {
                    }
                    else {
                        if (dist < 100)
                            nhits_mu_100++;
                        if (dist < min_dist_mu)
                            min_dist_mu = dist;
                        if (dist > max_dist_mu)
                            max_dist_mu = dist;
                    }
                    sum_ToT_mu += reco_hit_tot;
                    if (float(reco_floor) > 16.5) {
                        ToT_border_mu += reco_hit_tot;
                    }
                }

                if (!ismuon && iscasc) {
                    nhits_casc++;
                    if (std::isnan(dist)) {
                    }
                    else {
                        if (dist < 100)
                            nhits_casc_100++;
                        if (dist < min_dist_casc)
                            min_dist_casc = dist;
                        if (dist > max_dist_casc)
                            max_dist_casc = dist;
                    }
                    sum_ToT_casc += reco_hit_tot;
                    if (float(reco_floor) > 16.5) {
                        ToT_border_casc += reco_hit_tot;
                    }
                }

                // My new dist : current hit from lambda point not from vtx//
                //double mydist = (reco_current_hit.pos - d_track_min_vec).len();


                if (ismuon && iscasc) {
                    if (mydist < 50.)
                        nhits_casc_mu_50++;
                    if (mydist < 30.)
                        nhits_casc_mu_30++;
                }
                if (ismuon && !iscasc) {
                    if (mydist < 50.)
                        nhits_mu_50++;
                    if (mydist < 30.)
                        nhits_mu_30++;
                }
                if (!ismuon && iscasc) {
                    if (mydist < 50.)
                        nhits_casc_50++;
                    if (mydist < 30.)
                        nhits_casc_30++;
                }

                ToT_trig += reco_current_hit.tot;
                if (max_ToT_trig < reco_current_hit.tot) {
                    max_ToT_trig = reco_current_hit.tot;
                }

                /// ---------- In time hits from ka.Katerina --------- ///
                if (fabs(fabs(time_of_largest_pulse - reco_current_hit.t) - distance_from_LOM / c_light) < 500.0) {

                    float DistOMVTX = (d_track_min_vec - reco_current_hit.pos).len();

                    Ntrack_all2++;

                    if ((ML1lambda < -30.0)) {  // previously 30 or 10
                        Nallbehind2++;
                    }

                    if ((fabs(reco_tres) < 400.0)) {
                        if ((ML1lambda < -30.0)) {  // previously 30 or 10
                            Nallbehind3++;
                        }
                    }

                    if ((fabs(reco_tres) < 250.0)) {  // this means om behind vtx
                        Ntrack_all++;

                        if (fabs(reco_tres) <= 10.0) {
                            if (ML1lambda < strict_lambda_min) {
                                strict_lambda_min = ML1lambda;
                            }
                        }

                        if ((ML1lambda < -30.0) && (DistOMVTX < 165.0)) {  // previously 30 or 10
                            Nallbehind++;
                        }

                        if ((fabs(reco_tres) < 10.0) && (reco_d_closest < 120.0)) {

                            NtrackIT++;
                            ToT_IT += reco_current_hit.tot;
                            sum_pos_z_it += reco_current_hit.pos.z;

                            if (ML1lambda > lambdaIT_max) { lambdaIT_max = ML1lambda; }

                            if ((ML1lambda > 50.0)) { NtrackIT50++; }
                            if ((ML1lambda > 30.0)) { NtrackIT30++; }
                            if ((ML1lambda > 10.0)) { NtrackIT10++; }

                            if ((ML1lambda_2 > 50.0)) { NtrackIT50_2++; }
                            if ((ML1lambda_2 > 30.0)) { NtrackIT30_2++; }
                            if ((ML1lambda_2 > 10.0)) { NtrackIT10_2++; }

                            if ((fabs(ML1lambda_3) > 50.0)) { NtrackIT50_3++; }
                            if ((fabs(ML1lambda_3) > 30.0)) { NtrackIT30_3++; }
                            if ((fabs(ML1lambda_3) > 10.0)) { NtrackIT10_3++; }

                            //  if (  reco_d_track - reco_d_track_min > lambdaIT_max ){ lambdaIT_max = ML1lambda;}
                            if (reco_cos_angle < 0.0) {
                                if (reco_d_track - reco_d_track_min > lambdaCH_max) {
                                    lambdaCH_max = ML1lambda;
                                }
                            }
                        }
                        else {  // end of if ( (fabs(reco_tres) <= 30.0 ) && (d_closest < 120.0) )  {
                            if ((reco_tres < -10.0) && (reco_d_closest < 120.0)) {
                                NtrackEarly++;
                            }  // end of if ( (fabs(reco_tres) < -30.0) && (d_closest < 120
                            if ((reco_tres > 20.0) && (reco_d_closest < 120.0)) {
                                NtrackLate++;
                            }
                        }

                    }  // end of if ( (fabs(reco_tres) < 250.0) ){

                }  // end of 500.0
            }      // end of trig

            // calculate tot of all in time hits (also not trig)
            if (fabs(fabs(time_of_largest_pulse - reco_current_hit.t) - distance_from_LOM / c_light) < 500.0) {
                if ((fabs(reco_tres) < 10.0) && (reco_d_closest < 120.0)) {
                    ToT_allIT += reco_current_hit.tot;
                }
            }

        }  // end of loop over i_hit

        if (!(std::isnan(theta_max)) || (std::isnan(theta_min))) {
            diff_theta = fabs(theta_max - theta_min);
        }
        if (nhits_mu > 0) {
            myratio50_muon = nhits_mu_50 / nhits_mu;
            myratio30_muon = nhits_mu_30 / nhits_mu;
            ratio_closehits_muon = nhits_mu_100 / nhits_mu;
            redToT_muon = sum_ToT_mu / nhits_mu;
            // borderToT_mu = /sum_ToT_mu;
            //like myratio50_cascmuon but over muons instead
            myratio50_cascmuon_over_mu = nhits_casc_mu_50 / nhits_mu;
            myratio30_cascmuon_over_mu = nhits_casc_mu_30 / nhits_mu;
            ratio_closehits_cascmuon_over_mu = nhits_casc_mu_100 / nhits_mu;
            redToT_cascmuon_over_mu = sum_ToT_casc_mu / nhits_mu;
            //like myratio50_casc, but over muons instead
            myratio50_casc_over_mu = nhits_casc_50 / nhits_mu;
            myratio30_casc_over_mu = nhits_casc_30 / nhits_mu;
            ratio_closehits_casc_over_mu = nhits_casc_100 / nhits_mu;
            redToT_casc_over_mu = sum_ToT_casc / nhits_mu;
        }
        if (nhits_casc_mu > 0) {
            myratio50_cascmuon = nhits_casc_mu_50 / nhits_casc_mu;
            myratio30_cascmuon = nhits_casc_mu_30 / nhits_casc_mu;
            ratio_closehits_cascmuon = nhits_casc_mu_100 / nhits_casc_mu;
            redToT_cascmuon = sum_ToT_casc_mu / nhits_casc_mu;
        }
        if (nhits_casc > 0) {
            myratio50_casc = nhits_casc_50 / nhits_casc;
            myratio30_casc = nhits_casc_30 / nhits_casc;
            ratio_closehits_casc = nhits_casc_100 / nhits_casc;
            redToT_casc = sum_ToT_casc / nhits_casc;
        }

        diff_dist_mu = max_dist_mu - min_dist_mu;
        diff_dist_casc_mu = max_dist_casc_mu - min_dist_casc_mu;
        diff_dist_casc = max_dist_casc - min_dist_casc;

        print_mode() << " diff_dist_mu : " << diff_dist_mu << " diff_dist_casc_mu : " << diff_dist_casc_mu << "\n";
        print_mode() << " min_dist_mu : " << min_dist_mu << " max_dist_mu : " << max_dist_mu << "\n";
        print_mode() << " diff_theta : " << diff_theta << "\n";
        print_mode() << " nhits_casc : " << nhits_casc << " nhits_casc_50 : " << nhits_casc_50 << " nhits_casc_mu_50  : " << nhits_casc_mu_50 << "\n";

        if (num_triggered_hits > 0) {
            mean_tres_it = float(sum_tres) / float(num_triggered_hits);
        }

        // *******  Track Length Calc ******* //
        if (strict_lambda_min > 9999990.0) {
            strict_lambda_min = lambdaIT_max;
        }
        TrLengthIT = lambdaIT_max;

        ///////////////// Ntrack Analysis////////////////////
        // ******  Cuts concerning the in time stuff ****** //
        // **** Early hits, Shower like hits, Late hits etc **** //

        if (NtrackIT10 > 0) {
            ratio110 = float(NtrackEarly) / float(NtrackIT10);
            ratio410 = float(Nallbehind2) / float(NtrackIT10);
        }
        else {
            ratio110 = ratio410 = -0.1;
        }

        if (NtrackIT30 > 0) {
            ratio130 = float(NtrackEarly) / float(NtrackIT30);
            ratio430 = float(Nallbehind2) / float(NtrackIT30);
            ratiol = float(NtrackLate) / float(NtrackIT30);
        }
        else {
            ratio130 = ratio430 = ratiol = -0.1;
        }

        if (NtrackIT > 0) {
            ratio1 = float(NtrackEarly) / float(NtrackIT);
            ratio2 = float(Nallbehind) / float(NtrackIT);
            ratio4 = float(Nallbehind2) / float(NtrackIT);
            ratio5 = float(Nallbehind3) / float(NtrackIT);
            ratio6 = float(NtrackIT10) / float(NtrackIT);
        }
        else {
            ratio1 = ratio2 = ratio4 = ratio5 = ratio6 = -0.1;
        }

        if (Ntrack_all > 0) {
            ratio3 = float(NtrackIT) / float(Ntrack_all);
            ratio310 = float(NtrackIT10) / float(Ntrack_all);
            ratio330 = float(NtrackIT30) / float(Ntrack_all);
        }
        else {
            ratio3 = ratio310 = ratio330 = -0.1;
        }

        for (unsigned int i_hit = 0; i_hit < event.hits.size(); i_hit++) {  // loop over event hits

            double time = -999.9;
            double d_track = -999.9;
            double d_photon = -999.9;
            double d_closest = -999.9;
            double cos_angle = -999.9;
            Hit hit = event.hits[i_hit];
            int trig = hit.trig;

            if (trig > 0) {
                cherenkov_pars(best_trk, hit.pos, hit.dir, time, d_track, d_photon, d_closest, cos_angle);

                if (cos_angle < 0 && d_closest < 120 && fabs(hit.t - time) < 10.0) {  // cher.

                    if (d_track < minML1lamda) {
                        minML1lamda = d_track;
                    }

                    if (d_track > maxML1lamda) {
                        maxML1lamda = d_track;
                    }
                }
            }
        }

        TrLengthIT_2 = maxML1lamda;
        if (TrLengthIT_2 > 0.) {
            itoverlen = float(NtrackIT) / TrLengthIT_2;
        }
        else {
            itoverlen = -1;
        }

        if (minML1lamda < 999.9) {
            TrLengthIT_3 = maxML1lamda - minML1lamda;
        }
        else {
            TrLengthIT_3 = 0.;
        }


	//====================================================================================================
	//==================================== Lilia's variables =============================================
	//====================================================================================================
	double total_weight_per_track_MAXhits= 0.0;
	int count_pmts_nohit = 0;
	float detHeight = 650.0 - 120.0; //from stuff in reco_vertex_in definition (maybe needs update?)
	float Lmax = detHeight + ( 2*range_detector_max - detHeight ) * zenith;
	float potential_length = 1.0; //dummy


	for (unsigned i_hit = 0; i_hit < event.hits.size(); ++i_hit) {
	  Hit reco_current_hit = event.hits[i_hit];
	  int trig = reco_current_hit.trig;
	  int reco_current_dom_id = reco_current_hit.dom_id;
	  int reco_current_channel_id = reco_current_hit.channel_id;
	  Pmt reco_current_pmt = det.get_pmt(reco_current_dom_id, reco_current_channel_id);
	  Dom reco_current_dom = det.get_dom(reco_current_pmt);
	  int reco_floor = reco_current_dom.floor_id;
	  double reco_hit_tot = reco_current_hit.tot;
	  
            
	  cherenkov_pars(best_trk, reco_current_hit.pos, reco_current_hit.dir, reco_time, reco_d_track, reco_d_photon, reco_d_closest, reco_cos_angle);
	  reco_tres = reco_current_hit.t - reco_time;  // true time - expected by Chere time


	  //pmt whithin photon's reach but with no hit in the end
	  double dotProduct = reco_current_pmt.dir.dot(best_trk.dir); //if not hit, which is the current PMT?
	  double normPMT = sqrt( reco_current_pmt.dir.x*reco_current_pmt.dir.x + reco_current_pmt.dir.y*reco_current_pmt.dir.y + reco_current_pmt.dir.z*reco_current_pmt.dir.z );
	  double normTrk = sqrt( best_trk.dir.x*best_trk.dir.x + best_trk.dir.y*best_trk.dir.y + best_trk.dir.z*best_trk.dir.z );
	  double sup_angle = acos( dotProduct / (normPMT*normTrk) )*TMath::RadToDeg();
	  double current_ML1lambda_2 = reco_d_track - reco_d_track_min; //current emission point
	  if(reco_d_closest <= 200 && current_ML1lambda_2 < TrLengthIT_3){//Lilia's thesis (pg. 87)
	    if(sup_angle > 6 && sup_angle <= 86){
	      if( !(trig > 0) ){
		count_pmts_nohit++;
	      }//no trig. hit
	    }//supplementary angle condition (is it needed?)
	  }//if in cylinder

	  
	  if (trig > 0) {
	    
	    //if( fabs(reco_tres) < 10.0 ){//Lilia's L1 pulse condition?
	    
	    total_weight_per_track_MAXhits += reco_d_closest / reco_d_closest_max;

	    //}//end if L1 pulse condition

	  }//end if trig. hit
	}//end loop over hits



	//-------------
	//var1
	Log_No_PMTs_with_dist_weights_MAXhits = ( TMath::Log10(total_weight_per_track_MAXhits) + 1. ) / 4.;
	
	//var2
	No_OMs_length_potential_length_chosen = 360.0 * ( num_triggered_doms / potential_length ) / 90.0; //1.maybe num_cherenkov_doms or other computation? 2.length in definition cancel out? if not, use TrLengthIT_3?

	//var3
	//sum_ToT_mu AND sum_ToT_casc_mu // Widthtot = TMath::Log10(total_pmt_width) / 6.;
	
	//var4
	PMTs_hit_to_nohit= TMath::Log10(num_triggered_hits) / TMath::Log10(count_pmts_nohit);

	//var5
	MaxLenPos_OvrMaxDist = potential_length / Lmax; //max lentgh of track over maximum length in detector 
	
	//var6
	TrLenIT_3_OvrMaxDist = TrLengthIT_3 / Lmax; //1st-last photon over maximum length in detector

	//var7
	//E_mu_max
	
	//=====================================================================================================
	//=====================================================================================================

	


        /*
            ====================================
            |                                  |
            |           Event Cuts             |
            |                                  |
            ====================================
        */

        // L0 test (anti-bugs): if true, everything is fine
        //L0_test = found_reco_track && best_trk.rec_type == 4000 && jlik > 0 && jbeta0 > 0 && Elen > 0 && Slen > 0 && all_in && Ereco > 0 && TMath::Log10(best_trk.E) > 0;
        L0_test = found_reco_track && best_trk.rec_type == 4000 && jlik > 0 && jbeta0 > 0 && Elen > 0 && Slen > 0 && all_in && Ereco > 0;

        // Also check neutrino energy, if neutrino event
        if (nu_bool) L0_test = L0_test && nu.E > 0;

        // Print L0 test results
        print_mode() << "\033[1m[L0 Test]\033[22m : ";
        if (L0_test) {
            print_mode() << "\033[32mPassed\033[39m\n";
            number_of_L0_events++;
        }
        else {
            print_mode() << "\033[31mFailed\033[39m\n";

            // Show what went wrong
            if (!found_reco_track) print_mode() << " -> found_reco_track == false\n";
            if (!best_trk.rec_type == 4000) print_mode() << " -> best_trk.rec_type == " << best_trk.rec_type << "\n";
            if (!jlik > 0) print_mode() << " -> jlik == " << jlik << "\n";
            if (!jbeta0 > 0) print_mode() << " -> jbeta0 == " << jbeta0 << "\n";
            if (!Elen > 0) print_mode() << " -> Elen == " << Elen << "\n";
            if (!Slen > 0) print_mode() << " -> Slen == " << Slen << "\n";
            if (!all_in) print_mode() << " -> all_in == false\n";
            if (!Ereco > 0) print_mode() << " -> Ereco == " << Ereco << "\n";
        }

        // Anti-noise test:
	//  L0_antinoise = false;
	// if (jlik > 50 ){
	//   anti_noise=true;
	// }
        anti_noise = false;
        if (jlik > 50 && GNhit > 20){
            anti_noise=true;
	    number_of_anti_noise++;
        }

        // Print anti-noise results
        print_mode() << "\033[1m[Anti-noise]\033[22m : ";
        print_mode() << (anti_noise ? "\033[32mPassed\033[39m\n" : "\033[31mFailed\033[39m\n");

        // L1 test: if true, everything is fine Arca8
        L1_test = false;
        if (logEreco > 2.7 && jlik > 50.0 && TrLengthIT_2 > 100.0 && logbeta0 < -1.5) {
            L1_test = true;
            number_of_L1_events++;
        }

        // ARCA 6
        // logEreco > 2, lik > 40, logbeta < -1.5, 


        // Print L1 test results
        print_mode() << "\033[1m[L1 Test]\033[22m : ";
        print_mode() << (L1_test ? "\033[32mPassed\033[39m\n" : "\033[31mFailed\033[39m\n");


        // Guard Clause : L0 Test
        if (!L0_test) continue;
        
        // Guard Clause : anti-noise Test
	// if (!anti_noise) continue;

        // Fill processed events loop
        processedEventsTree->Fill();

        associated_n_events++;

    } // END of events loop

    qualityParametersTree->Fill(); // Per file variables

    /*
        ====================================
        |                                  |
        |          Final Output            |
        |                                  |
        ====================================
    */

    print_mode(0) << "\n ~ Statistics ~ \n # Events: " << number_of_events << "\n # Reco Tracks: " << number_of_reco_tracks << "\n # 3DMuon Trigevts: " << events_muons  << "\n # L0 events: " << number_of_L0_events << "\n # anti-noise events: " << number_of_anti_noise << "\n # L1 events: " << number_of_L1_events;
    print_mode.endl();


    print_mode(0) << "\nWriting file...\n";
    file_out->Write();

    print_mode(0) << "Closing file...\n";
    file_out->Close();

    // Calculate run-time
    clock.stop();
    unsigned minutes = (unsigned)(clock.getTime() / 60000);
    unsigned seconds = (unsigned)((clock.getTime() / 1000) - minutes * 60);

    std::cout << "Program ending. Run Time: " << minutes << " min : " << seconds << " sec" << std::endl;
    return 0;
}
//============================================================================================================================



/*
    ====================================
    |                                  |
    |           Functions              |
    |                                  |
    ====================================
*/



bool cherenkov_condition(const double& cos_angle, const double& d_closest, const double& tres) {
    return cos_angle < 0 && d_closest < 120 && fabs(tres) < 10;
}

bool cherenkov_condition(const Trk& track, const Hit& hit) {
    double time;
    double d_track;
    double d_photon;
    double d_closest;
    double cos_angle;

    cherenkov_pars(track, hit.pos, hit.dir, time, d_track, d_photon, d_closest, cos_angle);

    return cos_angle < 0 && d_closest < 120 && fabs(hit.t - time) < 10;
}

bool cherenkov_condition_cascade(const Trk& track, const Hit& hit) {
    /*
        parametrization of the shower maximum //number extracted from $JPP/software/JPhysics/JGeanz.hh
        const float a0 = 1.85;
        const float a1 = 0.62;
        const float b = 0.54;
    */
    const double dndl = 0.0298;
    const double water_index = 1.3499;
    const double c_light = 299792458 * 1e-9;
    const double v_light = c_light / (water_index + dndl);

    double d_shower = 0;
    if (track.rec_type == JPP_RECONSTRUCTION_TYPE && track.E > TMath::Exp(-1.85 / 0.62)) d_shower = (0.62 * TMath::Log(track.E) + 0.85) * 0.54;

    Vec v = hit.pos - (track.pos + track.dir * d_shower);
    double d_shower_photon = v.len();
    double cascade_time = track.t + d_shower / c_light + d_shower_photon / v_light;

    return d_shower_photon < 120 && fabs(hit.t - cascade_time) < 10;
}

// Fluxes
float random_flux(float E) {
    float f = 1.2 * std::pow(10, -8) * std::pow(E, -2);  // flux to compare to Lilia's
    return 0.5 * 10000.0 * f;                            // have flux in GeV^-1 s^-1 m^-2 sr^-1
}

float LoI_flux(float E) {
    float f = 1.2 * std::pow(10, -8) * std::pow(E, -2) * std::exp(-E / 3000000.0);  // flux from LOI in GeV^-1 s^-1 cm^-2 sr^-1
    return 0.5 * 10000.0 * f;                                                       // have flux in GeV^-1 s^-1 m^-2 sr^-1
}

float icecube_HESE_flux(float E) {
    float E_over_100_TeV = E / 100000.0;
    float f = 2.3 * std::pow(10, -18) * std::pow(E_over_100_TeV, -2.5);  // flux  from IceCube HES\E combined analysis in GeV^-1 s^-1 cm^-2 sr^-1
    return 0.5 * 10000.0 * f;                                            // have flux in GeV^-1 s^-1 m^-2 sr^-1
}

float icecube_flux_diff2019(float E) {
    float E_over_100_TeV = E / 100000.0;
    float f = 1.44 * std::pow(10, -18) * std::pow(E_over_100_TeV, -2.28);  // flux  from IceCube diffuse astrophysics numuCC analysis in GeV^-1 s^-1 cm^-2 sr^-1
    return 0.5 * 10000.0 * f;                                              // have flux in GeV^-1 s^-1 m^-2 sr^-1
}

float Antares_flux(float E) {
    float E_over_100_TeV = E / 100000.0;
    float f = pow(10, -18) * pow(E_over_100_TeV, -2.0);  // flux  from IceCube diffuse astrophysics numuCC analysis in GeV^-1 s^-1 cm^-2 sr^-1
    return 0.5 * 10000.0 * f;                            // have flux in GeV^-1 s^-1 m^-2 sr^-1
}

float new_flux(float E) {//for tests (or other)
  float E_over_100_TeV = E / 100000.0;
  float f = 1.44 * std::pow(10, -18) * std::pow(E_over_100_TeV, -2.37);  // flux  from IceCube diffuse astrophysics numuCC analysis in GeV^-1 s^-1 cm^-2 sr^-1
  return 0.5 * 10000.0 * f;                                              // have flux in GeV^-1 s^-1 m^-2 sr^-1
}
//--------------------------------------------------------------------------------------------------

// Organize Hits in maps
std::map<int, std::map<int, std::map<int, std::vector<Hit>>>> makeHitsMap(const std::vector<Hit>& hits, Det& det) {
    const Hit* hit;
    Pmt pmt;
    Dom dom;

    int dom_id;
    int channel_id;
    int line;
    int floor;

    std::map<int, std::map<int, std::map<int, std::vector<Hit>>>> result_map;

    for (unsigned i = 0; i < hits.size(); ++i) {
        hit = &hits[i];
        dom_id = hit->dom_id;
        channel_id = hit->channel_id;

        pmt = det.get_pmt(dom_id, channel_id);
        dom = det.get_dom(pmt);

        line = dom.line_id;
        floor = dom.floor_id;

        if (hit->trig > 0) {
            result_map[line][floor][channel_id].push_back(*hit);
        }
    }
    return result_map;
}

/*
---------
    Find Distance between Muon's Vertex and Photon's Production Point (lambda)
        !Caution!
            Will fail if not well reconstructed
--------
*/
/*
std::vector<float> find_lambda(const Vec& muPos, const Vec& muDir, const Vec& omPos, const float& theta_cherenkov) {
    // Initialize variables
    float theta_muTrackDir[2] = { 0, 0 };
    float cos_theta_cherenkov = std::cos(theta_cherenkov * TMath::DegToRad());
    const float c_water = 0.2222495796574987;

    float alpha;
    float beta;
    float gamma;
    float discriminant;
    Vec muPos_t[2];
    float solution[2];
    float track[2];
    float lambda;
    Vec muPos_t_final;
    float om_productionPoint_distance;
    float time;

    // Prepare equation coefficients
    alpha = (muDir.dot(muDir)).dot((muDir.dot(muDir - std::pow(cos_theta_cherenkov, 2))));
    beta = -2 * (muDir.dot((omPos - muPos))).dot((muDir * muDir - std::pow(cos_theta_cherenkov, 2)));
    gamma = muDir.dot((omPos - muPos).dot(muDir.dot(omPos - muPos))) - (omPos - muPos).dot((omPos - muPos)) * std::pow(cos_theta_cherenkov, 2);

    // Start solving
    discriminant = beta * beta - 4 * alpha * gamma;
    solution[0] = (-beta + std::sqrt(discriminant)) / (2 * alpha);
    solution[1] = (-beta - std::sqrt(discriminant)) / (2 * alpha);

    muPos_t[0] = muPos + muDir * solution[0];
    muPos_t[1] = muPos + muDir * solution[1];

    track[0] = (omPos - muPos_t[0]).len();
    track[1] = (omPos - muPos_t[1]).len();

    theta_muTrackDir[0] = std::acos(muDir.dot(omPos - muPos_t[0]) / track[0] * TMath::RadToDeg());
    theta_muTrackDir[1] = std::acos(muDir.dot(omPos - muPos_t[1]) / track[1] * TMath::RadToDeg());

    if ((solution[0] > 0 && solution[0] < 0.1) || (solution[1] > 0 && solution[1] < 0.1))
        lambda = 0;

    if (theta_muTrackDir[0] <= theta_muTrackDir[1]) {
        lambda = solution[0];
        muPos_t_final = muPos_t[0];
    }
    else {
        lambda = solution[1];
        muPos_t_final = muPos_t[1];
    }

    om_productionPoint_distance = (muPos_t_final - omPos).len();
    time = om_productionPoint_distance / c_water;

    // Return result
    std::vector<float> result;
    result.push_back(lambda);
    result.push_back(time);
    return result;
}
*/
/*
  // BDT-selection (KM3NeT ASTRO 2021 001) Alfonso Garcia, Rasa Muller, Aart Heijboer - June 2, 2021
  // Anna variables
  float downLikMax;       // Likelihood of best down-going solution from JGandalf’s prefit list.
  float upLikMax;         // Likelihood of up down-going solution from JGandalf’s prefit list.
  float cosZenith_min;    // Lowest reconstructed cosine of zenith angle from JGandalf’s prefit list.
  float cosZenith_max;    // Highest reconstructed cosine of zenith angle from JGandalf’s prefit list.
  float prefitN_1degree;  // Number of solutions from JGandalf’s prefit list with a < 1 with respect to best reconstructed track.
  float prefitN_down;     // Number of down-going solutions from JGandalf’s prefit list.
  float prefitN_up;       // Number of up-going solutions from JGandalf’s prefit list.

  float cascNHits;        // Number of hits that fulfill the Cerenkov(cascade) hypothesis using JGandalf’s best solution.
  float cascNHits_100m;   // Number of hits that fulfill the Cerenkov(cascade) hypothesis using JGandalf’s best solution less than 100 m away from the reconstructed vertex.
  float cascDistanceMin;  // Closest hit to the reconstructed vertex that fulfill the Cerenkov(cascade) hypothesis using JGandalf’s best solution.
  float cascDistanceMax;  // Closest hit to the reconstructed vertex that fulfill the Cerenkov(cascade) hypothesis using JGandalf’s best solution.

  float nHits;        // Number of hits that fulfill the Cerenkov(muon) hypothesis using JGandalf’s best solution.
  float nHits_100m;   // Number of hits that fulfill the Cerenkov(muon) hypothesis using JGandalf’s best solution less than 100 m away from the reconstructed vertex.
  float distanceMin;  // Closest hit to the reconstructed vertex that fulfill the Cerenkov(muon) hypothesis using JGandalf’s best solution.
  float distanceMax;  // Closest hit to the reconstructed vertex that fulfill the Cerenkov(muon) hypothesis using JGandalf’s best solution.

  float cascmuNHits;        // Number of hits that fulfill the Cerenkov(cascademu) hypothesis using JGandalf’s best solution.
  float cascmuNHits_100m;   // Number of hits that fulfill the Cerenkov(cascademu) hypothesis using JGandalf’s best solution less than 100 m away from the reconstructed vertex.
  float cascmuDistanceMin;  // Closest hit to the reconstructed vertex that fulfill the Cerenkov(cascademu) hypothesis using JGandalf’s best solution.
  float cascmuDistanceMax;  // Closest hit to the reconstructed vertex that fulfill the Cerenkov(cascademu) hypothesis using JGandalf’s best solution.
  float thetaDiff;          // Difference between thetamax and thetamin
  */
  //=============================================================================================================================
