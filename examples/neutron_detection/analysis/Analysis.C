/*
 * Analysis.C
 *
 *  Created on: Jun 15, 2019
 *      Author: mbenoit
 */
#pragma includepath "../../../src:../../../src/objects";

#include <Math/DisplacementVector2D.h>
#include <Math/Vector2D.h>
#include <Math/Vector3D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TH1D.h>
#include <TH1I.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TSpectrum.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstring>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
// FIXME: these includes should be absolute and provided with installation?
#include "../../../src/objects/MCParticle.hpp"
#include "../../../src/objects/Pixel.hpp"
#include "../../../src/objects/PixelCharge.hpp"
#include "../../../src/objects/PixelHit.hpp"
#include "../../../src/objects/PropagatedCharge.hpp"
#include "../../../src/modules/DetectorHistogrammer/Cluster.hpp"

#ifdef __MAKECINT__
#pragma link C++ class vector < allpix::PixelHit* > +;
#pragma link C++ class vector < allpix::MCTrack* > +;
#pragma link C++ class vector < allpix::MCParticle* > +;
#endif

std::string format(const char* formatchar, ...) {
    va_list args;
    va_start(args, formatchar);
    std::string format(formatchar);
    size_t len = std::vsnprintf(NULL, 0, format.c_str(), args);
    va_end(args);
    std::vector<char> vec(len + 1);
    va_start(args, formatchar);
    std::vsnprintf(&vec[0], len + 1, format.c_str(), args);
    va_end(args);
    return &vec[0];
}

std::vector<allpix::Cluster> doClustering(std::vector<allpix::PixelHit*> hits) {
    std::vector<allpix::Cluster> clusters;
    std::map<const allpix::PixelHit*, bool> usedPixel;

    if(hits.size() < 1) {
        return clusters;
    }

    auto pixel_it = hits.begin();
    for(; pixel_it != hits.end(); pixel_it++) {
        const allpix::PixelHit* pixel_hit = (*pixel_it);

        // Check if the pixel has been used:
        if(usedPixel[pixel_hit]) {
            continue;
        }

        // Create new cluster
        allpix::Cluster cluster(pixel_hit);
        usedPixel[pixel_hit] = true;
        // cout << "Creating new cluster with seed: " << pixel_hit->getPixel().getIndex().X()<< " " <<
        // pixel_hit->getPixel().getIndex().Y() <<std::endl;

        auto touching = [&](const allpix::PixelHit* pixel) {
            auto pxi1 = pixel->getIndex();
            for(auto& cluster_pixel : cluster.getPixelHits()) {

                auto distance = [](unsigned int lhs, unsigned int rhs) { return (lhs > rhs ? lhs - rhs : rhs - lhs); };

                auto pxi2 = cluster_pixel->getIndex();
                if(distance(pxi1.x(), pxi2.x()) <= 1 && distance(pxi1.y(), pxi2.y()) <= 1) {
                    return true;
                }
            }
            return false;
        };

        // Keep adding pixels to the cluster:
        for(auto other_pixel = pixel_it + 1; other_pixel != hits.end(); other_pixel++) {
            const allpix::PixelHit* neighbor = (*other_pixel);

            // Check if neighbor has been used or if it touches the current cluster:
            if(usedPixel[neighbor] || !touching(neighbor)) {
                continue;
            }

            cluster.addPixelHit(neighbor);
                        cout << "Adding pixel: " << neighbor->getPixel().getIndex().X() << " " <<
                        neighbor->getPixel().getIndex().Y() << " " << neighbor->getSignal()
                        << std::endl;
            usedPixel[neighbor] = true;
            other_pixel = pixel_it;
        }
        clusters.push_back(cluster);
    }
    return clusters;
}

const allpix::PixelHit* FindSeedPixel(allpix::Cluster& cluster) {

    int temp_tot = 0;
    int max_index = 0;
    int idx = 0;
    const allpix::PixelHit* seed_pixel = new allpix::PixelHit();

    for(auto pixel : cluster.getPixelHits()) {
        if(pixel->getSignal() > temp_tot) {
            temp_tot = pixel->getSignal();
            max_index = idx;
            seed_pixel = pixel;
        }
        idx++;
    }

    std::cout << format("SEED PIXEL X: %i Y: %i", seed_pixel->getIndex().x(), seed_pixel->getIndex().y()) << std::endl;
    return seed_pixel;
}

void TrimCluster(allpix::Cluster& cluster, int N = 10) {

    auto seed_pixel = FindSeedPixel(cluster);
    auto new_cluster = allpix::Cluster(seed_pixel);
    auto seed_pixel_idx = seed_pixel->getIndex();

    for(auto pixel : cluster.getPixelHits()) {
        auto pix_addr = pixel->getIndex();
        if((pix_addr.x() == seed_pixel_idx.x()) && pix_addr.y() == seed_pixel_idx.y()) {

        } else if(abs(int((pix_addr.x() - seed_pixel_idx.x()))) <= N && abs(int(pix_addr.y() - seed_pixel_idx.y())) <= N) {
            new_cluster.addPixelHit(pixel);
        } else {
            std::cout << format("Trimming pixel X: %i Y:%i", pix_addr.x(), pix_addr.y()) << std::endl;
        }
    }

    cluster = new_cluster;
}

void PrintCluster(allpix::Cluster& cluster) {

    std::cout << "--------- CLUSTER begin -----------" << std::endl;
    for(auto pixel : cluster.getPixelHits()) {
        auto pix_addr = pixel->getIndex();
        std::cout << format("X : %i Y : %i ", pix_addr.x(), pix_addr.y()) << std::endl;
        ;
    }
    std::cout << "--------- CLUSTER end   -----------" << std::endl;
}

void AnalysisB10(string input_file_folder,
              string input_file_name,
              string detector_name,
              double pitch = 0.055,
              int mat_size = 256,
              int maxTOT = 2e9) {

    // gROOT->ProcessLine(".L /eda/allpix2/allpix-squared/lib/libAllpixObjects.so");
    // gROOT->ProcessLine(".L /eda/allpix2/allpix-squared/lib/libAllpixModuleDetectorHistogrammer.so");

    // std::string input_file_folder = format("examples/UCN_Detection/test/");
    cout << "Opening file " << input_file_name << std::endl;
    auto file = TFile::Open((input_file_folder + "/" +input_file_name).c_str());
    cout << "file opened" << std::endl;

    TTree* pixel_hit_tree = static_cast<TTree*>(file->Get("PixelHit"));
    TTree* mctrack_tree = static_cast<TTree*>(file->Get("MCTrack"));
    TTree* mcparticles_tree = static_cast<TTree*>(file->Get("MCParticle"));
    int nevents = mcparticles_tree->GetEntries();
    cout << "got the trees" << std::endl;

    cout << "branches picked" << std::endl;

    std::vector<allpix::PixelHit*> input_hits_bot;
    std::vector<allpix::MCTrack*> input_tracks;
    std::vector<allpix::MCParticle*> input_particles_bot;
    cout << "created vectors" << std::endl;

    pixel_hit_tree->FindBranch(detector_name.c_str())->SetObject(&input_hits_bot);
    mctrack_tree->FindBranch("global")->SetObject(&input_tracks);
    mcparticles_tree->FindBranch(detector_name.c_str())->SetObject(&input_particles_bot);
    cout << "Assigned object vectors to branches" << std::endl;

    std::string output_file = input_file_folder +"/"+ detector_name + "_output_plots.root";
    TFile* outfile = new TFile(output_file.c_str(), "recreate");
    outfile->cd();

    TH1D* resx_bot = new TH1D("resx_bot", "Residual X ", 200, -0.05, 0.05);
    TH1D* resy_bot = new TH1D("resy_bot", "Residual Y ", 200, -0.05, 0.05);


    TH2D* hitmap_bot = new TH2D("hitmap_bot", "Cluster position", 500, 0, 2*7.04, 500, 0, 2*7.04);

    TH1D* hitmapx_bot = new TH1D("hitmapx_bot", "Cluster position in X", 2500, 0, 2*7.04);

    TH1I* clu_tot_bot = new TH1I("clu_tot_bot", "Cluster TOT", 1000, 0, maxTOT);

    int n_alpha = 0;
    int n_lithium = 0;
    int n_conversion = 0;
    int n_coinc = 0;
    int n_SiCapture = 0;

    double resx_bot_coinc = 0;
    double resy_bot_coinc = 0;

    for(int i = 0; i < pixel_hit_tree->GetEntries(); ++i) {
        int totalTOT = 0;

        cout << "----------------------Event " << i << "---------------------" << std::endl;
        cout << "################### MCTruth #######################" << std::endl;
        pixel_hit_tree->GetEvent(i);
        mctrack_tree->GetEvent(i);
        mcparticles_tree->GetEvent(i);
        bool has_int_bot = false;
        bool has_int_top = false;
        for(auto& p : input_tracks) {
            int id = p->getParticleID();
            auto start = p->getStartPoint();
            auto end = p->getEndPoint();
            auto energy = p->getKineticEnergyInitial() - p->getKineticEnergyFinal();
            cout << format("PID: %i ", id) << endl;
            cout << format("(%f %f %f ) (%f %f %f ) dedx : %f MeV",
                           start.X(),
                           start.Y(),
                           start.Z(),
                           end.X(),
                           end.Y(),
                           end.Z(),
                           energy)
                 << std::endl;
            if(id == 1000140290 || id == 1000140300) {
                n_SiCapture++;
            }
        }
        cout << "################### Bottom sensor ################" << std::endl;

        auto clusters = doClustering(input_hits_bot);
        cout << format("Found %i clusters, with %i pixels", clusters.size(), input_hits_bot.size()) << std::endl;

        for(auto& p : input_particles_bot) {
            cout << format("PID: %i ", p->getParticleID()) << endl;
        }

        auto comp_clu = [](allpix::Cluster& a, allpix::Cluster& b) { return a.getCharge() > b.getCharge(); };
        std::sort(clusters.begin(), clusters.end(), comp_clu);
        for(auto& clu : clusters) {

                   	std::cout << "old cluster" << std::endl;
                   	PrintCluster(clu);
                   	TrimCluster(clu);
            
                   	std::cout << "new cluster" << std::endl;
                   	PrintCluster(clu);

            double posX = clu.getPosition().X() ;
            double posY = clu.getPosition().Y();

            int cluTOT = clu.getCharge();
            cout << format(" X : %f Y : %f  %i, size %i", clu.getPosition().X(), clu.getPosition().Y(), cluTOT, clu.getSize()) << std::endl;

            for(auto p : clu.getMCParticles()) {

                if(p == NULL) {
                    continue;
                }
                int id = p->getParticleID();

                if(((id == 1000030070) || (id == 1000020040)) && has_int_bot == false) {
                    // double mcX = p->getGlobalStartPoint().X();
                    // double mcY = p->getGlobalStartPoint().Y();
                    double mcX = p->getLocalStartPoint().X();
                    double mcY = p->getLocalStartPoint().Y();                   
                    cout << "Neutron interaction detected" << std::endl;
                    if(cluTOT > 0 && clu.getSize() > 1) {
                        resx_bot->Fill(posX - mcX );
                        resy_bot->Fill(posY - mcY );
                        cout << "posX : " << posX << " mc truth X : " << mcX << endl;
                        cout << "posY : " << posY << " mc truth Y : " << mcY << endl;

                    }
                    resx_bot_coinc = posX - mcX + pitch / 2;
                    resy_bot_coinc = posY - mcY + pitch / 2;

                    if(clu.getSize() > 3)
                        clu_tot_bot->Fill(cluTOT);
                    totalTOT += cluTOT;
                    hitmap_bot->Fill(posX, posY);
                    hitmapx_bot->Fill(posX);
                    n_conversion++;
                    has_int_bot = true;
                }
                if(id == 1000030070) {
                    n_lithium++;
                    cout << "Lithium ion cnt " << n_lithium << std::endl;
                }
                if(id == 1000020040) {
                    n_alpha++;
                    cout << "Alpha cnt  " << n_alpha << std::endl;
                }

                if(has_int_bot) {
                    break;
                }
            }
        }
    }



    double efficiency = 100 * double(n_conversion) / nevents;
    cout << format("[Interaction report] %i Lithium, %i Alpha, %i coincidence, %i conversion detected, %i Si Capture, %f %% "
                   "efficiency ",
                   n_lithium,
                   n_alpha,
                   n_coinc,
                   n_conversion,
                   n_SiCapture,
                   efficiency)
         << std::endl;

    std::vector<TH1D*> resplots{resx_bot, resy_bot};
    for(auto& plot : resplots) {
        plot->SetLineWidth(2);
        plot->SetLineColor(kBlue);
        // plot->SetFillColor(kRed);
        // plot->SetFillStyle(3353);
        plot->GetXaxis()->SetTitle("X_{neutron} - X_{reco} [mm]");
        gStyle->SetOptFit(0011);
    }

    TCanvas* can = new TCanvas();
    can->Draw();
    can->SetWindowSize(1400, 800);
    can->Divide(2, 1);
    can->cd(1);
    resx_bot->Draw();
    can->cd(2);
    resy_bot->Draw();
    //    can->cd(3);
    //    resx_top->Draw();
    //    can->cd(4);
    //    resy_top->Draw();
    // auto fitx = resx_bot->Fit("gaus", "S", "", -0.01, 0.01);
    // auto fity = resy_bot->Fit("gaus", "S", "", -0.01, 0.01);

    // cout << format("resolution X : %f um Y: %f um", fitx->Parameter(2) * 1000, fity->Parameter(2) * 1000) << std::endl;

    can->Print((input_file_folder + "ResidualXY_Bottom.png").c_str());


    TCanvas* can2 = new TCanvas();
    can2->Draw();
    //can2->Divide(2);
    //can2->cd(1);
    hitmap_bot->Draw("colz");
    //can2->cd(2);
   // hitmap_top->Draw("colz");

    TCanvas* can3 = new TCanvas();
    can3->Draw();
    can3->SetWindowSize(1400, 800);

    // can3->Divide(2);
    // can3->cd(1);
    clu_tot_bot->SetFillColor(kRed);
    clu_tot_bot->SetLineColor(kRed);
    clu_tot_bot->SetFillStyle(3353);
    clu_tot_bot->SetLineWidth(2);
    clu_tot_bot->GetXaxis()->SetRangeUser(0, maxTOT);
    clu_tot_bot->GetXaxis()->SetTitle("TOT [A. U.]");
    clu_tot_bot->SetStats(0);
    clu_tot_bot->Draw();
    // can3->cd(2);
    // clu_tot_top->Draw();
    can3->Print((input_file_folder + "clusterTOT_Bottom.png").c_str());

    // TCanvas* can6 = new TCanvas();
    // can6->Draw();
    // can6->SetWindowSize(1400, 800);

    // // can3->Divide(2);
    // // can3->cd(1);
    // clu_tot_top->SetFillColor(kRed);
    // clu_tot_top->SetLineColor(kRed);
    // clu_tot_top->SetFillStyle(3353);
    // clu_tot_top->SetLineWidth(2);
    // clu_tot_top->GetXaxis()->SetRangeUser(0, maxTOT);
    // clu_tot_top->GetXaxis()->SetTitle("TOT [A. U.]");
    // clu_tot_top->SetStats(0);
    // clu_tot_top->Draw();
    // // can3->cd(2);
    // // clu_tot_top->Draw();
    // can6->Print((input_file_folder + "clusterTOT_Top.png").c_str());

    // TCanvas* can7 = new TCanvas();
    // can7->Draw();
    // can7->SetWindowSize(1400, 800);

    // // can3->Divide(2);
    // // can3->cd(1);
    // clu_tot_com->SetFillColor(kRed);
    // clu_tot_com->SetLineColor(kRed);
    // clu_tot_com->SetFillStyle(3353);
    // clu_tot_com->SetLineWidth(2);
    // clu_tot_com->GetXaxis()->SetRangeUser(0, 2 * maxTOT);
    // clu_tot_com->GetXaxis()->SetTitle("TOT [A. U.]");
    // clu_tot_com->SetStats(0);
    // clu_tot_com->Draw();
    // // can3->cd(2);
    // // clu_tot_top->Draw();
    // can7->Print((input_file_folder + "clusterTOT_Combined.png").c_str());

    TCanvas* can8 = new TCanvas();
    can8->Draw();
    can8->SetWindowSize(1400, 800);
    hitmapx_bot->Draw();
    can8->Print((input_file_folder + "hitmap_X_bottom.png").c_str());

    // TCanvas* can9 = new TCanvas();
    // can9->Draw();
    // can9->SetWindowSize(1400, 800);
    // hitmapx_top->Draw();
    // can9->Print((input_file_folder + "hitmap_X_top.png").c_str());

    resx_bot->Write();
    resy_bot->Write();
    //resx_top->Write();
    //resy_top->Write();
    //resx_com->Write();
    //resy_com->Write();
    hitmap_bot->Write();
    //hitmap_top->Write();
    hitmapx_bot->Write();
   // hitmapx_top->Write();
    //clu_tot_top->Write();
    clu_tot_bot->Write();
    //clu_tot_com->Write();

    // outfile->Close();
    //    file->Close();

    //    delete mcparticles_tree;
    //    delete mctrack_tree;
    //    delete file;
    //    delete outfile;
    //    delete pixel_hit_tree;
    //    delete resx_bot;
    //	delete resy_bot;
    //	delete resx_top;
    //	delete resy_top;
    //	delete resx_com;
    //	delete resy_com;
    //	delete clu_tot_top;
    //	delete clu_tot_bot;
    //	delete hitmap_bot;
    //	delete hitmap_top;
}


void AnalysisPE(string input_file_folder,
              string input_file_name,
              string detector_name,
              double pitch = 0.055,
              int mat_size = 256,
              int maxTOT = 2e9) {

    // gROOT->ProcessLine(".L /eda/allpix2/allpix-squared/lib/libAllpixObjects.so");
    // gROOT->ProcessLine(".L /eda/allpix2/allpix-squared/lib/libAllpixModuleDetectorHistogrammer.so");

    // std::string input_file_folder = format("examples/UCN_Detection/test/");
    cout << "Opening file " << input_file_name << std::endl;
    auto file = TFile::Open((input_file_folder + "/" +input_file_name).c_str());
    cout << "file opened" << std::endl;

    TTree* pixel_hit_tree = static_cast<TTree*>(file->Get("PixelHit"));
    TTree* mctrack_tree = static_cast<TTree*>(file->Get("MCTrack"));
    TTree* mcparticles_tree = static_cast<TTree*>(file->Get("MCParticle"));
    int nevents = mcparticles_tree->GetEntries();
    cout << "got the trees" << std::endl;

    cout << "branches picked" << std::endl;

    std::vector<allpix::PixelHit*> input_hits_bot;
    std::vector<allpix::MCTrack*> input_tracks;
    std::vector<allpix::MCParticle*> input_particles_bot;
    cout << "created vectors" << std::endl;

    pixel_hit_tree->FindBranch(detector_name.c_str())->SetObject(&input_hits_bot);
    mctrack_tree->FindBranch("global")->SetObject(&input_tracks);
    mcparticles_tree->FindBranch(detector_name.c_str())->SetObject(&input_particles_bot);
    cout << "Assigned object vectors to branches" << std::endl;

    std::string output_file = input_file_folder +"/"+ detector_name + "_output_plots.root";
    TFile* outfile = new TFile(output_file.c_str(), "recreate");
    outfile->cd();

    TH1D* resx_bot = new TH1D("resx_bot", "Residual X ", 200, -0.05, 0.05);
    TH1D* resy_bot = new TH1D("resy_bot", "Residual Y ", 200, -0.05, 0.05);


    TH2D* hitmap_bot = new TH2D("hitmap_bot", "Cluster position", 500, 0, 2*7.04, 500, 0, 2*7.04);

    TH1D* hitmapx_bot = new TH1D("hitmapx_bot", "Cluster position in X", 2500, 0, 2*7.04);

    TH1I* clu_tot_bot = new TH1I("clu_tot_bot", "Cluster TOT", 1000, 0, maxTOT);

    int n_alpha = 0;
    int n_lithium = 0;
    int n_conversion = 0;
    int n_coinc = 0;
    int n_SiCapture = 0;

    double resx_bot_coinc = 0;
    double resy_bot_coinc = 0;

    for(int i = 0; i < pixel_hit_tree->GetEntries(); ++i) {
        int totalTOT = 0;

        cout << "----------------------Event " << i << "---------------------" << std::endl;
        cout << "################### MCTruth #######################" << std::endl;
        pixel_hit_tree->GetEvent(i);
        mctrack_tree->GetEvent(i);
        mcparticles_tree->GetEvent(i);
        bool has_int_bot = false;
        bool has_int_top = false;
        for(auto& p : input_tracks) {
            int id = p->getParticleID();
            auto start = p->getStartPoint();
            auto end = p->getEndPoint();
            auto energy = p->getKineticEnergyInitial() - p->getKineticEnergyFinal();
            cout << format("PID: %i ", id) << endl;
            cout << format("(%f %f %f ) (%f %f %f ) dedx : %f MeV",
                           start.X(),
                           start.Y(),
                           start.Z(),
                           end.X(),
                           end.Y(),
                           end.Z(),
                           energy)
                 << std::endl;
            if(id == 1000140290 || id == 1000140300) {
                n_SiCapture++;
            }
        }
        cout << "################### Bottom sensor ################" << std::endl;

        auto clusters = doClustering(input_hits_bot);
        cout << format("Found %i clusters, with %i pixels", clusters.size(), input_hits_bot.size()) << std::endl;

        for(auto& p : input_particles_bot) {
            cout << format("PID: %i ", p->getParticleID()) << endl;
        }

        auto comp_clu = [](allpix::Cluster& a, allpix::Cluster& b) { return a.getCharge() > b.getCharge(); };
        std::sort(clusters.begin(), clusters.end(), comp_clu);
        for(auto& clu : clusters) {

                   	std::cout << "old cluster" << std::endl;
                   	PrintCluster(clu);
                   	TrimCluster(clu);
            
                   	std::cout << "new cluster" << std::endl;
                   	PrintCluster(clu);

            double posX = clu.getPosition().X() ;
            double posY = clu.getPosition().Y();

            int cluTOT = clu.getCharge();
            cout << format(" X : %f Y : %f  %i, size %i", clu.getPosition().X(), clu.getPosition().Y(), cluTOT, clu.getSize()) << std::endl;

            for(auto p : clu.getMCParticles()) {

                if(p == NULL) {
                    continue;
                }
                int id = p->getParticleID();
                cout << "Particle ID : " << id << endl;
                if(id == 2212) {
                    // double mcX = p->getGlobalStartPoint().X();
                    // double mcY = p->getGlobalStartPoint().Y();
                    double mcX = p->getLocalStartPoint().X();
                    double mcY = p->getLocalStartPoint().Y();                   
                    cout << "Neutron interaction detected" << std::endl;
                    if(cluTOT > 0 && clu.getSize() > 1) {
                        resx_bot->Fill(posX - mcX );
                        resy_bot->Fill(posY - mcY );
                        cout << "posX : " << posX << " mc truth X : " << mcX << endl;
                        cout << "posY : " << posY << " mc truth Y : " << mcY << endl;

                    }
                    resx_bot_coinc = posX - mcX + pitch / 2;
                    resy_bot_coinc = posY - mcY + pitch / 2;

                    if(clu.getSize() > 3)
                        clu_tot_bot->Fill(cluTOT);
                    totalTOT += cluTOT;
                    hitmap_bot->Fill(posX, posY);
                    hitmapx_bot->Fill(posX);
                    n_conversion++;
                    has_int_bot = true;
                }
                if(id == 2212) {
                    n_lithium++;
                    cout << "Lithium ion cnt " << n_lithium << std::endl;
                }
                if(id == 1000020040) {
                    n_alpha++;
                    cout << "Alpha cnt  " << n_alpha << std::endl;
                }

                if(has_int_bot) {
                    break;
                }
            }
        }
    }



    double efficiency = 100 * double(n_conversion) / nevents;
    cout << format("[Interaction report] %i Protons, %i Alpha, %i coincidence, %i conversion detected, %i Si Capture, %f %% "
                   "efficiency ",
                   n_lithium,
                   n_alpha,
                   n_coinc,
                   n_conversion,
                   n_SiCapture,
                   efficiency)
         << std::endl;

    std::vector<TH1D*> resplots{resx_bot, resy_bot};
    for(auto& plot : resplots) {
        plot->SetLineWidth(2);
        plot->SetLineColor(kBlue);
        // plot->SetFillColor(kRed);
        // plot->SetFillStyle(3353);
        plot->GetXaxis()->SetTitle("X_{neutron} - X_{reco} [mm]");
        gStyle->SetOptFit(0011);
    }

    TCanvas* can = new TCanvas();
    can->Draw();
    can->SetWindowSize(1400, 800);
    can->Divide(2, 1);
    can->cd(1);
    resx_bot->Draw();
    can->cd(2);
    resy_bot->Draw();
    //    can->cd(3);
    //    resx_top->Draw();
    //    can->cd(4);
    //    resy_top->Draw();
    // auto fitx = resx_bot->Fit("gaus", "S", "", -0.01, 0.01);
    // auto fity = resy_bot->Fit("gaus", "S", "", -0.01, 0.01);

    // cout << format("resolution X : %f um Y: %f um", fitx->Parameter(2) * 1000, fity->Parameter(2) * 1000) << std::endl;

    can->Print((input_file_folder + "ResidualXY_Bottom.png").c_str());


    TCanvas* can2 = new TCanvas();
    can2->Draw();
    //can2->Divide(2);
    //can2->cd(1);
    hitmap_bot->Draw("colz");
    //can2->cd(2);
   // hitmap_top->Draw("colz");

    TCanvas* can3 = new TCanvas();
    can3->Draw();
    can3->SetWindowSize(1400, 800);

    // can3->Divide(2);
    // can3->cd(1);
    clu_tot_bot->SetFillColor(kRed);
    clu_tot_bot->SetLineColor(kRed);
    clu_tot_bot->SetFillStyle(3353);
    clu_tot_bot->SetLineWidth(2);
    clu_tot_bot->GetXaxis()->SetRangeUser(0, maxTOT);
    clu_tot_bot->GetXaxis()->SetTitle("TOT [A. U.]");
    clu_tot_bot->SetStats(0);
    clu_tot_bot->Draw();
    // can3->cd(2);
    // clu_tot_top->Draw();
    can3->Print((input_file_folder + "clusterTOT_Bottom.png").c_str());

    // TCanvas* can6 = new TCanvas();
    // can6->Draw();
    // can6->SetWindowSize(1400, 800);

    // // can3->Divide(2);
    // // can3->cd(1);
    // clu_tot_top->SetFillColor(kRed);
    // clu_tot_top->SetLineColor(kRed);
    // clu_tot_top->SetFillStyle(3353);
    // clu_tot_top->SetLineWidth(2);
    // clu_tot_top->GetXaxis()->SetRangeUser(0, maxTOT);
    // clu_tot_top->GetXaxis()->SetTitle("TOT [A. U.]");
    // clu_tot_top->SetStats(0);
    // clu_tot_top->Draw();
    // // can3->cd(2);
    // // clu_tot_top->Draw();
    // can6->Print((input_file_folder + "clusterTOT_Top.png").c_str());

    // TCanvas* can7 = new TCanvas();
    // can7->Draw();
    // can7->SetWindowSize(1400, 800);

    // // can3->Divide(2);
    // // can3->cd(1);
    // clu_tot_com->SetFillColor(kRed);
    // clu_tot_com->SetLineColor(kRed);
    // clu_tot_com->SetFillStyle(3353);
    // clu_tot_com->SetLineWidth(2);
    // clu_tot_com->GetXaxis()->SetRangeUser(0, 2 * maxTOT);
    // clu_tot_com->GetXaxis()->SetTitle("TOT [A. U.]");
    // clu_tot_com->SetStats(0);
    // clu_tot_com->Draw();
    // // can3->cd(2);
    // // clu_tot_top->Draw();
    // can7->Print((input_file_folder + "clusterTOT_Combined.png").c_str());

    TCanvas* can8 = new TCanvas();
    can8->Draw();
    can8->SetWindowSize(1400, 800);
    hitmapx_bot->Draw();
    can8->Print((input_file_folder + "hitmap_X_bottom.png").c_str());

    // TCanvas* can9 = new TCanvas();
    // can9->Draw();
    // can9->SetWindowSize(1400, 800);
    // hitmapx_top->Draw();
    // can9->Print((input_file_folder + "hitmap_X_top.png").c_str());

    resx_bot->Write();
    resy_bot->Write();
    //resx_top->Write();
    //resy_top->Write();
    //resx_com->Write();
    //resy_com->Write();
    hitmap_bot->Write();
    //hitmap_top->Write();
    hitmapx_bot->Write();
   // hitmapx_top->Write();
    //clu_tot_top->Write();
    clu_tot_bot->Write();
    //clu_tot_com->Write();

    // outfile->Close();
    //    file->Close();

    //    delete mcparticles_tree;
    //    delete mctrack_tree;
    //    delete file;
    //    delete outfile;
    //    delete pixel_hit_tree;
    //    delete resx_bot;
    //	delete resy_bot;
    //	delete resx_top;
    //	delete resy_top;
    //	delete resx_com;
    //	delete resy_com;
    //	delete clu_tot_top;
    //	delete clu_tot_bot;
    //	delete hitmap_bot;
    //	delete hitmap_top;
}