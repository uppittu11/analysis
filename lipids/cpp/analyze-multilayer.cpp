#include "/Users/shamappk/group-code/common.hpp"
#include <ctime>
int main() {
    Trajectory traj("traj.dcd", "dcd");
    std::cout << "Loaded traj.dcd with " << traj.n_frames() << " frames" << std::endl;
    system_info info("system-info.txt"); 
    system("mkdir -p data");
    std::ofstream o_area("data/area.txt"); 
    o_area << "# frame area(A^2)" << std::endl; 
    std::ofstream o_apl("data/apl.txt");
    o_apl << "# frame apl_top(A^2) apl_bottom(A^2)" << std::endl; 
    std::ofstream o_s2("data/s2.txt");
    o_s2 << "# frame s2_top s2_bottom" << std::endl;
    std::ofstream o_tilt("data/tilt.txt"); 
    o_tilt << "# frame avg_tilt_top std_tilt_top avg_tilt_bottom std_tilt_bottom" << std::endl; 
    std::ofstream o_lipids_per_layer("data/lipids-per-layer.txt");
    o_lipids_per_layer << "# frame n_top n_bottom" << std::endl;
    std::ofstream o_cer_angles("data/cer-angles.txt"); 
    o_cer_angles << "# frame outer-acute outer-obtuse inner-acute inner-obtuse" << std::endl; 
    std::vector<float> area;
    float bin_width = 1.0/6; 
    int normdir = 2;
    int latdir1 = 0, latdir2 = 1; 
    histogram_1D profile_hist(-traj.box.back().period[normdir]/2, 
            traj.box.back().period[normdir]/2, bin_width);
    histogram_1D angle_hist(0, 180, 1.0);
    std::map<std::string, histogram_1D> profiles;
    std::map<std::string, std::map<std::string, histogram_1D> > layer_hists;
    std::vector<std::string> layers; layers.push_back("top"); layers.push_back("bottom");
    const char *init_components[] = {"ucer2-tails", "ucer2-headgroups",
        "ucer2-unequal-tail-part", "ucer2-fa-tail-full", "ffac24-tail",
        "ffac24-headgroups", "ffac24", "ffac24-tail-unequal", 
        "ffac24-tail-equal", "ucer2", "all", "chol"}; 
    std::vector<std::string> components(std::begin(init_components), 
                std::end(init_components)); 
    for(auto const& component: components) {
        for(auto const& layer: layers) {
            layer_hists[layer][component] = profile_hist; 
        }
    }
    profiles["water"] = profile_hist;
    coord_t tnormal; tnormal.push_back(0); tnormal.push_back(0); tnormal.push_back(0); 
    coord_t bnormal; bnormal.push_back(0); bnormal.push_back(0); bnormal.push_back(0); 
    tnormal[normdir] = 1.0; bnormal[normdir] = -1.0; 
    float inner_zone = 4.0;  // if com ±4*6 Å considered "inner" leaflets
    for(int frame = 0; frame < traj.n_frames(); frame++) {
        clock_t begin_time = std::clock();
        std::vector<gbb> molecules = info.coordlist_to_gbbs(traj.xyz[frame]);
        for(auto &molecule : molecules) {
            misc_cg::set_number(molecule); 
            misc_cg::set_groups(molecule);
            molecule.load_mass("/Users/shamappk/lipid-prototypes/cg/"
                    + molecule.smolecule_id + "/masses.txt"); 
            molecule.wrap_com(traj.box[frame]);
        }
        coord_t lipids_com = calc_total_com(molecules, "water");
        std::map<std::string, std::vector<gbb> > tails_by_layer; 
        std::map<std::string, int> lipids_per_layer; 
        std::map<std::string, int> cer_angles; 
        cer_angles["inner-acute"] = 0; cer_angles["outer-acute"] = 0;
        cer_angles["inner-obtuse"] = 0; cer_angles["outer-obtuse"] = 0;
        lipids_per_layer["top"] = 0; lipids_per_layer["bottom"] = 0; 
        for(auto &molecule : molecules) {
            if(molecule.v_com[normdir] > lipids_com[normdir]) {
                molecule.descriptors["layer"] = "top"; 
            }
            else {
                molecule.descriptors["layer"] = "bottom"; 
            }
            if(abs(molecule.v_com[normdir]) < inner_zone) {
                molecule.descriptors["inorout"] = "inner";
            }
            else {
                molecule.descriptors["inorout"] = "outer";
            }
            std::string layer = molecule.descriptors["layer"];
            if(molecule.smolecule_id == "ucer2-mapping2") {
                tails_by_layer[layer].push_back(
                        make_sub_gbb_from_group(molecule, "fa_tail_equal")); 
                tails_by_layer[layer].push_back(
                        make_sub_gbb_from_group(molecule, "sph_tail")); 
                lipids_per_layer[layer]++; 
                add_groups_to_zhist(molecule, 
                        layer_hists[layer]["ucer2-headgroups"],
                        "headgroup", lipids_com[normdir], normdir);
                add_groups_to_zhist(molecule, 
                        layer_hists[layer]["ucer2-unequal-tail-part"],
                        "fa_tail_unequal", lipids_com[normdir], normdir); 
                add_groups_to_zhist(molecule, 
                        layer_hists[layer]["ucer2-fa-tail-full"],
                        "fa_tail_full", lipids_com[normdir], normdir); 
                add_groups_to_zhist(molecule, 
                        layer_hists[layer]["ucer2-tails"],
                        "all_tail_beads", lipids_com[normdir], normdir); 
                add_groups_to_zhist(molecule,
                        layer_hists[layer]["ucer2"], "all", lipids_com[normdir], normdir);
                add_groups_to_zhist(molecule, 
                        layer_hists[layer]["all"], "all", lipids_com[normdir], normdir); 
                double tails_angle = calc_angle(molecule.v_coord[3], molecule.v_coord[8],
                        molecule.v_coord[14]);
                std::string ac_ob; 
                if(tails_angle < 90.0) {
                    ac_ob="acute";
                }
                else {
                    ac_ob="obtuse";
                }
                cer_angles[molecule.descriptors["inorout"]+"-"+ac_ob]++;
            }
            if(molecule.smolecule_id == "chol") {
                tails_by_layer[layer].push_back(
                        make_sub_gbb_from_group(molecule, "all"));
                add_groups_to_zhist(molecule, layer_hists[layer]["chol"],
                        "all", lipids_com[normdir], normdir);
                add_groups_to_zhist(molecule, layer_hists[layer]["all"],
                        "all", lipids_com[normdir], normdir); 
                lipids_per_layer[layer]++;
            }
            if(molecule.smolecule_id == "c24ffa") {
                tails_by_layer[layer].push_back(
                        make_sub_gbb_from_group(molecule, "tail"));
                add_groups_to_zhist(molecule, layer_hists[layer]["ffac24-headgroups"],
                        "headgroup", lipids_com[normdir], normdir);
                add_groups_to_zhist(molecule, layer_hists[layer]["ffac24"],
                        "all", lipids_com[normdir], normdir);
                add_groups_to_zhist(molecule, layer_hists[layer]["ffac24-tail"],
                        "tail", lipids_com[normdir], normdir); 
                add_groups_to_zhist(molecule, layer_hists[layer]["ffac24-tail-equal"],
                        "tail_equal", lipids_com[normdir], normdir); 
                add_groups_to_zhist(molecule, layer_hists[layer]["ffac24-tail-unequal"],
                        "tail_unequal", lipids_com[normdir], normdir); 
                add_groups_to_zhist(molecule, layer_hists[layer]["all"],
                        "all", lipids_com[normdir], normdir); 
                lipids_per_layer[layer]++;
            }
            if(molecule.smolecule_id == "water") {
                add_groups_to_zhist(molecule, profiles["water"], "all",
                        lipids_com[normdir], normdir); 
            }
        }
        float area_A2 = traj.box[frame].period[latdir1] * traj.box[frame].period[latdir2] * 36;
        area.push_back(area_A2);
        o_area << frame << " " << area_A2 << std::endl;
        o_s2 << frame << " " << nematic_order_gbbs(tails_by_layer["top"]) 
            << "  " << nematic_order_gbbs(tails_by_layer["bottom"]) << std::endl; 
        coord_t top_tilt = calc_tilt_angles(tails_by_layer["top"], tnormal); 
        coord_t bottom_tilt = calc_tilt_angles(tails_by_layer["bottom"], bnormal);
        angle_hist.insert_list(top_tilt); angle_hist.insert_list(bottom_tilt);
        o_tilt << frame << " " << average(top_tilt) << " " 
            << stdev(top_tilt)/sqrt(tails_by_layer["top"].size())
            << " " << average(bottom_tilt) << " " 
            << stdev(bottom_tilt)/sqrt(tails_by_layer["bottom"].size()) 
            << std::endl;
        o_apl << frame << " " << area_A2/lipids_per_layer["top"] << " "
            << area_A2/lipids_per_layer["bottom"] << std::endl;
        o_lipids_per_layer << frame << " " << lipids_per_layer["top"] << " " 
            << lipids_per_layer["bottom"] << std::endl;
        o_cer_angles << frame 
            << " " << cer_angles["outer-acute"]
            << " " << cer_angles["outer-obtuse"]
            << " " << cer_angles["inner-acute"] 
            << " " << cer_angles["inner-obtuse"] 
            << std::endl;
        if(frame % 10 == 0) { 
            std::cout << "Analyzed " << frame << " frames. "
                << float(std::clock() - begin_time) / CLOCKS_PER_SEC
                << " seconds per frame." << std::endl; 
            begin_time = std::clock(); 
        }
    }
    double avg_bin_vol = average(area) * bin_width;
    for(auto it = profiles.begin(); it != profiles.end(); ++it) {
        it->second.normalize(traj.n_frames()*avg_bin_vol);
        std::string filename = "data/profile_" + it->first + ".txt";
        it->second.print(filename); 
    }
    for(auto layer_it = layer_hists.begin(); layer_it != layer_hists.end(); ++layer_it) {
        for(auto tail_it = layer_it->second.begin();
                tail_it != layer_it->second.end(); ++tail_it) {
            tail_it->second.normalize(traj.n_frames()*avg_bin_vol);
            std::string fn = "data/profile_" + layer_it->first + "-"
                + tail_it->first + ".txt";
            tail_it->second.print(fn);
        }
    }
    angle_hist.normalize(traj.n_frames()); 
    angle_hist.print("data/hist_angles.txt");
    return 0; 
}
