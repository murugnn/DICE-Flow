#include "engine.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <unordered_map>
#include <string>

// ============================================================================
// GRAPH LOADER
// ============================================================================
template<typename WeightType>
class GraphLoader {
public:
    static Graph<WeightType> load_from_file(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }

        std::unordered_map<std::string, int> name_map;
        std::vector<std::tuple<std::string, std::string, WeightType>> edges;
        
        std::string line;
        // Skip header if present
        std::getline(file, line);
        // Check if first line is actually data (does not contain "source")
        if (line.find("source") == std::string::npos && 
            line.find("Source") == std::string::npos) {
            std::istringstream iss(line);
            std::string source, target;
            WeightType weight;
            if (iss >> source >> target >> weight) {
                edges.push_back({source, target, weight});
                if(name_map.find(source) == name_map.end()) name_map[source] = name_map.size();
                if(name_map.find(target) == name_map.end()) name_map[target] = name_map.size();
            }
        }

        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;
            
            std::istringstream iss(line);
            std::string source, target;
            WeightType weight;
            
            if (iss >> source >> target >> weight) {
                edges.push_back({source, target, weight});
                if(name_map.find(source) == name_map.end()) name_map[source] = name_map.size();
                if(name_map.find(target) == name_map.end()) name_map[target] = name_map.size();
            }
        }
        file.close();

        // Construct Graph
        Graph<WeightType> graph(name_map.size());
        graph.id_to_name.resize(name_map.size());
        
        for (auto& p : name_map) {
            graph.name_to_id[p.first] = p.second;
            graph.id_to_name[p.second] = p.first;
        }

        for (auto& [src, tgt, w] : edges) {
            graph.add_edge(name_map[src], name_map[tgt], w);
        }

        std::cout << "Loaded graph: " << graph.n_vertices << " nodes, " 
                  << edges.size() << " edges" << std::endl;

        return graph;
    }
};

// ============================================================================
// DICE ANALYZER (Bioinformatics Logic)
// ============================================================================
template<typename WeightType>
class DiCEAnalyzer {
private:
    Graph<WeightType> graph;
    std::string graph_name;

public:
    DiCEAnalyzer(const Graph<WeightType>& g, const std::string& name)
        : graph(g), graph_name(name) {}

    void compute_all_centralities(const std::string& output_file) {
        std::cout << "\n=== Computing Centralities for " << graph_name << " ===" << std::endl;
        
        // Uses the BaseRadixHeap inside CentralityCalculator for speed
        CentralityCalculator<WeightType> calc(graph);
        
        std::cout << "Computing betweenness centrality..." << std::endl;
        auto betweenness = calc.compute_betweenness();
        
        std::cout << "Computing eigenvector centrality..." << std::endl;
        auto eigenvector = calc.compute_eigenvector_weighted();

        // Write results
        std::ofstream out(output_file);
        out << "gene\tbetweenness\teigenvector\n";
        
        for (int i = 0; i < graph.n_vertices; i++) {
            out << graph.id_to_name[i] << "\t" 
                << std::fixed << std::setprecision(6)
                << betweenness[i] << "\t"
                << eigenvector[i] << "\n";
        }
        out.close();
        std::cout << "Results written to: " << output_file << std::endl;
    }

    void virtual_knockout_analysis(const std::string& gene_list_file,
                                   const std::string& output_file) {
        std::cout << "\n=== Virtual Knockout Analysis ===" << std::endl;

        std::vector<std::string> genes_to_test;
        std::ifstream gfile(gene_list_file);
        if (!gfile.is_open()) {
            std::cerr << "Warning: Cannot open knockout list " << gene_list_file << std::endl;
            return;
        }
        std::string gene;
        while (std::getline(gfile, gene)) {
            if (!gene.empty() && gene[0] != '#') {
                // Trim trailing whitespace
                gene.erase(gene.find_last_not_of(" \n\r\t")+1);
                genes_to_test.push_back(gene);
            }
        }
        gfile.close();

        std::cout << "Testing " << genes_to_test.size() << " genes" << std::endl;

        double baseline_aspl = compute_aspl(graph);
        std::cout << "Baseline ASPL: " << baseline_aspl << std::endl;

        std::ofstream out(output_file);
        out << "gene\tvitality_score\taspl_knockout\taspl_change\n";

        for (const auto& gene_name : genes_to_test) {
            if (graph.name_to_id.find(gene_name) == graph.name_to_id.end()) {
                continue;
            }

            int gene_id = graph.name_to_id[gene_name];
            Graph<WeightType> ko_graph = create_knockout_graph(gene_id);
            double ko_aspl = compute_aspl(ko_graph);
            double vitality = ko_aspl - baseline_aspl;

            out << gene_name << "\t" << vitality << "\t" << ko_aspl << "\t" << vitality << "\n";
            std::cout << "  " << gene_name << ": vitality = " << vitality << std::endl;
        }
        out.close();
    }

private:
    Graph<WeightType> create_knockout_graph(int knockout_id) {
        Graph<WeightType> ko_graph(graph.n_vertices);
        ko_graph.name_to_id = graph.name_to_id;
        ko_graph.id_to_name = graph.id_to_name;

        for (int u = 0; u < graph.n_vertices; u++) {
            if (u == knockout_id) continue;
            for (const auto& e : graph.adj[u]) {
                if (e.to != knockout_id) {
                    ko_graph.add_edge(u, e.to, e.weight);
                }
            }
        }
        return ko_graph;
    }

    double compute_aspl(const Graph<WeightType>& g) {
        ShortestPathEngine<WeightType> engine(g);
        double total_dist = 0.0;
        int count = 0;
        // Sample for speed
        int samples = std::min(100, g.n_vertices);
        
        for (int i = 0; i < samples; i++) {
            int source = (i * g.n_vertices) / samples;
            auto distances = engine.compute_sssp(source);
            for (int j = 0; j < g.n_vertices; j++) {
                if (j != source && distances[j] != std::numeric_limits<WeightType>::max()) {
                    total_dist += distances[j];
                    count++;
                }
            }
        }
        return count > 0 ? total_dist / count : 0.0;
    }
};

// ============================================================================
// MAIN ENTRY POINT
// ============================================================================
int main(int argc, char* argv[]) {
    // Usage: ./dice_analyzer <normal.txt> <tumor.txt> <output_dir> [mode] [gene_list]
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] 
                  << " <normal_net> <tumor_net> <out_dir> [mode: full|bench] [knockout_list]" << std::endl;
        return 1;
    }

    std::string normal_file = argv[1];
    std::string tumor_file = argv[2];
    std::string output_dir = argv[3];
    std::string mode = (argc > 4) ? argv[4] : "full";
    std::string gene_list = (argc > 5) ? argv[5] : "";

    try {
        // --------------------------------------------------------------------
        // BENCHMARK MODE (The "Thesis Proof")
        // --------------------------------------------------------------------
        if (mode == "bench") {
            std::cout << "Loading graph for Recursive Duan Benchmark..." << std::endl;
            // Load only one file for benchmarking speed
            auto graph = GraphLoader<long long>::load_from_file(normal_file);
            ShortestPathEngine<long long> engine(graph);

            std::cout << "Benchmarking: Duan 2025 (Recursive) vs Standard Dijkstra" << std::endl;
            
            double total_duan = 0;
            double total_dijk = 0;
            
            srand(time(0));
            // Run 10 Random SSSP Queries
            for(int i=0; i<10; i++) {
                int src = rand() % graph.n_vertices;
                int tgt = rand() % graph.n_vertices;

                // 1. Run Actual Recursive Duan
                auto m1 = engine.compute_duan_recursive(src, tgt);
                
                // 2. Run Standard Dijkstra
                auto m2 = engine.compute_dijkstra(src, tgt);
                
                // 3. Truth Check (Allow floating point epsilon for double conversion)
                if (std::abs(m1.path_cost - m2.path_cost) > 1e-6) {
                    std::cout << " [FAIL] MISMATCH! Src:" << src << " Tgt:" << tgt 
                              << " Duan:" << m1.path_cost 
                              << " Dijk:" << m2.path_cost << std::endl;
                } else {
                    std::cout << " [PASS] Match." << std::endl;
                }

                total_duan += m1.runtime_ms;
                total_dijk += m2.runtime_ms;
            }

            // Print special formatted string for Python script to catch
            std::cout << "BENCHMARK_RESULT:DuanRecursive:" << total_duan 
                      << ":Dijkstra:" << total_dijk << std::endl;
            return 0;
        }

        // --------------------------------------------------------------------
        // FULL ANALYSIS MODE (The "Biological Discovery")
        // --------------------------------------------------------------------
        std::cout << "Running Standard Pipeline (DiCE Analysis)..." << std::endl;
        
        std::cout << "Loading NORMAL network..." << std::endl;
        auto normal_graph = GraphLoader<long long>::load_from_file(normal_file);
        
        std::cout << "Loading TUMOR network..." << std::endl;
        auto tumor_graph = GraphLoader<long long>::load_from_file(tumor_file);

        // Analyze Normal
        DiCEAnalyzer<long long> normal_analyzer(normal_graph, "Normal");
        normal_analyzer.compute_all_centralities(output_dir + "/centrality_normal.txt");

        // Analyze Tumor
        DiCEAnalyzer<long long> tumor_analyzer(tumor_graph, "Tumor");
        tumor_analyzer.compute_all_centralities(output_dir + "/centrality_tumor.txt");

        // Optional: Knockout Analysis
        if (!gene_list.empty() && !tumor_graph.id_to_name.empty()) {
            tumor_analyzer.virtual_knockout_analysis(
                gene_list,
                output_dir + "/knockout_analysis.txt"
            );
        }

        std::cout << "\n=== Analysis Complete ===" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}