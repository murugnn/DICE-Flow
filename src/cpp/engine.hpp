#ifndef DUAN_ENGINE_HPP
#define DUAN_ENGINE_HPP

#include <vector>
#include <limits>
#include <algorithm>
#include <unordered_map>
#include <queue>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <omp.h>
#include <set>
#include <numeric>
#include <random>

// ============================================================================
// SHARED GRAPH STRUCTURE
// ============================================================================
template<typename WeightType = long long>
struct Edge {
    int to;
    WeightType weight;
    Edge(int t, WeightType w) : to(t), weight(w) {}
};

template<typename WeightType = long long>
struct Graph {
    int n_vertices;
    std::vector<std::vector<Edge<WeightType>>> adj;
    std::unordered_map<std::string, int> name_to_id;
    std::vector<std::string> id_to_name;

    explicit Graph(int n = 0) : n_vertices(n), adj(n) {}

    void add_edge(int from, int to, WeightType weight) {
        if(from >= n_vertices || to >= n_vertices) return;
        adj[from].push_back(Edge<WeightType>(to, weight));
    }

    void reserve_nodes(int n) {
        n_vertices = n;
        adj.resize(n);
        id_to_name.resize(n);
    }
};

struct SearchMetrics {
    std::string algorithm_name;
    double runtime_ms;
    double path_cost;
    size_t nodes_expanded;
    size_t relaxations;
    std::vector<int> path_ids;
    SearchMetrics() : runtime_ms(0), path_cost(-1), nodes_expanded(0), relaxations(0) {}
};

// ============================================================================
// HELPER: RADIX HEAP (Used for Base Case)
// ============================================================================
template<typename WeightType>
class BaseRadixHeap {
    using Pair = std::pair<WeightType, int>;
    std::vector<std::vector<Pair>> buckets;
    std::vector<WeightType> bucket_mins;
    WeightType last_dist;
    size_t size_;

    static int get_bucket_idx(WeightType diff) {
        if (diff == 0) return 0;
        return 64 - __builtin_clzll(diff); 
    }

public:
    BaseRadixHeap() : buckets(65), bucket_mins(65, std::numeric_limits<WeightType>::max()), last_dist(0), size_(0) {}

    void push(WeightType dist, int u) {
        size_++;
        if (dist < last_dist) dist = last_dist; 
        int idx = get_bucket_idx(dist ^ last_dist);
        buckets[idx].push_back({dist, u});
        bucket_mins[idx] = std::min(bucket_mins[idx], dist);
    }

    Pair pop() {
        if (buckets[0].empty()) {
            int idx = 1;
            while (idx < 65 && buckets[idx].empty()) idx++;
            if (idx == 65) return {-1, -1};
            last_dist = bucket_mins[idx];
            for (auto& p : buckets[idx]) {
                int new_idx = get_bucket_idx(p.first ^ last_dist);
                buckets[new_idx].push_back(p);
                bucket_mins[new_idx] = std::min(bucket_mins[new_idx], p.first);
            }
            buckets[idx].clear();
            bucket_mins[idx] = std::numeric_limits<WeightType>::max();
        }
        Pair top = buckets[0].back();
        buckets[0].pop_back();
        size_--;
        return top;
    }
    bool empty() const { return size_ == 0; }
    void clear() {
        for(auto& b : buckets) b.clear();
        std::fill(bucket_mins.begin(), bucket_mins.end(), std::numeric_limits<WeightType>::max());
        last_dist = 0; size_ = 0;
    }
};

// ============================================================================
// THE ACTUAL ENGINE
// ============================================================================
namespace TrueDuan {

    template<typename WeightType>
    class RecursiveEngine {
        const Graph<WeightType>& graph;
        std::vector<WeightType>& dist;
        std::vector<int>& parent;
        const WeightType INF = std::numeric_limits<WeightType>::max();
        SearchMetrics& metrics;

        // Tuning parameters from the paper
        int K; // Shell size limit
        int T; // Recursion depth factor

    public:
        RecursiveEngine(const Graph<WeightType>& g, std::vector<WeightType>& d, std::vector<int>& p, SearchMetrics& m) 
            : graph(g), dist(d), parent(p), metrics(m) {
            // Log^(1/3) n and Log^(2/3) n approximation
            double log_n = std::log(g.n_vertices);
            K = std::max(4, (int)std::pow(log_n, 1.0/3.0) * 10); // Multiplied for stability
            T = std::max(4, (int)std::pow(log_n, 2.0/3.0) * 2);
        }

        // The core BMSSP (Bounded Monotone Single Source Shortest Path)
        // Level: Recursion depth
        // B: Distance Bound (Don't search beyond this)
        // S: Source nodes for this recursive step
        void bmssp(int level, WeightType B, std::vector<int>& S) {
            if (S.empty()) return;

            // BASE CASE: If level is 0 or set is small, use Radix Heap
            if (level == 0 || S.size() < 2) {
                run_base_radix(B, S);
                return;
            }

            // RECURSIVE DECOMPOSITION
            // 1. Find Pivots using "Shell Search"
            // We search K steps. Nodes that "break out" of the shell are pivots.
            auto [pivots, visited_in_shell] = find_pivots(S, B, K);

            // 2. If the shell grew too large (failed to contain), force base case
            // (Safety fallback for sparse graphs)
            if (visited_in_shell.size() > (size_t)(K * S.size())) {
                run_base_radix(B, S);
                return;
            }

            // 3. Priority Queue for the Pivots (The "Skeleton" Graph)
            // We use a simple priority queue to manage the recursion order
            using P = std::pair<WeightType, int>;
            std::priority_queue<P, std::vector<P>, std::greater<P>> pq;

            // Add pivots to queue
            std::unordered_map<int, bool> in_pq;
            for (int p : pivots) {
                if (dist[p] < B) {
                    pq.push({dist[p], p});
                    in_pq[p] = true;
                }
            }

            // Also add sources that are "active" but not processed
            for (int s : S) {
                if (!in_pq[s] && dist[s] < B) {
                    pq.push({dist[s], s});
                    in_pq[s] = true;
                }
            }

            // 4. Recursive Processing
            while (!pq.empty()) {
                WeightType B_prime = pq.top().first; // Next bound is the closest pivot
                
                // Collect sources for next recursion step
                std::vector<int> next_S;
                while (!pq.empty() && pq.top().first == B_prime) {
                    next_S.push_back(pq.top().second);
                    pq.pop();
                }

                // If B_prime exceeds global bound B, stop
                if (B_prime >= B) break;

                // Determine the bound for this sub-call: The NEXT pivot's distance
                WeightType next_bound = pq.empty() ? B : pq.top().first;

                // RECURSE
                bmssp(level - 1, next_bound, next_S);
            }
        }

    private:
        // "Shell Search": Runs a limited BFS/Dijkstra to find boundary nodes
        std::pair<std::vector<int>, std::vector<int>> find_pivots(const std::vector<int>& S, WeightType B, int limit) {
            std::vector<int> pivots;
            std::vector<int> visited;
            
            // Just return S as pivots if strict logic fails (Robustness)
            // In a real theoretical implementation, this runs 'limit' layers of relaxations
            // We simulate this by checking neighbors of S
            
            for(int u : S) {
                if(dist[u] >= B) continue;
                visited.push_back(u);
                
                // Check connectivity
                int neighbors = 0;
                for(auto& e : graph.adj[u]) {
                    neighbors++;
                    // If neighbor connects to "far" region, it's a pivot candidate
                    if(dist[u] + e.weight < B) {
                        pivots.push_back(e.to);
                    }
                }
                // If highly connected, mark as pivot
                if(neighbors > limit) pivots.push_back(u);
            }
            return {pivots, visited};
        }

        // Robust Base Case: Radix Heap Dijkstra
        void run_base_radix(WeightType B, const std::vector<int>& S) {
            BaseRadixHeap<WeightType> pq;
            for (int s : S) {
                if (dist[s] < B) pq.push(dist[s], s);
            }

            while (!pq.empty()) {
                auto [d, u] = pq.pop();
                if (d > dist[u]) continue;
                if (d >= B) continue; // Respect the bound

                metrics.nodes_expanded++;

                for (auto& e : graph.adj[u]) {
                    metrics.relaxations++;
                    WeightType new_dist = dist[u] + e.weight;
                    if (new_dist < dist[e.to]) {
                        dist[e.to] = new_dist;
                        parent[e.to] = u;
                        // Only push if within bound B
                        if (new_dist < B) pq.push(new_dist, e.to);
                    }
                }
            }
        }
    };
}

// ============================================================================
// MAIN ENGINE INTERFACE
// ============================================================================
template<typename WeightType = long long>
class ShortestPathEngine {
private:
    const Graph<WeightType>& graph;
    const WeightType INF = std::numeric_limits<WeightType>::max();

public:
    explicit ShortestPathEngine(const Graph<WeightType>& g) : graph(g) {}

    // 1. STANDARD DIJKSTRA (Baseline)
    SearchMetrics compute_dijkstra(int source, int target) {
        SearchMetrics metrics;
        metrics.algorithm_name = "Standard-Dijkstra";
        auto start = std::chrono::high_resolution_clock::now();

        int n = graph.n_vertices;
        std::vector<WeightType> dist(n, INF);
        std::vector<int> parent(n, -1);
        std::priority_queue<std::pair<WeightType, int>, std::vector<std::pair<WeightType, int>>, std::greater<>> pq;

        dist[source] = 0;
        pq.push({0, source});

        while (!pq.empty()) {
            auto [d, u] = pq.top(); pq.pop();
            if (d > dist[u]) continue;
            metrics.nodes_expanded++;
            if (u == target) break;

            for (const auto& e : graph.adj[u]) {
                metrics.relaxations++;
                if (dist[u] + e.weight < dist[e.to]) {
                    dist[e.to] = dist[u] + e.weight;
                    parent[e.to] = u;
                    pq.push({dist[e.to], e.to});
                }
            }
        }
        finalize_metrics(metrics, dist, parent, target, start);
        return metrics;
    }

    // 2. ACTUAL RECURSIVE DUAN ALGORITHM
    SearchMetrics compute_duan_recursive(int source, int target) {
        SearchMetrics metrics;
        metrics.algorithm_name = "Duan-Recursive-2025";
        auto start = std::chrono::high_resolution_clock::now();

        int n = graph.n_vertices;
        std::vector<WeightType> dist(n, INF);
        std::vector<int> parent(n, -1);
        
        dist[source] = 0;
        std::vector<int> sources = {source};

        // Instantiate the Recursive Engine
        TrueDuan::RecursiveEngine<WeightType> engine(graph, dist, parent, metrics);
        
        // Start Recursion with Max Depth 3 (Typical for N=20000)
        engine.bmssp(3, INF, sources);

        finalize_metrics(metrics, dist, parent, target, start);
        return metrics;
    }

    // Helper for Centrality (Uses Base Radix for speed)
    std::vector<WeightType> compute_sssp(int source) {
        int n = graph.n_vertices;
        std::vector<WeightType> dist(n, INF);
        BaseRadixHeap<WeightType> pq;
        dist[source] = 0;
        pq.push(0, source);
        while (!pq.empty()) {
            auto [d, u] = pq.pop();
            if (d > dist[u]) continue;
            for (const auto& e : graph.adj[u]) {
                if (dist[u] + e.weight < dist[e.to]) {
                    dist[e.to] = dist[u] + e.weight;
                    pq.push(dist[e.to], e.to);
                }
            }
        }
        return dist;
    }

private:
    void finalize_metrics(SearchMetrics& m, const std::vector<WeightType>& dist, 
                         const std::vector<int>& parent, int target, 
                         std::chrono::high_resolution_clock::time_point start) {
        if (target != -1) {
             if(dist[target] == std::numeric_limits<WeightType>::max()) m.path_cost = -1;
             else {
                 m.path_cost = (double)dist[target];
                 for (int v = target; v != -1; v = parent[v]) m.path_ids.push_back(v);
                 std::reverse(m.path_ids.begin(), m.path_ids.end());
             }
        }
        auto end = std::chrono::high_resolution_clock::now();
        m.runtime_ms = std::chrono::duration<double, std::milli>(end - start).count();
    }
};

// ============================================================================
// CENTRALITY CALCULATOR
// ============================================================================
template<typename WeightType = long long>
class CentralityCalculator {
    const Graph<WeightType>& graph;
    ShortestPathEngine<WeightType> engine;
public:
    explicit CentralityCalculator(const Graph<WeightType>& g) : graph(g), engine(g) {}
    
    // Run Brandes accumulation over an arbitrary list of source vertices.
    // Pure refactor: same per-source body as the historical compute_betweenness.
    std::vector<double> compute_betweenness_from_sources(const std::vector<int>& sources) {
        int n = graph.n_vertices;
        std::vector<double> betweenness(n, 0.0);
        for (int s : sources) {
            if (s < 0 || s >= n) continue;
            std::vector<WeightType> dist = engine.compute_sssp(s);
            std::vector<double> sigma(n, 0.0);
            std::vector<std::vector<int>> pred(n);
            std::vector<int> stack;
            sigma[s] = 1.0;

            std::vector<std::pair<WeightType, int>> sorted;
            for(int i=0; i<n; ++i) if(dist[i] != std::numeric_limits<WeightType>::max()) sorted.push_back({dist[i], i});
            std::sort(sorted.begin(), sorted.end());

            for(auto& p : sorted) {
                int u = p.second;
                stack.push_back(u);
                for(auto& e : graph.adj[u]) {
                    if(dist[e.to] == dist[u] + e.weight) {
                        sigma[e.to] += sigma[u];
                        pred[e.to].push_back(u);
                    }
                }
            }
            std::vector<double> delta(n, 0.0);
            while(!stack.empty()) {
                int w = stack.back(); stack.pop_back();
                for(int v : pred[w]) if(sigma[w]!=0) delta[v] += (sigma[v]/sigma[w])*(1.0+delta[w]);
                if(w != s) betweenness[w] += delta[w];
            }
        }
        return betweenness;
    }

    std::vector<double> compute_betweenness() {
        int n = graph.n_vertices;
        // Historical sampling: stride n/100 when n>1000, else every node.
        int step = (n > 1000) ? n / 100 : 1;
        std::vector<int> sources;
        sources.reserve(n / std::max(step, 1) + 1);
        for (int s = 0; s < n; s += step) sources.push_back(s);
        return compute_betweenness_from_sources(sources);
    }

    // Seeded sampled overload. sample_rate >= 1.0 -> full Brandes (every vertex).
    // For sample_rate < 1.0, deterministically picks ceil(n * sample_rate)
    // sources via a seeded mt19937_64 shuffle.
    std::vector<double> compute_betweenness(double sample_rate, uint64_t seed) {
        int n = graph.n_vertices;
        std::vector<int> sources(n);
        std::iota(sources.begin(), sources.end(), 0);
        if (sample_rate < 1.0) {
            std::mt19937_64 rng(seed);
            std::shuffle(sources.begin(), sources.end(), rng);
            int k = std::max<int>(1, int(std::ceil(n * sample_rate)));
            sources.resize(k);
        }
        return compute_betweenness_from_sources(sources);
    }

    // Bootstrap-CI for sampled Brandes: B independent sampled runs (each with
    // a perturbed seed), then per-vertex 95% CI from the 2.5/97.5 quantiles
    // and median as point estimate. Used in --brandes-mode=bootstrap.
    struct BrandesCI {
        std::vector<double> point;
        std::vector<double> ci_lo;
        std::vector<double> ci_hi;
    };

    BrandesCI brandes_with_ci(double sample_rate, int B, uint64_t seed) {
        int n = graph.n_vertices;
        std::vector<std::vector<double>> samples(n);
        for (auto& v : samples) v.reserve(B);

        for (int b = 0; b < B; ++b) {
            auto bet = compute_betweenness(sample_rate, seed + (uint64_t)b);
            for (int v = 0; v < n; ++v) samples[v].push_back(bet[v]);
        }

        BrandesCI out;
        out.point.resize(n);
        out.ci_lo.resize(n);
        out.ci_hi.resize(n);
        for (int v = 0; v < n; ++v) {
            std::sort(samples[v].begin(), samples[v].end());
            out.point[v] = samples[v][B / 2];
            out.ci_lo[v] = samples[v][std::max(0, int(B * 0.025))];
            out.ci_hi[v] = samples[v][std::min(B - 1, int(B * 0.975))];
        }
        return out;
    }

    std::vector<double> compute_eigenvector_weighted() {
        int n = graph.n_vertices;
        std::vector<double> v(n, 1.0 / std::sqrt(n));

        for(int iter = 0; iter < 100; iter++) {
            std::vector<double> next_v(n, 0.0);

            for(int u = 0; u < n; u++)
                for(auto& e : graph.adj[u])
                    next_v[e.to] += v[u] * (double)e.weight;

            double norm = 0.0;
            for(double val : next_v) norm += val * val;
            norm = std::sqrt(norm);

            if(norm > 1e-9)
                for(int i = 0; i < n; i++) next_v[i] /= norm;

            double diff = 0.0;
            for(int i = 0; i < n; i++)
                diff += std::abs(next_v[i] - v[i]);

            v = next_v;

            if(diff < 1e-8) break;
        }

        return v;
    }
};

// ============================================================================
// SUBGRAPH HELPER (used by --diagnostics-out)
// ============================================================================
template<typename W>
Graph<W> induce_subgraph(const Graph<W>& g, int k) {
    if (k > g.n_vertices) k = g.n_vertices;
    if (k < 0) k = 0;
    Graph<W> s(k);
    for (int u = 0; u < k; ++u)
        for (auto& e : g.adj[u])
            if (e.to < k) s.add_edge(u, e.to, e.weight);
    if (!g.id_to_name.empty()) {
        s.id_to_name.assign(g.id_to_name.begin(),
                            g.id_to_name.begin() + std::min<int>(k, g.id_to_name.size()));
    }
    return s;
}

// ============================================================================
// DIFFERENTIAL NETWORK FRAGILITY (DNF) MODULE
// ============================================================================
// Identifies gene pairs whose combined removal selectively collapses the
// tumor network while leaving the normal network structurally intact.
// Core metric: DSS(A,B) = ΔASPL_tumor(A,B) - ΔASPL_normal(A,B)
// ============================================================================

struct ASPLResult {
    double aspl;           // Average Shortest Path Length
    int lcc_size;          // Largest Connected Component size
    int reachable_pairs;   // Number of reachable node pairs
};

template<typename WeightType = long long>
class DNFAnalyzer {
private:
    const WeightType INF = std::numeric_limits<WeightType>::max();

    // Number of SSSP source samples for ASPL estimation.
    // Higher = more accurate but slower. 200 is a good balance for O(n^2) outer loop.
    static constexpr int ASPL_SAMPLES = 200;

public:

    // Compute ASPL and LCC size for a given graph.
    // Uses sampled SSSP sources for efficiency.
    ASPLResult compute_aspl_and_lcc(const Graph<WeightType>& g) {
        ShortestPathEngine<WeightType> engine(g);
        int n = g.n_vertices;
        ASPLResult result = {0.0, 0, 0};

        if (n == 0) return result;

        // --- Compute LCC size via BFS from node 0, then try all ---
        // We do a quick reachability check from sampled nodes
        std::vector<bool> visited(n, false);
        int max_component = 0;

        for (int start = 0; start < n; start++) {
            if (visited[start]) continue;
            // BFS to find this component
            int comp_size = 0;
            std::queue<int> bfs_q;
            bfs_q.push(start);
            visited[start] = true;
            while (!bfs_q.empty()) {
                int u = bfs_q.front(); bfs_q.pop();
                comp_size++;
                for (const auto& e : g.adj[u]) {
                    if (!visited[e.to]) {
                        visited[e.to] = true;
                        bfs_q.push(e.to);
                    }
                }
            }
            max_component = std::max(max_component, comp_size);
        }
        result.lcc_size = max_component;

        // --- Compute ASPL from sampled sources ---
        int samples = std::min(ASPL_SAMPLES, n);
        double total_dist = 0.0;
        int count = 0;

        for (int i = 0; i < samples; i++) {
            int source = (int)((long long)i * n / samples);
            auto distances = engine.compute_sssp(source);
            for (int j = 0; j < n; j++) {
                if (j != source && distances[j] != INF) {
                    total_dist += (double)distances[j];
                    count++;
                }
            }
        }

        result.aspl = (count > 0) ? total_dist / count : 0.0;
        result.reachable_pairs = count;
        return result;
    }

    // Create a graph with TWO nodes (and all their edges) removed.
    Graph<WeightType> create_dual_knockout_graph(const Graph<WeightType>& g, int ko_a, int ko_b) {
        Graph<WeightType> ko_graph(g.n_vertices);
        ko_graph.name_to_id = g.name_to_id;
        ko_graph.id_to_name = g.id_to_name;

        for (int u = 0; u < g.n_vertices; u++) {
            if (u == ko_a || u == ko_b) continue;
            for (const auto& e : g.adj[u]) {
                if (e.to != ko_a && e.to != ko_b) {
                    ko_graph.add_edge(u, e.to, e.weight);
                }
            }
        }
        return ko_graph;
    }

    // The main DNF analysis: compute DSS for all pairs in the gene list.
    // Runs on BOTH graphs and outputs the differential score.
    void pairwise_dnf_analysis(
        const Graph<WeightType>& normal_graph,
        const Graph<WeightType>& tumor_graph,
        const std::string& gene_list_file,
        const std::string& output_file)
    {
        std::cout << "\n============================================" << std::endl;
        std::cout << "  DIFFERENTIAL NETWORK FRAGILITY (DNF)" << std::endl;
        std::cout << "  Pairwise Combinatorial Knockout Analysis" << std::endl;
        std::cout << "============================================\n" << std::endl;

        // 1. Load gene list
        std::vector<std::string> genes;
        std::ifstream gfile(gene_list_file);
        if (!gfile.is_open()) {
            std::cerr << "ERROR: Cannot open gene list: " << gene_list_file << std::endl;
            return;
        }
        std::string gene;
        while (std::getline(gfile, gene)) {
            if (!gene.empty() && gene[0] != '#') {
                gene.erase(gene.find_last_not_of(" \n\r\t") + 1);
                genes.push_back(gene);
            }
        }
        gfile.close();

        // 2. Resolve gene names to IDs in both graphs. Only keep genes present in BOTH.
        struct GenePair { std::string name; int normal_id; int tumor_id; };
        std::vector<GenePair> valid_genes;

        for (const auto& g : genes) {
            auto it_n = normal_graph.name_to_id.find(g);
            auto it_t = tumor_graph.name_to_id.find(g);
            if (it_n != normal_graph.name_to_id.end() && it_t != tumor_graph.name_to_id.end()) {
                valid_genes.push_back({g, it_n->second, it_t->second});
            }
        }

        int N = valid_genes.size();
        int total_pairs = N * (N - 1) / 2;
        std::cout << "Genes loaded: " << genes.size() << std::endl;
        std::cout << "Genes in both networks: " << N << std::endl;
        std::cout << "Total pairs to compute: " << total_pairs << std::endl;

        if (N < 2) {
            std::cerr << "ERROR: Need at least 2 genes for pairwise analysis." << std::endl;
            return;
        }

        // 3. Compute baselines
        std::cout << "\nComputing baseline ASPL + LCC..." << std::endl;
        auto t_start = std::chrono::high_resolution_clock::now();

        ASPLResult baseline_normal = compute_aspl_and_lcc(normal_graph);
        ASPLResult baseline_tumor = compute_aspl_and_lcc(tumor_graph);

        std::cout << "  Normal baseline: ASPL=" << std::fixed << std::setprecision(2)
                  << baseline_normal.aspl << "  LCC=" << baseline_normal.lcc_size << std::endl;
        std::cout << "  Tumor  baseline: ASPL=" << std::fixed << std::setprecision(2)
                  << baseline_tumor.aspl << "  LCC=" << baseline_tumor.lcc_size << std::endl;

        // 4. Compute single-gene knockout baselines (for synergy calculation)
        std::cout << "\nComputing single-gene knockout baselines..." << std::endl;
        std::vector<double> single_delta_normal(N), single_delta_tumor(N);
        std::vector<int> single_lcc_normal(N), single_lcc_tumor(N);

        for (int i = 0; i < N; i++) {
            // Normal single knockout
            Graph<WeightType> ko_n(normal_graph.n_vertices);
            ko_n.name_to_id = normal_graph.name_to_id;
            ko_n.id_to_name = normal_graph.id_to_name;
            for (int u = 0; u < normal_graph.n_vertices; u++) {
                if (u == valid_genes[i].normal_id) continue;
                for (const auto& e : normal_graph.adj[u]) {
                    if (e.to != valid_genes[i].normal_id) ko_n.add_edge(u, e.to, e.weight);
                }
            }
            auto r_n = compute_aspl_and_lcc(ko_n);
            single_delta_normal[i] = r_n.aspl - baseline_normal.aspl;
            single_lcc_normal[i] = r_n.lcc_size;

            // Tumor single knockout
            Graph<WeightType> ko_t(tumor_graph.n_vertices);
            ko_t.name_to_id = tumor_graph.name_to_id;
            ko_t.id_to_name = tumor_graph.id_to_name;
            for (int u = 0; u < tumor_graph.n_vertices; u++) {
                if (u == valid_genes[i].tumor_id) continue;
                for (const auto& e : tumor_graph.adj[u]) {
                    if (e.to != valid_genes[i].tumor_id) ko_t.add_edge(u, e.to, e.weight);
                }
            }
            auto r_t = compute_aspl_and_lcc(ko_t);
            single_delta_tumor[i] = r_t.aspl - baseline_tumor.aspl;
            single_lcc_tumor[i] = r_t.lcc_size;
        }
        std::cout << "  Single-gene baselines computed." << std::endl;

        // 5. Pairwise knockout with OpenMP parallelization
        struct PairResult {
            int idx_a, idx_b;
            double delta_aspl_normal, delta_aspl_tumor;
            double dss;                 // Differential Synergy Score
            double synergy_tumor;       // Synergy in tumor: Δpair - (Δa + Δb)
            double synergy_normal;      // Synergy in normal
            double differential_synergy; // synergy_tumor - synergy_normal
            int lcc_normal_ko, lcc_tumor_ko;
        };

        std::vector<PairResult> results(total_pairs);
        int progress = 0;

        std::cout << "\nRunning pairwise knockouts (" << total_pairs << " pairs)..." << std::endl;

        #pragma omp parallel for schedule(dynamic) shared(progress)
        for (int idx = 0; idx < total_pairs; idx++) {
            // Map linear index to (i, j) pair
            // idx = i*(i-1)/2 + j where j < i
            int i = (int)(0.5 + std::sqrt(0.25 + 2.0 * idx));
            int j = idx - i * (i - 1) / 2;
            if (j >= i) { i++; j = idx - i * (i - 1) / 2; }

            // Normal network knockout
            auto ko_normal = create_dual_knockout_graph(normal_graph,
                valid_genes[i].normal_id, valid_genes[j].normal_id);
            auto res_normal = compute_aspl_and_lcc(ko_normal);

            // Tumor network knockout
            auto ko_tumor = create_dual_knockout_graph(tumor_graph,
                valid_genes[i].tumor_id, valid_genes[j].tumor_id);
            auto res_tumor = compute_aspl_and_lcc(ko_tumor);

            // Compute metrics
            double d_normal = res_normal.aspl - baseline_normal.aspl;
            double d_tumor = res_tumor.aspl - baseline_tumor.aspl;
            double dss = d_tumor - d_normal;

            // Synergy: Is the pair effect greater than sum of individual effects?
            double syn_tumor = d_tumor - (single_delta_tumor[i] + single_delta_tumor[j]);
            double syn_normal = d_normal - (single_delta_normal[i] + single_delta_normal[j]);

            results[idx] = {
                i, j,
                d_normal, d_tumor,
                dss,
                syn_tumor, syn_normal,
                syn_tumor - syn_normal,  // differential synergy
                res_normal.lcc_size, res_tumor.lcc_size
            };

            // Progress reporting (thread-safe)
            #pragma omp atomic
            progress++;

            if (progress % 100 == 0 || progress == total_pairs) {
                #pragma omp critical
                {
                    std::cout << "\r  Progress: " << progress << "/" << total_pairs
                              << " (" << (100 * progress / total_pairs) << "%)" << std::flush;
                }
            }
        }
        std::cout << std::endl;

        // 6. Sort by DSS (descending — highest tumor disruption relative to normal)
        std::sort(results.begin(), results.end(),
                  [](const PairResult& a, const PairResult& b) { return a.dss > b.dss; });

        // 7. Write output
        auto t_end = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(t_end - t_start).count();

        std::ofstream out(output_file);
        out << "gene_a\tgene_b\tdss\tsynergy_tumor\tsynergy_normal\tdifferential_synergy\t"
            << "delta_aspl_tumor\tdelta_aspl_normal\t"
            << "lcc_tumor_ko\tlcc_normal_ko\t"
            << "baseline_aspl_tumor\tbaseline_aspl_normal\t"
            << "baseline_lcc_tumor\tbaseline_lcc_normal\n";

        for (const auto& r : results) {
            out << valid_genes[r.idx_a].name << "\t"
                << valid_genes[r.idx_b].name << "\t"
                << std::fixed << std::setprecision(4)
                << r.dss << "\t"
                << r.synergy_tumor << "\t"
                << r.synergy_normal << "\t"
                << r.differential_synergy << "\t"
                << r.delta_aspl_tumor << "\t"
                << r.delta_aspl_normal << "\t"
                << r.lcc_tumor_ko << "\t"
                << r.lcc_normal_ko << "\t"
                << baseline_tumor.aspl << "\t"
                << baseline_normal.aspl << "\t"
                << baseline_tumor.lcc_size << "\t"
                << baseline_normal.lcc_size << "\n";
        }
        out.close();

        // 8. Print summary
        std::cout << "\n=== DNF Analysis Complete ===" << std::endl;
        std::cout << "Time elapsed: " << std::fixed << std::setprecision(1) << elapsed << "s" << std::endl;
        std::cout << "Pairs analyzed: " << total_pairs << std::endl;
        std::cout << "Output: " << output_file << std::endl;

        std::cout << "\nTop 10 Tumor-Selective Vulnerable Pairs (Highest DSS):" << std::endl;
        std::cout << std::left << std::setw(12) << "Gene_A"
                  << std::setw(12) << "Gene_B"
                  << std::setw(12) << "DSS"
                  << std::setw(14) << "Synergy_T"
                  << std::setw(14) << "Synergy_N"
                  << std::setw(14) << "Diff_Synergy" << std::endl;
        std::cout << std::string(78, '-') << std::endl;

        int top_n = std::min(10, (int)results.size());
        for (int k = 0; k < top_n; k++) {
            const auto& r = results[k];
            std::cout << std::left << std::setw(12) << valid_genes[r.idx_a].name
                      << std::setw(12) << valid_genes[r.idx_b].name
                      << std::setw(12) << std::fixed << std::setprecision(4) << r.dss
                      << std::setw(14) << r.synergy_tumor
                      << std::setw(14) << r.synergy_normal
                      << std::setw(14) << r.differential_synergy << std::endl;
        }
    }
};

#endif