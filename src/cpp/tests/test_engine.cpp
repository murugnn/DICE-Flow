#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "../engine.hpp"

// Tiny directed graph: A -> B -> C, edges weight 1.
// Hand-computed Brandes betweenness (unnormalized): A=0, B=1, C=0.
TEST_CASE("Brandes betweenness on a 3-node path graph") {
    Graph<long long> g(3);
    g.id_to_name = {"A", "B", "C"};
    g.name_to_id = {{"A", 0}, {"B", 1}, {"C", 2}};
    g.add_edge(0, 1, 1);
    g.add_edge(1, 2, 1);

    CentralityCalculator<long long> calc(g);
    auto bet = calc.compute_betweenness();

    REQUIRE(bet.size() == 3);
    CHECK(bet[0] == doctest::Approx(0.0));
    CHECK(bet[1] == doctest::Approx(1.0));
    CHECK(bet[2] == doctest::Approx(0.0));
}

// Seeded overload: same 3-node graph; sample_rate=1.0 must equal full Brandes.
TEST_CASE("Sampled Brandes equals full Brandes when sample_rate=1.0") {
    Graph<long long> g(3);
    g.add_edge(0, 1, 1);
    g.add_edge(1, 2, 1);

    CentralityCalculator<long long> calc(g);
    auto bet_full = calc.compute_betweenness();
    auto bet_sampled = calc.compute_betweenness(1.0, /*seed=*/42);
    REQUIRE(bet_full.size() == bet_sampled.size());
    for (size_t i = 0; i < bet_full.size(); ++i)
        CHECK(bet_full[i] == doctest::Approx(bet_sampled[i]));
}

// Star 0 -> {1,2,3,4}: centre is a source-only, no node sits on a path
// between any other pair, so all betweenness values are zero.
TEST_CASE("Brandes on directed star with no transit yields zero centrality") {
    Graph<long long> g(5);
    g.add_edge(0, 1, 1);
    g.add_edge(0, 2, 1);
    g.add_edge(0, 3, 1);
    g.add_edge(0, 4, 1);
    CentralityCalculator<long long> calc(g);
    auto bet = calc.compute_betweenness();
    for (auto v : bet) CHECK(v == doctest::Approx(0.0));
}

// Bootstrap CI must produce per-vertex bounds that bracket the point estimate.
TEST_CASE("brandes_with_ci returns CI bounds spanning the point estimate") {
    Graph<long long> g(10);
    for (int i = 0; i < 9; ++i) g.add_edge(i, i + 1, 1);
    g.id_to_name.resize(10);

    CentralityCalculator<long long> calc(g);
    auto res = calc.brandes_with_ci(/*sample_rate=*/0.3, /*B=*/20, /*seed=*/42);
    CHECK(res.point.size() == 10);
    CHECK(res.ci_lo.size() == 10);
    CHECK(res.ci_hi.size() == 10);
    for (size_t i = 0; i < res.point.size(); ++i) {
        CHECK(res.ci_lo[i] <= res.point[i] + 1e-6);
        CHECK(res.ci_hi[i] >= res.point[i] - 1e-6);
    }
}
