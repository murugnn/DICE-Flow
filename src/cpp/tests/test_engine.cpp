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
