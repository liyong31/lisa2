// Author: Yong Li (liyong@ios.ac.cn)
#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <unordered_map>

#include <spot/tl/formula.hh>
#include <spot/tl/print.hh>
#include <spot/tl/parse.hh>
#include <spot/tl/ltlf.hh>
#include <spot/tl/relabel.hh>
#include <spot/tl/simplify.hh>
#include <spot/twaalgos/isdet.hh>
#include <spot/twaalgos/translate.hh>
#include <spot/twaalgos/hoa.hh>
#include <spot/twaalgos/product.hh>
#include <spot/twaalgos/remprop.hh>
#include <spot/twaalgos/dualize.hh>

#include <spot/misc/bddlt.hh>
#include <spot/misc/escape.hh>
// #include <spot/misc/game.hh>
#include <spot/twaalgos/game.hh>
#include <spot/twaalgos/synthesis.hh>
#include <spot/misc/timer.hh>
#include <spot/twa/twagraph.hh>
#include <spot/twaalgos/aiger.hh>
#include <spot/twaalgos/degen.hh>
#include <spot/twaalgos/determinize.hh>
#include <spot/twaalgos/parity.hh>
#include <spot/twaalgos/sbacc.hh>
#include <spot/twaalgos/totgba.hh>
#include <spot/twaalgos/translate.hh>
#include <spot/twaalgos/simulation.hh>
#include <spot/twaalgos/split.hh>
#include <spot/twaalgos/toparity.hh>

spot::twa_graph_ptr
to_dpa(const spot::twa_graph_ptr &split);

// solve games
int solve_game(spot::twa_graph_ptr nba, std::vector<std::string>& input, std::vector<std::string>& output);
