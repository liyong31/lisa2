// Author: Yong Li (liyong@ios.ac.cn)

#include "spotsynt.hh"

using namespace spot;
using namespace std;


spot::twa_graph_ptr
to_dpa(const spot::twa_graph_ptr &split)
{
  // if the input automaton is deterministic, degeneralize it to be sure to
  // end up with a parity automaton
  auto dpa = spot::tgba_determinize(split,
                                    false, true, true, false);
  dpa->merge_edges();
  if (false)
    dpa = spot::sbacc(dpa);
  spot::reduce_parity_here(dpa, true);
  spot::change_parity_here(dpa, spot::parity_kind_max,
                           spot::parity_style_odd);
  assert((
      [&dpa]() -> bool {
        bool max, odd;
        dpa->acc().is_parity(max, odd);
        return max && odd;
      }()));
  assert(spot::is_deterministic(dpa));
  return dpa;
}

int solve_game(spot::twa_graph_ptr nba, std::vector<string>& input, std::vector<string>& output)
{
  
  // now invoke ltlsynt
  spot::bdd_dict_ptr dict = nba->get_dict();

  cout << "number of input aps: " << input.size() << " number of output aps: " << output.size() << endl;
  bdd all_inputs = bddtrue;
  bdd all_outputs = bddtrue;
  for (string &in : input)
  {
    formula f = formula::ap(in);
    // cout << "in: " << f << endl;
    unsigned v = nba->register_ap(f);
    all_inputs &= bdd_ithvar(v);
  }
  //cout << "done with input propositions" << endl;
  for (string &out : output)
  {
    formula f = formula::ap(out);
    // cout << "out: " << f << endl;
    unsigned v = nba->register_ap(f);
    all_outputs &= bdd_ithvar(v);
  }
  spot::process_timer timer;
  timer.start();
  spot::stopwatch sw;
  spot::twa_graph_ptr dpa = nullptr;
  sw.start();
  /*
  sw.start();
  auto split = split_2step(nba, all_inputs, true);
  double split_time = sw.stop();
  std::cerr << "split inputs and outputs done in " << split_time
              << " seconds\nautomaton has "
              << split->num_states() << " states\n";
  sw.start();
  dpa = to_dpa(split);
  std::cerr << "determinization done\nDPA has "
              << dpa->num_states() << " states, "
              << dpa->num_sets() << " colors\n";
  dpa->merge_states();
  double paritize_time = sw.stop();
  std::cerr << "simplification done\nDPA has "
              << dpa->num_states() << " states\n"
              << "determinization and simplification took "
              << paritize_time << " seconds\n";
  */
  // this setting seems to be more efficient
  auto tmp = to_dpa(nba);
  // ofstream os1("dpa.hoa");
  // spot::print_hoa(os1, nba);
  std::cerr << "DPA has "
            << tmp->num_states() << " states, "
            << tmp->num_sets() << " colors\n";
  tmp->merge_states();
  double paritize_time = sw.stop();
  std::cerr << "simplification done\nDPA has "
            << tmp->num_states() << " states\n"
            << "simplification took: "
            << paritize_time*1000.0 << "ms...\n";
  sw.start();
  dpa = split_2step(tmp, all_outputs, true);
  spot::colorize_parity_here(dpa, true);
  double split_time = sw.stop();
  std::cerr << "split inputs and outputs done in " << split_time*1000.0
            << "ms...\nautomaton has "
            << tmp->num_states() << " states\n";
  
  unsigned nb_states_dpa = dpa->num_states();
  sw.start();
  
  // auto owner = complete_env(dpa);
  // auto pg = spot::parity_game(dpa, owner);
  bool win = solve_parity_game(dpa);
  double bgame_time = sw.stop();
  std::cerr << "solve game in: " << bgame_time *1000.0 << "ms...\n";

  if (win)
  {
    std::cout << "REALIZABLE\n";
    // if (!opt_real)
    // {
    //     if (want_time)
    //         sw.start();
    //     auto strat_aut =
    //         strat_to_aut(pg, strategy[1], dpa, all_outputs);
    //     if (want_time)
    //         strat2aut_time = sw.stop();

    //     // output the winning strategy
    //     if (opt_print_aiger)
    //         spot::print_aiger(std::cout, strat_aut);
    //     else
    //     {
    //         automaton_printer printer;
    //         printer.print(strat_aut, timer);
    //     }
    // }
    return 0;
  }
  else
  {
    std::cout << "UNREALIZABLE\n";
    return 1;
  }
}
