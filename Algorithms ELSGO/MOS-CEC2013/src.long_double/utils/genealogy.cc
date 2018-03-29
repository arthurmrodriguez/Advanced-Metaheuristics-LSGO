#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/fusion/include/adapt_struct.hpp>

#include <iostream>
#include <ostream>
#include <fstream>
#include <string>

using namespace boost::spirit;
using namespace boost::spirit::qi;

///////////////////////////////////////////////////////////////////////////////
//  Aux structs
///////////////////////////////////////////////////////////////////////////////

struct BasicIndiv {
  unsigned id_;
  unsigned islandid_;
  long double fitness_;
  unsigned technique_;
};

BOOST_FUSION_ADAPT_STRUCT(
			  BasicIndiv,
			  (unsigned, id_)
			  (unsigned, islandid_)
			  (long double, fitness_)
			  (unsigned, technique_)
			  )

struct IndivGen {
  unsigned act_;
  unsigned gen_;
  BasicIndiv ind_;
  BasicIndiv dad_;
};

BOOST_FUSION_ADAPT_STRUCT(
			  IndivGen,
			  (unsigned, act_)
			  (unsigned, gen_)
			  (BasicIndiv, ind_)
			  (BasicIndiv, dad_)
			  )

struct SexualCrossover {
  unsigned gen_;
  BasicIndiv ch1_;
  BasicIndiv ch2_;
  BasicIndiv dad_;
  BasicIndiv mom_;
};

BOOST_FUSION_ADAPT_STRUCT(
			  SexualCrossover,
			  (unsigned, gen_)
			  (BasicIndiv, ch1_)
			  (BasicIndiv, ch2_)
			  (BasicIndiv, dad_)
			  (BasicIndiv, mom_)
			  )

struct MultipleCrossover {
  unsigned gen_;
  std::vector<BasicIndiv> children;
  std::vector<BasicIndiv> parents;
};

BOOST_FUSION_ADAPT_STRUCT(
        MultipleCrossover,
        (unsigned, gen_)
        (std::vector<BasicIndiv>, children)
        (std::vector<BasicIndiv>, parents)
        )

class Individual {
public:

  Individual () : id_(0), islandid_(0), fitness_(0.0), technique_(0), nparents_(0), firstgen_(0), lastgen_(0) {}
  unsigned id_;
  unsigned islandid_;
  long double fitness_;
  unsigned technique_;
  unsigned parents_ids_[5];
  unsigned parents_islandids_[5];
  unsigned nparents_;
  unsigned firstgen_;
  unsigned lastgen_;

  std::vector< std::pair<unsigned, unsigned> > children;

  std::ostream& Print (std::ostream& os) const {
    os << "ID: " << id_ << std::endl;
    os << "ISLAND ID: " << islandid_ << std::endl;
    os << "FITNESS: " << fitness_ << std::endl;
    os << "TECHNIQUE: " << technique_ << std::endl;
    os << "FIRST GEN: " << firstgen_ << std::endl;
    os << "LAST GEN: " << lastgen_ << std::endl;
    os << "PARENTS: " << std::endl;

    for (unsigned i = 0; i < nparents_; i++)
      os << "  PARENT " << i << ": " << parents_ids_[i] << "  " << parents_islandids_[i] << std::endl;

    os << "CHILDREN: " << std::endl;

    for (unsigned i = 0; i < children.size(); i++)
      os << "  CHILD " << i << ": " << children[i].first << " in generation " << children[i].second << std::endl;

    return os;

  }

};

std::ostream& operator<< (std::ostream& os, const Individual& ind) {
  return ind.Print(os);
}

////////////////////////////////////////////////////////////////////////////
//  Main program
////////////////////////////////////////////////////////////////////////////
int main (int argc, char** argv) {

  // Type renaming, for easy writing
  typedef std::string::const_iterator iterator_type;

  // Rules for the parser
  rule<iterator_type, BasicIndiv(), space_type> mut_extra, best;
  rule<iterator_type, unsigned(), space_type> generation, unary_lit, stats;
  rule<iterator_type, std::vector<long double>(), space_type> stats_techs;
  rule<iterator_type, BasicIndiv(), space_type>individual;
  rule<iterator_type, IndivGen(), space_type> unary_op;
  rule<iterator_type, SexualCrossover(), space_type> sexual_cross;
  rule<iterator_type, MultipleCrossover(), space_type> mult_cross;

  stats = uint_ [_val = _1] >> *char_;
  stats_techs = uint_ >> lit("score") >> '=' >> '[' >> *(char_("0-9a-zA-Z,.=")) >> ']' >> lit("participation") >> '=' >> '[' >> (double_ % ',') [_val = _1] >> ']' >> *char_;
  generation %= '(' >> uint_ >> ')' >> ':';
  individual %= uint_ >> uint_ >> ',' >> double_ >> uint_;
  mut_extra %= '|' >> lit("dad") >> ':' >> individual;
  unary_lit = lit("Deceased") [_val = 2] | lit("New Genome") [_val = 1] | lit("Mutated") [_val = 3];
  unary_op %= unary_lit >> generation >> individual >> (mut_extra | attr(BasicIndiv()));
  sexual_cross %= lit("Children") >> generation >> individual >> ';' >> individual >> '|' >> lit("dad") >> ':' >> individual >> ';' >> lit("mom") >> ':' >> individual;
  mult_cross %= lit("Children") >> generation >> individual % ';' >> '|' >> (lit("parent") >> ':' >> individual) % ';';
  best = lit("Best") >> generation >> individual [_val = _1];

  // Check for the correct number of arguments
  if (argc < 3) {
    std::cerr << "Wrong number of arguments." << std::endl;
    std::cerr << "Usage: " << argv[0] << " log_file max_genomes" << std::endl;
    exit(-1);
  }

  std::string log_file = argv[1];

  // Open the file for reading
  std::ifstream input (log_file.c_str());

  if (!input.good()) {
    std::cerr << "Error: the file '" << log_file << "' does not exist or could not be open." << std::endl;
    exit(-2);
  }

  // Number of techniques (parsed from first statistics line)
  unsigned nTechs = 0;

  // Last generation of the experiment
  unsigned lastgen = 0;

  // Parse maximum number of genomes in the log file
  unsigned max_gens = strtol(argv[2], NULL, 10);

  // Vector for the genealogy of size max_gens
  std::vector<Individual*> genealogy (max_gens);

  // Variable to store the highest ID read
  unsigned highestID = 0;

  // Variable to store the ID of the best genome
  unsigned best_id = 0;

  // Buffer to read from file
  std::string str;

  // Output to a file, rather than to cout
  std::ofstream stdout_f ((log_file + "_stdout").c_str());

  // Main loop for processing lines
  while (!getline (input, str).eof()) {

#ifndef NDEBUG
    stdout_f << "===============================" << std::endl;
    stdout_f << "Read: " << str << std::endl;
    stdout_f << "===============================" << std::endl;
#endif

    // Variables to store the results of parsing
    IndivGen ind;
    SexualCrossover inds;
    MultipleCrossover mult_inds;

    // Need to pass the iterators of 'str' in named variables
    iterator_type iter = str.begin ();
    iterator_type end  = str.end   ();

    // Return value and parsed action type
    bool r;
    unsigned action;

    // Select the rule to apply
    if (str[0] == 'N' || str[0] == 'D' || str[0] == 'M') {

      // Parse a unary operator (New Genome, Deceased or Mutated)
      r = phrase_parse(iter, end, unary_op, qi::space, ind);
      action = ind.act_;

      // Special case: if 'Mutated' but not previously genome in the genealogy,
      // we should consider this as a creation of a new individual
      if (action == 3 && (genealogy[ind.ind_.id_]) == NULL) {
        action = 5;
      }

    }
    else if (str[0] == 'C') {

      size_t pos = str.find("dad");

      if (pos != std::string::npos) {
        // Parse a Sexual Crossover
        r = phrase_parse(iter, end, sexual_cross, qi::space, inds);
        action = 4;
      }
      else {
        // Parse a Multiple Crossover
        r = phrase_parse(iter, end, mult_cross, qi::space, mult_inds);
        action = 6;
      }

    }
    else if (str[0] == 'B') {

      BasicIndiv ind;

      // Parse a Multiple Crossover
      r = phrase_parse(iter, end, best, qi::space, ind);

      best_id = ind.id_;

    }
    else {

      static bool firstTime = false;

      if (!firstTime) {
        std::vector<long double> tmp;
        r = phrase_parse(iter, end, stats_techs, qi::space, tmp);
        nTechs = tmp.size();
        firstTime = true;
      }
      else
        // Parse statistics
        r = phrase_parse(iter, end, stats, qi::space, lastgen);

      continue;

    }

    // Aux variables to store the individuals to be added to the genealogy
    Individual *theInd1, *theInd2;

    if (r && iter == end) {

      switch (action) {
      case 1:
        theInd1 = new Individual();
        theInd1->firstgen_ = ind.gen_;
        theInd1->id_ = ind.ind_.id_;
        theInd1->islandid_ = ind.ind_.islandid_;
        theInd1->fitness_ = ind.ind_.fitness_;
        theInd1->technique_ = ind.ind_.technique_;

#ifndef NDEBUG
        stdout_f << "New Genome created: " << std::endl << *theInd1 << std::endl;
#endif

        genealogy[theInd1->id_] = theInd1;

        if (theInd1->id_ > highestID)
          highestID = theInd1->id_;

        break;
      case 2:
        (genealogy[ind.ind_.id_])->lastgen_ = ind.gen_;

#ifndef NDEBUG
        stdout_f << "Genome deceased: " << std::endl << *(genealogy[ind.ind_.id_]) << std::endl;
#endif

        break;
      case 3:
        (genealogy[ind.ind_.id_])->fitness_ = ind.ind_.fitness_;

#ifndef NDEBUG
        stdout_f << "Genome mutated: " << std::endl << *(genealogy[ind.ind_.id_]) << std::endl;
#endif

        break;
      case 4:
        theInd1 = new Individual();
        theInd1->firstgen_ = inds.gen_;
        theInd1->id_ = inds.ch1_.id_;
        theInd1->islandid_ = inds.ch1_.islandid_;
        theInd1->fitness_ = inds.ch1_.fitness_;
        theInd1->technique_ = ind.ind_.technique_;
        theInd1->nparents_ = 2;
        theInd1->parents_ids_[0] = inds.dad_.id_;
        theInd1->parents_islandids_[0] = inds.dad_.islandid_;
        theInd1->parents_ids_[1] = inds.mom_.id_;
        theInd1->parents_islandids_[1] = inds.mom_.islandid_;

#ifndef NDEBUG
        stdout_f << "New Genome created: " << std::endl << *theInd1 << std::endl;
#endif

        genealogy[theInd1->id_] = theInd1;

        if (theInd1->id_ > highestID)
          highestID = theInd1->id_;


        theInd2 = new Individual();
        theInd2->firstgen_ = inds.gen_;
        theInd2->id_ = inds.ch2_.id_;
        theInd2->islandid_ = inds.ch2_.islandid_;
        theInd2->fitness_ = inds.ch2_.fitness_;
        theInd2->technique_ = ind.ind_.technique_;
        theInd2->nparents_ = 2;
        theInd2->parents_ids_[0] = inds.dad_.id_;
        theInd2->parents_islandids_[0] = inds.dad_.islandid_;
        theInd2->parents_ids_[1] = inds.mom_.id_;
        theInd2->parents_islandids_[1] = inds.mom_.islandid_;

#ifndef NDEBUG
        stdout_f << "New Genome created: " << std::endl << *theInd2 << std::endl;
#endif

        genealogy[theInd2->id_] = theInd2;

        if (theInd2->id_ > highestID)
          highestID = theInd2->id_;

        break;
      case 5:

#ifndef NDEBUG
        stdout_f << "New only mutated individual" << std::endl;
#endif

        theInd1 = new Individual();
        theInd1->firstgen_ = ind.gen_;
        theInd1->id_ = ind.ind_.id_;
        theInd1->islandid_ = ind.ind_.islandid_;
        theInd1->fitness_ = ind.ind_.fitness_;
        theInd1->technique_ = ind.ind_.technique_;

        theInd1->nparents_ = 1;
        theInd1->parents_ids_[0] = ind.dad_.id_;
        theInd1->parents_islandids_[0] = ind.dad_.islandid_;

#ifndef NDEBUG
        stdout_f << "New Genome created: " << std::endl << *theInd1 << std::endl;
#endif

        genealogy[theInd1->id_] = theInd1;

        if (theInd1->id_ > highestID)
          highestID = theInd1->id_;

        break;

      case 6:

#ifndef NDEBUG
        stdout_f << "New individual(s) my multiple crossover" << std::endl;
#endif

        for (unsigned i = 0; i < mult_inds.children.size(); i++) {

          theInd1 = new Individual();
          theInd1->firstgen_ = mult_inds.gen_;
          theInd1->id_ = mult_inds.children[i].id_;
          theInd1->islandid_ = mult_inds.children[i].islandid_;
          theInd1->fitness_ = mult_inds.children[i].fitness_;
          theInd1->technique_ = ind.ind_.technique_;

          theInd1->nparents_ = mult_inds.parents.size();

          for (unsigned j = 0; j < mult_inds.parents.size(); j++) {
            theInd1->parents_ids_[j] = mult_inds.parents[j].id_;
            theInd1->parents_islandids_[j] = mult_inds.parents[j].islandid_;
          }

#ifndef NDEBUG
        stdout_f << "New Genome created: " << std::endl << *theInd1 << std::endl;
#endif

          genealogy[theInd1->id_] = theInd1;

          if (theInd1->id_ > highestID)
            highestID = theInd1->id_;

        }
        break;
      default:
        break;
      }

    }
    else {
      std::cerr << "Parsing of line: " << str << " FAILED!!!" << std::endl;
    }

  }

  std::vector< std::vector<long double> > participations (lastgen + 1, std::vector<long double> (nTechs, 0.0));
  std::vector< std::vector<long double> > inc_fit        (lastgen + 1, std::vector<long double> (nTechs, 0.0));
  std::vector<long double> entropy  (lastgen + 1, 0.0);
  std::vector<long double> endogamy (lastgen + 1, 0.0);
  std::vector<unsigned> new_genomes (lastgen + 1, 0);

  std::vector<unsigned> parents_ids;
  std::vector<bool> ancestors (highestID + 1, false);

  parents_ids.push_back(best_id);

  unsigned next_pos = 0;

  do {

    unsigned current_id = parents_ids[next_pos];
    Individual& current_indiv = *(genealogy[current_id]);

    stdout_f << "Next pos: " << next_pos << std::endl;
    stdout_f << "Current ID: " << current_id << std::endl;

    new_genomes[current_indiv.firstgen_]++;

    for (unsigned i = 0; i < current_indiv.nparents_; i++) {

      genealogy[current_indiv.parents_ids_[i]]->children.push_back(std::make_pair(current_id, current_indiv.firstgen_));
      participations[current_indiv.firstgen_][genealogy[current_indiv.parents_ids_[i]]->technique_]++;
      inc_fit[current_indiv.firstgen_][current_indiv.technique_] += (current_indiv.fitness_ - genealogy[current_indiv.parents_ids_[i]]->fitness_) / current_indiv.nparents_;
      endogamy[current_indiv.firstgen_] += (current_indiv.technique_ == genealogy[current_indiv.parents_ids_[i]]->technique_) ? (long double) 1.0 / (long double) current_indiv.nparents_ : 0.0;

      if (!ancestors[current_indiv.parents_ids_[i]]) {
        parents_ids.push_back(current_indiv.parents_ids_[i]);
        stdout_f << "New Genome added: " << std::endl << *(genealogy[current_indiv.parents_ids_[i]]) << std::endl << std::endl;
        ancestors[current_indiv.parents_ids_[i]] = true;
      }
      else {
        stdout_f << "Added children: " << std::endl << *(genealogy[current_indiv.parents_ids_[i]]) << std::endl << std::endl;
      }

    }

    next_pos++;

  } while (next_pos < parents_ids.size());

  std::ofstream entropy_f ((log_file + "_entropy").c_str());
  std::ofstream inc_fit_f ((log_file + "_inc_fit").c_str());
  std::ofstream endogamy_f ((log_file + "_endogamy").c_str());

  for (unsigned i = 1; i <= lastgen; i++) {

    entropy_f << i << ",";
    inc_fit_f << i << ",";
    endogamy_f << i << ",";

    unsigned total = 0;
    long double entropy = 0.0;

    for (unsigned j = 0; j < nTechs; j++)
      total += participations[i][j];

    for (unsigned j = 0; j < nTechs; j++) {
      long double prop = (total > 0) ? (long double) participations[i][j] / (long double) total : 0.0;

//      entropy_f << participations[i][j] << ",";
      entropy_f << prop << ",";

      if (prop > 0)
        entropy += prop * log2(prop);

      if (j == nTechs - 1)
         inc_fit_f << inc_fit[i][j];
      else
         inc_fit_f << inc_fit[i][j] << ",";
    }

    entropy *= -1;

//    entropy_f << "  => entropy = " << entropy << std::endl;
    entropy_f << entropy << std::endl;
    inc_fit_f << std::endl;

    endogamy_f << ((new_genomes[i] == 0) ? 0 : endogamy[i] / (long double) new_genomes[i]) << std::endl;

  }

  entropy_f.close();
  inc_fit_f.close();
  endogamy_f.close();

  stdout_f << "Number of techniques: " << nTechs << std::endl;
  stdout_f << "Last generation: " << lastgen << std::endl;
  stdout_f << "Highest ID: " << highestID << std::endl;

  stdout_f.close();

  return 0;

}
