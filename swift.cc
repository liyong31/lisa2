#include <algorithm>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <unordered_map>
#include <string>
#include <vector>
#include <queue>
#include <list>
#include <set>
#include <chrono> 

// spot 
#include <spot/twaalgos/hoa.hh>
#include <spot/twaalgos/translate.hh>
#include <spot/twaalgos/isdet.hh>
#include <spot/twaalgos/product.hh>
#include <spot/twaalgos/complement.hh>
#include <spot/twaalgos/sum.hh>
#include <spot/twaalgos/word.hh>
#include <spot/twaalgos/degen.hh>
#include <spot/twaalgos/remprop.hh>
#include <spot/twaalgos/minimize.hh>

#include <spot/twa/twagraph.hh>

#include <spot/tl/formula.hh>
#include <spot/tl/print.hh>
#include <spot/tl/parse.hh>
#include <spot/tl/relabel.hh>
#include <spot/tl/ltlf.hh>
#include <spot/tl/simplify.hh>
#include <spot/tl/nenoform.hh>

#include <spot/misc/optionmap.hh>
#include <spot/misc/timer.hh>

#include "mona.hh"
#include "spotutil.hh"
#include "debug.hh"
#include "ltlf2fol.hh"
// #include "spotsynt.hh"

//#include "dfwamin3.hh"

extern "C" {
#include <mona/bdd.h>
#include <mona/dfa.h>
#include <mona/mem.h>
}

#include <map>

using namespace spot;
using namespace std;

#define info(o) cout << "[INFO] " << ( o ) << endl
#define erro(o) cerr << "[ERRO] " << ( o ) << endl

// options

static struct opt_t
{
	const char* _ltlfile_name = nullptr;
	const char* _parfile_name = nullptr;
	const char* _outfile_name = nullptr;

	bool _dag = false;
	bool _mona = false;

	bool _synthesis = false;
	// bool _out_start = false;
	bool _env_first = false;

	bool _verbose = false;

}* opt;

twa_graph_ptr minimize_explicit(twa_graph_ptr A)
{
	twa_graph_ptr C = spot::minimize_wdba(A);
	//A = spot::minimize_obligation(A);
	// check equivalence of two automata
#ifdef DEBUG
	string word = is_twa_equivalent(A, C);
	if(word.size() == 0)
	{
		cout << "A: equivalent two automata" << endl;
	}
#endif
	return C;
}


void print_usage()
{
	cout << "Usage: swift [OPTION...] [FILENAME[/COL]...]" << endl;
	cout << "Read a formula file and output the number of states of the constructed DFA" << endl << endl;
	cout << " Input options:" << endl;
	cout << " -h  " << "                  show this help page" << endl;
	// cout << " -syn" << "                  synthesize after DFA construction (default false)" << endl;
	cout << " -dag" << "                  optimisations including demorgan and remove duplicates" << endl;
	cout << " -mona"<< "                  use MONA DFA for composition (Spot DFA used by default)" << endl;
	// cout << " -part" << " <file>          the file specifying the input and output propositions" << endl;
	cout << " -ltlf" << " <file>          the file specifying the input LTLf formula" << endl;
	cout << " -out" << " <file>           the file for the minimal DFA" << endl;
	cout << " -v  " << "                  verbose output" << endl;
	// cout << " -env" << "                  environment plays first" << endl;
}

void parse_opt(int argc, char** argv)
{
	// first one is lisa, separated by space
	if(argc == 1)
	{
		print_usage();
	}
	for(int i = 1; i < argc; i ++)
	{
		string s(argv[i]);
		//cout << argv[i] << endl;
		if(s.size() == 0)
		{
			continue;
		}
		if(s == "-dag")
		{
			opt->_dag = true;
			//cout << "hello" << s << endl;
			continue;
		}else 
		if(s == "-mona")
		{
			opt->_mona = true;
			//cout << "hello" << s << endl;
			continue;
		}else
		if(s == "-syn")
		{
			opt->_synthesis = true;
			continue;
		}else
		if(s == "-env")
		{
			opt->_env_first = true;
			continue;
		}
		if(s == "-v")
		{
			opt->_verbose = true;
			continue;
		}
		else if(s == "-ltlf" && i + 1 < argc)
		{
			opt->_ltlfile_name = argv[i + 1];
			//cout << "hello" << argv[i+1] << endl;
			i ++;
			continue;
		}else if(s == "-out" && i + 1 < argc)
		{
			opt->_outfile_name = argv[i + 1];
			//cout << "hello" << argv[i+1] << endl;
			i ++;
			continue;
		}else
		if(s == "-part" && i + 1 < argc)
		{
			opt->_parfile_name = argv[i + 1];
			i ++;
			continue;
		}else
		if(s == "-h")
		{
			print_usage();
			exit(0);
		}else
		{
			erro("wrong input options: " + s);
			print_usage();
			exit(-1);
		}
	}
	// validity checking
	if(opt->_ltlfile_name == nullptr )
	{
		erro( "missing LTLf file name");
		exit(-1);
	}
	if(opt->_synthesis && ( opt->_parfile_name == nullptr))
	{
		erro("missing proposition partition file name");
		exit(-1);
	}

}

void
read_from_part_file(const char *file_name, vector<string>& input, vector<string>& output)
{
    // const char * file_name
    ifstream part_file(file_name);
    if (part_file.is_open())
    {
        bool flag = false;
        string line;
        while(getline(part_file, line))
        {
            if(str_contain(line, "inputs"))
            {
                string delimiter = ":";
                line = line.substr(line.find(delimiter) + 1);
                str_split(line, input, ' ');
            }else
            if(str_contain(line, "outputs"))
            {
                string delimiter = ":";
                line = line.substr(line.find(delimiter) + 1);
                str_split(line, output, ' ');
            }else
            {
                cout << "read partfile error!" <<endl;
                cout << file_name <<endl;
                cout << line <<endl;
                exit(-1);
            }
        }
    }
}

std::string opToString(spot::op op)
{
    switch (op)
    {
        case spot::op::Not: return "!";
        case spot::op::And: return "&&";
        case spot::op::Or: return "||";
        default: return "";
    }
}

// formula syntax tree (with shared node for the same formula)
class Node;
typedef shared_ptr<Node> node_ptr;


class Node {

public:
	spot::op opName; // usually it is \/ or /\

	spot::formula formula;

	int refCount;

	std::set<node_ptr> children;
	std::set<node_ptr> parents;

	Node(const spot::op& op): opName(op), refCount(1) {}

	bool isLeaf() {
		return opName !=  op::And && opName != op::Or;
	}

	~Node() { children.clear(); parents.clear(); }

};


struct Stat {

	int unique_ast;
	int repeat_ast;
	int formula_size;
	// int repeat_i;

	Stat() : unique_ast(0), repeat_ast(0), formula_size(0) {}  

};

struct OpNums {
	int numDFAConversions;
	int numCompositions;

	OpNums(): numDFAConversions(0), numCompositions(0) {}

};

OpNums estimateOpNums(const spot::formula& f, bool  pushInside)
{
	OpNums currOpNums;

    if (f.kind() == op::And || f.kind() == op::Or)
    {
		int numChildren = f.size();
        for (formula child: f)
        {
            OpNums childOpNums = estimateOpNums(child, pushInside);
			currOpNums.numDFAConversions += childOpNums.numDFAConversions;
			currOpNums.numCompositions += childOpNums.numCompositions;
        }
		currOpNums.numCompositions += (numChildren - 1);
    } else if (pushInside && f.kind() == op::G && f[0].kind() == op::And)
    {
		int numChildren = f[0].size();
		for (auto child : f[0]) {
			formula gConjunct = formula::G(child);
			OpNums childOpNums = estimateOpNums(gConjunct, pushInside);
			currOpNums.numDFAConversions += childOpNums.numDFAConversions;
			currOpNums.numCompositions += childOpNums.numCompositions;
		}
		currOpNums.numCompositions += (numChildren - 1);
    } else if (pushInside && f.kind() == op::F && f[0].kind() == op::Or)
    {
		int numChildren = f[0].size();
		for (auto child : f[0]) {
			formula fDisjunct = formula::F(child);
			OpNums childOpNums = estimateOpNums(fDisjunct, pushInside);
			currOpNums.numDFAConversions += childOpNums.numDFAConversions;
			currOpNums.numCompositions += childOpNums.numCompositions;
		}
		currOpNums.numCompositions += (numChildren - 1);
    } else
    {
		currOpNums.numDFAConversions = 1;
    }

    return currOpNums;
}

int getASTSize(const spot::formula& f, bool pushInside)
{
	int size = 1;
    if (f.kind() == op::And || f.kind() == op::Or)
    {
		int numChildren = f.size();
        for (formula child: f)
        {
            int childSize = getASTSize(child, pushInside);
			size += childSize;
        }
    } else if (pushInside && f.kind() == op::G && f[0].kind() == op::And)
    {
		for (auto child : f[0]) {
			formula gConjunct = formula::G(child);
			int childSize = getASTSize(gConjunct, pushInside);
			size += childSize;
		}
    } else if (pushInside && f.kind() == op::F && f[0].kind() == op::Or)
    {
		for (auto child : f[0]) {
			formula fDisjunct = formula::F(child);
			int childSize = getASTSize(fDisjunct, pushInside);
			size += childSize;
		}
    } else
    {
		size = 1;
    }

    return size;
}


OpNums estimateOpNums(std::set<node_ptr>& computeTable, node_ptr root)
{
	OpNums currOpNums;

    if (computeTable.find(root) != computeTable.end()) {
		// no need to compute the child again
		// so get it for free
		return currOpNums;
    }

	if (root->isLeaf()) {
		currOpNums.numDFAConversions = 1;
	}else {
		// first do children nodes
		for (node_ptr child : root->children) {
			OpNums childOpNums = estimateOpNums(computeTable, child);
			currOpNums.numDFAConversions += childOpNums.numDFAConversions;
			currOpNums.numCompositions += childOpNums.numCompositions;
		}
		currOpNums.numCompositions += (root->children.size() - 1);
	}

	computeTable.insert(root);
    return currOpNums;
}


node_ptr createAST(std::map<formula, node_ptr>& nodeMap
, const spot::formula& f, Stat & stat)
{
	stat.formula_size ++;
	auto it = nodeMap.find(f);
    if (it != nodeMap.end()) {
		stat.repeat_ast ++;
		// cout << "one formula repeat" << endl;
		// ref count increase by 1
		it->second->refCount = it->second->refCount + 1;
        return it->second;
    }

	// cout << "formula: " << f << endl;
    node_ptr currNode = std::make_shared<Node>(f.kind());
	currNode->formula = f;
	// record this node
	nodeMap[f] = currNode;
    
    if (currNode->opName == op::And || currNode->opName == op::Or)
    {
        for (formula child: f)
        {
            node_ptr childNode = createAST(nodeMap, child, stat);
			childNode->parents.insert(currNode);
            currNode->children.insert(childNode);
        }
    } else if (currNode->opName == op::G && f[0].kind() == op::And)
    {
		currNode->opName = op::And;
		for (formula conjunct : f[0])
		{
			formula gConjunct = formula::G(conjunct);
			node_ptr childNode = createAST(nodeMap, gConjunct, stat);
			childNode->parents.insert(currNode);
			currNode->children.insert(childNode);
		}
    } else if (currNode->opName == op::F && f[0].kind() == op::Or)
    {
		currNode->opName = op::Or;
		for (formula disjunct : f[0])
		{
			formula fDisjunct = formula::F(disjunct);
			node_ptr childNode = createAST(nodeMap, fDisjunct, stat);
			childNode->parents.insert(currNode);
			currNode->children.insert(childNode);
		}
	// after negative normal form, no implies?
    } else
    {
		// we reach the leaf Node
		// There is nothing to do with this
		currNode->formula = f;
    }

    return currNode;
}

void printAST(node_ptr node, int indent = 0)
{
    if (node == nullptr)
        return;

    std::string indentStr(indent, ' ');

    if (node->opName == op::And || node->opName == op::Or )
    {
        std::cout << indentStr << "Operator: " << opToString(node->opName) << " id: " << node << std::endl;
    }
    else
    {
        std::cout << indentStr << "Formula: " << (node->formula) << " id: " << node << std::endl;
    }

    for (auto child : node->children)
        printAST(child, indent + 4);
}

void printAST(spot::formula node, int indent = 0)
{

    std::string indentStr(indent, ' ');

    if (node.kind() == op::And || node.kind() == op::Or )
    {
        std::cout << indentStr << "Operator: " << opToString(node.kind()) << std::endl;
		for (auto child : node)
        	printAST(child, indent + 4);
    }
    else if (node.kind() == op::F && node[0].kind() == op::Or)
    {
        std::cout << indentStr << "Operator: " << opToString(op::Or) << std::endl;
		for (auto child : node[0])
        	printAST(formula::F(child), indent + 4);
    }else if (node.kind() == op::G && node[0].kind() == op::And) {
		std::cout << indentStr << "Operator: " << opToString(op::And) << std::endl;
		for (auto child : node[0])
        	printAST(formula::G(child), indent + 4);
	}else {
		std::cout << indentStr << "Formula: " << str_psl(node, false) << std::endl;
	}
   
}


void reduceAST(std::map<formula, node_ptr>& nodeMap,
	node_ptr root, Stat& stat) {
	
	if (root == nullptr) {
		return;
	}

	// first do children nodes
	for (node_ptr child : root->children) {
		reduceAST(nodeMap, child, stat);
	}

	std::map<node_ptr, int> numOfOccurs;
	std::map<node_ptr, std::set<node_ptr>> parentSet;
	// obtain shared children, and compute the most shared grandChild 
	// among children
	for (node_ptr child : root->children) {
		if (! child->isLeaf()) {
			// only one step further
			for (node_ptr grandChild : child->children) {
				if (grandChild->refCount > 1) {
					numOfOccurs[grandChild] += 1;
					if (numOfOccurs[grandChild] == 1) {
						// first time
						std::set<node_ptr> parent;
						// add child to the parent set of grandChild
						parent.insert(child);
						parentSet[grandChild] = parent;
					}
					else {
						parentSet[grandChild].insert(child);
					}
				}
			}
		}
	}

	//  map iterate in ascending order of keys
	node_ptr maxGrandChild = nullptr;
	int maxOccurs = 1;
	auto it = numOfOccurs.begin();
	while (it != numOfOccurs.end()) {
		auto val = *it;
		if (val.second > maxOccurs) {
			maxGrandChild = val.first;
			maxOccurs = val.second;
		}
		it ++;
	}

	if (maxGrandChild != nullptr) {
		// all parents in the children of this child temp
		std::set<node_ptr> allChildParents = parentSet[maxGrandChild];

		// first parent of temp
		// mergeNode is the operator of the first parent
		node_ptr maxAndChildParentNode = std::make_shared<Node>((*allChildParents.begin())->opName);
		// operator of current Node
		node_ptr childrenParentNode = std::make_shared<Node>(root->opName);
		// there are equal number of children and parent?
		bool mergeToRoot = (allChildParents.size() == root->children.size());

		// create a new mergeNode and put this child on this node
		maxAndChildParentNode->children.insert(maxGrandChild);
		// also put a node for the operator of current node
		maxAndChildParentNode->children.insert(childrenParentNode);
		// add merged node to child of current node
		root->children.insert(maxAndChildParentNode);
		// put one of its child to child of mergeNode
		// even lower ?
		maxGrandChild->parents.insert(maxAndChildParentNode);
		// 
		childrenParentNode->parents.insert(maxAndChildParentNode);

		// stat.unique_ast -= (allChildParents.size() - 3);
		// all parents of the child, temp
		for (node_ptr childParent : allChildParents) {
			childParent->children.erase(maxGrandChild);
			childParent->parents.erase(root);
		
			root->children.erase(childParent);
			
			maxGrandChild->parents.erase(childParent);

			maxGrandChild->refCount--;

			if (childParent->children.size() > 1) {
				childrenParentNode->children.insert(childParent);
				childParent->parents.insert(childrenParentNode);
			} else {
				// there exists one child, lift the child up?
				(*(childParent->children.begin()))->parents.insert(childrenParentNode);
				childrenParentNode->children.insert(*(childParent->children.begin()));
			}
		}

		if (mergeToRoot) {
			root->opName = maxAndChildParentNode->opName;
			root->children.clear();
		    root->children.insert(maxAndChildParentNode->children.begin(), maxAndChildParentNode->children.end()) ;
			root->parents.clear();
			root->parents.insert(maxAndChildParentNode->parents.begin(), maxAndChildParentNode->parents.end()) ;
			// whether we should check we can lift this node up?
			// delete maxAndChildParentNode;
			root->refCount = maxAndChildParentNode->refCount;
		} else {
			maxAndChildParentNode->parents.insert(root);
		}
	}
	
}

// create a new AST
node_ptr removeDuplicates(std::map<formula, node_ptr>& nodeMap,
	node_ptr root, Stat& stat) {
	if (root == nullptr)
        return nullptr;
	spot::formula f;
	bool isLeaf = root->isLeaf();
	std::set<node_ptr> children;
	if (isLeaf) {
		// obtain the formula
		f = root->formula;
	}else {
		// traverse children first
		std::vector<formula> operands;
		for (node_ptr child: root->children) {
			node_ptr childNode = removeDuplicates(nodeMap, child, stat);
			operands.push_back(childNode->formula);
			children.insert(childNode);
		}
		f = formula::multop(root->opName, operands);
	}
	stat.formula_size ++;
	auto it = nodeMap.find(f);
    if (it != nodeMap.end()) {
		stat.repeat_ast ++;
		// cout << "one formula repeat" << endl;
		// ref count increase by 1
        return it->second;
    }else {
		// this is the first node ever created in the tree
		// we create a new node
		node_ptr newNode = std::make_shared<Node>(f.kind());
		newNode->formula = f;
		newNode->children.insert(children.begin(), children.end());
		// we also add parents to those nodes
		for (node_ptr childNode : children) {
			childNode->parents.insert(newNode);
		}
		nodeMap[f] = newNode;
		return newNode;
	}

}

struct GreaterThanBySize
{
  bool operator()(twa_graph_ptr& p1, twa_graph_ptr& p2) const
  {
    return p1->num_states() >= p2->num_states();
  }
};

bool isSinkSpot(twa_graph_ptr automaton, bool positive) {
  bool has_one_state = automaton->num_states() == 1;
  // only have sink loop
  // init - l -> init, init -!l -> loop 
  int init_state = automaton->get_init_state_number();

  bool is_initial_state_final = true;
  return has_one_state && (is_initial_state_final == positive);
}

// twa_graph_ptr composeOr(bdd_dict_ptr dict,
// 	std::map<node_ptr, twa_graph_ptr>& map, node_ptr node);
// twa_graph_ptr composeAnd(bdd_dict_ptr dict,
// 	std::map<node_ptr, twa_graph_ptr>& map, node_ptr node);

twa_graph_ptr composeSpot(bdd_dict_ptr dict, spot::op opType, 
	std::map<node_ptr, twa_graph_ptr>& map, node_ptr node);
twa_graph_ptr composeAST(bdd_dict_ptr dict,
	std::map<node_ptr, twa_graph_ptr>& map, node_ptr node);

twa_graph_ptr composeSpot(bdd_dict_ptr dict, spot::op opType, 
	std::map<node_ptr, twa_graph_ptr>& map, node_ptr node) {

	priority_queue<twa_graph_ptr, std::vector<twa_graph_ptr>, GreaterThanBySize> childrenSorted;
	for (node_ptr child : node->children) {
		twa_graph_ptr dfa = composeAST(dict, map, child);
		childrenSorted.push(dfa);
	}

	while (childrenSorted.size() > 1) {
		twa_graph_ptr first = childrenSorted.top();
        childrenSorted.pop();
        twa_graph_ptr second = childrenSorted.top();
        childrenSorted.pop();
		// cout << "A: " << first->num_states() << endl;
		// cout << "B: " << second->num_states() << endl;
		twa_graph_ptr P ;
		if (opType == spot::op::And) {
			P = spot::product(first, second);
		}else {
		  P = spot::product_or(first, second);
		}
		if (opt->_verbose) {
          		cout << "A: " << first->num_states() << " B: " << second->num_states() << " product: " << P->num_states();
          	}
		P = minimize_explicit(P);	
		// cout << "P: " << P->num_states() << endl;
		if (opt->_verbose) {
			 cout << " min_product: " << P->num_states() << endl;
		}
	
		childrenSorted.push(P);
	}

	return childrenSorted.top();
}

twa_graph_ptr composeAST(bdd_dict_ptr dict,
	std::map<node_ptr, twa_graph_ptr>& map, node_ptr node)
{
	// cout << "entering : " << node->formula << endl;
	auto it = map.find(node);
    if (it != map.end()) {
        return it->second;
    }
	twa_graph_ptr res = nullptr;
	if (node->opName == spot::op::And
	|| node->opName == spot::op::Or)
	{
		res = composeSpot(dict, node->opName, map, node);
	}else if (node->isLeaf()){
		// cout << "Start translating... " << str_psl(node->formula, true) << endl;
		res = trans_formula((node->formula), dict, 7);
		// cout << "States node: " << res->num_states() << endl;
		res = minimize_explicit(res);
		// dict->register_all_variables_of(res, node);
		// cout << "Leaf node: " << res->num_states() << endl;
	}else {
		cout << "Should not enter here" << endl;
		exit(-1);
	}
	map[node] = res;
	return res;
} 

class MonaDFA {
public:
	std::vector<string> props;
	DFA* dfa;

	MonaDFA(std::vector<string> prop, DFA* d) {
		props = prop; // copy
		dfa = d;
	}

	~MonaDFA() {
		dfaFree(dfa);
	}

	DFA* get() {
		return dfa;
	}

	int numStates() {
		return dfa->ns;
	}
};

typedef shared_ptr<MonaDFA> dfa_ptr;

bool isSink(dfa_ptr automaton, bool positive) {
  bool has_one_state = automaton->numStates() == 1;
  bool is_initial_state_final = automaton->get()->f[0] == 1;
  return has_one_state && (is_initial_state_final == positive);
}

// dfa_ptr composeOrMona(char** vars,
	// std::map<node_ptr, dfa_ptr>& map, node_ptr node);
// dfa_ptr composeMona(std::map<dfa_ptr, std::vector<string>> namesMap,  dfaProductType opType,
// 	std::map<node_ptr, dfa_ptr>& map, node_ptr node);
dfa_ptr composeASTMona(std::map<string, int>& namesMap, std::map<node_ptr, dfa_ptr>& map, node_ptr node);

struct GreaterThanByMonaSize
{
  bool operator()(dfa_ptr& d1, dfa_ptr& d2) const
  {
    return d1->numStates() >= d2->numStates();
  }
};

dfa_ptr composeMona(dfaProductType opType, std::map<string, int>& namesMap
, std::map<node_ptr, dfa_ptr>& map, node_ptr node) {
	// cout << "Compose node or: " << endl;
	// auto cmp = [](const dfa_ptr d1, const dfa_ptr d2) {return d1->ns > d2->ns; };
    std::priority_queue<dfa_ptr, std::vector<dfa_ptr>, GreaterThanByMonaSize> 
                queue;
	for (node_ptr child : node->children) {
		dfa_ptr dfa = composeASTMona(namesMap, map, child);
		queue.push(dfa);
	}

	bool positive = (opType == dfaProductType::dfaOR);

        while (queue.size() > 1) {
			// cout << "queue size: " << queue.size() << endl;
          dfa_ptr lhs = queue.top();
          queue.pop();
          dfa_ptr rhs = queue.top();
          queue.pop();
		//   cout << "A: " << lhs->numStates() << " B: " << rhs->numStates() << endl;
		  DFA* prod = dfaProduct(lhs->get(), rhs->get(), opType);
		// cout << "prod: " << prod->ns << endl;

		  DFA* min = dfaMinimize(prod);
		  if (opt->_verbose) {
          		cout << "A: " << lhs->numStates() << " B: " << rhs->numStates() << " product: " << prod->ns << " min_product: " << min->ns << endl;
          	}
		// cout << "min: " << min->ns << endl;
		  dfaFree(prod);
		  std::vector<string> str;
          dfa_ptr dfa = std::make_shared<MonaDFA>(str, min);
		   
		  if (isSink(dfa, positive)) {
      		return dfa;
    	  }
          queue.push(dfa);
        }
    // we get the final one
	// map[node] = queue.top();
    return queue.top();
}

dfa_ptr composeASTMona(std::map<string, int>& namesMap, 
	std::map<node_ptr, dfa_ptr>& map, node_ptr node)
{
	// cout << "entering : " << node->formula << endl;
	auto it = map.find(node);
    if (it != map.end()) {
        return it->second;
    }
	dfa_ptr res = nullptr;
	if (node->opName == spot::op::And)
	{
		res = composeMona(dfaProductType::dfaAND, namesMap, map, node);
	}
	else if (node->opName == spot::op::Or)
	{
		res = composeMona(dfaProductType::dfaOR, namesMap, map, node);
	}else if (node->isLeaf()){
		formula f = node->formula;
		string mona_file_name = "./ltlf.mona" + get_current_time_string();
    	ofstream ofs(mona_file_name, ofstream::out);
    	ofs << "#LTLf formula" << endl;
    	ofs << "#" << str_psl(f, true) << endl;
    	formula bnf = get_bnf(f);
    	ofs << "# Backus normal form" << endl;
	    ofs << "#" << str_psl(bnf, true) << endl;
    	// the BNF form, and then convert it to fol formula
    	ltlf_to_fol(ofs, bnf);
    	ofs.close();
    	string dfa_file_name = "./mona.dfa" + get_current_time_string();
    	string command = "mona -u -xw " + mona_file_name+ " >" + dfa_file_name;
    	int r = system(command.c_str());
		int numAps = get_size_formula_ap(f);
		char ** vars;
		std::vector<char> filename_cstr(dfa_file_name.c_str(),
                                  dfa_file_name.c_str() + dfa_file_name.size() + 1);
		// char* monaFile = "./mona.dfa";
		DFA* dfa = dfaImport(filename_cstr.data(), &vars, 0);
		// cout << "dfa size: " << dfa->ns << endl;
		// dict->register_all_variables_of(res, node);
		int len = 0;
		std::vector<string> propNames;
		while (len < numAps) {
			propNames.push_back(*(vars + len));
			// cout << "length " << len << endl;
			len ++;
		}
		// cout << "prop size: " << propNames.size() << endl;
		int map[namesMap.size()];
		for(int i = 0; i < numAps; i ++) {
			map[i] = namesMap[propNames[i]];
			// cout << " name " << propNames[i] << endl;
			// cout << " old " << i << " new " <<  namesMap[propNames[i]] << endl;
        	free(*(vars + i));
    	}
		dfaReplaceIndices(dfa, map);
		dfa_ptr sdfa ( new MonaDFA(propNames, dfa));
		// cout << "success free: " << propNames.size() << endl;
		res = sdfa;
		std::vector<char> mona_cstr(mona_file_name.c_str(),
                                  mona_file_name.c_str() + mona_file_name.size() + 1);
		std::remove(mona_cstr.data());
		std::remove(filename_cstr.data());
		// we know length 
		// cout << "Leaf node: " << res->numStates() << endl;
		// cout << "Prop num: " << len << endl;
	}else {
		cout << "Should not enter here" << endl;
		exit(-1);
	}
	map[node] = res;
	return res;
} 

twa_graph_ptr fromMonaDFA(dfa_ptr dfa, const std::vector<string>& names) {
	//DFA *a, char *filename, int num, char *names[], char orders[]
	// dfaExport(dfa->get(), "mona.dfa", )
  spot::bdd_dict_ptr dict = spot::make_bdd_dict();

  twa_graph_ptr aut = spot::make_twa_graph(dict);
  std::map<int, bdd> bdd_props;
  int numStates = dfa->numStates();
  
    // get bdd repr for propositions
    for(int i = 0; i < names.size(); i ++)
    {
		int index = aut->register_ap(names[i]);
        bdd p = bdd_ithvar(index);
        bdd_props[index] = p;
    }
    // add another state
	int index = aut->register_ap(ALIVE_AP);
    bdd alive = bdd_ithvar(index);
    bdd_props[index] = alive;
    
    aut->set_buchi();
    aut->prop_state_acc();
    // add one extra state for accepting state
    aut->new_states(numStates + 1);

  paths state_paths, pp;
  trace_descr tp;
  int i, j, any = 0;

	DFA* a = dfa->get();
//   printf("DFA for formula with free variables: ");

//   for (i = 0; i < no_free_vars; i++)
//     printf ("%s ", free_variables[i]);
  aut->set_init_state(a->s);
//   printf("\nInitial state: %d\n", a->s);
//   printf("Accepting states: ");
//   for (i = 0; i < a->ns; i++)
    // if (a->f[i] == 1)
    //   printf ("%d ", i);

//   printf("\n");

//   printf ("Transitions:\n");

  for (i = 0; i < a->ns; i++) {
    state_paths = pp = make_paths(a->bddm, a->q[i]);
    if (a->f[i] == 1) {
      aut->new_edge(i, numStates, !alive);
    }
    while (pp) {
      //   printf ("State %d: ", i);
      bdd tr = bddtrue;
      for (tp = pp->trace; tp; tp = tp->next) {

        if (tp->value) {
          tr &= bdd_props[tp->index];
          // printf("1");
        } else {
          tr &= !bdd_props[tp->index];
          // printf("0");
        }
      }

      // printf (" -> state %d\n", pp->to);
      aut->new_edge(i, pp->to, tr & alive);
      pp = pp->next;
    }
	kill_paths(state_paths);
  }
	// now set accepting states
    aut->new_edge(numStates, numStates, !alive, {0});
    aut->merge_edges();
	return aut;
}

void outputMonaDFA(string filename, dfa_ptr dfa, const std::vector<string>& names) {
	std::vector<char> filename_cstr(filename.c_str(),
                                  filename.c_str() + filename.size() + 1);
  char** arr = new char*[names.size()];
  for (size_t i = 0; i < names.size(); i++) {
    arr[i] = new char[names[i].size() + 1];
    strcpy(arr[i], names[i].c_str());
  }

  char* orders = new char[names.size()];
  memset(orders, 2, names.size());
//   cout << "Hekk" << endl;
  dfaExport(dfa->get(), filename_cstr.data(), names.size(), arr, orders);
//   cout << "Hekk" << endl;
  for (size_t i = 0; i < names.size(); i++) {
    delete[] arr[i];
  }
  delete[] arr;
}

void synthesisMonaDFA(dfa_ptr dfa, const std::vector<string>& names, string partfile) {
	string filename = "mona.dfa";


    // // remember the nullptr terminator
    // result.reserve(names.size());

    // std::transform(begin(names), end(names),
    //                std::back_inserter(result),
    //                [](std::string &s) { return s.data(); }
    //               );
	// dfaExport(dfa->get(), "mona.dfa", names.size(), vars, 0);
	outputMonaDFA(filename, dfa, names);
	cout << "Hekk" << endl;
	clock_t syn_start = clock();
	// syft(partfile, filename, "1");
	clock_t syn_end = clock();
	cout << "Time spent by Syft: " << 1000.0 * (syn_end - syn_start)/CLOCKS_PER_SEC << "ms ..." << endl;

    // return read_from_mona_file("mona.dfa", dict);
}

int main(int argc, char** argv)
{
    opt_t o;
    opt = &o;
    parse_opt(argc, argv);
	
    ifstream ltlfile;
	ltlfile.open(opt->_ltlfile_name);
    string line;
    clock_t c_start = clock();
    formula input_f;
	// cout << ltlfile << endl;
	cout << "Starting the decomposition phase" << endl;

    getline(ltlfile, line);
    //cout << "formula: " << line << endl;
    auto pf1 = spot::parse_infix_psl(line.c_str());
    if (pf1.format_errors(std::cerr))
    {
        std::cerr << "Error: " << line << std::endl;
        return -1;
    }
    // formula 
    input_f = pf1.f;


	// cout << "formula: " << input_f << endl;
	input_f = get_nnf(input_f);
	// cout << "formula: " <<  input_f << endl;
	std::map<formula, node_ptr> nodeMap;
	Stat stat;
	cout << "Creating DAG for the formula" << endl;
    node_ptr dag = createAST(nodeMap, input_f, stat);
	cout << "Nodes in syntax tree: " << getASTSize(input_f, false) << endl;
	OpNums inputNums = estimateOpNums(input_f, false);
	//currOpNums.numDFAConversions += numChildren;
	//	currOpNums.numCompositions += (numChildren - 1);
	cout << "Estimated number of DFA conversion operaions in syntax tree: " << inputNums.numDFAConversions << endl;
	cout << "Estimated number of composition operaions in syntax tree: " << inputNums.numCompositions << endl;
	cout << "Nodes in pushed syntax tree: " << getASTSize(input_f, true) << endl;
	inputNums = estimateOpNums(input_f, true);
	cout << "Estimated number of DFA conversion operaions in pushed syntax tree: " << inputNums.numDFAConversions << endl;
	cout << "Estimated number of composition operaions in pushed syntax tree: " << inputNums.numCompositions << endl;
	cout << "Nodes in DAG: " << (stat.formula_size - stat.repeat_ast) << endl;
	std::set<node_ptr> computedTable;
	inputNums = estimateOpNums(computedTable, dag);
	cout << "Estimated number of DFA conversion operaions in DAG: " << inputNums.numDFAConversions << endl;
	cout << "Estimated number of composition operaions in DAG: " << inputNums.numCompositions << endl;

	if (opt->_dag) {
		cout << "Reducing DAG for the formula" << endl;
		reduceAST(nodeMap, dag, stat);
		std::map<formula, node_ptr> newNodeMap;
		Stat newStat;
		dag = removeDuplicates(newNodeMap, dag, newStat);
		cout << "Nodes in the reduced DAG: " << (newStat.formula_size - newStat.repeat_ast) << endl;
		computedTable.clear();
		inputNums = estimateOpNums(computedTable, dag);
		cout << "Estimated number of DFA conversion operaions in reduced DAG: " << inputNums.numDFAConversions << endl;
		cout << "Estimated number of composition operaions in reduced DAG: " << inputNums.numCompositions << endl;
	}

	if (opt->_verbose) {
		cout << "Pushed syntax tree:" << endl;
		printAST(input_f, 0);

		cout << "Reduced DAG:" << endl;
		printAST(dag, 0);
	}
    
	// printAST(newTree);
	clock_t c_end = clock();
	cout << "Time spent in decomposition: " << 1000.0 * (c_end - c_start)/CLOCKS_PER_SEC << "ms ..." << endl;

	int numStates = -1;
	clock_t comp_start = clock();
	spot::bdd_dict_ptr dict = spot::make_bdd_dict();

	twa_graph_ptr dfa;
	dfa_ptr monaDfa;
	std::vector<string> names;
	if (!opt->_mona) {
		std::map<node_ptr, twa_graph_ptr> dfaTable;
		dfa = composeAST(dict, dfaTable, dag);
		numStates = dfa->num_states();
	}else {
		std::map<string, int> namesMap;
        std::map<node_ptr, dfa_ptr> monaMap;
        std::set<formula> aps;
        get_formula_aps(input_f, aps);
        // std::vector<std::string> names;
        names.reserve(aps.size());
        for (auto ap : aps) {
            names.push_back(ap.ap_name());
            // cout << " " << str_psl(ap, false) << endl;
        }
        std::sort(names.begin(), names.end());
        int i = 0;
        for (auto name : names) {
            // cout << " name: " << name << " index: " << i << endl;
            namesMap[name] = i++;
        }
        monaDfa = composeASTMona(namesMap, monaMap, dag);
		numStates = monaDfa->numStates();
	}
	clock_t comp_end = clock();
	cout << "Time spent in composition phase: " << 1000.0 * (comp_end - comp_start)/CLOCKS_PER_SEC << "ms ..." << endl;
	cout << "States in minmal DFA: " << numStates << endl;
	c_end = clock();
	// cout << "Breakdown: " << sizes << endl;
	if (opt->_outfile_name != nullptr) {
		if (opt->_mona) {
			outputMonaDFA(opt->_outfile_name, monaDfa, names);
		}else {
			// output
    		ofstream outfile(opt->_outfile_name);
    		print_hoa(outfile, dfa);
		}
	}
	cout << "Total runtime for DFA construction: " << 1000.0 * (c_end - c_start)/CLOCKS_PER_SEC << "ms ..." << endl;
	
	// synthesis 
	// synthesisW
	// if (!opt->_synthesis) {
	// 	return 0;
	// }
	// if (opt->_mona) {
	// 	// dfa = fromMonaDFA(monaDfa, names);
	// 	comp_start = clock();
	// 	string filename = "mona.dfa";
	// 	outputMonaDFA(filename, monaDfa, names);
	// 	dfa = read_from_mona_file(filename.c_str(), dict);
	// 	comp_end = clock();
	// 	cout << "Time spent in DFA conversion: " << 1000.0 * (comp_end - comp_start)/CLOCKS_PER_SEC << "ms ..." << endl;
	// }
	// if (opt->_mona) {
	// 	// convert to spot DFA
	// 	cout << "Synthesis by Syft..." << endl;
	// 	comp_start = clock();
	// 	synthesisMonaDFA(monaDfa, names, opt->_parfile_name);
	// 	comp_end = clock();
	// 	cout << "Time spent in synthesis: " << 1000.0 * (comp_end - comp_start)/CLOCKS_PER_SEC << "ms ..." << endl;
	// }else {
		// vector<string> input;
		// vector<string> output;
		// // cout << "read part file " << endl;
		// if(opt->_parfile_name != nullptr)
		// {
		// 	read_from_part_file(opt->_parfile_name, input, output);
		// }else
		// {
		// 	cerr << "Please input the file name for inputs and outputs" << endl;
		// 	exit(-1);
		// }

		// output.push_back(ALIVE_AP);
		// solve_game(dfa, input, output);
	// }
	
	opt = nullptr;
	// c_end = clock();
	// cout << "Breakdown: " << sizes << endl;
	// cout << "Total runtime for synthesis: " << 1000.0 * (c_end - c_start)/CLOCKS_PER_SEC << "ms ..." << endl;
}
