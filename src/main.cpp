/*
 Arjun
 Copyright (c) 2019, Mate Soos and Kuldeep S. Meel. All rights reserved.
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 */
#include <fstream>
#include <boost/program_options.hpp>

using boost::lexical_cast;
namespace po = boost::program_options;
using std::string;
using std::vector;

#if defined(__GNUC__) && defined(__linux__)
#include <fenv.h>
#endif

#include <iostream>
#include <iomanip>
#include <random>
#include <algorithm>
#include <map>
#include <set>
#include <vector>
#include <atomic>
#include <signal.h>


#include "time_mem.h"
//#include "GitSHA1.h"
#include "MersenneTwister.h"

#include <cryptominisat5/cryptominisat.h>
#include "cryptominisat5/dimacsparser.h"
#include "cryptominisat5/streambuffer.h"

using namespace CMSat;
using std::cout;
using std::cerr;
using std::endl;
using std::map;
using std::set;
using std::fstream;
using std::ios;


typedef vector<int> intvec;

po::options_description mis_options = po::options_description("MIS options");
po::options_description help_options;
po::variables_map vm;
po::positional_options_description p;
string command_line;
CMSat::SATSolver* solver = NULL;
double startTime;
vector<Lit> tmp;
vector<char> seen;
uint32_t orig_num_vars;
int32_t varout;  
enum ModeType {one_mode, many_mode};

//assert indic[var] to FASLE to force var==var+orig_num_vars
vector<uint32_t> var_to_indic; //maps an ORIG VAR to an INDICATOR VAR
vector<uint32_t> indic_to_var; //maps an INDICATOR VAR to ORIG VAR
map<int32_t, int32_t>  Yvar_to_Ypvar; //maps for Y VAR to Y' VAR in F(X,Y) \and \lnot F(X,Y')
map<int32_t, int32_t>  Yvar_to_Cvar; // maps for Y var to C var for c-> y xor y'
map<int32_t, int32_t>  tempvar_to_temp2var; // maps for allvariable except X and Y var to temp2  var for ~F(X,Y)
vector<vector<int32_t>> clauses_list; // clause lits
vector<uint32_t> sampling_set_tmp1;
vector<uint32_t>* sampling_set = NULL;
vector<uint32_t> x_vars;
vector<uint32_t> this_indic2;
vector <Lit> assumptions;
vector<uint32_t> postive_unate;
vector<uint32_t> negative_unate;
vector<uint32_t> zvariables;
string posunate_details;
string negunate_details;


struct Config {
    int verb = 0;
    int seed = 0;
    int bva = 0;
    int bve = 1;
    int guess = 1;
    int force_by_one = 1;
    int simp_at_start = 1;
    int always_one_by_one = 1;
    int recompute_sampling_set = 0;
    int backward_only = 0;
};
Config conf;
MTRand mtrand;

void print_indep_set()
{
    cout << "c ind ";
    for(const uint32_t s: *sampling_set) {
        cout << s+1 << " ";
    }
    cout << "0" << endl;

//     cout << "c inc ";
//     for(const uint32_t s: *sampling_set) {
//         cout << incidence[s] << " ";
//     }
//     cout << "0" << endl;

    cout << "c set size: " << std::setw(8)
    << sampling_set->size()
    << " fraction of original: "
    <<  std::setw(6) << std::setprecision(4)
    << (double)sampling_set->size()/(double)orig_num_vars
    << endl << std::flush;
}


void add_mis_options()
{
    std::ostringstream my_epsilon;
    std::ostringstream my_delta;
    std::ostringstream my_kappa;

    mis_options.add_options()
    ("help,h", "Prints help")
    ("version", "Print version info")
    ("input", po::value<string>(), "file to read")
    ("verb,v", po::value(&conf.verb)->default_value(conf.verb), "verbosity")
    ("seed,s", po::value(&conf.seed)->default_value(conf.seed), "Seed")
    ("bva", po::value(&conf.bva)->default_value(conf.bva), "bva")
    ("bve", po::value(&conf.bve)->default_value(conf.bve), "bve")
    ("guess", po::value(&conf.guess)->default_value(conf.guess), "Guess small set")
    ("one", po::value(&conf.always_one_by_one)->default_value(conf.always_one_by_one), "always one-by-one mode")
    ("simpstart", po::value(&conf.simp_at_start)->default_value(conf.simp_at_start), "simp at startup")
    ("recomp", po::value(&conf.recompute_sampling_set)->default_value(conf.recompute_sampling_set), "Recompute sampling set even if it's part of the CNF")
    ("byforce", po::value(&conf.force_by_one)->default_value(conf.force_by_one), "Force 1-by-1 query")
    ("backwardonly", po::value(&conf.backward_only)->default_value(conf.backward_only), "Only do backwards query")

    ;

    help_options.add(mis_options);
}

void add_supported_options(int argc, char** argv)
{
    add_mis_options();
    p.add("input", 1);

    try {
        po::store(po::command_line_parser(argc, argv).options(help_options).positional(p).run(), vm);
        if (vm.count("help"))
        {
            cout
            << "Probably Approximate counter" << endl;

            cout
            << "approxmc [options] inputfile" << endl << endl;

            cout << help_options << endl;
            std::exit(0);
        }

        if (vm.count("version")) {
            //cout << "[mis] Version: " << get_version_sha1() << endl;
            std::exit(0);
        }

        po::notify(vm);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::unknown_option> >& c
    ) {
        cerr
        << "ERROR: Some option you gave was wrong. Please give '--help' to get help" << endl
        << "       Unkown option: " << c.what() << endl;
        std::exit(-1);
    } catch (boost::bad_any_cast &e) {
        std::cerr
        << "ERROR! You probably gave a wrong argument type" << endl
        << "       Bad cast: " << e.what()
        << endl;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::invalid_option_value> >& what
    ) {
        cerr
        << "ERROR: Invalid value '" << what.what() << "'" << endl
        << "       given to option '" << what.get_option_name() << "'"
        << endl;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::multiple_occurrences> >& what
    ) {
        cerr
        << "ERROR: " << what.what() << " of option '"
        << what.get_option_name() << "'"
        << endl;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::required_option> >& what
    ) {
        cerr
        << "ERROR: You forgot to give a required option '"
        << what.get_option_name() << "'"
        << endl;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::too_many_positional_options_error> >& what
    ) {
        cerr
        << "ERROR: You gave too many positional arguments. Only the input CNF can be given as a positional option." << endl;
        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::ambiguous_option> >& what
    ) {
        cerr
        << "ERROR: The option you gave was not fully written and matches" << endl
        << "       more than one option. Please give the full option name." << endl
        << "       The option you gave: '" << what.get_option_name() << "'" <<endl
        << "       The alternatives are: ";
        for(size_t i = 0; i < what.alternatives().size(); i++) {
            cout << what.alternatives()[i];
            if (i+1 < what.alternatives().size()) {
                cout << ", ";
            }
        }
        cout << endl;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::invalid_command_line_syntax> >& what
    ) {
        cerr
        << "ERROR: The option you gave is missing the argument or the" << endl
        << "       argument is given with space between the equal sign." << endl
        << "       detailed error message: " << what.what() << endl
        ;
        std::exit(-1);
    }
}

void readInAFile(const string& filename, uint32_t var_offset, bool get_sampling_set)
{
    #ifndef USE_ZLIB
    FILE * in = fopen(filename.c_str(), "rb");
    DimacsParser<StreamBuffer<FILE*, FN> > parser(solver, NULL, 0);
    #else
    gzFile in = gzopen(filename.c_str(), "rb");
    DimacsParser<StreamBuffer<gzFile, GZ> > parser(solver, NULL, 0);
    #endif

    if (in == NULL) {
        std::cerr
        << "ERROR! Could not open file '"
        << filename
        << "' for reading: " << strerror(errno) << endl;

        std::exit(-1);
    }

    if (!parser.parse_DIMACS(in, false, var_offset)) {
        exit(-1);
    }
    if (get_sampling_set) {
        *sampling_set = parser.sampling_vars;
    }
    clauses_list = parser.clauses_list;
    x_vars = parser.x_vars;
    #ifndef USE_ZLIB
        fclose(in);
    #else
        gzclose(in);
    #endif
}

void init_samping_set(bool recompute)
{
  
    if (sampling_set->size() > 50) {
        cout
        << "[mis] Sampling var set contains over 100 variables, not displaying"
        << endl;
    } else {
      
        cout << "[mis] Y set : ";
        for (auto v: *sampling_set) {
            cout << v+1 << ", ";
        }
        cout << endl;
    }
    
    cout << "[mis] Orig size Y vars   : " << sampling_set->size() << endl;
    
    if (x_vars.size() > 100) {
        cout
        << "[mis] x var set contains over 100 variables, not displaying"
        << endl;
    } else {
      
        cout << "[mis] x var set: ";
        for (auto v: x_vars) {
            cout << v+1 << ", ";
        }
        cout << endl;
    }
    
    cout << "[mis] Orig size X vars  : " << x_vars.size() << endl;
    
    return;
}

void add_fixed_clauses()
{
    vector<uint32_t>::iterator find_var;
    var_to_indic.resize(orig_num_vars, var_Undef);
    indic_to_var.resize(solver->nVars(), var_Undef);
    int32_t var;
    uint32_t this_indic;
    uint32_t this_indic_2;
    vector <Lit> vec_dnf_var;
    solver->new_var();
    this_indic=solver-> nVars()-1;
    solver->new_vars(2*sampling_set->size());
    //for(var: *sampling_set) {
    for (auto v = sampling_set->begin(); v != sampling_set->end() ; ++v)
    {
      var = *v;
	  Yvar_to_Ypvar[var]=this_indic;
	  this_indic_2=this_indic+sampling_set->size();
	  Yvar_to_Cvar[var]=this_indic_2;
	  if (conf.verb){
	    cout<<"for Yvar to Ypvar "<<var+1<<" : "<<this_indic+1<<endl;}
	  //c variable for c-> y xor y'
	   if (conf.verb){
	      cout<<"for Yvar to Cvar "<<var+1<<" : "<<this_indic_2+1<<endl;}
	    tmp.clear();
	    tmp.push_back(Lit(Yvar_to_Cvar[var],true));
	    tmp.push_back(Lit(var,false));
	    tmp.push_back(Lit(Yvar_to_Ypvar[var],true));
	    solver->add_clause(tmp);
	    if (conf.verb){
	      cout<<"adding clause : "<<tmp << endl;}
	    tmp.clear();
	    tmp.push_back(Lit( Yvar_to_Cvar[var],true));
	    tmp.push_back(Lit(var,true));
	    tmp.push_back(Lit( Yvar_to_Ypvar[var],false));
	    solver->add_clause(tmp);
	    if (conf.verb){
	      cout<<"adding clause : "<<tmp << endl;} 
	    this_indic = this_indic+1;
	    
    }
    if (conf.verb){cout<<"adding clauses for ~F(X,Y)"<< endl;}
    for(vector<intvec>::iterator it = clauses_list.begin(); it != clauses_list.end(); ++it)
    {
        //it is now a pointer to a intvec
    solver->new_var();
    this_indic=solver-> nVars()-1;
	tmp.clear();
    tmp.push_back(Lit(this_indic,false));
    zvariables.push_back(this_indic);
        for(intvec::iterator jt = it->begin(); jt != it->end(); ++jt)
        {
	    
	    
	    // jt is now a pointer to an integer.
	    
	    var=std::abs(*jt)-1;
	    if (Yvar_to_Ypvar.count(var)>0) {
		      if (*jt>0){
			tmp.push_back(Lit(Yvar_to_Ypvar[var],false));}
		      else {
			tmp.push_back(Lit(Yvar_to_Ypvar[var],true));} }
	    else {
	      find_var = find(x_vars.begin(), x_vars.end(), var);
	      if(find_var!=x_vars.end())
	      {
		    if (*jt>0){
		    tmp.push_back(Lit(var,false));}
		  else {
		    tmp.push_back(Lit(var,true));}
	      }
		
	    
        }}
        solver->add_clause(tmp); 
	if (conf.verb){cout<<"adding clause : "<<tmp <<endl;}
    }
    tmp.clear();
    for(uint32_t i = 0; i < zvariables.size() ; i++) {
        tmp.push_back(Lit(zvariables[i],true));

    }
    solver->add_clause(tmp); 
    if (conf.verb){cout<<"adding clause : "<<tmp <<endl;}
//exit(1);
    
}








void find_unate()
{
 lbool ret = l_Undef;
   for( auto testvar= sampling_set->rbegin(); testvar != sampling_set->rend(); ++testvar)
    {
      auto var=*testvar;
      assumptions.push_back(Lit(Yvar_to_Cvar[var], false));
    }
    for (auto var: *sampling_set) 
    {
     
      if (conf.verb){
      cout<<"test var "<<var+1<< endl;}
      assumptions.pop_back();
      assumptions.push_back(Lit(Yvar_to_Ypvar[var], false));
      assumptions.push_back(Lit(var, true));
      if (conf.verb){cout<<"assumptions : "<<assumptions<<endl;}
      solver->set_max_confl(50);
      solver->set_no_confl_needed();
      ret = solver->solve(&assumptions);
      if (conf.verb==1){
	cout<<"conflit details:"<<solver->get_last_conflicts()<<endl;
	cout<<"last decision:"<<solver->get_last_decisions()<<endl;
	cout<<"get_last_propagations:"<<solver->get_last_propagations()<<endl;
	
      }
      if(conf.verb==2){ 
	cout <<" v "; 
	for (uint32_t var = 0; var < solver->nVars(); var++) {
	    cout<<((solver->get_model()[var] == l_False)? "-": "" )<<var+1<<" ";}
	cout << " 0"<<endl;}
      assumptions.pop_back();
      assumptions.pop_back();
      if (ret == l_False) {
	  postive_unate.push_back(var);
	  tmp.clear();
	  tmp.push_back(Lit(var,false));
	  solver->add_clause(tmp);
	  if (conf.verb){cout<<"positive unate : added clause : "<<tmp<<endl;}
	  tmp.clear();
	  tmp.push_back(Lit(Yvar_to_Ypvar[var],false));
	  if (conf.verb){cout<<"positive unate : added clause : "<<tmp<<endl;}
	  solver->add_clause(tmp);
      }
      else {
	assumptions.push_back(Lit(Yvar_to_Ypvar[var], true));
	assumptions.push_back(Lit(var, false));
	if (conf.verb){cout<<" Not positive unate .. Searching for negative unate"<<endl;
	cout<<"assumptions : "<<assumptions<<endl;}
	solver->set_max_confl(50);
	solver->set_no_confl_needed();
	lbool ret = solver->solve(&assumptions);
	assumptions.pop_back();
	assumptions.pop_back();
	if(conf.verb==2){
	  cout<<ret<<endl;
	  cout <<" v "; 
	  for (uint32_t var = 0; var < solver->nVars(); var++) {
	      cout<<((solver->get_model()[var] == l_False)? "-": "" )<<var+1<<" ";}
	  cout << " 0"<<endl;}
	if (ret == l_False) { 
	  negative_unate.push_back(var);
	  tmp.clear();
	  tmp.push_back(Lit(var,true));
	  solver->add_clause(tmp);
	  if (conf.verb){cout<<"negative unate : added clause : "<<tmp<<endl;}
	  tmp.clear();
	  tmp.push_back(Lit(Yvar_to_Ypvar[var],true));
	  if (conf.verb){cout<<"negative unate : added clause : "<<tmp<<endl;}
	  solver->add_clause(tmp);
	}
	else{
	  tmp.clear();
	  tmp.push_back(Lit(Yvar_to_Cvar[var],false));
	  if (conf.verb){cout<<"no unate : added clause : "<<tmp<<endl;}
	  solver->add_clause(tmp);} 
      }
    }
   cout<<"no positive unates "<<postive_unate.size()<<endl;
   cout<<"no. negative_unate "<<negative_unate.size()<<endl;
   for(uint32_t i = 0; i < postive_unate.size() ; i++) {
     // cout<<"positive unate "<<postive_unate[i]+1<<endl;
     posunate_details += std::to_string(postive_unate[i]+1)+" ";
     }
    for(uint32_t i = 0; i < negative_unate.size() ; i++) {
     // cout<<"negative_unate "<<negative_unate[i]+1<<endl;
       negunate_details += std::to_string(negative_unate[i]+1)+" ";
    }
}

void init_solver_setup()
{
 
  double myTime = cpuTime();
    solver = new SATSolver();
    if (conf.verb > 2) {
        solver->set_verbosity(conf.verb-2);
    }

    //parsing the input
    if (vm.count("input") == 0) {
        cout << "ERROR: you must pass a file" << endl;
    }
    const string inp = vm["input"].as<string>();

    //Read in file and set sampling_set in case we are starting with empty
    readInAFile(inp.c_str(), 0, sampling_set->empty());
    init_samping_set(conf.recompute_sampling_set);
    orig_num_vars = solver->nVars();
    
    add_fixed_clauses();
   
    cout << "[mis] Setup time: " << (cpuTime()-myTime) << endl;
    
    find_unate();
    
    cout << "[mis] Total time: " << (cpuTime()-myTime) << endl;

    
    
}




int main(int argc, char** argv)
{
    #if defined(__GNUC__) && defined(__linux__)
    feenableexcept(FE_INVALID   |
                   FE_DIVBYZERO |
                   FE_OVERFLOW
                  );
    #endif
    sampling_set = &sampling_set_tmp1;
    //Reconstruct the command line so we can emit it later if needed
    for(int i = 0; i < argc; i++) {
        command_line += string(argv[i]);
        if (i+1 < argc) {
            command_line += " ";
        }
    }

    add_supported_options(argc, argv);
    //cout << "[mis] Arjun Version: " << get_version_sha1() << endl;
    cout
    << "c executed with command line: "
    << command_line
    << endl;
    cout << "[mis] using seed: " << conf.seed << endl;
    mtrand.seed(conf.seed);
    init_solver_setup();


    std::string s = vm["input"].as<string>();
    std::string delimiter = ".qdimacs";
    std ::string filename = s.substr(0, s.find(delimiter));
    string tempfile=filename+"_vardetails";
    fstream new_file; 
    new_file.open(tempfile,ios::out);  
    if(!new_file) 
      {
        cout<<"File creation failed";
      }
    else
      {
        
        new_file<<"Posunate :"<<posunate_details<<endl;
        new_file<<"Negunate :"<<negunate_details<<endl;
        new_file.close();
      }
    //cout << solver->get_text_version_info();
    
    
    delete solver;
    return 0;
}
