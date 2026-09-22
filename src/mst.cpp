#include <cstdlib>
#include <R.h>
#include "constants.h"
#include "linkage_group_DH.h"
#include "genetic_map_DH.h"
#include "genetic_map_RIL.h"
#include <Rdefines.h>


SEXP elem(SEXP list, const char *str);

int trace;
extern "C" SEXP mst(SEXP Plist, SEXP data)
{
  SEXP map = R_NilValue;
  string pop_type;
  genetic_map* barley;

  trace = INTEGER(elem(Plist, "trace"))[0];

  pop_type = CHAR(STRING_ELT(elem(Plist, "population_type"),0));
  if (pop_type == "DH")
    barley = new genetic_map_DH();
  else
    barley = new genetic_map_RIL();

  barley->read_raw_mapping_data(Plist, data);
  // Reporting strategy
  // 1. generate_map
  //    --> alloc return list 'map' (number-of-connected-components)
  //    --> for each lg allocate a list of length 2 ('dist','imputed_values')
  //    --> catch current_linkage_group raw_data in 'imputed_values' component
  //        from dump().
  // 2. write_output
  //    --> catch distances in 'dist' for each lg
  //    --> name the distances
  //    --> set dimnames(imputed_values)
  barley->generate_map(map);
  barley->write_output(map);

  // Release the working set (raw genotype data and the n x n distance matrix),
  // which was previously leaked on every call. Deleting a C++ object performs
  // no R allocation, so this cannot trigger a garbage collection of 'map'.
  // Note that an Rf_error raised above longjmps past this point; making the
  // error paths leak free requires R_UnwindProtect and is left for later.
  delete barley;

  return(map);
}

void print_vector(vector<int> tmp)
{
  for (unsigned int ii = 0 ; ii < tmp.size(); ii++)
    {
      Rprintf("%d,", ii);
    }
  Rprintf("\n");
}
