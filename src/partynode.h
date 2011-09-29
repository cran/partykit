
#define ID_NODE                 0               /* id_node */
#define SPLIT_NODE              1               /* split_node */
#define KIDS_NODE               2               /* kids_node */
#define SURROGATES_NODE         3               /* surrogates_node */
#define INFO_NODE               4               /* info_node */
#define LENGTH_NODE             INFO_NODE + 1   /* length node object */

int fitted_node(SEXP node, SEXP data, SEXP vmatch, SEXP perm, SEXP perms, int obs);
