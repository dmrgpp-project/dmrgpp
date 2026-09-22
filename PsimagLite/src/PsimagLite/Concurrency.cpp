#include "Concurrency.h"

namespace PsimagLite {

SizeType          Concurrency::mode = 0;
LabelDisabled     Concurrency::mpi_disabled;
CodeSectionParams Concurrency::code_section_params(1);
} // namespace PsimagLite
