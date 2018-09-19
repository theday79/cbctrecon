// For testing CbctRecon
#if USE_TINYREFL
#include <tinyrefl/api.hpp>
#include "cbctrecon.h"
#include "cbctrecon.h.tinyrefl"
#else
#include "cbctrecon.h"
#endif

int main()
{
  std::cout << "Running cbctrecon_test!" << std::endl;
  return 0;
}
