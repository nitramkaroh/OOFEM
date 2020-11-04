#pragma once

namespace oofem
{
  // Helper function to so that the generic permute() only fails if
  // actually instantiated
  template <typename T> constexpr bool error_when_instantiated()
  {
    return false;
  }
}
