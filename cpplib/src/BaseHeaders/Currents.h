#pragma once

#include "BaseTypes.h"
#include "Concepts.h"

namespace cpplib {

    // Forward declarations
    class SimpleAtom;
    class CompositeAtom;

    namespace currents {

        using AtomTypeRequest = CompositeAtom;
        using AtomTypeData = SimpleAtom;

    }
}
