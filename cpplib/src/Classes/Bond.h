#include <string>
#include "../BaseHeaders/BaseTypes.h"
#include "../BaseHeaders/Concepts.h"
namespace cpplib {
	struct Bond {
	public:
		// Declarations
		using AtomIndex = basic_types::AtomIndex;

		// Data
		AtomIndex first = 0;
		AtomIndex second = 0;

		// Constructors
		constexpr Bond() noexcept = default;
		constexpr Bond(const AtomIndex a1, const AtomIndex a2) noexcept : first(a1), second(a2) {
		};

		// Operators
		constexpr auto operator<=>(const Bond& other) const noexcept = default;

		// Functions
		constexpr void validate() noexcept {
			if (first > second) ::std::swap(first, second);
		}
		::std::string ToStr() const {
			::std::string res("(");
			res += ::std::to_string(this->first);
			res += ", ";
			res += ::std::to_string(this->second);
			res += ")";
			// (1, 2)
			return res;
		}
	};
	struct BondEx : public Bond {
	public:
		// Declarations
		using LengthType = basic_types::FloatingPointType; // Float or Double
		using base = Bond; // Bond

		// Data
		LengthType length = 0.0; // Bond Length 

		// Constructors
		constexpr BondEx() = default;

		// Always uses Bond::validate function
		constexpr BondEx(const base& bond, float len) noexcept : base(bond), length(len) {
			base::validate();
		}
		constexpr BondEx(AtomIndex a1, AtomIndex a2, LengthType l) noexcept
			: length(l) {
			if (a1 < a2) {
				base::first = a1;
				base::second = a2;
			} else {
				base::first = a2;
				base::second = a1;
			}
		}

		constexpr BondEx(AtomIndex a1, AtomIndex a2) noexcept : Bond(a1, a2) {
		}

		// Compares only "base". Ignores length.
		constexpr bool operator==(const BondEx& other) const noexcept {
			return base::operator==(other);
		}

		// Compares in next order:
		// 1. "base" 
		// 2. length
		constexpr auto operator<=>(const BondEx& other) const noexcept = default;

		// to_string constant function
		::std::string ToStr() const {
			::std::string res("(");
			res += ::std::to_string(first);
			res += ", ";
			res += ::std::to_string(second);
			res += ", {\"distance\": ";
			res += ::std::to_string(length);
			res += "})";
			return res;
			// (1, 2, {"distance": 1.0})
		}
	};

	static_assert(BondConcept<Bond>, "Bond must satisfy BondConcept");
	static_assert(BondConcept<BondEx>, "BondEx must satisfy BondConcept");
}