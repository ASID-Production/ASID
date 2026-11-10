#include "gtest/gtest.h"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <string>

#include "../src/Classes/Distances.h"

namespace {

class TemporaryData {
public:
	TemporaryData() : directory_(std::filesystem::temp_directory_path() /
		("cpplib-singleton-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()))) {
		std::filesystem::create_directory(directory_);
	}
	~TemporaryData() {
		std::error_code ignored;
		std::filesystem::remove_all(directory_, ignored);
	}
	std::filesystem::path path(const char* name) const { return directory_ / name; }
private:
	std::filesystem::path directory_;
};

class TestData : public cpplib::SingletonBase<TestData> {
	friend class cpplib::SingletonBase<TestData>;
	TestData(const std::filesystem::path& file, cpplib::STATE& state) : SingletonBase(file) {
		std::ifstream input(file);
		state = input.is_open() ? cpplib::STATE::OK : cpplib::STATE::FILE_NOT_FOUND;
	}
public:
	TestData() : SingletonBase(std::filesystem::path{}) {}
};

TEST(SingletonLoadingTest, RetriesFailureAndRejectsDifferentPaths) {
	using cpplib::STATE;
	TemporaryData data;
	const auto first = data.path("first.ini");
	const auto second = data.path("second.ini");

	EXPECT_FALSE(TestData::exists());
	EXPECT_EQ(TestData::try_call({}).second, STATE::INVALID_ARGUMENT);
	EXPECT_EQ(TestData::try_call(first).second, STATE::FILE_NOT_FOUND);
	EXPECT_THROW(TestData::instance(first), std::runtime_error);
	EXPECT_FALSE(TestData::exists());

	std::ofstream(first) << "data";
	const TestData* loaded = nullptr;
	ASSERT_NO_THROW(loaded = &TestData::instance(first));
	EXPECT_TRUE(TestData::exists());
	EXPECT_EQ(TestData::state(), STATE::OK);
	EXPECT_EQ(&TestData::instance(first.parent_path() / "." / first.filename()), loaded);
	EXPECT_EQ(TestData::try_call(second).second, STATE::INVALID_ARGUMENT);
	EXPECT_THROW(TestData::instance(second), std::runtime_error);
	EXPECT_EQ(&TestData::instance(first), loaded);
	EXPECT_EQ(TestData::reload(second), STATE::FILE_NOT_FOUND);
	EXPECT_EQ(&TestData::instance(first), loaded);
}

TEST(DistancesLoadingTest, DefaultPathUsesLoadedModuleDirectory) {
	EXPECT_EQ(cpplib::default_relpath<cpplib::Distances>, "module/BondLength.ini");
	EXPECT_TRUE(cpplib::get_module_directory().is_absolute());
	EXPECT_EQ(cpplib::Distances::instance().filename(),
		std::filesystem::weakly_canonical(cpplib::get_module_directory() / "BondLength.ini"));
}

TEST(DistancesLoadingTest, RejectsTruncatedReloadAndKeepsLoadedData) {
	TemporaryData data;
	const auto invalid = data.path("invalid.ini");
	const auto& loaded = cpplib::Distances::instance();
	const auto limit = loaded.maxDistance(1, 1);
	std::ofstream(invalid) << "1\n1 1 0.2";

	EXPECT_EQ(cpplib::Distances::reload(invalid), cpplib::STATE::BAD_SIZE);
	EXPECT_EQ(&cpplib::Distances::instance(), &loaded);
	EXPECT_EQ(loaded.maxDistance(1, 1), limit);
	EXPECT_TRUE(cpplib::Distances::exists());
}

} // namespace
