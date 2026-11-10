#include "gtest/gtest.h"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <string>
#include <thread>
#include <atomic>
#include <vector>

#include "../../src/Classes/Distances.h"

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
	explicit TestData(const std::filesystem::path& file) : SingletonBase(file) {
		std::ifstream input(file);
        std::getline(input, value);
		if (!input.is_open()) throw cpplib::DataLoadError(cpplib::STATE::FILE_NOT_FOUND, file);
	}
public:
    std::string value;
	TestData() : SingletonBase(std::filesystem::path{}) {}
};

void CheckSingletonLoading() {
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
    const auto original_cwd = std::filesystem::current_path();
    std::filesystem::current_path(first.parent_path());
    EXPECT_EQ(&TestData::instance(), loaded);
    std::filesystem::current_path(original_cwd);
	EXPECT_EQ(TestData::state(), STATE::OK);
	EXPECT_EQ(&TestData::instance(first.parent_path() / "." / first.filename()), loaded);
	EXPECT_EQ(TestData::try_call(second).second, STATE::INVALID_ARGUMENT);
	EXPECT_THROW(TestData::instance(second), std::runtime_error);
	EXPECT_EQ(&TestData::instance(first), loaded);
    std::ofstream(second) << "different data";
    EXPECT_THROW(TestData::instance(second), cpplib::DataLoadError);
    EXPECT_EQ(&TestData::instance(), loaded);
    EXPECT_EQ(loaded->value, "data");
    std::atomic<bool> bad_read{false};
    std::vector<std::thread> readers;
    for (int i = 0; i < 4; ++i) readers.emplace_back([&] {
        for (int j = 0; j < 1000; ++j) {
            const auto& object = TestData::instance();
            if (&object != loaded || object.value != "data") bad_read.store(true);
        }
    });
    for (auto& reader : readers) reader.join();
    EXPECT_FALSE(bad_read.load());




}

TEST(SingletonLoadingDeathTest, RetriesFailureAndRejectsDifferentPaths) {
    ASSERT_EXIT({ CheckSingletonLoading(); std::exit(::testing::Test::HasFailure() ? 1 : 0); },
                ::testing::ExitedWithCode(0), "");
}

TEST(DistancesLoadingTest, DefaultPathUsesLoadedModuleDirectory) {
	EXPECT_EQ(cpplib::default_relpath<cpplib::Distances>, "module/BondLength.ini");
	EXPECT_TRUE(cpplib::get_module_directory().is_absolute());
	EXPECT_EQ(cpplib::Distances::instance().filename(),
		std::filesystem::weakly_canonical(cpplib::get_module_directory() / "BondLength.ini"));
}


} // namespace
