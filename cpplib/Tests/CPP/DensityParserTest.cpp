#include "gtest/gtest.h"
#include "../../src/Classes/Reader.h"
#include <chrono>
TEST(DensityParserTest, RejectsMalformedFilesInRelease) {
    const auto path=std::filesystem::temp_directory_path()/
        ("cpplib-density-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count())+".txt");
    struct Cleanup { std::filesystem::path file; ~Cleanup() {std::error_code ec; std::filesystem::remove(file,ec);} } cleanup{path};
    for(const char* input:{"", "-1\n", "119\n", "1\n1 .1 1 0\n", "1\n1 .1 1 577\n",
        "1\n1 0 1 1\n2\n", "1\n1 .1 1 1\nnan\n", "1\n1 .1 1 2\n2\n",
        "1\n1 .1 1 1\n2junk\n", "1\n1 .1 1 1\n2\n1 .1 1 1\n2\n"}) {
        std::ofstream(path)<<input;
        EXPECT_THROW(cpplib::DensityParser::parse_file(path),std::exception)<<input;
    }
    std::ofstream(path)<<"1\n1 .1 1 1\n2\n";
    const auto parsed=cpplib::DensityParser::parse_file(path);
    ASSERT_EQ(parsed.size(),2);
    EXPECT_DOUBLE_EQ(parsed[1].spline.radius(),1);
}
