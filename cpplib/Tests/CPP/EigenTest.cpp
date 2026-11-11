#include "gtest/gtest.h"
#include "../../src/Classes/Geometry.h"
#include <random>
using cpplib::geometry::Matrix;
namespace {
void checkEigen(const Matrix<double>& a,double scale) {
    const auto e=a.EigenvaluesAndVectors();
    EXPECT_GE(e.values[0],e.values[1]); EXPECT_GE(e.values[1],e.values[2]);
    for(int i=0;i<3;++i) {
        for(int k=0;k<3;++k) {
            double av=0; for(int j=0;j<3;++j) av+=(a.El(k,j)/scale)*e.vectors[i][j];
            EXPECT_NEAR(av,(e.values[i]/scale)*e.vectors[i][k],1e-12);
        }
        for(int j=0;j<3;++j) {
            double dot=0; for(int k=0;k<3;++k) dot+=e.vectors[i][k]*e.vectors[j][k];
            EXPECT_NEAR(dot,i==j?1:0,1e-12);
        }
    }
}
}
TEST(EigenTest, DegeneracyAndScale) {
    for(const auto& values:std::vector<std::array<double,9>>{
        {3,0,0,0,2,0,0,0,1},{2,0,0,0,2,0,0,0,-1},
        {1,0,0,0,1,0,0,0,1},{-1,0,0,0,0,0,0,0,1},
        {2,1,1,1,2,1,1,1,2},{1,1e-12,0,1e-12,1,0,0,0,1+1e-12}}) {
        for(double scale:{1e-200,1e-7,1.0,1e200}) {
            auto a=values; for(auto& x:a) x*=scale;
            checkEigen(Matrix<double>(a),scale);
        }
    }
    checkEigen(Matrix<double>{},1);
}
TEST(EigenTest, RandomSymmetricResiduals) {
    std::mt19937 gen(138); std::uniform_real_distribution<double> dist(-1,1);
    for(int n=0;n<1000;++n) {
        Matrix<double> a;
        for(int i=0;i<3;++i) for(int j=i;j<3;++j) a.El(i,j)=a.El(j,i)=dist(gen);
        checkEigen(a,1);
    }
}
TEST(EigenTest, NonfiniteIsExplicit) {
    for(double invalid:{std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity()}) {
        Matrix<double> a(1.0); a.El(0,0)=invalid;
        for(double value:a.EigenvaluesAndVectors().values) EXPECT_TRUE(std::isnan(value));
    }
}
