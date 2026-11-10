#include "gtest/gtest.h"
#include "../../src/Classes/Splines.h"
#include "../../src/Classes/Reader.h"
#include <limits>

using cpplib::CubicSpline;
namespace {
double cubic(double r) { const auto t=r-1; return t*t*t-3*t+2; }
void checkCubic(const CubicSpline<>& spline, double r) {
    double f,g,h; spline.eval(r,f,g,h);
    EXPECT_NEAR(f,cubic(r),1e-11);
    EXPECT_NEAR(g,3*(r-1)*(r-1)-3,1e-10);
    EXPECT_NEAR(h,6*(r-1),1e-10);
}
}
TEST(SplinesTest, RuntimeValidation) {
    EXPECT_THROW((CubicSpline<>(1,2,0,{})),std::invalid_argument);
    EXPECT_THROW((CubicSpline<>(1,2,577,{})),std::invalid_argument);
    EXPECT_THROW((CubicSpline<>(1,2,1,{1})),std::invalid_argument);
    for (double r : {0.0,-1.0,std::numeric_limits<double>::quiet_NaN()})
        EXPECT_THROW((CubicSpline<>(r,2,1,{1,0})),std::invalid_argument);
    EXPECT_THROW((CubicSpline<>(2,1,1,{1,0})),std::invalid_argument);
    EXPECT_THROW((CubicSpline<>({0,0},{1,0})),std::invalid_argument);
    EXPECT_THROW((CubicSpline<>({0,2,1},{2,1,0})),std::invalid_argument);
    EXPECT_THROW((CubicSpline<>({-1,2},{1,0})),std::invalid_argument);
    EXPECT_THROW((CubicSpline<>({0,1},{1,std::numeric_limits<double>::infinity()})),std::invalid_argument);
}
TEST(SplinesTest, RejectsUnrepresentableCoefficients) {
    EXPECT_THROW((CubicSpline<>({0,1e-310},{1,0})),std::invalid_argument);
}
TEST(SplinesTest, ArbitraryGridsReproduceCubicAndDerivatives) {
    for (const auto& knots : std::vector<std::vector<double>>{{1,2},{1,1.1,1.7,2},{1,1.25,1.5,1.75,2}}) {
        std::vector<double> values; for (double r:knots) values.push_back(cubic(r));
        CubicSpline<> spline(knots,values);
        for (double r:knots) checkCubic(spline,r);
        for (double r:{0.0,0.5,1.01,1.3,1.9,1.999999}) checkCubic(spline,r);
    }
}
TEST(SplinesTest, LogarithmicGridAndFinalIntervalSearch) {
    std::vector<double> values; for (int i=0;i<=4;++i) values.push_back(cubic(std::pow(2.0,i/4.0)));
    CubicSpline<> spline(1,2,4,values);
    checkCubic(spline,1.8);
    const double target=cubic(1.99);
    EXPECT_NEAR(spline.find_value(1,target),1.99,1e-12);
    double f,g,h; spline.eval(3,f,g,h); EXPECT_EQ(f,0); EXPECT_EQ(g,0); EXPECT_EQ(h,0);
    EXPECT_THROW(spline.eval(-1,f,g,h),std::invalid_argument);
    EXPECT_THROW(spline.eval(std::numeric_limits<double>::quiet_NaN(),f,g,h),std::invalid_argument);
}
TEST(SplinesTest, RadialHessianAndCenterLimit) {
    cpplib::RadialSpline spline(1,2,1,{2,0});
    auto center=spline.evaluate({0,0,0});
    EXPECT_DOUBLE_EQ(center.val,4);
    for (int i=0;i<3;++i) {
        EXPECT_EQ(center.grad[i],0);
        for (int j=0;j<3;++j) EXPECT_DOUBLE_EQ(center.hess.El(i,j),i==j ? -6 : 0);
    }
    const cpplib::RadialSpline::PointType p(0.3,0.4,0.2);
    const auto value=spline.evaluate(p);
    const double step=1e-5;
    for (int i=0;i<3;++i) {
        auto plus=p,minus=p; plus[i]+=step; minus[i]-=step;
        const auto vp=spline.evaluate(plus), vm=spline.evaluate(minus);
        EXPECT_NEAR(value.grad[i],(vp.val-vm.val)/(2*step),1e-8);
        for (int j=0;j<3;++j) {
            EXPECT_NEAR(value.hess.El(j,i),(vp.grad[j]-vm.grad[j])/(2*step),1e-8);
            EXPECT_NEAR(value.hess.El(i,j),value.hess.El(j,i),1e-12);
        }
    }
}
TEST(SplinesTest, ExplicitZeroNodeHasSmoothRadialCenter) {
    cpplib::RadialSpline radial;
    radial.spline.create({0,1,2},{2,1,0});
    const auto center=radial.evaluate({0,0,0});
    EXPECT_EQ(center.val,2);
    EXPECT_LT(center.hess.El(0,0),0);
    EXPECT_EQ(center.hess.El(0,0),center.hess.El(1,1));
    EXPECT_EQ(center.hess.El(1,1),center.hess.El(2,2));
}

TEST(SplinesTest, TinyRadiusDoesNotUnderflowToCenter) {
    cpplib::RadialSpline radial;
    radial.spline.create({0,1,2},{2,1,0});
    const auto value=radial.evaluate({1e-200,0,0});
    EXPECT_TRUE(std::isfinite(value.val));
    EXPECT_TRUE(std::isfinite(value.hess.El(0,0)));
    const auto center=radial.evaluate({0,0,0});
    for(int i=0;i<3;++i) EXPECT_NEAR(value.hess.El(i,i),center.hess.El(i,i),1e-12);
}
