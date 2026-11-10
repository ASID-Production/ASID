#include "gtest/gtest.h"
#include "../../src/Classes/BaderOperator.h"
using namespace cpplib;
TEST(BaderRegressionTest, NewtonIsScaleInvariantAndHandlesSingularity) {
    for(double scale:{1e-100,1.0,1e100}) {
        CriticalPoint cp(CriticalPoint::TYPE::B,{0.01,0.02,0.03});
        TripleDouble value;
        value.hess.El(0,0)=2*scale; value.hess.El(1,1)=-3*scale; value.hess.El(2,2)=-4*scale;
        value.grad={0.02*scale,-0.06*scale,-0.12*scale}; cp.UpdateVal(value);
        EXPECT_NEAR(cp.NewtonRaphsonPredict().r(),0,1e-15);
        value.hess.El(2,2)=0; cp.UpdateVal(value);
        EXPECT_EQ(cp.NewtonRaphsonPredict().a,cp.pos_.a);
    }
}
TEST(BaderRegressionTest, RfoFollowsCriticalPointSignature) {
    for(auto type:{CriticalPoint::TYPE::B,CriticalPoint::TYPE::R,CriticalPoint::TYPE::C}) {
        CriticalPoint cp(type,{0,0,0}); TripleDouble value;
        const int positives=static_cast<int>(type);
        for(int i=0;i<3;++i) {value.grad[i]=0.01; value.hess.El(i,i)=i<positives?2:-2;}
        cp.UpdateVal(value); const auto step=cp.FindNextPosition();
        for(int i=0;i<3;++i) EXPECT_LT(step[i]*value.grad[i]*value.hess.El(i,i),0);
        EXPECT_LE(step.r(),CriticalPoint::MAX_STEP);
    }
    CriticalPoint nuclear(CriticalPoint::TYPE::N,{1,2,3});
    EXPECT_EQ(nuclear.FindNextPosition().a,nuclear.pos_.a);
}
TEST(BaderRegressionTest, RefinementConvergesWithIndependentNewtonBudget) {
    std::vector<RadialSpline> splines; splines.emplace_back(1,2,1,std::vector<double>{2,0});
    PointsSoA<double> psoa;
    psoa.cartesians={-3,-3,-3,3,3,3};
    psoa.addPoint({-1.5,0,0},0,0); psoa.addPoint({1.5,0,0},1,0);
    psoa.calculateSpatialIndexes(1); psoa.computeRuntimeFlatShifts();
    // A large RFO stopping threshold forces immediate exit; Newton must still converge.
    for(double stop:{1.0,1e-20}) {
        CriticalPoint cp(CriticalPoint::TYPE::B,{0.02,0.01,0.01});
        cp.UpdateVal(CriticalPoint::CalculateEDinPoint(cp.pos_,psoa,splines));
        BaderOperator::OptimizePoint(cp,psoa,splines,stop);
        EXPECT_LT(cp.value_.getGradSq(),1e-12);
        EXPECT_LT(cp.pos_.r(),1e-6);
        EXPECT_EQ(cp.type_,CriticalPoint::TYPE::B);
    }
}
