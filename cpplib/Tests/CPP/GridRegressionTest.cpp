#include "gtest/gtest.h"
#include "../../src/Classes/Geometry.h"
using namespace cpplib::geometry;
TEST(GridRegressionTest, NonperiodicVisitsAdjacentBucketsWithoutBoundaryContacts) {
    for(double cutoff:{1.0,2.0,6.0,12.0}) {
        Cell<double> cell(10,10,10,90,90,90);
        std::vector<Point<double>> points{{0.01,0.01,0.01},{0.99,0.99,0.99},{0.49,0.5,0.5},{0.51,0.5,0.5}};
        SG<double> grid(cell,false,cutoff); grid.updateOnlyPoints(points);
        const auto contacts=grid.findAllContacts();
        size_t expected=0;
        for(size_t i=0;i<points.size();++i) for(size_t j=i+1;j<points.size();++j)
            if((cell.fracToCart()*(points[i]-points[j])).rSq()<=cutoff*cutoff) ++expected;
        EXPECT_EQ(contacts.size(),expected);
        for(const auto& bond:contacts) EXPECT_EQ(bond.shiftcode.get_code(),13);
    }
}
TEST(GridRegressionTest, PbcBoundaryContactIsOptional) {
    Cell<double> cell(10,10,10,90,90,90);
    for(bool periodic:{false,true}) {
        SG<double> grid(cell,periodic,1);
        grid.updateOnlyPoints({{0.01,0.5,0.5},{0.99,0.5,0.5}});
        EXPECT_EQ(grid.findAllContacts().size(),periodic?1:0);
    }
}
