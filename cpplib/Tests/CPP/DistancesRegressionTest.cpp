#include "gtest/gtest.h"
#include "../../src/Classes/Distances.h"
using cpplib::Distances;
TEST(DistancesRegressionTest, FiltersOverlapsAndNonbonds) {
    const auto& d=Distances::instance();
    const auto lo=d.minDistance(1,1), hi=d.maxDistance(1,1);
    ASSERT_GT(lo,0); ASSERT_GT(hi,lo);
    using Grid=cpplib::geometry::SG<Distances::FloatingPointType>;
    std::vector<Distances::AtomTypeBase> types{1,1};
    for(auto length:{lo/2,(lo+hi)/2,hi*2}) {
        const bool bonded=length>=lo && length<=hi;
        const bool contact=length<=hi;
        EXPECT_EQ(static_cast<signed char>(d.isBond(1,1,length)),length<lo?-1:bonded?1:0);
        std::vector<Distances::PointType> points{{0,0,0},{length,0,0}};
        std::vector<Grid::BondWithShift> bonds{{0,1,cpplib::geometry::ShiftCode(13),length*length}};
        EXPECT_EQ(d.filter_bond_list(bonds,types,points).size(),bonded?1:0);
        EXPECT_EQ(d.filter_bond_list(bonds,types,points,Distances::BondFilter::Contacts).size(),contact?1:0);
        using OldGrid=cpplib::geometry::SpatialGrid<Distances::FloatingPointType>;
        std::vector<OldGrid::BondWithShift> old{{0,1,cpplib::geometry::ShiftCode(13)}};
        EXPECT_EQ(d.filter_bond_list(old,types,points,Distances::PointType::distance).size(),bonded?1:0);
        EXPECT_EQ(d.filter_bond_list(old,types,points,Distances::PointType::distance,
                                    Distances::BondFilter::Contacts).size(),contact?1:0);
    }
}

TEST(DistancesRegressionTest, GeometricNeighborsRetainOverlaps) {
    const auto& d=Distances::instance();
    const auto length=d.minDistance(1,1)/2;
    cpplib::geometry::Cell<double> cell(10,10,10,90,90,90);
    cpplib::geometry::SG<double> grid(cell,false,2);
    grid.updateOnlyPoints({{.5,.5,.5},{.5+length/10,.5,.5}});
    EXPECT_EQ(grid.findAllContacts().size(),1);
}
