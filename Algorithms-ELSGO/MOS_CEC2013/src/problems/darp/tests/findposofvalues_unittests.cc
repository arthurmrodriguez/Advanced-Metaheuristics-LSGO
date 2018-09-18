#include "gtest/gtest.h"
#include "../aux.h"
#include <vector>

using namespace std;

TEST(findPosOfValue,Empty) {
  vector<int> values;
  EXPECT_EQ(-1,findPosOfValue(0,values));
}

TEST(findPosOfValue,OnlyOne) {
  vector<int> values;
  values.push_back(1);
  EXPECT_EQ(0,findPosOfValue(1,values));
  EXPECT_EQ(-1,findPosOfValue(0,values));
}

TEST(findPosOfValue,Two) {
  vector<int> values;
  values.push_back(1);
  values.push_back(2);
  EXPECT_EQ(0,findPosOfValue(1,values));
  EXPECT_EQ(1,findPosOfValue(2,values));
  EXPECT_EQ(-1,findPosOfValue(0,values));
  EXPECT_EQ(-1,findPosOfValue(-1,values));
  EXPECT_EQ(-1,findPosOfValue(-2,values));
}

TEST(findPosOfValue,TwoInverted) {
  vector<int> values;
  values.push_back(2);
  values.push_back(1);
  EXPECT_EQ(1,findPosOfValue(1,values));
  EXPECT_EQ(0,findPosOfValue(2,values));
  EXPECT_EQ(-1,findPosOfValue(0,values));
  EXPECT_EQ(-1,findPosOfValue(-1,values));
  EXPECT_EQ(-1,findPosOfValue(-2,values));
}

TEST(findPosOfValue,TenSorted) {
  vector<int> values;
  for (int i=0; i<=10; i++) values.push_back(i);

  EXPECT_EQ(0,findPosOfValue(0,values));
  EXPECT_EQ(1,findPosOfValue(1,values));
  EXPECT_EQ(0,findPosOfValue(0,values));

  EXPECT_EQ(10,findPosOfValue(10,values));

  EXPECT_EQ(-1,findPosOfValue(11,values));
  EXPECT_EQ(-1,findPosOfValue(-1,values));
  EXPECT_EQ(-1,findPosOfValue(-2,values));
}

TEST(findPosOfValue,TenInverted) {
  vector<int> values;
  int max_value = 10;
  for (int i=0; i<=max_value; i++) values.push_back(max_value-i);

  EXPECT_EQ(max_value,findPosOfValue(0,values));
  EXPECT_EQ(max_value-1,findPosOfValue(1,values));
  EXPECT_EQ(max_value,findPosOfValue(0,values));

  EXPECT_EQ(0,findPosOfValue(max_value,values));

  EXPECT_EQ(-1,findPosOfValue(max_value+1,values));
  EXPECT_EQ(-1,findPosOfValue(-1,values));
  EXPECT_EQ(-1,findPosOfValue(-2,values));
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
