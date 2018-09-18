#include "../daropstests_aux.h"

// Hack for testing private members
#define protected public
#define private public
#include "../../KMedoidsReqDistributor.cc"

#undef protected
#undef private

#include <time.h>

using namespace std;


void fillClusterIdsOnlyFirstCluster(int clusterids[], int size, int nodes) {
  for (int i=0;i<size; i++) clusterids[i] = nodes-1;
}

void fillClusterIdsOnlyLastCluster(int clusterids[], int size, int nodes) {
  for (int i=0;i<size; i++) clusterids[i] = nodes-1;
}

void fillClusterIdsUniformly(int clusterids[], int size, int nodes) {
  for (int i=0;i<size; i++) clusterids[i] = i % nodes;
}

void fillClusterIdsRandomly(int clusterids[], int size, int nodes) {
  for (int i=0;i<size; i++) clusterids[i] = GARandomInt(0,nodes-1);
}


void checkCreateDistReqs(void (*fillIds)(int [], int , int ) ) {
  int id,pos,load;
  bool crit;
  long fb,fe;
  Vertex::VertexType t;

  id = 1; pos = 1; crit = true; load = 1; fb = 5; fe = 10; t = Vertex::PICKUP;
  Vertex v1(id,pos,crit,load,fb,fe,t);

  id = 2; pos = 2; crit = true; load = 1; fb = 5; fe = 10; t = Vertex::PICKUP;
  Vertex v2(id,pos,crit,load,fb,fe,t);

  id = -2; pos = 3; crit = true; load = 1; fb = 5; fe = 10; t = Vertex::DELIVERY;
  Vertex v3(id,pos,crit,load,fb,fe,t);

  id = -1; pos = 4; crit = true; load = 1; fb = 5; fe = 10; t = Vertex::DELIVERY;
  Vertex v4(id,pos,crit,load,fb,fe,t);

  id = 3; pos = 5; crit = true; load = 1; fb = 5; fe = 10; t = Vertex::PICKUP;
  Vertex v5(id,pos,crit,load,fb,fe,t);


  Request r1(&v1,&v2);  // Dummy requests, only the delivery vertex is being used for the computation of the matrix
  Request r2(&v2,&v1);
  Request r3(&v3,&v3);
  Request r4(&v4,&v3);

  vector<Request> reqs; reqs.push_back(r1); reqs.push_back(r2); reqs.push_back(r3); reqs.push_back(r4);

  int nodes = 3;

  int clusterids[reqs.size()];
  fillIds(clusterids,reqs.size(),nodes);

  vector< vector<Request> > dist_reqs = createDistReqs(reqs,nodes,clusterids);
  ASSERT_EQ(nodes,dist_reqs.size());

  int sum_reqs = 0;
  for (int nodepos=0; nodepos<nodes; nodepos++) sum_reqs += dist_reqs[nodepos].size();
  ASSERT_EQ(reqs.size(),sum_reqs);

  for (int i=0; i<reqs.size(); i++) {
    Request& req           = reqs[i];
    int      assigned_node = clusterids[i];
    bool     found         = find( dist_reqs[assigned_node].begin(), dist_reqs[assigned_node].end(), req ) != dist_reqs[assigned_node].end();
    ASSERT_EQ( true, found );
  }

}

TEST(KMedoidsReqDistributorTest,checkcreateDistReqs) {
  checkCreateDistReqs(fillClusterIdsOnlyFirstCluster);
  checkCreateDistReqs(fillClusterIdsOnlyLastCluster);
  checkCreateDistReqs(fillClusterIdsUniformly);
  checkCreateDistReqs(fillClusterIdsUniformly);
}


int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);

  GAGenome::optCriterion(GAGenome::MINIMIZATION);

  setUpTestingEnvironment();

  return RUN_ALL_TESTS();
}
