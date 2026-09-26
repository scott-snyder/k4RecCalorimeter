/**
 * @file RecCalorimeter/src/components/PairCaloClustersPi0.cpp
 * @author Zhibo Wu, scott snyder <snyder@bnl.gov>
 * @date Rewritten Sep, 2026
 * @brief Make pi0 candidates from cluster pairs.
 */


#include "PairCaloClustersPi0.h"
// Key4HEP
#include "k4FWCore/MetadataUtils.h"
#include "k4FWCore/GaudiChecks.h"

#include "TLorentzVector.h"
#include "TVector3.h"

#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/maximum_weighted_matching.hpp"

// Include the <cmath> header for sqrt, pow
#include <cmath>
#include <fstream>
#include <ranges>


DECLARE_COMPONENT(PairCaloClustersPi0)


namespace {


/// Graph data structures.  This is the recommended structure for a
/// non-sparse graph not being modified dynamically.
/// The maximum_weighted_matching algorithm requires the edge weight
/// as an internal property.
using EdgeProps = boost::property<boost::edge_weight_t, double>;
using Graph = boost::adjacency_list<boost::vecS,
                                    boost::vecS,
                                    boost::undirectedS,
                                    boost::no_property,
                                    EdgeProps>;
using Vertex = boost::graph_traits<Graph>::vertex_descriptor;
using Edge = boost::graph_traits<Graph>::edge_descriptor;


/// Helper: Turn a pair of iterators into a range.
template <class IT>
auto make_range (const std::pair<IT, IT>& p)
{
  return std::ranges::subrange (p.first, p.second);
}


/// Helper: Extract a 4-vector from a cluster (assuming it comes
/// from the origin).
TLorentzVector getTLV (const edm4hep::Cluster& cl)
{
  double e = cl.getEnergy();
  TVector3 disp(cl.getPosition().x, cl.getPosition().y, cl.getPosition().z);
  return TLorentzVector(disp * (e / disp.Mag()), e);
}


/**
 * @brief Create graph corresponding to a set of clusters.
 * @param inClusters 
 */
Graph makeGraph(const edm4hep::ClusterCollection& inClusters,
                std::vector<TLorentzVector>& pairs,
                double minClusterEnergy,
                double maxDR,
                double masspeak,
                double masslow,
                double masshigh)
{
  Graph g(inClusters.size());

  double wsum = 0;

  for (size_t i = 0; i < inClusters.size(); ++i) {
    const auto& cl_i = inClusters.at(i);
    if (cl_i.getEnergy() < minClusterEnergy)
      continue;
    TLorentzVector tlv_i = getTLV (cl_i);

    for (size_t j = i + 1; j < inClusters.size(); j++) {
      const auto& cl_j = inClusters.at(j);
      if (cl_j.getEnergy() < minClusterEnergy)
        continue;
      TLorentzVector tlv_j = getTLV (cl_j);

      TLorentzVector vpair = tlv_i + tlv_j;
      double invM = vpair.M();
      if (invM > masslow && invM < masshigh && tlv_i.DeltaR(tlv_j) < maxDR) {
        double w = std::pow (invM - masspeak, 2);
        wsum += w;
        boost::add_edge (i, j, w, g);
        pairs.push_back (vpair);
      }      
    }
  }

  for (auto e : make_range (boost::edges (g))) {
    double w = boost::get (boost::edge_weight, g, e);
    boost::put (boost::edge_weight, g, e, 2*wsum - w);
  }

  return g;
}


std::vector<Edge> findEdges (const Graph& g)
{
  size_t nv = boost::num_vertices(g);
  std::vector<Vertex> mate (nv);
  boost::maximum_weighted_matching (g, mate.data());

  std::vector<Edge> out;
  for (Vertex v1 = 0; v1 < nv; ++v1) {
    Vertex v2 = mate[v1];
    if (v2 != boost::graph_traits<Graph>::null_vertex() && v1 < v2) {
      out.push_back (boost::edge (v1, v2, g).first);
    }
  }

  return out;
}


} // anonymous namespace


PairCaloClustersPi0::PairCaloClustersPi0(const std::string& name, ISvcLocator* svcLoc)
    : Gaudi::Algorithm(name, svcLoc) {
  declareProperty("inClusters", m_inClusters, "Input cluster collection");
  declareProperty("reconstructedPi0", m_reconstructedPi0, "Output1: Reconstructed pi0 collection");
  declareProperty("unpairedClusters", m_unpairedClusters, "Output2: Unpaired cluster collection");
  declareProperty("pairedClusters", m_pairedClusters, "Output3: Paired cluster collection");
}

StatusCode PairCaloClustersPi0::initialize() {
  K4_GAUDI_CHECK( Gaudi::Algorithm::initialize() );

  // If there are shapeParameters metadata in the input cluster collection, ship them to the output cluster collections
  auto shapeParameterNames = k4FWCore::getCollectionParameter<std::vector<std::string>>(
                                 m_inClusters.objKey(), edm4hep::labels::ShapeParameterNames, this)
                                 .value_or(std::vector<std::string>{});
  if (shapeParameterNames.size() > 0) {
    k4FWCore::putCollectionParameter(m_pairedClusters.objKey(), edm4hep::labels::ShapeParameterNames,
                                     shapeParameterNames, this);
    k4FWCore::putCollectionParameter(m_unpairedClusters.objKey(), edm4hep::labels::ShapeParameterNames,
                                     shapeParameterNames, this);
  }
  // print pi0 mass window
  info() << "pi0 mass window for cluster pairing: [" << m_massLow << "," << m_massHigh << "] GeV, peak= " << m_massPeak
         << " GeV" << endmsg;

  return StatusCode::SUCCESS;
}

StatusCode PairCaloClustersPi0::execute(const EventContext&) const {
  verbose() << "-------------------------------------------" << endmsg;

  // Get the input collection with clusters
  const edm4hep::ClusterCollection* inClusters = m_inClusters.get();

  // Initialize output clusters
  edm4hep::ReconstructedParticleCollection* reconstructedPi0 = m_reconstructedPi0.createAndPut();
  edm4hep::ClusterCollection* unpairedClusters = m_unpairedClusters.createAndPut();
  edm4hep::ClusterCollection* pairedClusters = m_pairedClusters.createAndPut();

  K4_GAUDI_CHECK( doPairing (*inClusters,
                             *reconstructedPi0,
                             *pairedClusters,
                             *unpairedClusters) );
  std::cout << std::format("aaa {} {} {} {}\n", name(),
                       reconstructedPi0->size(),
                       pairedClusters->size(),
                       unpairedClusters->size());
  std::cout.flush();

  return StatusCode::SUCCESS;
}

/// Cluster pairing
StatusCode PairCaloClustersPi0::doPairing(const edm4hep::ClusterCollection& inClusters,
                                          edm4hep::ReconstructedParticleCollection& reconstructedPi0s,
                                          edm4hep::ClusterCollection& pairedClusters,
                                          edm4hep::ClusterCollection& unpairedClusters) const
{
  size_t nclust = inClusters.size();
  std::vector<TLorentzVector> pairs;
  Graph g = makeGraph (inClusters,
                       pairs,
                       m_minClusterEnergy,
                       m_maxDR,
                       m_massPeak,
                       m_massLow,
                       m_massHigh);
  std::vector<Edge> edges = findEdges (g);

  std::vector<bool> used_clusts (nclust);
  for (const Edge& e : edges) {
    const auto& cl1 = inClusters.at(boost::source(e, g));
    const auto& cl2 = inClusters.at(boost::target(e, g));
    TLorentzVector tlv1 = getTLV (cl1);
    TLorentzVector tlv2 = getTLV (cl2);
    TLorentzVector tlv_pi = tlv1 + tlv2;
    edm4hep::MutableReconstructedParticle this_pi0
      (111, tlv_pi.E(),
       edm4hep::Vector3f(tlv_pi.Px(), tlv_pi.Py(), tlv_pi.Pz()),
       edm4hep::Vector3f(0, 0, 0), 0., tlv_pi.M(), 0.,
       edm4hep::CovMatrix4f());
    this_pi0.addToClusters(cl1);
    this_pi0.addToClusters(cl2);
    reconstructedPi0s.push_back(this_pi0);
    pairedClusters.push_back (cl1.clone());
    pairedClusters.push_back (cl2.clone());
    used_clusts[boost::source(e, g)] = true;
    used_clusts[boost::target(e, g)] = true;
  }

  for (size_t i = 0; i < nclust; ++i) {
    if (!used_clusts[i]) {
      unpairedClusters.push_back (inClusters.at(i).clone());
    }
  }

  return StatusCode::SUCCESS;
}
