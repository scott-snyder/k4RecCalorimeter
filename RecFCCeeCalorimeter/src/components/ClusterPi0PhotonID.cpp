// k4FWCore / k4Interface
#include "k4FWCore/DataHandle.h"
#include "k4FWCore/MetadataUtils.h"
#include "k4FWCore/Transformer.h"

// Gaudi
#include "GaudiKernel/ToolHandle.h"

// EDM4hep
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/ParticleIDCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/Vector3f.h"
#include "edm4hep/utils/ParticleIDUtils.h"

// STL
#include <optional>
#include <string>
#include <vector>

/** @struct ClusterPi0PhotonID
 *
 * Gaudi MultiTransformer that identifies photons and pi0 candidates from a
 * collection of clusters with associated classification scores (e.g. produced by
 * TRAPPISTPi0PhotonInference).
 * For each input cluster, the algorithm retrieves the classification score from the cluster
 * shape parameters using a configurable shape parameter name. If the score is greater than
 * a user-defined threshold, the cluster is identified as a pi0 candidate; otherwise,
 * it is identified as a photon candidate. A corresponding ReconstructedParticle is created
 * and added to the corresponding output collection together with a ParticleID object that contains
 * the identification information.
 *
 * input: edm4hep::ClusterCollection
 * output: edm4hep::ReconstructedParticleCollection, edm4hep::ParticleIDCollection
 *
 *  @author Jacopo Fanini
 *  @date   2026-07
 *
 */

struct ClusterPi0PhotonID final
    : k4FWCore::MultiTransformer<std::tuple<edm4hep::ReconstructedParticleCollection, edm4hep::ParticleIDCollection>(
          const edm4hep::ClusterCollection&)> {
public:
  ClusterPi0PhotonID(const std::string& name, ISvcLocator* svcLoc)
      : MultiTransformer(name, svcLoc, {KeyValues("inClusters", {"clustersWithScore"})},
                         {KeyValues("outParticles", {"IdentifiedParticles"}),
                          KeyValues("outParticleIDs", {"IdentifiedParticleIDs"})}) {}

  StatusCode initialize() override {
    // Retrieve shape parameter names for input cluster collection
    auto inputKey = this->inputLocations("inClusters")[0];
    auto shapeParameterNames =
        k4FWCore::getCollectionParameter<std::vector<std::string>>(inputKey, edm4hep::labels::ShapeParameterNames, this)
            .value_or(std::vector<std::string>{});
    debug() << "Input cluster has " << shapeParameterNames.size() << " names in metadata" << endmsg;
    // Find and save the index of the shape parameter to be used
    auto it = std::find(shapeParameterNames.begin(), shapeParameterNames.end(), m_shapeParameterScoreName.value());
    if (it != shapeParameterNames.end()) {
      m_scoreIndex = std::distance(shapeParameterNames.begin(), it);
      debug() << "Feature " << m_shapeParameterScoreName.value() << " found in position " << m_scoreIndex.value()
              << " of collection metadata" << endmsg;
    } else {
      throw std::runtime_error("Feature " + m_shapeParameterScoreName.value() + " not found in collection metadata");
    }
    // Add the ParticleID algorithm name to the output collection metadata
    auto outputKey = this->outputLocations("outParticleIDs")[0];
    k4FWCore::putCollectionParameter(outputKey, edm4hep::labels::PIDAlgoName, m_particleIDAlgorithmName.value(), this);
    m_pidMeta = m_particleIDAlgorithmName.value();
    debug() << "ParticleID algorithm name set to " << m_pidMeta.algoName << endmsg;
    debug() << "ParticleID algorithm type set to " << m_pidMeta.algoType() << endmsg;

    return StatusCode::SUCCESS;
  }

  std::tuple<edm4hep::ReconstructedParticleCollection, edm4hep::ParticleIDCollection>
  operator()(const edm4hep::ClusterCollection& clustersWithScore) const override {
    edm4hep::ReconstructedParticleCollection outputParticles;
    edm4hep::ParticleIDCollection outputParticleIDs;
    for (const auto& cluster : clustersWithScore) {
      // Retrieve the score from the shape parameters and if it exceeds the threshold
      // create and fill a new ReconstructedParticle
      float clusterScore = cluster.getShapeParameters(m_scoreIndex.value());
      debug() << "Cluster score: " << clusterScore << ", threshold: " << m_threshold.value() << endmsg;
      if (clusterScore > m_threshold.value()) {
        buildParticleAndID(outputParticles, outputParticleIDs, cluster, 111, m_pi0Mass.value(), clusterScore);
      } else {
        buildParticleAndID(outputParticles, outputParticleIDs, cluster, 22, 0.0, (1. - clusterScore));
      }
    }
    return std::make_tuple(std::move(outputParticles), std::move(outputParticleIDs));
  }

  StatusCode finalize() override { return StatusCode::SUCCESS; }

private:
  std::optional<std::size_t> m_scoreIndex = std::nullopt;
  edm4hep::utils::ParticleIDMeta m_pidMeta;

  Gaudi::Property<std::string> m_shapeParameterScoreName{
      this, "ClusterScoreShapeParameterName", "ClusterTRAPPISTScore",
      "Name of the shape parameter to be used for particle identification"};

  Gaudi::Property<float> m_threshold{this, "Threshold", 0.5,
                                     "Threshold for particle identification based on the shape parameter value"};

  Gaudi::Property<float> m_pi0Mass{this, "Pi0Mass", 0.1349768, "Mass of the Pi0 particle"};

  Gaudi::Property<std::string> m_particleIDAlgorithmName{this, "ParticleIDAlgorithmName", "TRAPPIST",
                                                         "Name of the algorithm used for particle identification"};

  edm4hep::Vector3f calculateMomentum(double energy, const edm4hep::Vector3f& position,
                                      const edm4hep::Vector3f& origin) const {
    double dirx = position.x - origin.x;
    double diry = position.y - origin.y;
    double dirz = position.z - origin.z;
    double quadsum_dir = std::sqrt(dirx * dirx + diry * diry + dirz * dirz);
    double px = energy * dirx / quadsum_dir;
    double py = energy * diry / quadsum_dir;
    double pz = energy * dirz / quadsum_dir;
    return edm4hep::Vector3f(px, py, pz);
  }

  void buildParticleAndID(edm4hep::ReconstructedParticleCollection& partCollection,
                          edm4hep::ParticleIDCollection& idCollection, const edm4hep::Cluster& cluster, int pdg,
                          float mass, float pidScore) const {

    auto particle = partCollection.create();
    auto particleID = idCollection.create();

    particleID.setParticle(particle);
    particleID.setPDG(pdg);
    particleID.setAlgorithmType(m_pidMeta.algoType());
    particleID.setLikelihood(pidScore);

    particle.addToClusters(cluster);

    const double energy = cluster.getEnergy();
    particle.setEnergy(energy);

    particle.setPDG(pdg);
    particle.setCharge(0);
    particle.setMass(mass);
    particle.setGoodnessOfPID(pidScore);

    edm4hep::Vector3f position{cluster.getPosition().x, cluster.getPosition().y, cluster.getPosition().z};
    // Assuming particle is coming from the IP, will change this once the tool to retrieve cluster direction is
    // implemented
    edm4hep::Vector3f origin{0, 0, 0};
    particle.setMomentum(calculateMomentum(energy, position, origin));
  }
};

DECLARE_COMPONENT(ClusterPi0PhotonID)