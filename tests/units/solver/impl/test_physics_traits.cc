#include <gtest/gtest.h>

#include "physics_traits.h"
#include "physics_traits_acoustic.h"
#include "physics_traits_elastic.h"
#include "sem_enums.h"

namespace solver {
namespace fe {
namespace test {

using utils::enums::physicType;

TEST(PhysicsTraitsTest, AcousticSpecializationName) {
  EXPECT_STREQ(PhysicsTraits<physicType::kAcoustic>::kName, "Acoustic");
}

TEST(PhysicsTraitsTest, ElasticSpecializationName) {
  EXPECT_STREQ(PhysicsTraits<physicType::kElastic>::kName, "Elastic");
}

// kAcoustoElastic has no explicit specialization → falls through to primary template.
TEST(PhysicsTraitsTest, PrimaryTemplateFallbackName) {
  EXPECT_STREQ(PhysicsTraits<physicType::kAcoustoElastic>::kName, "");
}

TEST(PhysicsTraitsTest, AcousticWavefieldTypeIsConcrete) {
  static_assert(!std::is_void_v<PhysicsTraits<physicType::kAcoustic>::WavefieldType>);
}

TEST(PhysicsTraitsTest, ElasticWavefieldTypeIsConcrete) {
  static_assert(!std::is_void_v<PhysicsTraits<physicType::kElastic>::WavefieldType>);
}

TEST(PhysicsTraitsTest, PrimaryTemplateWavefieldTypeIsVoid) {
  static_assert(std::is_void_v<PhysicsTraits<physicType::kAcoustoElastic>::WavefieldType>);
}

}  // namespace test
}  // namespace fe
}  // namespace solver
