#pragma once

#include <climits>
#include <cstdint>

// Mix the run seed and original record number into a stable, positive seed.
inline long ReplaySeed(std::uint32_t runSeed, std::uint32_t record)
{
  static_assert(sizeof(long) >= 8, "Replay seeds require 64-bit long");
  std::uint64_t value = (static_cast<std::uint64_t>(runSeed) << 32) | record;
  value += UINT64_C(0x9e3779b97f4a7c15);
  value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
  value ^= value >> 31;
  const long seed = static_cast<long>(value & LONG_MAX);
  return seed == 0 ? 1 : seed;
}
