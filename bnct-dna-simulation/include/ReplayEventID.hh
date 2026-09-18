#pragma once

#include "CommandLineParser.hh"
#include "G4Types.hh"
#include <cstdlib>

inline G4int ReplayEventID(G4int localID)
{
  static thread_local const G4int start = [] {
    auto *command = G4DNAPARSER::CommandLineParser::GetParser()->GetCommandIfActive("-start-record");
    return command ? static_cast<G4int>(std::strtol(command->GetOption().c_str(), nullptr, 10)) : 0;
  }();
  return localID + start;
}
