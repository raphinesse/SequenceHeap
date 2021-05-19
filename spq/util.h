// This file replaces the original util.h by Peter Sanders
// Copyright 2020–2021 Raphael von der Grün

#pragma once

#include <cassert>

// Delegate all custom Assert* macros to assert
#define Assert(C) assert(C)
#define Assert2(C) assert(C)

// Disable Debug4 unconditionally
#define Debug4(C)
