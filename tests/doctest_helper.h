#pragma once

#define REQUIRE_EQUAL(left, right) REQUIRE((left) == (right))
#define REQUIRE_SMALL(value, eps) REQUIRE((value) < (eps))
