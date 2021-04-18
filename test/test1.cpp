#include "cg_coeff.h"
#define BOOST_TEST_MODULE SqrTests
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(FailTest) {
  BOOST_CHECK_EQUAL(1.0, cg_coeff(1.0, 1.0, 1.0, 1., 2.0, 2.0));
}
