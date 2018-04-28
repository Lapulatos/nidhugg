#include <config.h>
#ifdef HAVE_BOOST_UNIT_TEST_FRAMEWORK

#include "2SAT.h"
#include "Event.h"

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(twoSAT_test)

VClock<IPid> clock;
DCEvent e1(IID<IPid>(0, 0));
DCEvent e2(IID<IPid>(0, 1));
DCEvent e3(IID<IPid>(0, 2));
DCEvent e4(IID<IPid>(0, 3));

BOOST_AUTO_TEST_CASE(satisfiable_1)
{
  Events2SAT sat;

  BOOST_CHECK(sat.addAssignment(&e1, &e2, true));
  BOOST_CHECK(sat.addAssignment(&e1, &e2, true));
  BOOST_CHECK(sat.addAssignment(&e3, &e4, true));
  BOOST_CHECK(sat.addAssignment(&e2, &e3, true));
}

BOOST_AUTO_TEST_CASE(satisfiable_2)
{
  Events2SAT sat;

  BOOST_CHECK(sat.addImplication(&e1, &e2, &e2, &e3));
  BOOST_CHECK(sat.addImplication(&e2, &e3, &e1, &e2));

  BOOST_CHECK(sat.solve());
  BOOST_CHECK(sat.getValuation(&e1, &e2) == true);
  BOOST_CHECK(sat.getValuation(&e2, &e3) == true);
}

BOOST_AUTO_TEST_CASE(satisfiable_3)
{
  Events2SAT sat;

  BOOST_CHECK(sat.addImplication(&e1, &e2, &e3, &e4));
  BOOST_CHECK(sat.addImplication(&e1, &e2, &e2, &e3));
  BOOST_CHECK(sat.addImplication(&e2, &e3, &e3, &e4));

  BOOST_CHECK(sat.solve());
}

BOOST_AUTO_TEST_CASE(satisfiable_4)
{
  Events2SAT sat;

  BOOST_CHECK(sat.addImplication(&e1, &e2, &e3, &e4));
  BOOST_CHECK(sat.addImplication(&e1, &e2, &e2, &e3));
  BOOST_CHECK(sat.addImplication(&e2, &e3, &e3, &e4));
  BOOST_CHECK(sat.addImplication(&e3, &e4, &e2, &e3));

  BOOST_CHECK(sat.solve());
}

BOOST_AUTO_TEST_CASE(satisfiable_5)
{
  Events2SAT sat;

  BOOST_CHECK(sat.addAssignment(&e2, &e3, true));
  BOOST_CHECK(sat.addImplication(&e1, &e2, &e3, &e4));
  BOOST_CHECK(sat.addImplication(&e1, &e2, &e2, &e3));
  BOOST_CHECK(sat.addImplication(&e2, &e3, &e3, &e4));
  BOOST_CHECK(sat.addImplication(&e3, &e4, &e2, &e3));

  BOOST_CHECK(sat.solve());
}

BOOST_AUTO_TEST_CASE(satisfiable_6)
{
  Events2SAT sat;

  BOOST_CHECK(sat.addAssignment(&e2, &e3, true));
  BOOST_CHECK(sat.addImplication(&e1, &e2, &e2, &e3));
  BOOST_CHECK(sat.addImplication(&e2, &e3, &e3, &e4));
  BOOST_CHECK(sat.addImplication(&e3, &e4, &e1, &e2));

  BOOST_CHECK(sat.solve());
  BOOST_CHECK(sat.getValuation(&e1, &e2) == true);
  BOOST_CHECK(sat.getValuation(&e2, &e3) == true);
  BOOST_CHECK(sat.getValuation(&e3, &e4) == true);
}

BOOST_AUTO_TEST_CASE(satisfiable_7)
{
  Events2SAT sat;

  /* suppose:
        e1
       /  \
      e2  e3
       \  /
        e4
  */

  BOOST_CHECK(sat.addAssignment(&e1, &e2, true));
  BOOST_CHECK(sat.addAssignment(&e1, &e3, true));
  BOOST_CHECK(sat.addAssignment(&e2, &e4, true));
  BOOST_CHECK(sat.addAssignment(&e3, &e4, true));

  // antisymmetry
  for (auto ev : {&e1, &e2, &e3, &e4})
    for (auto ev2 : {&e1, &e2, &e3, &e4})
        if (ev != ev2)
            BOOST_CHECK(sat.addImplication(ev, ev2, ev2, ev, true, false));

  BOOST_CHECK(sat.solve());
}

// same as satisfiable_7, but we add assignments after formulas
BOOST_AUTO_TEST_CASE(satisfiable_8)
{
  Events2SAT sat;

  /* suppose:
        e1
       /  \
      e2  e3
       \  /
        e4
  */

  // antisymmetry
  for (auto ev : {&e1, &e2, &e3, &e4})
    for (auto ev2 : {&e1, &e2, &e3, &e4})
        if (ev != ev2)
            BOOST_CHECK(sat.addImplication(ev, ev2, ev2, ev, true, false));

  BOOST_CHECK(sat.addAssignment(&e1, &e2, true));
  BOOST_CHECK(sat.addAssignment(&e1, &e3, true));
  BOOST_CHECK(sat.addAssignment(&e2, &e4, true));
  BOOST_CHECK(sat.addAssignment(&e3, &e4, true));

  BOOST_CHECK(sat.solve());
}

BOOST_AUTO_TEST_CASE(satisfiable_9)
{
  Events2SAT sat;

  /* suppose:
      e2-->e3
  */

  BOOST_CHECK(sat.addAssignment(&e2, &e3, true));

  // acyclic transitivity
  BOOST_CHECK(sat.addImplication(&e1, &e2, &e1, &e3));
  BOOST_CHECK(sat.addImplication(&e3, &e4, &e2, &e4));
  BOOST_CHECK(sat.addAssignment(&e1, &e2, true));

  BOOST_CHECK(sat.solve());
  BOOST_CHECK(sat.getValuation(&e1, &e3) == true);
}

BOOST_AUTO_TEST_CASE(satisfiable_10)
{
  Events2SAT sat;

  /* suppose:
         e1
         |
     e2  e3
         |
         e4
  */

  BOOST_CHECK(sat.addAssignment(&e1, &e3, true));
  BOOST_CHECK(sat.addAssignment(&e3, &e4, true));
  // transitive edge
  BOOST_CHECK(sat.addAssignment(&e1, &e4, true));

  // antisymmetry
  for (auto ev : {&e1, &e2, &e3, &e4}) {
    for (auto ev2 : {&e1, &e2, &e3, &e4}) {
        if (ev != ev2) {
            BOOST_CHECK(sat.addImplication(ev, ev2, ev2, ev, true, false));
            BOOST_CHECK(sat.addImplication(ev2, ev, ev, ev2, false, true));
        }
    }
  }

  // some transitivity
  BOOST_CHECK(sat.addImplication(&e2, &e1, &e2, &e3));
  BOOST_CHECK(sat.addImplication(&e2, &e3, &e2, &e4));
  BOOST_CHECK(sat.addImplication(&e3, &e2, &e1, &e2));
  BOOST_CHECK(sat.addImplication(&e4, &e2, &e3, &e2));
  // this should lead to adding an edge e4->e2 (from antisymmetry)
  // and from the transitivity to {e1,e3}->e2
  //BOOST_CHECK(sat.addAssignment(&e2, &e4, false));
  BOOST_CHECK(sat.addAssignment(&e4, &e2, true));

  BOOST_CHECK(sat.solve());
  BOOST_CHECK(sat.getValuation(&e2, &e4) == false);
  BOOST_CHECK(sat.getValuation(&e4, &e2) == true);
  BOOST_CHECK(sat.getValuation(&e1, &e2) == true);
  BOOST_CHECK(sat.getValuation(&e2, &e1) == false);
  BOOST_CHECK(sat.getValuation(&e3, &e2) == true);
  BOOST_CHECK(sat.getValuation(&e2, &e3) == false);
}

// simplified satisfiable_10 (due to a bug finding)
BOOST_AUTO_TEST_CASE(satisfiable_11)
{
  Events2SAT sat;

  /* suppose:
         e1
         |
     e2  e3
  */

  BOOST_CHECK(sat.addAssignment(&e1, &e3, true));

  // some transitivity
  BOOST_CHECK(sat.addImplication(&e2, &e1, &e2, &e3));
  BOOST_CHECK(sat.addImplication(&e3, &e2, &e1, &e2));
  // this should lead to adding e1->e2
  BOOST_CHECK(sat.addAssignment(&e3, &e2, true));

  BOOST_CHECK(sat.solve());
  BOOST_CHECK(sat.getValuation(&e3, &e2) == true);
  BOOST_CHECK(sat.getValuation(&e1, &e2) == true);
}

BOOST_AUTO_TEST_CASE(tautology_1)
{
  Events2SAT sat;

  BOOST_CHECK(sat.addImplication(&e1, &e2, &e1, &e2));
  BOOST_CHECK(sat.solve());
}

BOOST_AUTO_TEST_CASE(contradiction_1)
{
  Events2SAT sat;

  BOOST_CHECK(sat.addAssignment(&e1, &e2, true));
  BOOST_CHECK(!sat.addImplication(&e1, &e2, &e1, &e2, true, false));
  //BOOST_CHECK(!sat.solve());
}

BOOST_AUTO_TEST_CASE(contradiction_2)
{
  Events2SAT sat;

  BOOST_CHECK(sat.addAssignment(&e1, &e2, true));
  BOOST_CHECK(!sat.addAssignment(&e1, &e2, false));
  //BOOST_CHECK(!sat.solve());
}

BOOST_AUTO_TEST_CASE(unsatisfiable_2)
{
  Events2SAT sat;

  BOOST_CHECK(sat.addImplication(&e1, &e2, &e3, &e4));
  BOOST_CHECK(sat.addImplication(&e1, &e2, &e2, &e3));
  BOOST_CHECK(sat.addImplication(&e2, &e3, &e3, &e4));
  BOOST_CHECK(sat.addImplication(&e3, &e4, &e2, &e3));

  // make it unsat
  BOOST_CHECK(sat.addAssignment(&e2, &e3, true));
  BOOST_CHECK(!sat.addImplication(&e2, &e3, &e2, &e3, true, false));
  //BOOST_CHECK(!sat.solve());
}

BOOST_AUTO_TEST_CASE(unsatisfiable_3)
{
  Events2SAT sat;

  // X23 and ~X23 are on the same SCC component,
  // so this formula must be unsatisfiable
  BOOST_CHECK(sat.addImplication(&e3, &e4, &e2, &e3));
  BOOST_CHECK(sat.addImplication(&e2, &e3, &e2, &e3, true, false));
  BOOST_CHECK(sat.addImplication(&e2, &e3, &e3, &e4, false, true));

  BOOST_CHECK(!sat.solve());
}

BOOST_AUTO_TEST_CASE(sat_1)
{
  Events2SAT sat;

  // create a chain of implications
  BOOST_CHECK(sat.addImplication(&e1, &e2, &e2, &e3));
  BOOST_CHECK(sat.addImplication(&e2, &e3, &e3, &e4));
  BOOST_CHECK(sat.addImplication(&e3, &e4, &e4, &e1));

  // make the first variable in the chain 'true'
  BOOST_CHECK(sat.addAssignment(&e1, &e2, true));
  BOOST_CHECK(sat.solve());

  // check that the reset of variables is also true
  BOOST_CHECK(sat.getValuation(&e1, &e2) == true);
  BOOST_CHECK(sat.getValuation(&e2, &e3) == true);
  BOOST_CHECK(sat.getValuation(&e3, &e4) == true);
  BOOST_CHECK(sat.getValuation(&e4, &e1) == true);
}

BOOST_AUTO_TEST_SUITE_END()

#endif // HAVE_BOOST

