/*
 * border.hpp
 *
 *  Created on: 30 Nov 2020
 *      Author: Holger Schmitz
 */

#ifndef HUERTO_BOUNDARY_BORDER_HPP_
#define HUERTO_BOUNDARY_BORDER_HPP_

#include "../../types.hpp"
#include "../../simulation/simulation_context.hpp"

/**
 * Get the border extent for a source current
 */
bool getBorderExtent(Direction dir,
                     int thickness,
                     int distance,
                     Index &blow,
                     Index &bhigh,
                     bool isH,
                     SimulationContext &context,
                     bool restricted = true);


#endif /* HUERTO_BOUNDARY_BORDER_HPP_ */
