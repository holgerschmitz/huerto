/*
 * task.cpp
 *
 *  Created on: 29 Dec 2020
 *  Author: Holger Schmitz (holger@notjustphysics.com)
 */

#include "task.hpp"

void SimulationTaskRunner::init(const schnek::BlockList& blocks) {
  for (schnek::pBlock child : blocks) {
    SimulationTask *task = dynamic_cast<SimulationTask*>(child.get());
    if (task != NULL) {
      tasks[task->getPhase()].push_back(task);
    }
  }
}

void SimulationTaskRunner::executeTasks(std::string phase) {
  if (tasks.count(phase) > 0) {
    for (SimulationTask *t : tasks[phase]) {
      t->execute();
    }
  }
}
