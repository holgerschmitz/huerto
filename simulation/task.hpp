/*
 * task.hpp
 *
 *  Created on: 29 Dec 2020
 *  Author: Holger Schmitz (holger@notjustphysics.com)
 */

#ifndef HUERTO_SIMULATION_TASK_HPP_
#define HUERTO_SIMULATION_TASK_HPP_

#include <schnek/variables/block.hpp>

#include <string>
#include <map>
#include <list>

/**
 * Interface for a class that runs tasks at a specified phase in the simulation cycle
 */
class SimulationTask {
  public:
    virtual ~SimulationTask() {}
    virtual std::string getPhase() = 0;
    virtual void execute() = 0;
};

class SimulationTaskRunner {
  private:
    typedef std::list<SimulationTask*> TaskList;
    std::map<std::string, TaskList> tasks;
  public:
    void init(const schnek::BlockList& blocks);
    void executeTasks(std::string phase);
};

#endif /* HUERTO_SIMULATION_TASK_HPP_ */
