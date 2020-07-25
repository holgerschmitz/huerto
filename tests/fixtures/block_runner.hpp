/*
 * block_runner.hpp
 *
 *  Created on: 22 Jan 2020
 *      Author: Holger Schmitz
 */

#ifndef HUERTO_TESTS_FIXTURES_BLOCK_RUNNER_HPP_
#define HUERTO_TESTS_FIXTURES_BLOCK_RUNNER_HPP_

#include <mpi.h>

#include <schnek/variables/block.hpp>
#include <schnek/variables/blockclasses.hpp>
#include <schnek/parser/parser.hpp>
#include <schnek/parser/parsertoken.hpp>

#include <boost/test/unit_test.hpp>

#include <fstream>

class SingleBlockRunner
{
  protected:
    schnek::pBlock block;

    template<class BlockType>
    boost::shared_ptr<BlockType> createBlock(std::string blockName, std::string fileName)
    {
        boost::shared_ptr<BlockType> block;
        try
        {
          schnek::BlockClasses blocks;

          blocks.registerBlock(blockName).setClass<BlockType>();
          std::ifstream in(fileName);
          BOOST_REQUIRE(in);

          schnek::Parser P(blockName, blockName, blocks);

          schnek::pBlock application = P.parse(in);
          application->initAll();

          block = boost::dynamic_pointer_cast<BlockType>(application);
          BOOST_REQUIRE(!!block);
        }
        catch (schnek::ParserError &e)
        {
        BOOST_FAIL(
            "SingleBlockRunner: Parse error in " + e.getFilename() + " at line "
                + boost::lexical_cast<std::string>(e.getLine()) + ": " + e.message);
        }
        catch (schnek::VariableNotInitialisedException &e)
        {
          BOOST_FAIL("SingleBlockRunner: Variable was not initialised: " + e.getVarName());
        }
        catch (schnek::EvaluationException &e)
        {
          BOOST_FAIL("SingleBlockRunner: Error in evaluation: " + e.getMessage());
        }
        catch (schnek::VariableNotFoundException &e)
        {
          BOOST_FAIL("SingleBlockRunner: Error: " + e.getMessage());
        }
        catch (SchnekException &e)
        {
          BOOST_FAIL("SingleBlockRunner: An error occured");
        }
        catch (std::string &err)
        {
          BOOST_FAIL("SingleBlockRunner: FATAL ERROR: " + err);
        }

        return block;
    }
};

// ============================================================================

//template<class BlockType>
//std::shared_ptr<BlockType> SingleBlockRunner::createBlock<BlockType>(std::string blockName, std::string fileName)


#endif /* HUERTO_TESTS_FIXTURES_BLOCK_RUNNER_HPP_ */
