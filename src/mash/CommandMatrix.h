// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef INCLUDED_CommandMatrix
#define INCLUDED_CommandMatrix

#include "Command.h"
#include "Sketch.h"

namespace mash {

class CommandMatrix : public Command
{
public:
    
    struct CompareInput
    {
        CompareInput(const Sketch & sketchRefNew, const Sketch & sketchQueryNew, uint64_t indexRefNew, uint64_t indexQueryNew, uint64_t pairCountNew, const Sketch::Parameters & parametersNew, double maxPValueNew)
            :
            sketchRef(sketchRefNew),
            sketchQuery(sketchQueryNew),
            indexRef(indexRefNew),
            indexQuery(indexQueryNew),
            pairCount(pairCountNew),
            parameters(parametersNew),
            maxPValue(maxPValueNew)
            {}
        
        const Sketch & sketchRef;
        const Sketch & sketchQuery;
        
        uint64_t indexRef;
        uint64_t indexQuery;
        uint64_t pairCount;
        
        const Sketch::Parameters & parameters;
        double maxPValue;
    };
    
    struct CompareOutput
    {
        CompareOutput(const Sketch & sketchRefNew, const Sketch & sketchQueryNew, uint64_t indexRefNew, uint64_t indexQueryNew, uint64_t pairCountNew)
            :
            sketchRef(sketchRefNew),
            sketchQuery(sketchQueryNew),
            indexRef(indexRefNew),
            indexQuery(indexQueryNew),
            pairCount(pairCountNew)
        {
            pairs = new PairOutput[pairCount];
        }
        
        ~CompareOutput()
        {
            delete [] pairs;
        }
        
        struct PairOutput
        {
            uint64_t numer;
            uint64_t denom;
            double distance;
            double pValue;
            bool pass;
        };
        
        const Sketch & sketchRef;
        const Sketch & sketchQuery;
        
        uint64_t indexRef;
        uint64_t indexQuery;
        uint64_t pairCount;
        
        PairOutput * pairs;
    };
    
    CommandMatrix();
    
    int run() const; // override
    
private:
    
    void writeOutput(CompareOutput * output, bool table) const;
};

CommandMatrix::CompareOutput * compare(CommandMatrix::CompareInput * input);
void compareSketches(CommandMatrix::CompareOutput::PairOutput * output, const Sketch::Reference & refRef, const Sketch::Reference & refQry, uint64_t sketchSize, int kmerSize, double kmerSpace, double maxPValue);
double pValue(uint64_t x, uint64_t lengthRef, uint64_t lengthQuery, double kmerSpace, uint64_t sketchSize);

} // namespace mash

#endif