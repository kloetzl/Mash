// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef INCLUDED_CommandViz
#define INCLUDED_CommandViz

#include "Command.h"
#include "Sketch.h"

namespace mash {

class CommandViz : public Command
{
public:
    
    CommandViz();
    
    int run() const; // override
    
private:
	
	int printCounts(const Sketch & sketch) const;
	int writeJson(const Sketch & sketch) const;
};

} // namespace mash

#endif
