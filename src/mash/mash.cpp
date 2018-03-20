// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include "CommandBounds.h"
#include "CommandList.h"
#include "CommandSketch.h"
#include "CommandFind.h"
#include "CommandDistance.h"
#include "CommandMatrix.h"
#include "CommandScreen.h"
#include "CommandContain.h"
#include "CommandInfo.h"
#include "CommandViz.h"
#include "CommandPaste.h"

int main(int argc, const char ** argv)
{
    mash::CommandList commandList("mash");
    
    commandList.addCommand(new mash::CommandSketch());
    //commandList.addCommand(new CommandFind());
    commandList.addCommand(new mash::CommandDistance());
    commandList.addCommand(new mash::CommandMatrix());
    commandList.addCommand(new mash::CommandScreen());
#ifdef COMMAND_WITHIN
    commandList.addCommand(new mash::CommandContain());
#endif
#ifdef COMMAND_FIND
	commandList.addCommand(new mash::CommandFind());
#endif
    commandList.addCommand(new mash::CommandInfo());
    commandList.addCommand(new mash::CommandViz());
    commandList.addCommand(new mash::CommandPaste());
    commandList.addCommand(new mash::CommandBounds());
    
    return commandList.run(argc, argv);
}
