// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include "CommandViz.h"
#include "Sketch.h"
#include <iostream>
#include <vector>
#include <unistd.h>

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

namespace mash {

#ifdef ARCH_32
	#define HASH "MurmurHash3_x86_32"
#else
	#define HASH "MurmurHash3_x64_128"
#endif

std::vector<size_t> bishop(const Sketch & sketch, const Sketch::Reference &ref, ssize_t width, ssize_t height, size_t limit);
void print_matrix(const std::vector<size_t> &matrix, size_t width, size_t height);


CommandViz::CommandViz()
: Command()
{
    name = "viz";
    summary = "Display information about sketch files.";
    description = "Display information about sketch files.";
    argumentString = "<sketch>";
    
    useOption("help");
    // addOption("header", Option(Option::Boolean, "H", "", "Only show header info. Do not list each sketch. Incompatible with -d, -t and -c.", ""));
}

int CommandViz::run() const
{
    
    bool header = false; //options.at("header").active;
    bool tabular = false; //options.at("tabular").active;
    bool counts = false; //options.at("counts").active;
    bool dump = false; //options.at("dump").active;
    
    const string & file = arguments[0];
    
    if ( ! hasSuffix(file, suffixSketch) )
    {
        cerr << "ERROR: The file \"" << file << "\" does not look like a sketch." << endl;
        return 1;
    }
    
	Sketch sketch;
	Sketch::Parameters params;
	params.parallelism = 1;
	
	uint64_t referenceCount;
	
	if ( header )
	{
		referenceCount = sketch.initParametersFromCapnp(arguments[0].c_str());
	}
	else
	{
	    sketch.initFromFiles(arguments, params);
	    referenceCount = sketch.getReferenceCount();
	}

    if ( tabular )
    {
    	cout << "#Hashes\tLength\tID\tComment" << endl;
    }
    else
    {
	}
	
	string alphabet;
	sketch.getAlphabetAsString(alphabet);
	
	cout << "Header:" << endl;
	cout << "  Hash function (seed):          " << HASH << " (" << sketch.getHashSeed() << ")" << endl;
	cout << "  K-mer size:                    " << sketch.getKmerSize() << " (" << (sketch.getUse64() ? "64" : "32") << "-bit hashes)" << endl;
	cout << "  Alphabet:                      " << alphabet << (sketch.getNoncanonical() ? "" : " (canonical)") << (sketch.getPreserveCase() ? " (case-sensitive)" : "") << endl;
	cout << "  Target min-hashes per sketch:  " << sketch.getMinHashesPerWindow() << endl;
	cout << "  Sketches:                      " << referenceCount << endl;

	
	for (size_t j = 0; j < sketch.getReferenceCount(); j++) {
        const Sketch::Reference & ref = sketch.getReference(j);

        size_t height = 40;
        size_t width  = 61;

        for (int i = 0; i < 1000; i += 10)
        {
			auto matrix = bishop(sketch, ref, width, height, i);
			print_matrix(matrix, width, height);
			int check = usleep(10000);
        	printf("\033[%dA", height / 2);
        }

        for (int i = 0; i < height / 2; ++i)
        {
        	printf("\n");
        }
	}

    if ( ! header )
    {
        vector<vector<string>> columns(4);
        
        if ( ! tabular )
        {
			cout << endl;
			cout << "Sketches:" << endl;
		
			columns[0].push_back("[Hashes]");
			columns[1].push_back("[Length]");
			columns[2].push_back("[ID]");
			columns[3].push_back("[Comment]");
		}
        
        for ( uint64_t i = 0; i < sketch.getReferenceCount(); i++ )
        {
            const Sketch::Reference & ref = sketch.getReference(i);
            
            if ( tabular )
            {
            	cout
            		<< ref.hashesSorted.size() << '\t'
            		<< ref.length << '\t'
            		<< ref.name << '\t'
            		<< ref.comment << endl;
            }
            else
            {
				columns[0].push_back(std::to_string(ref.hashesSorted.size()));
				columns[1].push_back(std::to_string(ref.length));
				columns[2].push_back(ref.name);
				columns[3].push_back(ref.comment);
			}
        }
        
        if ( ! tabular )
        {
	        printColumns(columns, 2, 2, "-", 0);
	    }
    }
    
    return 0;
}

std::vector<size_t> bishop(const Sketch & sketch, const Sketch::Reference &ref, ssize_t width, ssize_t height, size_t limit) {

	bool use64 = sketch.getUse64();

	ssize_t pos_x = width - 3;
	ssize_t pos_y = height - 3;
	size_t K = sketch.getKmerSize();

	auto matrix = std::vector<size_t>(width * height);

	bool torus = true;

	auto walk = [&](int direction) {
		int x_diff[4] = {-1,1,-1,1};
		int y_diff[4] = {-1,-1,1,1};
		pos_x += x_diff[direction & 3];
		pos_y += y_diff[direction & 3];

		if (torus){
			if (pos_x < 0) pos_x = width - 1;
			if (pos_x >= width) pos_x = 0;
			if (pos_y < 0) pos_y = height - 1;
			if (pos_y >= height) pos_y = 0;
		} else {
			pos_x = std::max((ssize_t)0, std::min(pos_x, width - 1));
			pos_y = std::max((ssize_t)0, std::min(pos_y, height - 1));
		}
	};

	limit = std::min(limit, (size_t)ref.hashesSorted.size());
	if (limit == 0) limit = ref.hashesSorted.size();
	for ( size_t j = 0; j < limit; j++ )
	{
		auto hash_obj = ref.hashesSorted.at(j);
		size_t hash = use64 ? hash_obj.hash64 : hash_obj.hash32;
		// std::cerr << hash << std::endl;
		size_t bits = use64 ? 64 : 32;

		for ( int k = 0; k < bits / 2; ++k)
		{
			matrix[pos_y * width + pos_x]++;
			walk(hash);
			hash >>= 2;
		}
	}

	return matrix;
}

size_t ilog(size_t num) {
	size_t ret = 0;
	while (num) {
		ret++;
		num >>= 1;
	}
	return ret;
}

void print_matrix(const std::vector<size_t> & matrix, size_t width, size_t height) {
	const char symbols[] = {
		' ', '.', '_', '-',
		'o', '+',
		'=', '*', 'B', 'O',
		'X', '@', '%', '&',
		'#', '/', '^', 'S', 'E', '?'
	};

	const char *upper_half = "▀";

	size_t i;
    size_t j;
    size_t temp;
    int num_symbols = 20;

    for (i = 0; i < height; i+=2) {
        // printf("|");
        for (j = 0; j < width; j++) {
            // std::cerr << temptemp << std::endl;
            printf("%c[38;5;%dm", '\x1b', (matrix[j + width * i]) + 16);
            printf("%c[48;5;%dm", '\x1b', (matrix[j + width * (i + 1)]) + 16);
            // printf("%c", symbols[logo]);
            printf("%s", upper_half);
            printf("%c[0m", '\x1b');
        }
        printf("\n");
    }
}

} // namespace mash