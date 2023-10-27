/** @file gnuplotTest.cpp

    @brief Plots data using the gsGnuplot class

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Scholz
*/

#include <gismo.h>
#include <gismo_dev.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    bool plot = false;
    gsCmdLine cmd("Creates an example plot with gnuplot");
    cmd.addSwitch("show", "Open plot in okular", plot);

    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    //Create some example data
    gsVector<real_t> xData1(5), yData1(5), xData2(8), yData2(8);
    xData1 << 1, 2, 3, 4, 5;
    yData1 << 1, 2, 3, 4, 5;
    xData2 << 1, 2, 3, 4, 5, 6, 7, 8;
    yData2 << 5, 4, 3, 2, 1, 0, -1, -2;

    //Make a gsGnuplot object
    gsGnuplot<real_t> plotter;
    //Push the data to the plot object and set some optional style options, see Gnuplot documentation for all options
    plotter.addData(xData1, yData1, "first plot", "lt 1 lw 2 pt 7 lc 'blue' ", "linespoints");
    plotter.addData(xData2, yData2, "second plot", "lt 1 lw 2 pt 7 lc 'red' ", "linespoints");

    //Make a simple plot using the standard gnuplot parameters
    plotter.plot("plot1");

    //Put some options to make the plot look nicer.
    plotter.plot("plot2", "set xlabel 'x'\n"
            "set ylabel 'y'\n"
            "set title 'Title of the plot' \n"
            "set key center right box\n"
            "set grid back lc rgb '#808080' lt 0 lw 1\n"
            "set border 3 back lc rgb '#808080' lt 1\n"
            "set tics font ', 16' ");

    if (plot) // Show the resulting plots ?
    {
        gsFileManager::open("plot1.pdf");
        gsFileManager::open("plot2.pdf");
    }
    
    return EXIT_SUCCESS;
}
