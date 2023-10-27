/** @file gsGnuplot.h

    @brief Plots data using gnuplot

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Scholz
*/

#pragma once

//#include <fstream>

namespace gismo
{

/** @brief Collects datasets in the form of gsVectors and creates a script for plotting using gnuplot.
 * Every given dataset is plotted into the same grid. All options of gnuplot can be used. Gnuplot's documentation
 * can be found in http://www.gnuplot.info/documentation.html
 *
 * @tparam T The arithmetic type
 */
template<class T = real_t>
class gsGnuplot
{
public:
    typedef std::deque<gsVector<T> > Container;
    typedef typename Container::const_iterator citer;
    typedef std::deque<std::string> stringContainer;

public:
    /** @brief Add a dataset to be plotted. A title can be specified.
     *
     * @param xData
     * @param yData
     * @param title Title for the plot key
     * @param linestyle  Gnuplot line style options
     * @param with Type of data representations, e.g. "lines", "points" or "linespoints"
     */
    void addData(gsVector <T> xData, gsVector <T> yData, std::string title = "", std::string linestyle = "",
                 std::string with = "lp")
    {
        GISMO_ASSERT(xData.size() == yData.size(), "Data sets should have the same length");

        //Add title
        m_titles.push_back(give(title));
        m_linestyles.push_back(give(linestyle));
        m_with.push_back(give(with));
        m_xData.push_back(give(xData));
        m_yData.push_back(give(yData));
    }

    /** @brief Creates a script for gnuplot creates plot in pdf format.
     *
     * @param filename The name used for the datafile, gnuplot script and output.
     * @param options Gnuplot options to be set before the "plot" command
     */
    void plot(const std::string &filename, const std::string &options = "", std::streamsize precision = 8)
    {
        //Make the filenames
        std::string datafileName, plotscriptName, pdfName, plotCommand;
        datafileName = filename + ".csv";
        plotscriptName = filename + "Plotscript";
        pdfName = filename + ".pdf";
        plotCommand = "gnuplot " + plotscriptName;

        //gsInfo << datafileName.str() << "\n" << plotscriptName.str() << "\n" << pdfName.str() << "\n" << plotCommand.str() << "\n";

        //Write the data to a file
        std::ofstream file;
        file.open(datafileName.c_str());
        file.precision(precision);

        //Find the maximum data length
        index_t maxSize = 0;
        for (citer itX = m_xData.begin(); itX != m_xData.end(); itX++)
        {
            maxSize = std::max<index_t>(maxSize, itX->size());
        }

        for (index_t i = 0; i < maxSize; i++)
        {
            //Iterate over all datasets
            citer itX = m_xData.begin();
            citer itY = m_yData.begin();
            for (; itX != m_xData.end(); itX++)
            {
                if (itX->size() > i)
                {
                    file << (*itX)(i) << ",\t" << (*itY)(i);
                }
                else
                {
                    file << ",\t";
                }
                itY++;

                if (itY != m_yData.end()) file << ",\t";
                else file << "\n";
            }
        }
        file.close();

        //Create the plot script
        file.open(plotscriptName.c_str());
        file << "#!/usr/bin/gnuplot\n"
             << "set terminal pdf enhanced color size 11.7, 8.3\n"
             << "set datafile separator \",\"\n"
             << "set datafile commentschars '#'\n"
             << "set output \"" << pdfName << "\"\n\n"
             << "#User-defined options: \n" << options << "\n#***\n\n"
             << "plot ";
        //Create a plot command for each data set
        for (size_t k = 0; k != m_xData.size(); ++k)
        {
            file << "\"" << datafileName << "\" u " << 1 + 2 * k << ":" << 2 + 2 * k << " w "
                 << m_with[k] << " t \"" << m_titles[k] << "\" " << m_linestyles[k]
                 << (k + 1 == m_xData.size() ? "\n" : ",");
        }
        file.close();

        //Make the plot
        //GISMO_ENSURE(0==std::system("gnuplot --version"), "Gnuplot is not available");
        const int result = std::system(plotCommand.c_str());
        if (0!=result)
        {
            gsWarn<<"std::system: Plotting command `"
                  <<plotCommand<<"` has failed.\n";
        }
    }

private:
    /// @brief The titles of each dataset
    stringContainer m_titles;
    stringContainer m_linestyles;
    stringContainer m_with;
    Container m_xData, m_yData;
};
}
