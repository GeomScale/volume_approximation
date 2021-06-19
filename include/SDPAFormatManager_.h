// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLESTI_SDPA_FORMAT_MANAGER_H
#define VOLESTI_SDPA_FORMAT_MANAGER_H


#include "convex_bodies/spectrahedra/spectrahedron.h"

#include <string>
#include <sstream>


/// Reads/writes files according to the SDPA format for sdps.
/// Currently supported Format:
///
///
/// <objective function vector>
/// <i-th row of j-th matrix>
///
/// For example:
/// 2
/// 3
/// 1 1
///
///
/// \tparam NT Numerical Type
template <typename NT>
class SdpaFormatManager {
private:
    typedef std::string::iterator string_it;
    typedef std::list<NT> listVector;

    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;

    /// Return the first non white space/tab character and advance the iterator one position
    /// @param[in, out] it Current position
    /// @param[in] end End of string
    /// @return First non white space/tab character
    char consumeSymbol(string_it &at, string_it & end) {
        while (at != end) {
            if (*at != ' ' && *at != '\t') {
                char c = *at;
                at++;
                return c;
            }

            at++;
        }

        return '\0';
    }


    /// Determine if current line is a comment
    /// @param[in] line The current line
    /// @return true if line is a comment, false otherwise
    bool isCommentLine(std::string & line) {
        string_it at = line.begin();
        string_it end = line.end();
        char c = consumeSymbol(at, end);
        return c == '"' || c == '*';
    }


    /// Get an integer from the string
    /// \param[in] string
    /// \return an integer
    int fetchNumber(std::string &string) {
        std::stringstream stream(string);
        int num;
        stream >> num;
        return num;
    }


    /// Read a vector of the form {val1, val2, ..., valn}
    /// @param string Contains the vector
    /// @return a list with the n numbers
    listVector readVector(std::string &string) {
        std::stringstream stream(string);
        listVector vector;
        NT value;

        while (stream >> value) {
            vector.push_back(value);
        }

        return vector;
    }

public:

    /// Reads an SDPA format file
    /// \param[in] is An open stram pointing to the file
    /// \param[out] matrices the matrices A0, A1, A2, ..., An
    /// \param[out] objectiveFunction The objective function of the sdp
    void loadSDPAFormatFile(std::ifstream &is, std::vector<MT> &matrices, VT &objectiveFunction) {
        std::string line;
        std::string::size_type sz;

        std::getline(is, line, '\n');

        //skip comments
        while (isCommentLine(line)) {
            std::getline(is, line, '\n');
        }

        //read variables number
        int variablesNum = fetchNumber(line);

        if (std::getline(is, line, '\n').eof())
            throw std::runtime_error("Unexpected end of file");

        //read number of blocks
        int blocksNum = fetchNumber(line);

        if (std::getline(is, line, '\n').eof())
            throw std::runtime_error("Unexpected end of file");

        //read block structure vector
        listVector blockStructure = readVector(line);

        if (blockStructure.size() != blocksNum)
            throw std::runtime_error("Wrong number of blocks");

        if (std::getline(is, line, '\n').eof())
            throw std::runtime_error("Unexpected end of file");

        //read constant vector
        listVector constantVector = readVector(line);

        while (constantVector.size() < variablesNum) {
            if (std::getline(is, line, '\n').eof())
                throw std::runtime_error("Unexpected end of file");

            listVector t = readVector(line);
            constantVector.insert(std::end(constantVector), std::begin(t), std::end(t));
        }

        matrices = std::vector<MT>(variablesNum + 1);
        int matrixDim = 0;
        for (auto x : blockStructure)
            matrixDim += std::abs((int) x);

        //read constraint matrices
        for (int atMatrix = 0; atMatrix < matrices.size(); atMatrix++) {
            MT matrix;
            matrix.setZero(matrixDim, matrixDim);

            int offset = 0;

            for (auto blockSize : blockStructure) {

                if (blockSize > 0) { //read a block blockSize x blockSize
                    int at = 0;
                    int i = 0, j = 0;

                    while (at < blockSize * blockSize) {
                        if (std::getline(is, line, '\n').eof())
                            throw 1;

                        listVector vec = readVector(line);

                        for (double value : vec) {
                            matrix(offset + i, offset + j) = value;
                            at++;
                            if (at % (int) blockSize == 0) { // new row
                                i++;
                                j = 0;
                            } else { //new column
                                j++;
                            }
                        }
                    } /* while (at<blockSize*blockSize) */

                } else { //read diagonal block
                    blockSize = std::abs(blockSize);
                    int at = 0;

                    while (at < blockSize) {
                        if (std::getline(is, line, '\n').eof())
                            throw 1;

                        listVector vec = readVector(line);

                        for (double value : vec) {
                            matrix(offset + at, offset + at) = value;
                            at++;
                        }
                    } /* while (at<blockSize) */
                }

                offset += std::abs(blockSize);
            } /* for (auto blockSize : blockStructure) */

            //the LMI in SDPA format is >0, I want it <0
            if (atMatrix == 0) //F0 has - before it in SDPA format, the rest have +
                matrices[atMatrix] = matrix;
            else
                matrices[atMatrix] = -1 * matrix;
        }

        // return lmi and objective function
        objectiveFunction.setZero(variablesNum);
        int at = 0;

        for (auto value : constantVector)
            objectiveFunction(at++) = value;
    }

    /// Create a SDPA format file
    /// \param[in] os Open stream to file
    /// \param[in] matrices The matrices A0, ..., An
    /// \param[in] objectiveFunction The objective function of the sdp
    void writeSDPAFormatFile(std::ostream &os, std::vector<MT> const & matrices, VT const & objectiveFunction) {
        int dim = matrices.size() - 1;
        MT A0 = matrices[0];

        os << dim << "\n";
        os << 1 << "\n";
        os << A0.rows() << "\n";

        os << objectiveFunction.transpose() << "\n";

        for (int i = 0; i < A0.rows(); i++)
            os << A0.row(i) << "\n";

        for (int at=1 ; at<matrices.size() ; ++at)
            for (int i = 0; i < matrices[at].rows(); i++)
                os << -1 * matrices[at].row(i) << "\n";
    }


    /// Read a spectrahedron and a vector (objective function) from a SDPA format input file
    /// \tparam Point
    /// \param[in] is opened stream to input file
    /// \param[out] spectrahedron
    /// \param[out] objectiveFunction
    template <typename Point>
    void loadSDPAFormatFile(std::ifstream &is, Spectrahedron<NT, MT, VT> &spectrahedron, Point &objectiveFunction) {
        std::vector<MT> matrices;
        VT coeffs;
        loadSDPAFormatFile(is, matrices, coeffs);
        LMI<NT, MT, VT> lmi(matrices);
        spectrahedron = Spectrahedron<NT, MT, VT>(lmi);
        objectiveFunction = Point(coeffs);
    }


    /// Write a spectrahedron and a vector (objective function) to a SDPA format output file
    /// \tparam Point
    /// \param[in] is opened stream to output file
    /// \param[in] spectrahedron
    /// \param[in] objectiveFunction
    template <typename Point>
    void writeSDPAFormatFile(std::ostream &os, Spectrahedron<NT, MT, VT> const & spectrahedron, Point const & objectiveFunction) {
        writeSDPAFormatFile(os, spectrahedron.getLMI().getMatrices(), objectiveFunction.getCoefficients());
    }
};


#endif //VOLESTI_SDPA_FORMAT_MANAGER_H

