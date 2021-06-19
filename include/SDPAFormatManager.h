// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLESTI_SDPA_FORMAT_MANAGER_NEW_H
#define VOLESTI_SDPA_FORMAT_MANAGER_NEW_H


#include "convex_bodies/spectrahedra_new/newSpectrahedron.h"

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
public:
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

//public:

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
            if (atMatrix == 0) { //F0 has - before it in SDPA format, the rest have +
                matrices[atMatrix] = matrix;
            }else {
                matrices[atMatrix] = -1 * matrix;
            }
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

    void loadSDPAFormatFile(std::ifstream &is, std::vector< Eigen::SparseMatrix<NT> > &matrices, VT &objectiveFunction) {
    }

    template <typename Point>
    void loadSDPAFormatFile(std::ifstream &is, Spectrahedron<NT, Eigen::SparseMatrix<NT>, VT> &spectrahedron, Point &objectiveFunction) {
        std::vector<Eigen::SparseMatrix<NT>> matrices;
        VT coeffs;
        loadSDPAFormatFile(is, matrices, coeffs);
        LMI<NT, Eigen::SparseMatrix<NT>, VT> lmi(matrices);
        spectrahedron = Spectrahedron<NT, Eigen::SparseMatrix<NT>, VT>(lmi);
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

/// Reads an SDPA format file
    /// \param[in] is An open stram pointing to the file
    /// \param[out] matrices the matrices A0, A1, A2, ..., An
    /// \param[out] objectiveFunction The objective function of the sdp
template <typename NT, typename SMT, typename DMT, typename VT>
void ole(std::ifstream &is, std::vector<SMT> &matrices_sparse, std::vector<DMT> &matrices_dense, VT &objectiveFunction) {
        
        typedef std::string::iterator string_it;
        typedef std::list<NT> listVector;
        SdpaFormatManager<NT> spm;

        typedef Eigen::Triplet<NT> T;
        std::vector<T> tripletList;
        
        std::string line;
        std::string::size_type sz;

        std::getline(is, line, '\n');
        double ref = 0.00000000000001;

        //skip comments
        while (spm.isCommentLine(line)) {
            std::getline(is, line, '\n');
        }

        //read variables number
        int variablesNum = spm.fetchNumber(line);

        if (std::getline(is, line, '\n').eof())
            throw std::runtime_error("Unexpected end of file");

        //read number of blocks
        int blocksNum = spm.fetchNumber(line);

        if (std::getline(is, line, '\n').eof())
            throw std::runtime_error("Unexpected end of file");

        //read block structure vector
        listVector blockStructure = spm.readVector(line);

        if (blockStructure.size() != blocksNum)
            throw std::runtime_error("Wrong number of blocks");

        if (std::getline(is, line, '\n').eof())
            throw std::runtime_error("Unexpected end of file");

        //read constant vector
        listVector constantVector = spm.readVector(line);

        while (constantVector.size() < variablesNum) {
            if (std::getline(is, line, '\n').eof())
                throw std::runtime_error("Unexpected end of file");

            listVector t = spm.readVector(line);
            constantVector.insert(std::end(constantVector), std::begin(t), std::end(t));
        }

        matrices_dense = std::vector<DMT>(variablesNum + 1);
        matrices_sparse = std::vector<SMT>(variablesNum + 1);
        int matrixDim = 0;
        for (auto x : blockStructure)
            matrixDim += std::abs((int) x);

        //read constraint matrices
        for (int atMatrix = 0; atMatrix < matrices_dense.size(); atMatrix++) {
            DMT matrix;
            SMT smatrix;
            matrix.setZero(matrixDim, matrixDim);

            int offset = 0;

            for (auto blockSize : blockStructure) {

                if (blockSize > 0) { //read a block blockSize x blockSize
                    int at = 0;
                    int i = 0, j = 0;

                    while (at < blockSize * blockSize) {
                        if (std::getline(is, line, '\n').eof())
                            throw 1;

                        listVector vec = spm.readVector(line);

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

                        listVector vec = spm.readVector(line);

                        for (double value : vec) {
                            matrix(offset + at, offset + at) = value;
                            at++;
                        }
                    } /* while (at<blockSize) */
                }

                offset += std::abs(blockSize);
            } /* for (auto blockSize : blockStructure) */

            //the LMI in SDPA format is >0, I want it <0
            tripletList.clear();
            if (atMatrix == 0) { //F0 has - before it in SDPA format, the rest have +
                matrices_dense[atMatrix] = matrix;
                for (int i=0; i<matrixDim; i++){
                    for (int j=0; j<matrixDim; j++){
                        if (i>=j) {
                            tripletList.push_back(T(i, j, matrix(i,j)));
                        }
                    }
                }
                smatrix.resize(matrixDim, matrixDim);
                smatrix.setFromTriplets(tripletList.begin(), tripletList.end());
                //smatrix = smatrix.pruned(ref);
                matrices_sparse[atMatrix] = smatrix.pruned(ref);
                //std::cout<<Eigen::MatrixXd(smatrix.pruned(ref))<<"\n"<<std::endl;
            }else {
                matrices_dense[atMatrix] = -1 * matrix;
                for (int i=0; i<matrixDim; i++){
                    for (int j=0; j<matrixDim; j++){
                        if (i>=j) {
                            tripletList.push_back(T(i, j, -matrix(i,j)));
                        }
                    }
                }
                smatrix.resize(matrixDim, matrixDim);
                smatrix.setFromTriplets(tripletList.begin(), tripletList.end());
                //smatrix = smatrix.pruned(ref);
                //smatrix = (-1.0*matrix).sparseView();
                //smatrix = smatrix.pruned(ref);
                matrices_sparse[atMatrix] = smatrix.pruned(ref);
                //std::cout<<Eigen::MatrixXd(smatrix.pruned(ref))<<"\n"<<std::endl;
            }
        }

        // return lmi and objective function
        objectiveFunction.setZero(variablesNum);
        int at = 0;

        for (auto value : constantVector)
            objectiveFunction(at++) = value;
}

template <typename DMT, typename NT, typename MT, typename VT, typename Point>
void loadSDPASparseFormatFile(std::ifstream &is, Spectrahedron<NT, MT, VT> &spectrahedron, Point &objectiveFunction) {
        std::vector<MT> matrices_sparse;
        std::vector<DMT> matrices_dense;
        VT coeffs;
        ole<NT>(is, matrices_sparse, matrices_dense, coeffs);
        LMI<NT, MT, VT> lmi(matrices_sparse);
        spectrahedron = Spectrahedron<NT, MT, VT>(lmi);
        objectiveFunction = Point(coeffs);
}

std::vector<double> readVector3(std::string &string) {
    std::stringstream stream(string);
    std::vector<double> vector;
    double value;

    while (stream >> value) {
        vector.push_back(value);
    }

    return vector;
}


/// Reads an SDPA format file
/// \param[in] is An open stram pointing to the file
/// \param[out] matrices the matrices A0, A1, A2, ..., An
/// \param[out] objectiveFunction The objective function of the sdp
template <typename MT, typename NT, typename LMII, typename VT>
void loadSparseSDPAFormatFile(std::istream &is, LMII &lmi, VT &objectiveFunction) {
    std::string line;
    std::string::size_type sz;
    SdpaFormatManager<NT> spm;
    typedef std::list<NT> listVector;

    typedef Eigen::Triplet<NT> T;
    std::vector<T> tripletList;

    std::getline(is, line, '\n');

    //skip comments
    while (spm.isCommentLine(line)) {
        std::getline(is, line, '\n');
    }

    //read variables number
    int variablesNum = spm.fetchNumber(line);

    if (std::getline(is, line, '\n').eof())
        throw std::runtime_error("Unexpected end of file");

    //read number of blocks
    int blocksNum = spm.fetchNumber(line);

    if (std::getline(is, line, '\n').eof())
        throw std::runtime_error("Unexpected end of file");

    //read block structure vector
    std::vector<NT> blockSizes = readVector3(line);

    if (blockSizes.size() != blocksNum)
        throw std::runtime_error("Wrong number of blocks");

    if (std::getline(is, line, '\n').eof())
        throw std::runtime_error("Unexpected end of file");

    //read objective function
    listVector constantVector = spm.readVector(line);

    while (constantVector.size() < variablesNum) {
        if (std::getline(is, line, '\n').eof())
            throw std::runtime_error("Unexpected end of file");

        std::vector<NT> t = readVector3(line);
        constantVector.insert(std::end(constantVector), std::begin(t), std::end(t));
    }

    std::vector<MT> matrices = std::vector<MT>(variablesNum + 1);
    int matrixDim = 0;
    for (auto x : blockSizes)
        matrixDim += std::abs((int) x);

    for (int i=0 ; i<matrices.size() ; ++i)
        matrices[i].resize(matrixDim, matrixDim);
    // read constraint matrices
    // entries are of the form
    // <matno> <blkno> <i> <j> <entry>
    tripletList.clear();
    MT smatrix(matrixDim, matrixDim);
    int current_mat = 0;
    while (!std::getline(is, line, '\n').eof()) {
        std::vector<NT> t = readVector3(line);
        std::cout<<"t = "<<std::endl;
        for (int k = 0; k<5; k++){
            std::cout<<t[k]<<" ";
        }
        std::cout<<"\n";
        std::cout<<"current_mat = "<<current_mat<<std::endl;
        if (t[0] > current_mat) {
            smatrix.setFromTriplets(tripletList.begin(), tripletList.end());
            //smatrix = smatrix.pruned(ref);
            matrices[t[0] - 1] = smatrix;
            std::cout<<"number of matrix = "<<t[0]<<std::endl;
            //std::cout<<Eigen::MatrixXd(smatrix)<<"\n"<<std::endl;
            tripletList.clear();
            smatrix.resize(matrixDim, matrixDim);
            current_mat++;
        }
//            std::cout << line << "\n";
        int blockOffset = 0;
        for (int i=1; i<t[1] ; ++i)
            blockOffset += std::abs(blockSizes[i-1]);

        int i = t[2] + blockOffset-1;
        int j = t[3] + blockOffset-1;

        if (i <= j) {
            if (t[0] > 0){
                tripletList.push_back(T(j, i, -t[4]));
            } else {
                tripletList.push_back(T(j, i, t[4]));
            }
        }

        //matrices[t[0]](i,j) = t[4];
//            std::cout << i << " " << j << "\n";
        // matrix is symmetric
        // only upper triangular is provided
        // fill lower triangular
        //if (i!=j)
        //    matrices[t[0]](j,i) = t[4];
    }
    std::cout<<"current_mat = "<<current_mat<<std::endl;
    smatrix.setFromTriplets(tripletList.begin(), tripletList.end());
    std::cout<<"variablesNum = "<<variablesNum<<std::endl;
            //smatrix = smatrix.pruned(ref);
    matrices[variablesNum] = smatrix;
    
    //std::cout<<Eigen::MatrixXd(smatrix)<<"\n"<<std::endl;


    //for (int atMatrix=1 ; atMatrix<matrices.size() ; atMatrix++) {
        //the LMI in SDPA format is >0, I want it <0
        //F0 has - before it in SDPA format, the rest have +
      //  matrices[atMatrix] *= -1;
    //}
    // return lmi and objective function
    objectiveFunction.setZero(variablesNum);
    int at = 0;
    std::cout<<"constantVector.size() = "<<constantVector.size()<<std::endl;

    for (auto value : constantVector)
        objectiveFunction(at++) = value;
    lmi = LMII(matrices);
}




/// Reads an SDPA format file
/// \param[in] is An open stram pointing to the file
/// \param[out] matrices the matrices A0, A1, A2, ..., An
/// \param[out] objectiveFunction The objective function of the sdp
template <typename MT, typename SMT, typename NT, typename LMII, typename VT>
void loadBlockSparseSDPAFormatFile(std::istream &is, LMII &lmi, VT &objectiveFunction) {
    std::string line;
    std::string::size_type sz;
    SdpaFormatManager<NT> spm;
    typedef std::list<NT> listVector;

    typedef Eigen::Triplet<NT> T;
    std::vector<T> tripletList;

    std::getline(is, line, '\n');

    //skip comments
    while (spm.isCommentLine(line)) {
        std::getline(is, line, '\n');
    }

    //read variables number
    int variablesNum = spm.fetchNumber(line);

    if (std::getline(is, line, '\n').eof())
        throw std::runtime_error("Unexpected end of file");

    //read number of blocks
    int blocksNum = spm.fetchNumber(line);

    if (std::getline(is, line, '\n').eof())
        throw std::runtime_error("Unexpected end of file");

    //read block structure vector
    std::vector<NT> blockSizes = readVector3(line);

    if (blockSizes.size() != blocksNum)
        throw std::runtime_error("Wrong number of blocks");

    if (std::getline(is, line, '\n').eof())
        throw std::runtime_error("Unexpected end of file");

    //read objective function
    listVector constantVector = spm.readVector(line);

    while (constantVector.size() < variablesNum) {
        if (std::getline(is, line, '\n').eof())
            throw std::runtime_error("Unexpected end of file");

        std::vector<NT> t = readVector3(line);
        constantVector.insert(std::end(constantVector), std::begin(t), std::end(t));
    }

    std::vector<MT> matrices;// = std::vector<MT>(variablesNum + 1);
    std::vector<SMT> sparse_matrices;// = std::vector<SMT>(blocksNum);
    int matrixDim = 0;
    for (auto x : blockSizes)
        matrixDim += std::abs((int) x);

    //for (int i=0 ; i<matrices.size() ; ++i)
    //    matrices[i].resize(matrixDim, matrixDim);
    // read constraint matrices
    // entries are of the form
    // <matno> <blkno> <i> <j> <entry>
    tripletList.clear();
    
    int current_mat = 0, current_block = 1, num_of_blocks = blockSizes.size();
    matrixDim = blockSizes[0];
    SMT smatrix(matrixDim, matrixDim);
    MT A;
    std::cout<<"matrixDim = "<<matrixDim<<std::endl;
    while (!std::getline(is, line, '\n').eof()) {
        std::vector<NT> t = readVector3(line);
        //std::cout<<"t = "<<std::endl;
        //for (int k = 0; k<5; k++){
        //    std::cout<<t[k]<<" ";
        //}
        //std::cout<<"\n";
        std::cout<<"current_mat = "<<current_mat<<std::endl;
        std::cout<<"current_block = "<<current_block<<std::endl;
        if (t[1] != current_block || t[0] != current_mat) {
            smatrix.setFromTriplets(tripletList.begin(), tripletList.end());
            sparse_matrices.push_back(smatrix);
            //smatrix = smatrix.pruned(ref);
            
            std::cout<<"matrix to add = "<<std::endl;
            std::cout<<Eigen::MatrixXd(smatrix)<<"\n"<<std::endl;
            tripletList.clear();
            
            if (t[0] > current_mat) {
                std::cout<<"block matrix is ready"<<std::endl;
                if (current_block < num_of_blocks) {
                    std::cout<<"add the rest zero blocks"<<std::endl;
                    for (int ii = current_block; ii<num_of_blocks; ii++){
                        matrixDim = blockSizes[ii];
                        smatrix = SMT(matrixDim, matrixDim);
                        std::cout<<"matrix to add = "<<std::endl;
                        std::cout<<Eigen::MatrixXd(smatrix)<<"\n"<<std::endl;
                        sparse_matrices.push_back(smatrix);
                    }
                }
                std::cout<<"initialize block matrix"<<std::endl;
                A = MT(sparse_matrices);
                matrices.push_back(A);
                sparse_matrices.clear();
                current_mat++;
                current_block = t[1];
                if (current_block > 1) {
                    std::cout<<"add first zero blocks"<<std::endl;
                    for(int jj = 1; jj<current_block; jj++){
                        std::cout<<"block = "<<jj<<std::endl;
                        matrixDim = blockSizes[jj - 1];
                        smatrix = SMT(matrixDim, matrixDim);
                        std::cout<<"matrix to add = "<<std::endl;
                        std::cout<<Eigen::MatrixXd(smatrix)<<"\n"<<std::endl;
                        sparse_matrices.push_back(smatrix);
                    }
                }
            } else {
                if (t[1] - current_block > 1) {
                    std::cout<<"add intermediate zero blocks"<<std::endl;
                    for(int jj = current_block+1; jj<t[1]; jj++){
                        std::cout<<"block = "<<jj<<std::endl;
                        matrixDim = blockSizes[jj - 1];
                        smatrix = SMT(matrixDim, matrixDim);
                        std::cout<<"matrix to add = "<<std::endl;
                        std::cout<<Eigen::MatrixXd(smatrix)<<"\n"<<std::endl;
                        sparse_matrices.push_back(smatrix);
                    }
                }
                current_block = t[1];
            }
            matrixDim = blockSizes[current_block - 1];
            smatrix = SMT(matrixDim, matrixDim);
            std::cout<<"matrixDim = "<<matrixDim<<std::endl;
        }
//            std::cout << line << "\n";
        //int blockOffset = 0;
        //for (int i=1; i<t[1] ; ++i)
        //    blockOffset += std::abs(blockSizes[i-1]);

        int i = t[2] - 1;// + blockOffset-1;
        int j = t[3] - 1;// + blockOffset-1;

        if (i <= j) {
            if (t[0] > 0){
                tripletList.push_back(T(int(t[3]) - 1, int(t[2]) - 1, -t[4]));
            } else {
                tripletList.push_back(T(int(t[3]) - 1, int(t[2]) - 1, t[4]));
            }
        }

        //matrices[t[0]](i,j) = t[4];
//            std::cout << i << " " << j << "\n";
        // matrix is symmetric
        // only upper triangular is provided
        // fill lower triangular
        //if (i!=j)
        //    matrices[t[0]](j,i) = t[4];
    }
    std::cout<<"current_mat = "<<current_mat<<std::endl;
    smatrix.setFromTriplets(tripletList.begin(), tripletList.end());
    sparse_matrices.push_back(smatrix);
    
    if (current_block < num_of_blocks) {
        std::cout<<"add the rest zero blocks"<<std::endl;
        for (int ii = current_block; ii<num_of_blocks; ii++){
            matrixDim = blockSizes[ii];
            smatrix = SMT(matrixDim, matrixDim);
            std::cout<<"matrix to add = "<<std::endl;
            std::cout<<Eigen::MatrixXd(smatrix)<<"\n"<<std::endl;
            sparse_matrices.push_back(smatrix);
        }
    }
    A = MT(sparse_matrices);
    matrices.push_back(A);
    std::cout<<"variablesNum = "<<variablesNum<<std::endl;
            //smatrix = smatrix.pruned(ref);
    //matrices[variablesNum] = smatrix;
    
    //std::cout<<Eigen::MatrixXd(smatrix)<<"\n"<<std::endl;


    //for (int atMatrix=1 ; atMatrix<matrices.size() ; atMatrix++) {
        //the LMI in SDPA format is >0, I want it <0
        //F0 has - before it in SDPA format, the rest have +
      //  matrices[atMatrix] *= -1;
    //}
    // return lmi and objective function
    objectiveFunction.setZero(variablesNum);
    int at = 0;
    std::cout<<"constantVector.size() = "<<constantVector.size()<<std::endl;

    for (auto value : constantVector)
        objectiveFunction(at++) = value;
    lmi = LMII(matrices);
}



#endif //VOLESTI_SDPA_FORMAT_MANAGER_H

