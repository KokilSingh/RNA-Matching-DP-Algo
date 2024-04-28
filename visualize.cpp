/** @mainpage Design and Analysis of Algorithms - Assignment 2
* @authors Kokil Singh  <br>
Sarvesh Sutaone <br>
Sanjana Adapala <br>
Yashraj Mehta
* @section intro About the Assignment
* We implement a Dynamic Programming Algorithm to predict the secondary structure of a RNA.
* 
* @section report About the Report
* The report includes the following :
* - Background
* - Problem Formulation
* - DP Solution
* - Algorithm Discussion
* - Timing Analysis
* - Issues in Coding
* - Visualization
* - References
*
* @section documentation Code Documentation
* The code is documented using Doxygen. The documentation can be viewed in the Files section.
*/


/**
 * @file visualize.cpp
 * @brief This file contains RNA Folding algorithm visualization
 * 
 * This program visualizes the RNA folding problem by determining the maximum number
 * of possible pairings and writing the pairings to an output file.
 */

// Include relevant header files
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <cstdlib>
#include <queue>
#include <utility>
using namespace std;

//! Function to check if pairing between bases is possible or not.
// A can pair with U
// C can pair with G

bool complementary(char a, char b) {
/**
 * @brief Check if pairing between bases is possible.
 * 
 * A can pair with U, C can pair with G.
 * 
 * @param a First base
 * @param b Second base
 * @return True if pairing is possible, False otherwise
 */
    if ((a == 'A' && b == 'U') || (b == 'A' && a == 'U') || (a == 'C' && b == 'G') || (b == 'C' && a == 'G')) {
        return true;
    } else {
        return false;
    }
}

void countPairs(vector<vector<int>> &matching, vector<vector<int>> &mapped, string &B) {
/**
 * @brief Function to count the maximum number of pairings possible with the given rules
 * 
 * This function writes the input RNA sequence and the identified pairs to an output file named "pairings.txt".
 * 
 * @param matching Matrix storing the maximum number of pairings
 * @param mapped Matrix storing the mapping information for pairings
 * @param B Input RNA sequence
 */
    //! Starting from 5 to incorporate Rule 1 (No sharp turns)
    for (int k = 5; k < B.size(); k++) {
        for (int i = 0; i < B.size() - k; i++) {
            int j = i + k;
            int part_a = 0;
            int part_b = 0;
            for (int t = i; t < j - 4; t++) {
                if (t > 0 && complementary(B[t], B[j]) && part_a < 1 + matching[i][t - 1] + matching[t + 1][j - 1]) {
                    part_a = 1 + matching[i][t - 1] + matching[t + 1][j - 1];
                    part_b = t + 1;   
                } else if (complementary(B[t], B[j]) && t == 0 && part_a < 1 + matching[t + 1][j - 1]) {
                    part_a = 1 + matching[t + 1][j - 1];
                    part_b = t + 1;
                }
            }
            matching[i][j] = max(matching[i][j - 1], part_a);
            if (matching[i][j - 1] >= part_a) {
                mapped[i][j] = 0;
            } else {
                mapped[i][j] = part_b;
            }
        }
    }
}

 
void printPairs(vector<vector<int>> &matching, vector<vector<int>> &mapped, string B) {
/**
 * @brief Function to print the input sequence and the pairs to a file "pairings.txt"
 * 
 * This function handles file input, initializes data structures, calls functions
 * to count pairings and write them to a file, and executes a Python script for visualization.
 * 
 * @param argc Number of command-line arguments
 * @param argv Command-line arguments
 * @return 0 on success, 1 on failure
 */
    ofstream outputFile("pairings.txt");
    if (!outputFile.is_open()) {
        cout << "Error opening output file!" << endl;
        return;
    }
    //! Writing Input Sequence to the file
    outputFile << B << endl;
    int no_of_matches = matching[0][B.size() - 1];
    int row = 0;
    int col = B.size() - 1;
    bool isMatched = false;
    int left = 0;
    int right = B.size()-1;
    //! While loop to identify the base pairs
    queue<pair<int,int>> subproblems;
    subproblems.push({left,right});
    //cout<<"Pushed: 0 "<<left<<" "<<right<<"\n";
    while(no_of_matches>0 && !subproblems.empty())
    {
        
        left=subproblems.front().first;
        right=subproblems.front().second;
        //cout<<"In: "<<left<<" "<<left<<" "<<right<<"\n";
        subproblems.pop();
        if(mapped[left][right]>=1)
        {
            outputFile << mapped[left][right] << " " << right+1 << "\n";
            //cout << "Pairing Found: "<<mapped[left][right] << " " << right+1 << "\n";
            no_of_matches--;
            if(left+1<B.size() && right-1>mapped[left][right])
            {
                //! Internal subproblem
                subproblems.push({mapped[left][right],right-1});
                //cout<<"\tPushed 1: "<<left+1<<" "<<mapped[left][right]<<" "<<right-1<<"\n";
            }
            if(mapped[left][right]-2>left)
            {
                //! Left subproblem
                subproblems.push({left,mapped[left][right]-2});
                //cout<<"\tPushed 2: "<<left<<" "<<left<<" "<<mapped[left][right]-2<<"\n";
            }
            
            isMatched=true;
        }
        else if(right-1>left)
        {
            // If zero
            subproblems.push({left,right-1});
            //cout<<"\tPushed 3: "<<left<<" "<<left<<" "<<right-1<<"\n";
        }
    }
    
    if (isMatched == false) {
        outputFile << 0 << " " << 0;
        outputFile << endl;
    }
    outputFile.close();
}

int main(int argc, char *argv[]) {
    //! Check if the filename is provided as command line argument
    if (argc != 2) {
        cout << "Usage: " << argv[0] << " <filename>" << endl;
        return 1;
    }
    //! Open the input file
    string filename = argv[1];
    ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        cout << "Error opening file: " << filename << endl;
        return 1;
    }
    string B;
    inputFile >> B;
    inputFile.close();
    //! Initialize tables
    vector<vector<int>> matching(B.size(), vector<int>(B.size(), 0));
    vector<vector<int>> mapped(B.size(), vector<int>(B.size(), 0));
    //! Count the number of pairs
    countPairs(matching, mapped, B);
    cout << "\nMaximum Number of Pairs Possible: " << matching[0][B.size() - 1] << "\n"; 
    //! Write the number of pairs to a file
    printPairs(matching, mapped, B);
    //! Execute the Python file named "1.py"
    system("python 1.py");
    return 0;
}
