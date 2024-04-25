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

// Function to check if pairing between bases is possible or not.
// A can pair with U
// C can pair with G
bool complementary(char a, char b) {
    if ((a == 'A' && b == 'U') || (b == 'A' && a == 'U') || (a == 'C' && b == 'G') || (b == 'C' && a == 'G')) {
        return true;
    } else {
        return false;
    }
}

// Function to count the maximum number of pairings possible with the given rules 
void countPairs(vector<vector<int>> &matching, vector<vector<int>> &mapped, string &B) {
    // Starting from 5 to incorporate Rule 1 (No sharp turns)
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


// Function to print the input sequence and the pairs to a file "pairings.txt" 
void printPairs(vector<vector<int>> &matching, vector<vector<int>> &mapped, string B) {
    ofstream outputFile("pairings.txt");
    if (!outputFile.is_open()) {
        cout << "Error opening output file!" << endl;
        return;
    }
    // Writing Input Sequence to the file
    outputFile << B << endl;
    int no_of_matches = matching[0][B.size() - 1];
    int row = 0;
    int col = B.size() - 1;
    bool isMatched = false;
    int left = 0;
    int right = B.size()-1;
    // While loop to identify the base pairs
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
                // Internal subproblem
                subproblems.push({mapped[left][right],right-1});
                //cout<<"\tPushed 1: "<<left+1<<" "<<mapped[left][right]<<" "<<right-1<<"\n";
            }
            if(mapped[left][right]-2>left)
            {
                // Left subproblem
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
    // Check if the filename is provided as command line argument
    if (argc != 2) {
        cout << "Usage: " << argv[0] << " <filename>" << endl;
        return 1;
    }
    // Open the input file
    string filename = argv[1];
    ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        cout << "Error opening file: " << filename << endl;
        return 1;
    }
    string B;
    inputFile >> B;
    inputFile.close();
    // Initialize tables
    vector<vector<int>> matching(B.size(), vector<int>(B.size(), 0));
    vector<vector<int>> mapped(B.size(), vector<int>(B.size(), 0));
    // Count the number of pairs
    countPairs(matching, mapped, B);
    cout << "\nMaximum Number of Pairs Possible: " << matching[0][B.size() - 1] << "\n"; 
    // Write the number of pairs to a file
    printPairs(matching, mapped, B);
    // Execute the Python file named "1.py"
    system("python 1.py");
    return 0;
}
