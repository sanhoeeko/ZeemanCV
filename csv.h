#pragma once

#include <fstream>
#include <string>
#include <vector>
using namespace std;

template<typename ty>
void array_to_csv(string filename, vector<ty> data){
    int n = data.size();
    ofstream file(filename); //create output file stream with given filename
    for(int i=0; i<n; i++){ //iterate through data array
        file << data[i] << endl; //write each element to file followed by newline character
    }
    file.close(); //close file stream after writing is complete
}