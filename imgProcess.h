#pragma once

#include <opencv2/opencv.hpp> 
#include <string>
#include <vector>
#include <iostream>

using namespace cv;
using namespace std;

struct Shape {
    int x, y;
};

Shape read(string filename, uchar*& buffer){ //Note "&", buffer can be changed
    Mat image = imread(filename, IMREAD_COLOR);
    int num_lines = image.rows;
    int num_cols = image.cols;
    vector<Mat> channels;
    split(image, channels);
    Mat green_channel = channels[1]; // 0 = blue, 1 = green, 2 = red
    
    //copy data
    buffer = new uchar[num_lines * num_cols];
    uchar* ptr = buffer;
    uchar* cvptr = green_channel.ptr<uchar>();
    for (int i = 0; i < num_lines; i++) {
        for (int j = 0; j < num_cols; j++) {
            *ptr++ = *cvptr++;
        }
    }
    return Shape{ num_lines, num_cols };
}