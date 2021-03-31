#include <stdio.h>
#include "cs1300bmp.h"
#include <iostream>
#include <fstream>
#include "Filter.h"
#include <omp.h>

using namespace std;

#include "rdtsc.h"

//
// Forward declare the functions
//
Filter * readFilter(string filename);
double applyFilter(Filter *filter, cs1300bmp *input, cs1300bmp *output);

int
main(int argc, char **argv)
{

  if ( argc < 2) {
    fprintf(stderr,"Usage: %s filter inputfile1 inputfile2 .... \n", argv[0]);
  }

  //
  // Convert to C++ strings to simplify manipulation
  //
  string filtername = argv[1];

  //
  // remove any ".filter" in the filtername
  //
  string filterOutputName = filtername;
  string::size_type loc = filterOutputName.find(".filter");
  if (loc != string::npos) {
    //
    // Remove the ".filter" name, which should occur on all the provided filters
    //
    filterOutputName = filtername.substr(0, loc);
  }

  Filter *filter = readFilter(filtername);

  double sum = 0.0;
  int samples = 0;

  for (int inNum = 2; inNum < argc; inNum++) {
    string inputFilename = argv[inNum];
    string outputFilename = "filtered-" + filterOutputName + "-" + inputFilename;
    struct cs1300bmp *input = new struct cs1300bmp;
    struct cs1300bmp *output = new struct cs1300bmp;
    int ok = cs1300bmp_readfile( (char *) inputFilename.c_str(), input);

    if ( ok ) {
      double sample = applyFilter(filter, input, output);
      sum += sample;
      samples++;
      cs1300bmp_writefile((char *) outputFilename.c_str(), output);
    }
    delete input;
    delete output;
  }
  fprintf(stdout, "Average cycles per sample is %f\n", sum / samples);

}

struct Filter *
readFilter(string filename)
{
  ifstream input(filename.c_str());

  if ( ! input.bad() ) {
    int size = 0;
    input >> size;
    Filter *filter = new Filter(size);
    int div;
    input >> div;
    filter -> setDivisor(div);
    for (int i=0; i < size; i++) {
      for (int j=0; j < size; j++) {
        int value;
        input >> value;
        filter -> set(i,j,value);
      }
    }
    return filter;
  } else {
    cerr << "Bad input in readFilter:" << filename << endl;
    exit(-1);
  }
}


double
applyFilter(struct Filter *filter, cs1300bmp *input, cs1300bmp *output)
{
    long long cycStart, cycStop;
    cycStart = rdtscll();

    int columns = (input->width) -  1;
    int rows = (input->height) - 1;
    output -> width = columns + 1;
    output -> height = rows + 1;
    
    int size = filter -> getSize();
    int red1, red2, red3, green4, green5, green6, blue7, blue8, blue9 = 0;//this refers to the colors of the three planes intially at 0
    float divisor = filter -> getDivisor();
    int *FilterStart; //this array will be in charge or applying the filter to the array of pixels
    
    int i = 0;
    FilterStart = &i;//this will get the acc value of the plane and multiply filter to it
    
    #pragma omp parallel for
    for(int row = 1; row < rows; row++) {
        for(int col = 1; col < columns; col++) {
            red1 = 0;
            red2 = 0;
            red3 = 0;
            green4 = 0;
            green5 = 0;
            green6 = 0;
            blue7 = 0;
            blue8 = 0;
            blue9 = 0;
            
            red1 += (input->color[0][row+i-1][col-1]*FilterStart[i*size]);//the first plane will be intialized row col at 0
            red2 += (input->color[0][row+i][col-1]*FilterStart[(i+1)*size]);
            red3 += (input->color[0][row+i+1][col-1]*FilterStart[(i+2)*size]);
            green4 += (input->color[1][row+i-1][col-1]*FilterStart[i*size]);//this will start adding the values of the filter and plane at [1][0],[1][1][1][2]
            green5 += (input->color[1][row+i][col-1]*FilterStart[(i+1)*size]);
            green6 += (input->color[1][row+i+1][col-1]*FilterStart[(i+2)*size]);
            blue7 += (input->color[2][row+i-1][col-1]*FilterStart[i*size]);
            blue8 += (input->color[2][row+i][col-1]*FilterStart[(i+1)*size]);
            blue9 += (input->color[2][row+i+1][col-1]*FilterStart[(i+2)*size]);
            //this will be in charge of the values at [2][0] [2][1] [2][2]      
            red1 += (input->color[0][row+i-1][col]*FilterStart[i*size+1]);
            red2 += (input->color[0][row+i][col]*FilterStart[(i+1)*size+1]);//eventually all of these addtions and multiplciations will be added to eachother for the total val
            red3 += (input->color[0][row+i+1][col]*FilterStart[(i+2)*size+1]);
            green4 += (input->color[1][row+i-1][col]*FilterStart[i*size+1]);
            green5 += (input->color[1][row+i][col]*FilterStart[(i+1)*size+1]);
            green6 += (input->color[1][row+i+1][col]*FilterStart[(i+2)*size+1]);
            blue7 += (input->color[2][row+i-1][col]*FilterStart[i*size+1]);
            blue8 += (input->color[2][row+i][col]*FilterStart[(i+1)*size+1]);
            blue9 += (input->color[2][row+i+1][col]*FilterStart[(i+2)*size+1]);
                    
            red1 += (input->color[0][row+i-1][col+1]*FilterStart[i*size+2]);
            red2 += (input->color[0][row+i][col+1]*FilterStart[(i+1)*size+2]);
            red3 += (input->color[0][row+i+1][col+1]*FilterStart[(i+2)*size+2]);
            green4 += (input->color[1][row+i-1][col+1]*FilterStart[i*size+2]);
            green5 += (input->color[1][row+i][col+1]*FilterStart[(i+1)*size+2]);
            green6 += (input->color[1][row+i+1][col+1]*FilterStart[(i+2)*size+2]);
            blue7 += (input->color[2][row+i-1][col+1]*FilterStart[i*size+2]);
            blue8 += (input->color[2][row+i][col+1]*FilterStart[(i+1)*size+2]);
            blue9 += (input->color[2][row+i+1][col+1]*FilterStart[(i+2)*size+2]);             
            
            red1 = (red1+red2+red3) / divisor; // the divisor for each plane is 3 in total 9
            green4 = (green4+green5+green6) / divisor;
            blue7 = (blue7+blue8+blue9) / divisor;
            
            (red1 < 0) ? red1=0 : red1=red1; //ternaray will check if red1 val < 0 will be initalized at 0
            (red1 > 255) ? red1=255 : red1=red1; //same with val> 255 make 255
            
            (green4 < 0) ? green4=0 : green4=green4;
            (green4 > 255) ? green4=255 : green4 = green4;
            
            (blue7 < 0) ? blue7=0 : blue7=blue7;
            (blue7 > 255) ? blue7=255 : blue7=blue7;
            
            output -> color[0][row][col] = red1; //filter is applied to all 9 pixels!
            output -> color[1][row][col] = green4;
            output -> color[2][row][col] = blue7;
        }
    }
  cycStop = rdtscll();
  double diff = cycStop - cycStart;
  double diffPerPixel = diff / (output -> width * output -> height);
  fprintf(stderr, "Took %f cycles to process, or %f cycles per pixel\n",
	  diff, diff / (output -> width * output -> height));
  return diffPerPixel;
}