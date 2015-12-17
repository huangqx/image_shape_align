/*---
function [edgeVIds, faceEdgeIds] = facetop(faceVIds)
Input: 

Output

--*/

#include "mex.h"
#include "lattice2grid.h"
#include <vector>
#include <algorithm>
using namespace std;


void mexFunction(
     int nargout,
     mxArray *output[],
     int nargin,
     const mxArray *input[]) {
  /* check argument */
  if (nargin != 2) {
    mexErrMsgTxt("One input arguments required.");
  }
  if (nargout != 1) {
    mexErrMsgTxt("Incorrect number of output arguments."); 
  }
  
  double *data = (double*)mxGetData(input[0]);
  int dim = static_cast<int> (mxGetM(input[0]));
  int numOfVertices = static_cast<int> (mxGetN(input[0]));
  
  vector<ColorPT> lattice;
  lattice.resize(numOfVertices);

  for (int vId = 0; vId < numOfVertices; ++vId) {
    lattice[vId].pos.x = data[5*vId];
    lattice[vId].pos.y = data[5*vId+1];
    lattice[vId].color.color[0] = data[5*vId+2];
    lattice[vId].color.color[1] = data[5*vId+3];
    lattice[vId].color.color[2] = data[5*vId+4];
  }

  data = (double*)mxGetData(input[1]);
  ParaFI para;
  para.dimLatticeX = int(data[0]);
  para.dimLatticeY = int(data[1]);
  para.lowerCorner.x = data[2];
  para.lowerCorner.y = data[3];
  para.nCols = int(data[4]);
  para.nRows = int(data[5]);
  para.pixelGap = data[6];

  vector<Color3f> image;
  image.resize(para.nCols*para.nRows);

  FFD2Image(lattice, para, &image);

  output[0] = mxCreateDoubleMatrix(3, para.nCols*para.nRows, mxREAL);
  data = mxGetPr(output[0]);
  for (unsigned id = 0; id < para.nCols*para.nRows; ++id) {
    data[3*id] = image[id].color[0];
    data[3*id+1] = image[id].color[1];
    data[3*id+2] = image[id].color[2];
  }
}


