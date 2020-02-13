
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "mex.hpp"
#include "mexAdapter.hpp"


#ifndef _FNIRS1_MATLAB_ERROR_
#define _FNIRS1_MATLAB_ERROR_



void matlabError(const std::string &msg) {
  static matlab::mex::Function __func;
  static std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = 
    __func.getEngine();
  static matlab::data::ArrayFactory _factory_;
  static std::ostringstream _oss_;
  // Pass stream content to MATLAB fprintf function
  _oss_ << msg;
  matlabPtr->feval(u"error", 0,
		   std::vector<matlab::data::Array>({ 
		       _factory_.createScalar(_oss_.str()) }));
  // Clear stream buffer
  _oss_.str("");
  _oss_.clear();
};


#endif  // _FNIRS1_MATLAB_ERROR_

