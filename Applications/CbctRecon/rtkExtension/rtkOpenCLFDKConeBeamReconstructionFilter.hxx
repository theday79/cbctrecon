/*=========================================================================
 *
 *  Copyright RTK Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#ifndef RTKOPENCLFDKCONEBEAMRECONSTRUCTIONFILTER_HXX
#define RTKOPENCLFDKCONEBEAMRECONSTRUCTIONFILTER_HXX

#include "rtkOpenCLFDKConeBeamReconstructionFilter.h"

namespace rtk {

OpenCLFDKConeBeamReconstructionFilter ::
    OpenCLFDKConeBeamReconstructionFilter() {
  // Create each filter which are specific for OpenCL
  m_BackProjectionFilter = BackProjectionFilterType::New();
#if USE_CLFFT
  m_WeightFilter = WeightFilterType::New(); // Not yet implemented in OpenCL due
                                            // to heavy texture use in cuda
                                            // version
  m_RampFilter = RampFilterType::New();
  m_WeightFilter->SetInput(m_ExtractFilter->GetOutput());
  m_RampFilter->SetInput(m_WeightFilter->GetOutput());
#endif
  std::cout << "before getoutput" << std::endl;
  // Permanent internal connections
  m_BackProjectionFilter->SetInput(
      1, m_RampFilter
             ->GetOutput()); // m_rampFilter is not performed yet at this point
  std::cout << "after getoutput" << std::endl;

  // Default parameters
  m_BackProjectionFilter->InPlaceOn();
  m_BackProjectionFilter->SetTranspose(false);
}

void OpenCLFDKConeBeamReconstructionFilter ::GenerateData() {
  auto *openclbp = dynamic_cast<BackProjectionFilterType *>(
      m_BackProjectionFilter.GetPointer());

  // Init GPU memory
  openclbp->InitDevice();
#if USE_CLFFT
  // Because the Cuda inplace filter takes care of inheritance in
  // CUDreconstruction, we have to do it manually here
  this->Superclass::m_RampFilter = GetRampFilter();
#endif
  // Run reconstruction
  this->Superclass::GenerateData();

  // Transfer result to CPU image
  openclbp->CleanUpDevice();
}

} // end namespace rtk

#endif // RTKOPENCLFDKCONEBEAMRECONSTRUCTIONFILTER_HXX
