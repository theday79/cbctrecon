// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com

/* Do reconstruction algorithms */

// std
#include <array>
#include <cmath> // for asin, cos
#include <cstdlib>
#include <iostream> // for operator<<, endl, basic_ostream, cout
#include <string>   // for operator<<

// Local
#include "cbctrecon_types.h"

// PLM
#ifdef USE_OPENCL_PLM
#undef TIMEOUT
#undef CUDA_FOUND
#include <fdk_opencl.h>
#include <opencl_util.h>
#include <plm_image.h>
#include <proj_image.h>
#include <proj_image_filter.h>
#include <proj_matrix.h>

FloatImageType::Pointer PlastimatchOpenCLFDK(
    const FloatImageType::Pointer &spCurImg,
    const rtk::ThreeDCircularProjectionGeometry::Pointer &m_spCustomGeometry,
    std::array<float, 3> spacing, std::array<plm_long, 3> sizeOutput) {
  spCurImg->Update();

  //////////////////////// start PLM fdk /////////////////////////////
  // Copy of opencl_reconstruct_conebeam from fdk_opencl.cxx in plastimatch
  // (because we already have loaded projections into memory)
  Opencl_device ocl_dev{};
  // Proj_image *proj;
  // int image_num;
  // float scale;
  std::cout << "Loading plmopencl..." << std::endl;
  LOAD_LIBRARY_SAFE(libplmopencl);
  LOAD_SYMBOL(opencl_open_device, libplmopencl);
  LOAD_SYMBOL(opencl_load_programs, libplmopencl);
  LOAD_SYMBOL(opencl_kernel_create, libplmopencl);
  LOAD_SYMBOL(opencl_buf_create, libplmopencl);
  LOAD_SYMBOL(opencl_buf_write, libplmopencl);
  LOAD_SYMBOL(opencl_set_kernel_args, libplmopencl);
  LOAD_SYMBOL(opencl_kernel_enqueue, libplmopencl);
  LOAD_SYMBOL(opencl_buf_read, libplmopencl);

  std::cout << "Start opencl environment and enqueue kernels..." << std::endl;
  /* Set up devices and kernels */
  opencl_open_device(&ocl_dev);
  opencl_load_programs(&ocl_dev, "fdk_opencl.cl");
  opencl_kernel_create(&ocl_dev, "fdk_kernel_nn");

  /* Retrieve 2D image to get dimensions */
  const auto proj_dim = spCurImg->GetLargestPossibleRegion().GetSize();
  // proj = proj_dir->load_image(0);

  // Generate image sources for cone beam CT reconstruction
  std::array<float, 3> offset = {
      {static_cast<float>(-0.5 * sizeOutput[0] * spacing[0]),
       static_cast<float>(-0.5 * sizeOutput[1] * spacing[1]),
       static_cast<float>(-0.5 * sizeOutput[2] * spacing[2])}};

  auto *vol =
      new Volume(&sizeOutput[0], &offset[0], &spacing[0], nullptr, PT_FLOAT, 1);
  std::cout << "Buffer vol: " << vol->pix_size * vol->npix << std::endl;
  std::cout << "Buffer proj: " << proj_dim[1] * proj_dim[0] * sizeof(float)
            << std::endl;
  /* Set up device memory */
  auto ocl_buf_vol = opencl_buf_create(
      &ocl_dev, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
      static_cast<size_t>(vol->pix_size * vol->npix), vol->img);

  auto ocl_buf_img =
      opencl_buf_create(&ocl_dev, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR,
                        proj_dim[1] * proj_dim[0] * sizeof(float), // is 0
                        nullptr);

  auto ocl_buf_matrix =
      opencl_buf_create(&ocl_dev, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR,
                        12 * sizeof(float), nullptr);

  /* Copy volume dim (convert from size_t to int) */
  auto ocl_vol_dim = cl_int4{{static_cast<cl_int>(vol->dim[0]),
                              static_cast<cl_int>(vol->dim[1]),
                              static_cast<cl_int>(vol->dim[2])}};

  auto ocl_buf_vol_origin =
      opencl_buf_create(&ocl_dev, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                        3 * sizeof(float), &vol->origin[0]);

  auto ocl_buf_vol_spacing =
      opencl_buf_create(&ocl_dev, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                        3 * sizeof(float), &vol->spacing[0]);

  /* Copy projection image dim (convert from size_t to int) */
  auto ocl_proj_dim = cl_int2{
      {static_cast<cl_int>(proj_dim[0]), static_cast<cl_int>(proj_dim[1])}};

  auto ocl_buf_nrm =
      opencl_buf_create(&ocl_dev, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR,
                        3 * sizeof(float), nullptr);

  auto ocl_buf_ic =
      opencl_buf_create(&ocl_dev, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR,
                        2 * sizeof(float), nullptr);

  /* Calculate the scale */
  const auto scale =
      static_cast<float>(sqrt(3.0) / static_cast<double>(1 + proj_dim[2])) *
      10.0f;

  auto itShiftX = m_spCustomGeometry->GetProjectionOffsetsX().begin();
  auto itShiftY = m_spCustomGeometry->GetProjectionOffsetsY().begin();
  std::cout << "Proj. # " << std::endl;
  /* Project each image into the volume one at a time */
  for (auto image_num = 0; image_num < static_cast<int>(proj_dim[2]);
       image_num++) {
    /* Translate image and properties to PLM */
    auto *proj = new Proj_image;
    proj->init();
    proj->dim[0] = proj_dim[0];
    proj->dim[1] = proj_dim[1];
    proj->xy_offset[0] = *itShiftX;
    proj->xy_offset[1] = *itShiftY;
    itk::Index<3U> cur_idx = {{0, 0, image_num}};
    proj->img = &spCurImg->GetPixel(cur_idx);

    proj->pmat = new Proj_matrix;
    proj->pmat->ic[0] = 0.5 * proj_dim[0] - 0.5 + proj->xy_offset[0];
    proj->pmat->ic[1] = 0.5 * proj_dim[1] - 0.5 + proj->xy_offset[1];
    proj->pmat->sad =
        m_spCustomGeometry->GetSourceToIsocenterDistances()[image_num];

    const auto sid =
        m_spCustomGeometry->GetSourceToDetectorDistances()[image_num];
    const auto angle = m_spCustomGeometry->GetGantryAngles()[image_num];
    const double tgt[] = {0.0, 0.0, 0.0};
    const double cam[] = {tgt[0] + proj->pmat->sad * cos(angle),
                          tgt[1] - proj->pmat->sad * sin(angle), tgt[2]};
    const double vup[] = {0, 0, 1}; // R U sure about that?
    const double ps[] = {spCurImg->GetSpacing()[0], spCurImg->GetSpacing()[1]};

    proj->pmat->set(&cam[0], &tgt[0], &vup[0], sid, &proj->pmat->ic[0],
                    &ps[0]); // , proj->dim);

    if (image_num % (proj_dim[2] / 10) == 0) {
      std::cout << image_num << "::" << angle * 180 * itk::Math::one_over_pi
                << std::endl;
      proj->stats();
      std::cout << "Cam::" << proj->pmat->cam[0] << ", " << proj->pmat->cam[1]
                << ", " << proj->pmat->cam[2] << std::endl;
    }

    /* Apply ramp filter */
    proj_image_filter(proj);

    /* Copy image bytes to device */
    opencl_buf_write(&ocl_dev, ocl_buf_img,
                     proj_dim[1] * proj_dim[0] * sizeof(float), proj->img);

    /* Copy matrix to device (convert from double to float) */
    std::array<float, 12> matrix{};
    std::memcpy(&matrix.at(0), &proj->pmat->matrix[0], 12 * sizeof(float));

    opencl_buf_write(&ocl_dev, ocl_buf_matrix, 12 * sizeof(float), &matrix[0]);

    /* Copy ic to device (convert from double to float) */
    float ic[] = {static_cast<float>(proj->pmat->ic[0]),
                  static_cast<float>(proj->pmat->ic[1])};
    opencl_buf_write(&ocl_dev, ocl_buf_ic, 2 * sizeof(float), &ic[0]);

    /* Copy nrm to device (convert from double to float) */
    float nrm[] = {static_cast<float>(proj->pmat->nrm[0]),
                   static_cast<float>(proj->pmat->nrm[1]),
                   static_cast<float>(proj->pmat->nrm[2])};
    opencl_buf_write(&ocl_dev, ocl_buf_nrm, 3 * sizeof(float), &nrm[0]);

    /* Convert sad from double to float */
    const cl_float sad = proj->pmat->sad;

    /* Set fdk kernel arguments */
    opencl_set_kernel_args( // Marked as "dangerous" c-style function by
                            // cppcoreguidelines-pro-type-vararg
        &ocl_dev, sizeof(cl_mem), &ocl_buf_vol[0], sizeof(cl_mem),
        &ocl_buf_img[0], sizeof(cl_mem), &ocl_buf_matrix[0], sizeof(cl_int4),
        &ocl_vol_dim, sizeof(cl_mem), &ocl_buf_vol_origin[0], sizeof(cl_mem),
        &ocl_buf_vol_spacing[0], sizeof(cl_int2), &ocl_proj_dim, sizeof(cl_mem),
        &ocl_buf_nrm[0], sizeof(cl_mem), &ocl_buf_ic[0], sizeof(cl_float), &sad,
        sizeof(cl_float), &scale, static_cast<size_t>(0));

    /* Compute workgroup size */
    /* (Max local_work_size for (Greg's?) ATI RV710 is 128) */
    const size_t local_work_size = 128;
    const auto global_work_size = static_cast<size_t>(vol->npix);

    /* Invoke kernel */
    opencl_kernel_enqueue(&ocl_dev, global_work_size, local_work_size);
    ++itShiftX;
    ++itShiftY;

    delete proj->pmat;
    delete proj;
  }
  std::cout << "\nReading results..." << std::endl;
  /* Read back results */
  opencl_buf_read(&ocl_dev, ocl_buf_vol, vol->pix_size * vol->npix, vol->img);

  UNLOAD_LIBRARY(libplmopencl);
  return Plm_image(vol).itk_float();
  // End of opencl_reconstruct_conebeam
  //////////////////////// end PLM fdk /////////////////////////////
}
#endif
