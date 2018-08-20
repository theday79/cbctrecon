/*Utility functions for cbctrecon*/

template <typename ImageType>
double
GetFOVRadius(const rtk::ThreeDCircularProjectionGeometry::Pointer &geometry,
             const typename ImageType::Pointer &ProjStack) {

  using FOVfilterType = rtk::FieldOfViewImageFilter<ImageType, ImageType>;
  typename FOVfilterType::Pointer FOVfilter = FOVfilterType::New();
  FOVfilter->SetGeometry(geometry);
  FOVfilter->SetProjectionsStack(ProjStack.GetPointer());
  double x, z;
  double r_inf = -1.0;
  bool hasOverlap = FOVfilter->ComputeFOVRadius(
      FOVfilterType::FOVRadiusType::RADIUSINF, x, z, r_inf);
  // halffan gives r(BOTH)~25, r(SUP)~25, r(INF)~232 -> RADIUSINF also seems to
  // work for fullfan, so we'll use that.

  if (hasOverlap) {
    std::cout << "FOV (inf) radius was found: r=" << r_inf << std::endl;
  }
  double r_sup = -1.0;
  hasOverlap = FOVfilter->ComputeFOVRadius(
      FOVfilterType::FOVRadiusType::RADIUSSUP, x, z, r_sup);
  if (hasOverlap) {
    std::cout << "FOV (sup) radius was found: r=" << r_sup << std::endl;
  }
  return std::max(r_inf, r_sup);
}

template <typename T, typename ImageType>
bool GetOutputResolutionFromFOV(
    typename T::SizeType &sizeOutput, typename T::SpacingType &spacing,
    const rtk::ThreeDCircularProjectionGeometry::Pointer &geometry,
    const typename ImageType::Pointer &ProjStack,
    const QString &outputFilePath) {

  QFileInfo outFileInfo(outputFilePath);
  QDir outFileDir = outFileInfo.absoluteDir();

  if (outputFilePath.length() < 2 || !outFileDir.exists()) {
    double radius = GetFOVRadius<ImageType>(geometry, ProjStack);
    if (radius != -1.0) {
      sizeOutput[0] = 512; // AP
      sizeOutput[1] = 200; // SI
      sizeOutput[2] = 512; // LR
      spacing[0] = 2.0 * radius / sizeOutput[0];
      spacing[1] = 1.0;
      spacing[2] = 2.0 * radius / sizeOutput[2];
      return true;
    }
  }

  return false;
}

#if USE_OPENCL_PLM
FloatImageType::Pointer PlastimatchOpenCLFDK(
    const FloatImageType::Pointer &spCurImg,
    const rtk::ThreeDCircularProjectionGeometry::Pointer &m_spCustomGeometry,
    std::array<float, 3> spacing, std::array<plm_long, 3> sizeOutput) {
  spCurImg->Update();

  //////////////////////// start PLM fdk /////////////////////////////
  // Copy of opencl_reconstruct_conebeam from fdk_opencl.cxx in plastimatch
  // (because we already have loaded projections into memory)
  Opencl_device ocl_dev{};
  Opencl_buf *ocl_buf_vol;
  Opencl_buf *ocl_buf_img;
  Opencl_buf *ocl_buf_matrix;
  cl_int4 ocl_vol_dim{};
  Opencl_buf *ocl_buf_vol_origin;
  Opencl_buf *ocl_buf_vol_spacing;
  cl_int2 ocl_proj_dim{};
  Opencl_buf *ocl_buf_nrm;
  Opencl_buf *ocl_buf_ic;
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
  const itk::Size<3U> proj_dim = spCurImg->GetLargestPossibleRegion().GetSize();
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
  ocl_buf_vol =
      opencl_buf_create(&ocl_dev, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
                        vol->pix_size * vol->npix, vol->img);

  ocl_buf_img =
      opencl_buf_create(&ocl_dev, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR,
                        proj_dim[1] * proj_dim[0] * sizeof(float), // is 0
                        nullptr);

  ocl_buf_matrix =
      opencl_buf_create(&ocl_dev, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR,
                        12 * sizeof(float), nullptr);

  /* Copy volume dim (convert from size_t to int) */
  ocl_vol_dim.x = vol->dim[0];
  ocl_vol_dim.y = vol->dim[1];
  ocl_vol_dim.z = vol->dim[2];

  ocl_buf_vol_origin =
      opencl_buf_create(&ocl_dev, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                        3 * sizeof(float), &vol->origin[0]);

  ocl_buf_vol_spacing =
      opencl_buf_create(&ocl_dev, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                        3 * sizeof(float), &vol->spacing[0]);

  /* Copy projection image dim (convert from size_t to int) */
  ocl_proj_dim.x = proj_dim[0];
  ocl_proj_dim.y = proj_dim[1];

  ocl_buf_nrm =
      opencl_buf_create(&ocl_dev, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR,
                        3 * sizeof(float), nullptr);

  ocl_buf_ic =
      opencl_buf_create(&ocl_dev, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR,
                        2 * sizeof(float), nullptr);

  /* Calculate the scale */
  const float scale =
      static_cast<float>(sqrt(3.0) / static_cast<double>(1 + proj_dim[2])) *
      10.0f;

  std::vector<double>::const_iterator itShiftX, itShiftY;
  itShiftX = m_spCustomGeometry->GetProjectionOffsetsX().begin();
  itShiftY = m_spCustomGeometry->GetProjectionOffsetsY().begin();
  std::cout << "Proj. # " << std::endl;
  /* Project each image into the volume one at a time */
  for (int image_num = 0; image_num < static_cast<int>(proj_dim[2]);
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
    if (proj->img == ITK_NULLPTR) {
      std::cout << "No image!" << std::endl;
    }

    proj->pmat = new Proj_matrix;
    proj->pmat->ic[0] = 0.5 * proj_dim[0] - 0.5 + proj->xy_offset[0];
    proj->pmat->ic[1] = 0.5 * proj_dim[1] - 0.5 + proj->xy_offset[1];
    proj->pmat->sad =
        m_spCustomGeometry->GetSourceToIsocenterDistances()[image_num];

    const double sid =
        m_spCustomGeometry->GetSourceToDetectorDistances()[image_num];
    const double angle = m_spCustomGeometry->GetGantryAngles()[image_num];
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
    std::memcpy(matrix.data(), &proj->pmat->matrix[0], 12);

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
    size_t local_work_size = 128;
    size_t global_work_size = static_cast<float>(vol->npix);

    /* Invoke kernel */
    opencl_kernel_enqueue(&ocl_dev, global_work_size, local_work_size);
    itShiftX++;
    itShiftY++;

    delete proj->pmat;
    delete proj;
  }
  std::cout << std::endl;
  std::cout << "Reading results..." << std::endl;
  /* Read back results */
  opencl_buf_read(&ocl_dev, ocl_buf_vol, vol->pix_size * vol->npix, vol->img);

  UNLOAD_LIBRARY(libplmopencl);
  return Plm_image(vol).itk_float();
  // End of opencl_reconstruct_conebeam
  //////////////////////// end PLM fdk /////////////////////////////
}
#endif

#if USE_OPENCL_RTK
FloatImageType::Pointer RTKOpenCLFDK(
    const FloatImageType::Pointer &spCurImg,
    const rtk::ThreeDCircularProjectionGeometry::Pointer &m_spCustomGeometry,
    FloatImageType::SpacingType spacing, FloatImageType::SizeType sizeOutput,
    std::array<const double, 5> fdk_options) {

  try {
    spCurImg->Update();
  } catch (const std::exception &err) {
    std::cerr << "Couldn't update spCurImg: " << err.what() << std::endl;
  }
  // Generate image sources for cone beam CT reconstruction
  using ConstantImageSourceType = rtk::ConstantImageSource<FloatImageType>;
  ConstantImageSourceType::PointType origin;

  origin[0] = -0.5 * sizeOutput[0] * spacing[0]; // Y in DCM?
  origin[1] = -0.5 * sizeOutput[1] * spacing[1]; // Z in DCM?
  origin[2] = -0.5 * sizeOutput[2] * spacing[2]; // X in DCM?

  ConstantImageSourceType::Pointer constantImageSource =
      ConstantImageSourceType::New();
  constantImageSource->SetOrigin(origin);
  constantImageSource->SetSpacing(spacing);
  constantImageSource->SetSize(sizeOutput);
  constantImageSource->SetConstant(0.0); // initial value
  // constantImageSource->Update();

  // just for expressiveness
  const double fTruncCorFactor = fdk_options.at(0);
  const double fHannCut = fdk_options.at(1);
  const double fCosineCut = fdk_options.at(2);
  const double fHamming = fdk_options.at(3);
  const double fHannCutY = fdk_options.at(4);

  // FDK reconstruction filtering
  using FDKOPENCLType = rtk::OpenCLFDKConeBeamReconstructionFilter;
  FDKOPENCLType::Pointer feldkampOCL;
  feldkampOCL = FDKOPENCLType::New();

  feldkampOCL->SetInput(0, constantImageSource->GetOutput());
  feldkampOCL->SetInput(1, spCurImg);
  feldkampOCL->SetGeometry(m_spCustomGeometry);
  feldkampOCL->GetRampFilter()->SetTruncationCorrection(fTruncCorFactor);
  feldkampOCL->GetRampFilter()->SetHannCutFrequency(fHannCut);
  feldkampOCL->GetRampFilter()->SetHannCutFrequencyY(fHannCutY);
  feldkampOCL->GetRampFilter()->SetCosineCutFrequency(fCosineCut);
  feldkampOCL->GetRampFilter()->SetHammingFrequency(fHamming);

  feldkampOCL->Update();
  // feldkampOCL->PrintTiming(std::cout); Deprecated in rtk 1.4

  return feldkampOCL->GetOutput();
}
#endif

double
CbctRecon::GetValueFrom3DImageFloat(int reqX, int reqY, int reqZ,
                                    FloatImageType::Pointer &sp3DFloatImage) {
  if (sp3DFloatImage == nullptr) {
    return -1.0;
  }

  itk::ImageSliceConstIteratorWithIndex<FloatImageType> it(
      sp3DFloatImage, sp3DFloatImage->GetBufferedRegion());

  it.SetFirstDirection(0);  // x?
  it.SetSecondDirection(1); // y?
  it.GoToBegin();

  int idxX, idxY, idxZ;
  idxZ = 0;

  while (!it.IsAtEnd()) {
    if (idxZ == reqZ) {
      idxY = 0;
      while (!it.IsAtEndOfSlice()) {
        if (idxY == reqY) {
          idxX = 0;
          while (!it.IsAtEndOfLine()) {
            if (idxX == reqX) {
              double tmpVal = it.Get();
              return tmpVal;
            }
            ++it;
            idxX++;
          }
          break;
        }
        it.NextLine();
        idxY++;
      }
      break;
    }
    it.NextSlice();
    idxZ++;
  }

  return -2.0;
}

double CbctRecon::GetValueFrom3DImageUshort(
    int reqX, int reqY, int reqZ, UShortImageType::Pointer &sp3DUshortImage) {
  if (sp3DUshortImage == nullptr) {
    return 0;
  }

  itk::ImageSliceConstIteratorWithIndex<UShortImageType> it(
      sp3DUshortImage, sp3DUshortImage->GetBufferedRegion());

  it.SetFirstDirection(0);  // x?
  it.SetSecondDirection(1); // y?
  it.GoToBegin();

  int idxX, idxY, idxZ;
  idxZ = 0;

  while (!it.IsAtEnd()) {
    if (idxZ == reqZ) {
      idxY = 0;
      while (!it.IsAtEndOfSlice()) {
        if (idxY == reqY) {

          idxX = 0;
          while (!it.IsAtEndOfLine()) {
            if (idxX == reqX) {
              unsigned short tmpVal = it.Get();
              return tmpVal;
            }
            ++it;
            idxX++;
          }
          break;
        }
        it.NextLine();
        idxY++;
      }
      break;
    }
    it.NextSlice();
    idxZ++;
  }
  return 65535;
}

double WEPL_from_point(const std::array<size_t, 3> cur_point_id,
                       const std::array<double, 3> vec_basis,
                       const std::array<double, 3> vec_cubesize,
                       const std::array<size_t, 3> cubedim,
                       const FloatImageType::Pointer &vec_wepl_cube) {
  const double step_length = 0.1;
  const std::array<double, 3> step = {{vec_basis.at(0) * step_length,
                                       vec_basis.at(1) * step_length,
                                       vec_basis.at(2) * step_length}};

  const std::array<double, 3> inv_cubesize = {{1.0 / vec_cubesize.at(0),
                                               1.0 / vec_cubesize.at(1),
                                               1.0 / vec_cubesize.at(2)}};

  std::array<double, 3> point = {
      {static_cast<double>(cur_point_id.at(0)) * vec_cubesize.at(0),
       static_cast<double>(cur_point_id.at(1)) * vec_cubesize.at(1),
       static_cast<double>(cur_point_id.at(2)) * vec_cubesize.at(2)}};

  double out = 0.0;

  while (true) {
    // point_id = point / cube_size
    const std::array<int, 3> point_id = {
        {static_cast<int>(round(point.at(0) * inv_cubesize.at(0))),
         static_cast<int>(round(point.at(1) * inv_cubesize.at(1))),
         static_cast<int>(round(point.at(2) * inv_cubesize.at(2)))}};

    if (point_id.at(0) < 0.0 ||
        point_id.at(0) >= static_cast<int>(cubedim.at(0)) ||
        point_id.at(1) < 0.0 ||
        point_id.at(1) >= static_cast<int>(cubedim.at(1)) ||
        point_id.at(2) < 0.0 ||
        point_id.at(2) >= static_cast<int>(cubedim.at(2))) {
      break;
    }

    // get nearest neighbors:
    const std::array<double, 3> point_id_pos = {
        {point.at(0) * inv_cubesize.at(0), point.at(1) * inv_cubesize.at(1),
         point.at(2) * inv_cubesize.at(2)}};

    int idx_2 = -1;
    if (point_id.at(0) < (point_id_pos.at(0))) {
      idx_2 = 1;
    }

    int idy_2 = -1;
    if (point_id.at(1) < (point_id_pos.at(1))) {
      idx_2 = 1;
    }

    int idz_2 = -1;
    if (point_id.at(2) < (point_id_pos.at(2))) {
      idz_2 = 1;
    }

    if ((point_id.at(0) + idx_2) < 0.0 ||
        (point_id.at(0) + idx_2) >= static_cast<int>(cubedim.at(0)) ||
        (point_id.at(1) + idy_2) < 0.0 ||
        (point_id.at(1) + idy_2) >= static_cast<int>(cubedim.at(1)) ||
        (point_id.at(2) + idz_2) < 0.0 ||
        (point_id.at(2) + idz_2) >= static_cast<int>(cubedim.at(2))) {
      break;
    }

    std::array<double, 8> weights{};
    // x                    xyz
    weights.at(0) = sqrt( // 000 =
        pow(point_id.at(0) - point_id_pos.at(0), 2) +
        pow(point_id.at(1) - point_id_pos.at(1), 2) +
        pow(point_id.at(2) - point_id_pos.at(2), 2));
    weights.at(1) = sqrt( // 100 =
        pow(point_id.at(0) + idx_2 - point_id_pos.at(0), 2) +
        pow(point_id.at(1) - point_id_pos.at(1), 2) +
        pow(point_id.at(2) - point_id_pos.at(2), 2));
    // y
    weights.at(2) = sqrt( // 010 =
        pow(point_id.at(0) - point_id_pos.at(0), 2) +
        pow(point_id.at(1) + idy_2 - point_id_pos.at(1), 2) +
        pow(point_id.at(2) - point_id_pos.at(2), 2));
    weights.at(3) = sqrt( // 110 =
        pow(point_id.at(0) + idx_2 - point_id_pos.at(0), 2) +
        pow(point_id.at(1) + idy_2 - point_id_pos.at(1), 2) +
        pow(point_id.at(2) - point_id_pos.at(2), 2));
    // z
    weights.at(4) = sqrt( // 001 =
        pow(point_id.at(0) - point_id_pos.at(0), 2) +
        pow(point_id.at(1) - point_id_pos.at(1), 2) +
        pow(point_id.at(2) + idz_2 - point_id_pos.at(2), 2));
    weights.at(5) = sqrt( // 101 =
        pow(point_id.at(0) + idx_2 - point_id_pos.at(0), 2) +
        pow(point_id.at(1) - point_id_pos.at(1), 2) +
        pow(point_id.at(2) + idz_2 - point_id_pos.at(2), 2));
    weights.at(6) = sqrt( // 011 =
        pow(point_id.at(0) - point_id_pos.at(0), 2) +
        pow(point_id.at(1) + idy_2 - point_id_pos.at(1), 2) +
        pow(point_id.at(2) + idz_2 - point_id_pos.at(2), 2));
    weights.at(7) = sqrt( // 111 =
        pow(point_id.at(0) + idx_2 - point_id_pos.at(0), 2) +
        pow(point_id.at(1) + idy_2 - point_id_pos.at(1), 2) +
        pow(point_id.at(2) + idz_2 - point_id_pos.at(2), 2));

    const double sum_weights =
        std::accumulate(weights.begin(), weights.end(), 0.0);

    std::transform(
        weights.begin(), weights.end(), weights.begin(),
        [&sum_weights](double val) -> double { return val / sum_weights; });

    for (int i = 0; i < 8; i++) {
      // convert point_id to cube_ids
      FloatImageType::IndexType cube_id{};
      cube_id[0] = point_id.at(0) + (idx_2 * (i % 2)); // x= 0,1,0,1,0,1,0,1
      cube_id[1] =
          point_id.at(1) + (idy_2 * ((i / 2) % 2));    // y= 0,0,1,1,0,0,1,1
      cube_id[2] = point_id.at(2) + (idz_2 * (i / 4)); // z= 0,0,0,0,1,1,1,1

      out += vec_wepl_cube->GetPixel(cube_id) * weights.at(i);
    }

    // point = point - step
    point.at(0) -= step.at(0);
    point.at(1) -= step.at(1);
    point.at(2) -= step.at(2);
  }

  return out * step_length;
}

std::array<double, 3> get_basis_from_angles(double gantry, double couch) {
  gantry += 180.0;
  gantry *= M_PI / 180.0;
  couch *= M_PI / 180.0;

  std::array<double, 3> basis = {
      {sin(gantry) * cos(couch), -cos(gantry), sin(couch) * sin(gantry)}};
  return basis;
}

std::vector<double>
WEPL_trace_from_point(const std::array<size_t, 3> cur_point_id,
                      const std::array<double, 3> vec_basis,
                      const std::array<double, 3> vec_cubesize,
                      const std::array<size_t, 3> cubedim,
                      const FloatImageType::Pointer &vec_wepl_cube) {

  std::vector<double> cumWEPL; // cumulative WEPL

  const double step_length = 0.1;
  const std::array<double, 3> step = {{vec_basis.at(0) * step_length,
                                       vec_basis.at(1) * step_length,
                                       vec_basis.at(2) * step_length}};

  const std::array<double, 3> inv_cubesize = {{1.0 / vec_cubesize.at(0),
                                               1.0 / vec_cubesize.at(1),
                                               1.0 / vec_cubesize.at(2)}};

  std::array<double, 3> point = {
      {static_cast<double>(cur_point_id.at(0)) * vec_cubesize.at(0),
       static_cast<double>(cur_point_id.at(1)) * vec_cubesize.at(1),
       static_cast<double>(cur_point_id.at(2)) * vec_cubesize.at(2)}};

  double out_point = 0.0;
  while (true) {
    // point_id = point / cube_size
    const std::array<int, 3> point_id = {
        {static_cast<int>(round(point.at(0) * inv_cubesize.at(0))),
         static_cast<int>(round(point.at(1) * inv_cubesize.at(1))),
         static_cast<int>(round(point.at(2) * inv_cubesize.at(2)))}};

    if (point_id.at(0) < 0.0 ||
        point_id.at(0) >= static_cast<int>(cubedim.at(0)) ||
        point_id.at(1) < 0.0 ||
        point_id.at(1) >= static_cast<int>(cubedim.at(1)) ||
        point_id.at(2) < 0.0 ||
        point_id.at(2) >= static_cast<int>(cubedim.at(2))) {
      break;
    }

    // get nearest neighbors:
    const std::array<double, 3> point_id_pos = {
        {point.at(0) * inv_cubesize.at(0), point.at(1) * inv_cubesize.at(1),
         point.at(2) * inv_cubesize.at(2)}};

    int idx_2 = -1;
    if (point_id.at(0) < (point_id_pos.at(0))) {
      idx_2 = 1;
    }

    int idy_2 = -1;
    if (point_id.at(1) < (point_id_pos.at(1))) {
      idx_2 = 1;
    }

    int idz_2 = -1;
    if (point_id.at(2) < (point_id_pos.at(2))) {
      idz_2 = 1;
    }

    if ((point_id.at(0) + idx_2) < 0.0 ||
        (point_id.at(0) + idx_2) >= static_cast<int>(cubedim.at(0)) ||
        (point_id.at(1) + idy_2) < 0.0 ||
        (point_id.at(1) + idy_2) >= static_cast<int>(cubedim.at(1)) ||
        (point_id.at(2) + idz_2) < 0.0 ||
        (point_id.at(2) + idz_2) >= static_cast<int>(cubedim.at(2))) {
      break;
    }

    std::array<double, 8> weights{};
    // x                    xyz
    weights.at(0) = sqrt( // 000 =
        pow(point_id.at(0) - point_id_pos.at(0), 2) +
        pow(point_id.at(1) - point_id_pos.at(1), 2) +
        pow(point_id.at(2) - point_id_pos.at(2), 2));
    weights.at(1) = sqrt( // 100 =
        pow(point_id.at(0) + idx_2 - point_id_pos.at(0), 2) +
        pow(point_id.at(1) - point_id_pos.at(1), 2) +
        pow(point_id.at(2) - point_id_pos.at(2), 2));
    // y
    weights.at(2) = sqrt( // 010 =
        pow(point_id.at(0) - point_id_pos.at(0), 2) +
        pow(point_id.at(1) + idy_2 - point_id_pos.at(1), 2) +
        pow(point_id.at(2) - point_id_pos.at(2), 2));
    weights.at(3) = sqrt( // 110 =
        pow(point_id.at(0) + idx_2 - point_id_pos.at(0), 2) +
        pow(point_id.at(1) + idy_2 - point_id_pos.at(1), 2) +
        pow(point_id.at(2) - point_id_pos.at(2), 2));
    // z
    weights.at(4) = sqrt( // 001 =
        pow(point_id.at(0) - point_id_pos.at(0), 2) +
        pow(point_id.at(1) - point_id_pos.at(1), 2) +
        pow(point_id.at(2) + idz_2 - point_id_pos.at(2), 2));
    weights.at(5) = sqrt( // 101 =
        pow(point_id.at(0) + idx_2 - point_id_pos.at(0), 2) +
        pow(point_id.at(1) - point_id_pos.at(1), 2) +
        pow(point_id.at(2) + idz_2 - point_id_pos.at(2), 2));
    weights.at(6) = sqrt( // 011 =
        pow(point_id.at(0) - point_id_pos.at(0), 2) +
        pow(point_id.at(1) + idy_2 - point_id_pos.at(1), 2) +
        pow(point_id.at(2) + idz_2 - point_id_pos.at(2), 2));
    weights.at(7) = sqrt( // 111 =
        pow(point_id.at(0) + idx_2 - point_id_pos.at(0), 2) +
        pow(point_id.at(1) + idy_2 - point_id_pos.at(1), 2) +
        pow(point_id.at(2) + idz_2 - point_id_pos.at(2), 2));

    const double sum_weights =
        std::accumulate(weights.begin(), weights.end(), 0.0);

    std::transform(
        weights.begin(), weights.end(), weights.begin(),
        [&sum_weights](double val) -> double { return val / sum_weights; });

    for (int i = 0; i < 8; i++) {
      // convert point_id to cube_ids
      FloatImageType::IndexType cube_id{};
      cube_id[0] = point_id.at(0) + (idx_2 * (i % 2)); // x= 0,1,0,1,0,1,0,1
      cube_id[1] =
          point_id.at(1) + (idy_2 * ((i / 2) % 2));    // y= 0,0,1,1,0,0,1,1
      cube_id[2] = point_id.at(2) + (idz_2 * (i / 4)); // z= 0,0,0,0,1,1,1,1

      out_point += vec_wepl_cube->GetPixel(cube_id) * weights.at(i);
    }

    cumWEPL.push_back(out_point);

    // point = point - step
    point.at(0) -= step.at(0);
    point.at(1) -= step.at(1);
    point.at(2) -= step.at(2);
  }
  std::valarray<double> vdiff(cumWEPL.size());
  std::adjacent_difference(cumWEPL.begin(), cumWEPL.end(), &vdiff[0]);

  std::vector<double> out_vec;
  double revWEPL = cumWEPL.at(cumWEPL.size() - 1);
  size_t idx = 0;
  for (auto val : vdiff) {
    revWEPL -= val;
    out_vec.push_back(revWEPL * step_length);
  }

  return out_vec;
}

struct WEPLVector {
  double WEPL;
  FloatVector point;
};

std::vector<WEPLVector>
WEPLContourFromRtssContour(Rtss_contour_modern rt_contour,
                           const std::array<double, 3> vec_basis,
                           const FloatImageType::Pointer &vec_wepl_cube) {

  const std::array<double, 3> pixel_size = {{vec_wepl_cube->GetSpacing()[0],
                                             vec_wepl_cube->GetSpacing()[1],
                                             vec_wepl_cube->GetSpacing()[2]}};

  const std::array<size_t, 3> cubedim = {
      {vec_wepl_cube->GetLargestPossibleRegion().GetSize()[0],
       vec_wepl_cube->GetLargestPossibleRegion().GetSize()[1],
       vec_wepl_cube->GetLargestPossibleRegion().GetSize()[2]}};

  std::vector<WEPLVector> WEPL_contour;

  for (auto point : rt_contour.coordinates) {
    FloatImageType::PointType p;
    p.SetElement(0, point.x);
    p.SetElement(1, point.y);
    p.SetElement(2, point.z);
    FloatImageType::IndexType cur_idx{};
    vec_wepl_cube->TransformPhysicalPointToIndex(p, cur_idx);
    const std::array<size_t, 3> point_id = {{static_cast<size_t>(cur_idx[0]),
                                             static_cast<size_t>(cur_idx[1]),
                                             static_cast<size_t>(cur_idx[2])}};
    double wepl = WEPL_from_point(point_id, vec_basis, pixel_size, cubedim,
                                  vec_wepl_cube);
    WEPL_contour.emplace_back(wepl, point);
  }

  return WEPL_contour;
}
