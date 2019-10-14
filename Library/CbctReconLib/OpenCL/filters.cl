// Kernels for image filters, image support is not required for this file

__kernel void padding_kernel(__global float *input, __global float *output,
                             const int4 paddingIdx, const uint4 paddingDim,
                             const uint4 inputDim,
                             __constant float *truncationWeights,
                             const uint sizeWeights) {
  const unsigned int idx = get_global_id(0);
  int k = idx / (paddingDim.x * paddingDim.y);

  const unsigned int w = idx % (paddingDim.x * paddingDim.y);
  int j = w / paddingDim.x;
  int i = w % paddingDim.x;

  /*if (i >= paddingDim.x || j >= paddingDim.y || k >= paddingDim.z) {
    printf("Is it possible to get here? (%d, %d, %d)", i, j, k);
    return;
  }*/

  const unsigned long int out_idx = i + (j + k * paddingDim.y) * paddingDim.x;
  i -= paddingIdx.x;
  j -= paddingIdx.y;
  k -= paddingIdx.z;

  // out of input y/z dimensions
  if (j < 0 || j >= inputDim.y || k < 0 || k >= inputDim.z)
    output[out_idx] = 0.0f;
  // central part in CPU code
  else if (i >= 0 && i < inputDim.x)
    output[out_idx] = input[i + (j + k * inputDim.y) * inputDim.x];
  // left mirroring (equation 3a in [Ohnesorge et al, Med Phys, 2000])
  else if (i < 0 && -i < sizeWeights) {
    const int begRow = (j + k * inputDim.y) * inputDim.x;
    output[out_idx] =
        (2 * input[begRow + 1] - input[-i + begRow]) * truncationWeights[-i];
  }
  // right mirroring (equation 3b in [Ohnesorge et al, Med Phys, 2000])
  else if ((i >= inputDim.x) && (i - inputDim.x + 1) < sizeWeights) {
    const unsigned int borderDist = i - inputDim.x + 1;
    const int endRow = inputDim.x - 1 + (j + k * inputDim.y) * inputDim.x;
    output[out_idx] = (2 * input[endRow] - input[endRow - borderDist]) *
                      truncationWeights[borderDist];
  }
  // zero padding
  else
    output[out_idx] = 0.0f;
}

__kernel void multiply_kernel(__global float2 *projFFT, int4 fftDimension,
                              __global float2 *kernelFFT) {
  const unsigned int idx = get_global_id(0);
  // int k = idx / (fftDimension.x * fftDimension.y);

  const unsigned int w = idx % (fftDimension.x * fftDimension.y);
  int j = w / fftDimension.x;
  int i = w % fftDimension.x;

  // if (i >= fftDimension.x || j >= fftDimension.y || k >= fftDimension.z) {
  //  printf("Is it possible to get here? (%d, %d, %d)", i, j, k);
  //  return;
  //}

  // long int proj_idx = i + (j + k * fftDimension.y) * fftDimension.x;
  long int proj_idx = i + j * fftDimension.x;

  float2 result;
  result.x = projFFT[proj_idx].x * kernelFFT[i].x -
             projFFT[proj_idx].y * kernelFFT[i].y;
  result.y = projFFT[proj_idx].y * kernelFFT[i].x +
             projFFT[proj_idx].x * kernelFFT[i].y;
  projFFT[proj_idx] = result;
}

__kernel void multiply_kernel2D(__global float2 *projFFT, int4 fftDimension,
                                __global float2 *kernelFFT) {
  const unsigned int idx = get_global_id(0);
  // int k = idx / (fftDimension.x * fftDimension.y);

  const unsigned int w = idx % (fftDimension.x * fftDimension.y);
  int j = w / fftDimension.x;
  int i = w % fftDimension.x;

  // if (i >= fftDimension.x || j >= fftDimension.y || k >= fftDimension.z) {
  //  printf("Is it possible to get here? (%d, %d, %d)", i, j, k);
  //  return;
  //}

  long int kernel_idx = i + j * fftDimension.x;
  // long int proj_idx = kernel_idx + k * fftDimension.y * fftDimension.x;
  long int proj_idx = kernel_idx;

  float2 result;
  result.x = projFFT[proj_idx].x * kernelFFT[kernel_idx].x -
             projFFT[proj_idx].y * kernelFFT[kernel_idx].y;
  result.y = projFFT[proj_idx].y * kernelFFT[kernel_idx].x +
             projFFT[proj_idx].x * kernelFFT[kernel_idx].y;
  projFFT[proj_idx] = result;
}

__kernel void subtract_kernel2D(__global float *input, uint4 inputDimension,
                                __global float *subImg) {
  const unsigned int idx = get_global_id(0);
  // int k = idx / (inputDimension.x * inputDimension.y);

  const unsigned int w = idx % (inputDimension.x * inputDimension.y);
  // int j = w / inputDimension.x;
  // int i = w % inputDimension.x;

  // long int sub_idx = i + j * inputDimension.x;
  const float out_val = input[idx] - subImg[w];

  input[idx] = out_val;
}

// Actually is divide ln(65535/X) by ln(65535/Y)
__kernel void divide_kernel3D_loginv_Ushort(__global ushort *input,
                                            __global ushort *divImg,
                                            __global float *output) {
  const unsigned int idx = get_global_id(0);
  // mu_t = ln(65535/I)
  const float short_lim = log((float)USHRT_MAX);
  const float mu_in = short_lim - log((float)input[idx]);
  const float mu_div = short_lim - log((float)divImg[idx]);
  output[idx] = mu_in / mu_div;
}

// Adds value to input
__kernel void add_const_kernel(__global float *input, float value) {
  const unsigned int idx = get_global_id(0);
  const float out_val = input[idx] + value;

  input[idx] = out_val;
}

// Adds value to input or assigns to 0 or ln(65535) if outside bounds.
__kernel void add_const_with_thresh_kernel(__global float *input, float value) {
  const float short_lim = log((float)USHRT_MAX);
  const unsigned int idx = get_global_id(0);
  const float out_val = input[idx] + value;

  if (out_val < 0.0f) {
    input[idx] = 0.0f;
  } else if (out_val > short_lim) {
    input[idx] = short_lim;
  } else {
    input[idx] = out_val;
  }
}

// Adds add_value and then multiplies mul_value to input
__kernel void add_mul_const_kernel(__global float *input, const float add_value,
                                   const float mul_value) {
  const unsigned int idx = get_global_id(0);
  const float out_val = (input[idx] + add_value) * mul_value;
  input[idx] = out_val;
}

// Adds add_value and then multiplies mul_value to input or assigns to 0 or
// ln(65535) if result outside bounds.
__kernel void add_mul_const_with_thresh_kernel(__global float *input,
                                               const float add_value,
                                               const float mul_value) {
  const float short_lim = log((float)USHRT_MAX);
  const unsigned int idx = get_global_id(0);
  const float out_val = (input[idx] + add_value) * mul_value;
  if (out_val < 0.0f) {
    input[idx] = 0.0f;
  } else if (out_val > short_lim) {
    input[idx] = short_lim;
  } else {
    input[idx] = out_val;
  }
}

__kernel void min_max_kernel(__global const float *input,
                             __global float2 *subImg, const uint divider) {
  const uint idx = get_global_id(0);

  const uint first_idx_in = idx * divider;

  float local_arr[128]; // 128 is the max divider size

  int j = 0;
  for (int i = first_idx_in; i < (first_idx_in + divider); i++) {
    local_arr[j++] = input[i];
  }

  float min_val = FLT_MAX;
  float max_val = FLT_MIN;

  for (int i = 0; i < divider; i++) {
    if (min_val > local_arr[i])
      min_val = local_arr[i];
    if (max_val < local_arr[i])
      max_val = local_arr[i];
  }

  subImg[idx].x = min_val; // min
  subImg[idx].y = max_val; // max
}

__kernel void min_max_kernel2(__global const float2 *input,
                              __global float2 *subImg, uint divider,
                              const uint end_idx) {
  const uint idx = get_global_id(0);

  const uint first_idx_in = idx * divider;

  if (first_idx_in > end_idx) {
    subImg[idx].x = FLT_MAX; // min
    subImg[idx].y = FLT_MIN; // max
    return;
  }

  if ((first_idx_in + divider) > end_idx)
    divider = end_idx - first_idx_in;

  float2 local_arr[128]; // 128 is the max divider size

  int j = 0;
  for (int i = first_idx_in; i < (first_idx_in + divider); i++) {
    local_arr[j++] = input[i];
  }

  float min_val = FLT_MAX;
  float max_val = FLT_MIN;

  for (uint i = 0; i < divider; i++) {
    if (min_val > local_arr[i].x)
      min_val = local_arr[i].x;
    if (max_val < local_arr[i].y)
      max_val = local_arr[i].y;
  }

  subImg[idx].x = min_val; // min
  subImg[idx].y = max_val; // max
}

//// From http://geomalgorithms.com/a03-_inclusion.html
// Copyright 2000 softSurfer, 2012 Dan Sunday
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.

// a Point is defined by its coordinates {int x, y;}
//===================================================================

// isLeft(): tests if a point is Left|On|Right of an infinite line.
//    Input:  three points P0, P1, and P2
//    Return: >0 for P2 left of the line through P0 and P1
//            =0 for P2  on the line
//            <0 for P2  right of the line
//    See: Algorithm 1 "Area of Triangles and Polygons"
inline int isLeft(const float2 P0, const float2 P1, const float2 P2) {
  return ((P1.x - P0.x) * (P2.y - P0.y) - (P2.x - P0.x) * (P1.y - P0.y));
}
//===================================================================

// wn_PnPoly(): winding number test for a point in a polygon
//      Input:   P = a point,
//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//      Return:  wn = the winding number (=0 only when P is outside)
inline int wn_PnPoly(const float2 P, __constant float2 *V, const ulong n) {
  int wn = 0; // the  winding number counter

  // loop through all edges of the polygon
  for (ulong i = 0; i < n; i++) {          // edge from V[i] to  V[i+1]
    if (V[i].y <= P.y) {                   // start y <= P.y
      if (V[i + 1].y > P.y)                // an upward crossing
        if (isLeft(V[i], V[i + 1], P) > 0) // P left of  edge
          ++wn;                            // have  a valid up intersect
    } else {                               // start y > P.y (no test needed)
      if (V[i + 1].y <= P.y)               // a downward crossing
        if (isLeft(V[i], V[i + 1], P) < 0) // P right of  edge
          --wn;                            // have  a valid down intersect
    }
  }
  return wn;
}
//===================================================================

// Crop everything outside structure
__kernel void
crop_by_struct_kernel(__global ushort *dev_vol, __constant float2 *structure,
                      const ulong number_of_verti, const ulong4 vol_dim,
                      const float2 vol_offset, const float2 vol_spacing) {
  const ulong id = get_global_id(0);

  if (id >= vol_dim.x * vol_dim.y) {
    return;
  }

  if (number_of_verti < 3) {
    // If no struct, nothing is inside.
    // And a line or a point is not a struct
    dev_vol[id] = 0;
    return;
  }

  const ulong j = id / vol_dim.x;
  const ulong i = id - j * vol_dim.x;

  // Get (x,y,z) coordinates
  float2 vp;
  vp.x = vol_offset.x + (i * vol_spacing.x);
  vp.y = vol_offset.y + (j * vol_spacing.y);

  if (wn_PnPoly(vp, structure, number_of_verti) == 0) {
    // the point, vp, is outside the structure
    dev_vol[id] = 0;
  }
}

float median_by_bubble_sort(float list[], const uint n) {
  uint n_i = n;
  while (n_i >= n / 2) {
    uint new_n = 0;
    for (uint i = 1; i < n_i; ++i) {
      if (list[i - 1] > list[i]) {
        // Swap
        float t = list[i];
        list[i] = list[i - 1];
        list[i - 1] = t;
        new_n = i;
      }
    }
    n_i = new_n;
  }
  return list[n / 2];
}

// Convert input to intensity, then take the difference and apply a median
// filter then convert back to line integral
__kernel void log_i_to_i_subtract_median_i_to_log_i(
    __global const float *proj_raw, __global const float *proj_sca,
    __global float *proj_corr, const ulong2 size, const uint median_radius) {
  const ulong id = get_global_id(0);

  if (id >= size.x * size.y) {
    return;
  }
  const long j = id / size.x;
  const long i = id - j * size.x;

  // Optimized for small median radius, like 3
  float unsorted_vector[64];
  uint i_uv = 0;

  for (int y = -median_radius; y <= (int)median_radius; ++y) {
    const int jy = j + y;
    if (jy < 0 || jy > size.y) {
      continue;
    }
    const ulong id_jy = jy * size.x;
    for (int x = -median_radius; x <= (int)median_radius; ++x) {
      const int ix = i + x;
      if (ix < 0 || ix > size.x) {
        continue;
      }
      const ulong r_id = ix + id_jy;
      const float raw_val = proj_raw[r_id];

      // assuming *_val is not < 0, then *_i is in [0; 1]
      const float raw_i = exp(-raw_val);

      const float sca_val = proj_sca[r_id];
      const float sca_i = exp(-sca_val);

      unsorted_vector[i_uv++] = raw_i - sca_i;
    }
  }

  // lowest possible value of i_uv should be median_radius^2
  const float median = median_by_bubble_sort(unsorted_vector, i_uv);

  if (median > 0.0) {
    proj_corr[id] = -log(median);
  } else {
    proj_corr[id] = FLT_MAX;
  }
}

// Convert input to intensity, then take the difference and apply a median
// filter first in y direction, then x (the other way around is harder to
// imagine)
__kernel void log_i_to_i_subtract_median_y_x(__global const float *proj_raw,
                                             __global const float *proj_prim,
                                             __global float *proj_sca,
                                             const ulong2 size,
                                             const uint median_radius) {
  const ulong local_id = get_local_id(0);
  const ulong group_id = get_group_id(0);
  const uint actual_local_size = (size.x * size.y) / get_num_groups(0);

  const long first_id_in_group = actual_local_size * group_id - median_radius;
  const long id = first_id_in_group + local_id;

  if (id >= size.x * size.y || id < 0) {
    return;
  }
  const long j = id / size.x;
  const long i = id - j * size.x;

  __local float
      localbuffer[128];      // must be >= actual_local_size + 2*median_radius
  float unsorted_vector[64]; // must be larger than 2*median_radius + 1
  uint i_uv = 0;

  for (int y = -median_radius; y <= (int)median_radius; ++y) {
    const int jy = j + y;
    if (jy < 0 || jy >= size.y) {
      continue;
    }
    const ulong r_id = i + jy * size.x;

    const float raw_val = proj_raw[r_id];
    // assuming *_val is not < 0, then *_i is in [0; 1]
    const float raw_i = exp(-raw_val);

    const float prim_val = proj_prim[r_id];
    const float prim_i = exp(-prim_val);

    unsorted_vector[i_uv++] = raw_i - prim_i;
  }

  // lowest possible value of i_uv should be median_radius^2
  localbuffer[local_id] = median_by_bubble_sort(unsorted_vector, i_uv);

  barrier(CLK_LOCAL_MEM_FENCE);

  if (local_id < median_radius ||
      local_id >= (median_radius + actual_local_size)) {
    return;
  }

  i_uv = 0;
  for (int x = -median_radius; x <= (int)median_radius; ++x) {
    const int ix = local_id + x;
    if (ix < 0 || ix >= 128) {
      continue;
    }
    unsorted_vector[i_uv++] = localbuffer[ix];
  }

  // lowest possible value of i_uv should be median_radius^2
  const float median = median_by_bubble_sort(unsorted_vector, i_uv);

  // Do not convert back, we need the data in intensity for median_y and the
  // gaussian
  proj_sca[id] = median;
}

__kernel void i_to_log_i_kernel(__global float *proj_sca, const ulong2 size) {
  const ulong id = get_global_id(0);
  if (id >= size.x * size.y) {
    return;
  }

  const float i_val = proj_sca[id];
  if (i_val > 0.0) {
    proj_sca[id] = -log(i_val);
  } else {
    proj_sca[id] = FLT_MAX;
  }
}
