/* -----------------------------------------------------------------------
   See COPYRIGHT.TXT and LICENSE.TXT for copyright and license information
   ----------------------------------------------------------------------- */

// Sweep through all voxel values and renormalize them to the Hounsfield (hu)
// scale.

__kernel void fdk_kernel_nn(__global float *dev_vol, __global float *dev_img,
                            __constant float *dev_matrix, int4 vol_dim,
                            __constant float *vol_offset,
                            __constant float *vol_spacing, int2 proj_dim,
                            __constant float *nrm, __constant float *ic,
                            const float sad, const float scale) {
  uint id = get_global_id(0);
  uint id_l = get_local_id(0);

  int k = id / vol_dim.x / vol_dim.y;
  int j = (id - (k * vol_dim.x * vol_dim.y)) / vol_dim.x;
  int i = id - k * vol_dim.x * vol_dim.y - j * vol_dim.x;

  if (k >= vol_dim.z) {
    return;
  }

  // Get volume value from global memory
  float dev_vol_value = dev_vol[id];

  // Get (x,y,z) coordinates
  float4 vp;
  vp.x = vol_offset[0] + (i * vol_spacing[0]);
  vp.y = vol_offset[1] + (j * vol_spacing[1]);
  vp.z = vol_offset[2] + (k * vol_spacing[2]);

  // Matrix multiplication
  float4 ip;
  ip.x = (dev_matrix[0] * vp.x) + (dev_matrix[1] * vp.y) +
         (dev_matrix[2] * vp.z) + dev_matrix[3];
  ip.y = (dev_matrix[4] * vp.x) + (dev_matrix[5] * vp.y) +
         (dev_matrix[6] * vp.z) + dev_matrix[7];
  ip.z = (dev_matrix[8] * vp.x) + (dev_matrix[9] * vp.y) +
         (dev_matrix[10] * vp.z) + dev_matrix[11];

  // Change coordinate systems
  ip.x = ic[0] + ip.x / ip.z;
  ip.y = ic[1] + ip.y / ip.z;

  // Retrieve pixel location from 2D image
  int2 pos;
  pos.y = convert_int_rtn(ip.x);
  pos.x = convert_int_rtn(ip.y);

  // Clip against image dimensions
  if (pos.x < 0 || pos.x >= proj_dim.x || pos.y < 0 || pos.y >= proj_dim.y) {
    return;
  }

  // Get pixel from image
  // float4 voxel_data = read_imagef(dev_img, dev_img_sampler, pos);
  float pix_val = dev_img[pos.y * proj_dim.x + pos.x];

  // Get distance to voxel projected to panel normal
  float s = (nrm[0] * vp.x) + (nrm[1] * vp.y) + (nrm[2] * vp.z);

  // Conebeam weighting factor
  s = sad - s;
  s = (sad * sad) / (s * s);

  // Weight pixel and backproject into volume
  dev_vol[id] = dev_vol_value + (scale * s * pix_val);
}

__constant sampler_t projectionSampler =
    CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP | CLK_FILTER_LINEAR;

__kernel void OpenCLFDKBackProjectionImageFilterKernel(
    __global float *volume, __constant float *matrix,
    __read_only image2d_t projection, uint4 volumeDim) {
  uint volumeIndex = get_global_id(0);

  uint i = volumeDim.x * volumeDim.y;
  uint k = volumeIndex / i;
  uint j = (volumeIndex - (k * i)) / volumeDim.x;
  i = volumeIndex - k * i - j * volumeDim.x;

  if (k >= volumeDim.z)
    return;

  float2 ip;
  float ipz;

  // matrix multiply
  ip.x = matrix[0] * i + matrix[1] * j + matrix[2] * k + matrix[3];
  ip.y = matrix[4] * i + matrix[5] * j + matrix[6] * k + matrix[7];
  ipz = matrix[8] * i + matrix[9] * j + matrix[10] * k + matrix[11];
  ipz = 1 / ipz;
  ip.x = ip.x * ipz;
  ip.y = ip.y * ipz;

  // Get projection value and add
  float4 projectionData = read_imagef(projection, projectionSampler, ip);

  volume[volumeIndex] += ipz * ipz * projectionData.x;
}
