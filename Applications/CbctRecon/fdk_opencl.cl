/* -----------------------------------------------------------------------
   See COPYRIGHT.TXT and LICENSE.TXT for copyright and license information
   ----------------------------------------------------------------------- */
// Sweep through all voxel values and renormalize them to the Hounsfield (hu) scale.

__kernel void 
fdk_kernel_nn (
    __global float *dev_vol,
    __global float *dev_img,
    __constant float *dev_matrix,
    int4 vol_dim,
    __constant float *vol_offset,
    __constant float *vol_spacing,
    int2 proj_dim,
    __constant float *nrm,
    __constant float *ic,
    const float sad,
    const float scale
)
{
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
    ip.x = (dev_matrix[0] * vp.x) + (dev_matrix[1] * vp.y) + (dev_matrix[2] * vp.z) + dev_matrix[3];
    ip.y = (dev_matrix[4] * vp.x) + (dev_matrix[5] * vp.y) + (dev_matrix[6] * vp.z) + dev_matrix[7];
    ip.z = (dev_matrix[8] * vp.x) + (dev_matrix[9] * vp.y) + (dev_matrix[10] * vp.z) + dev_matrix[11];

    // Change coordinate systems
    ip.x = ic[0] + ip.x / ip.z;
    ip.y = ic[1] + ip.y / ip.z;

    // Retrieve pixel location from 2D image
    int2 pos;
    pos.y = convert_int_rtn (ip.x);
    pos.x = convert_int_rtn (ip.y);

    // Clip against image dimensions
    if (pos.x < 0 || pos.x >= proj_dim.x || pos.y < 0 || pos.y >= proj_dim.y)
    {
        return;
    }

    // Get pixel from image
    //float4 voxel_data = read_imagef(dev_img, dev_img_sampler, pos);
    float pix_val = dev_img[pos.y * proj_dim.x + pos.x];

    // Get distance to voxel projected to panel normal 
    float s = (nrm[0] * vp.x) + (nrm[1] * vp.y) + (nrm[2] * vp.z);

    // Conebeam weighting factor
    s = sad - s;
    s = (sad * sad) / (s * s);

    // Weight pixel and backproject into volume
    dev_vol[id] = dev_vol_value + (scale * s * pix_val);
}

__constant sampler_t projectionSampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP | CLK_FILTER_LINEAR;

__kernel void
OpenCLFDKBackProjectionImageFilterKernel(__global float *volume,
                                         __constant float *matrix,
                                         __read_only image2d_t projection,
                                         uint4 volumeDim)
{
  uint volumeIndex = get_global_id(0);

  uint i = volumeDim.x * volumeDim.y;
  uint k = volumeIndex / i;
  uint j = ( volumeIndex - (k * i) ) / volumeDim.x;
  i = volumeIndex - k * i - j * volumeDim.x;

  if (k >= volumeDim.z)
    return;

  float2 ip;
  float  ipz;

  // matrix multiply
  ip.x = matrix[0]*i + matrix[1]*j + matrix[ 2]*k + matrix[ 3];
  ip.y = matrix[4]*i + matrix[5]*j + matrix[ 6]*k + matrix[ 7];
  ipz  = matrix[8]*i + matrix[9]*j + matrix[10]*k + matrix[11];
  ipz = 1 / ipz;
  ip.x = ip.x * ipz;
  ip.y = ip.y * ipz;

  // Get projection value and add
  float4 projectionData = read_imagef(projection, projectionSampler, ip);

  volume[volumeIndex] += ipz * ipz * projectionData.x;
}


__kernel void
padding_kernel(__global float *input,
	__global float *output,
	const int4 paddingIdx,
	const uint4 paddingDim,
	const uint4 inputDim,
	__constant float *truncationWeights,
	const uint sizeWeights)
{
	const unsigned int idx = get_global_id(0);
	int k = idx / (paddingDim.x * paddingDim.y);

	const unsigned int w = idx % (paddingDim.x * paddingDim.y);
	int j = w / paddingDim.x;
	int i = w % paddingDim.x;

	//if (i >= paddingDim.x || j >= paddingDim.y || k >= paddingDim.z) {
	//	printf("Is it possible to get here? (%d, %d, %d)", i, j, k);
	//	return;
	//}

	const unsigned long int out_idx = i + (j + k*paddingDim.y) * paddingDim.x;
	i -= paddingIdx.x;
	j -= paddingIdx.y;
	k -= paddingIdx.z;

	// out of input y/z dimensions
	if (j < 0 || j >= inputDim.y || k < 0 || k >= inputDim.z)
		output[out_idx] = 0.0f;
	// central part in CPU code
	else if (i >= 0 && i < inputDim.x)
		output[out_idx] = input[i + (j + k*inputDim.y) * inputDim.x];
	// left mirroring (equation 3a in [Ohnesorge et al, Med Phys, 2000])
	else if (i < 0 && -i < sizeWeights)
	{
		const int begRow = (j + k*inputDim.y) * inputDim.x;
		output[out_idx] = (2 * input[begRow + 1] - input[-i + begRow]) * truncationWeights[-i];
	}
	// right mirroring (equation 3b in [Ohnesorge et al, Med Phys, 2000])
	else if ((i >= inputDim.x) && (i - inputDim.x + 1) < sizeWeights)
	{
		const unsigned int borderDist = i - inputDim.x + 1;
		const int endRow = inputDim.x - 1 + (j + k*inputDim.y) * inputDim.x;
		output[out_idx] = (2 * input[endRow] - input[endRow - borderDist]) * truncationWeights[borderDist];
	}
	// zero padding
	else
		output[out_idx] = 0.0f;
}


__kernel void multiply_kernel(
	__global float2 *projFFT,
	int4 fftDimension,
	__global float2 *kernelFFT)
{
	const unsigned int idx = get_global_id(0);
	int k = idx / (fftDimension.x * fftDimension.y);

	const unsigned int w = idx % (fftDimension.x * fftDimension.y);
	int j = w / fftDimension.x;
	int i = w % fftDimension.x;

	//if (i >= fftDimension.x || j >= fftDimension.y || k >= fftDimension.z) {
	//	printf("Is it possible to get here? (%d, %d, %d)", i, j, k);
	//	return;
	//}

	long int proj_idx = i + (j + k * fftDimension.y) * fftDimension.x;

	float2 result;
	result.x = projFFT[proj_idx].x * kernelFFT[i].x - projFFT[proj_idx].y * kernelFFT[i].y;
	result.y = projFFT[proj_idx].y * kernelFFT[i].x + projFFT[proj_idx].x * kernelFFT[i].y;
	projFFT[proj_idx] = result;
}

__kernel void multiply_kernel2D(
	__global float2 *projFFT,
	int4 fftDimension,
	__global float2 *kernelFFT)
{
	const unsigned int idx = get_global_id(0);
	int k = idx / (fftDimension.x * fftDimension.y);

	const unsigned int w = idx % (fftDimension.x * fftDimension.y);
	int j = w / fftDimension.x;
	int i = w % fftDimension.x;

	//if (i >= fftDimension.x || j >= fftDimension.y || k >= fftDimension.z) {
	//	printf("Is it possible to get here? (%d, %d, %d)", i, j, k);
	//	return;
	//}

	long int kernel_idx = i + j * fftDimension.x;
	long int proj_idx = kernel_idx + k * fftDimension.y * fftDimension.x;

	float2 result;
	result.x = projFFT[proj_idx].x * kernelFFT[kernel_idx].x - projFFT[proj_idx].y * kernelFFT[kernel_idx].y;
	result.y = projFFT[proj_idx].y * kernelFFT[kernel_idx].x + projFFT[proj_idx].x * kernelFFT[kernel_idx].y;
	projFFT[proj_idx] = result;
}

__kernel void subtract_kernel2D(
	__global float *input,
	uint4 inputDimension,
	__global float *subImg)
{
	const unsigned int idx = get_global_id(0);
	// int k = idx / (inputDimension.x * inputDimension.y);

	const unsigned int w = idx % (inputDimension.x * inputDimension.y);
	//int j = w / inputDimension.x;
	//int i = w % inputDimension.x;

	//long int sub_idx = i + j * inputDimension.x;

	input[idx] = input[idx] - subImg[w];;
}

// Actually is divide ln(65535/X) by ln(65535/Y) 
__kernel void divide_kernel3D_Ushort(
	__global ushort *input,
	__global ushort *divImg,
	__global float *output)
{
	const unsigned int idx = get_global_id(0);
	//mu_t = ln(65535/I)
	const float short_lim = log(65535.0f);
	const float mu_in = short_lim - log((float)input[idx]);
	const float mu_div = short_lim - log((float)divImg[idx]);
	output[idx] = mu_in / mu_div;
}

__kernel void add_const_kernel(
	__global float *input,
	float value)
{
	const unsigned int idx = get_global_id(0);

	input[idx] += value;
}

__kernel void min_max_kernel(
	__global const float *input,
	__global float2 *subImg,
	const uint divider)
{
	const uint idx = get_global_id(0);
	subImg[idx].x = 65535.0f; //min
	subImg[idx].y = -9999.0f; //max

	const uint first_idx_in = idx * divider;

	for(uint i = first_idx_in; i<(first_idx_in+divider); i++){ 
		if (subImg[idx].x > input[i])
			subImg[idx].x = input[i];
		if (subImg[idx].y < input[i])
			subImg[idx].y = input[i];
	}
}

__kernel void min_max_kernel2(
	__global const float2 *input,
	__global float2 *subImg,
	uint divider,
	const uint end_idx)
{
	const uint idx = get_global_id(0);
	subImg[idx].x = 65535.0f; //min
	subImg[idx].y = -9999.0f; //max

	const uint first_idx_in = idx * divider;
	
	if ((first_idx_in + divider) > end_idx)
			divider = end_idx - first_idx_in;

	for (uint i = first_idx_in; i < (first_idx_in + divider); i++) {
		if (subImg[idx].x > input[i].x)
			subImg[idx].x = input[i].x;
		if (subImg[idx].y < input[i].y)
			subImg[idx].y = input[i].y;
	}
}