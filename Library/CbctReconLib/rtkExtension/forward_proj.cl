
struct Ray {
  float3 o; // origin
  float3 d; // direction
};

// Intersection function of a ray with a box, followed "slabs" method
// http://www.siggraph.org/education/materials/HyperGraph/raytrace/rtinter3.htm
// The function has to be explicitely inlined, otherwise "multiple definition"
// errors will pop during linking. See
// http://choorucode.com/2011/03/15/cuda-device-function-in-header-file/ for
// more information

int intersectBox(struct Ray r, float *tnear, float *tfar, float3 boxMin,
                 float3 boxMax) {
  // Compute intersection of ray with all six bbox planes
  float3 invR;
  invR.x = 1.f / r.d.x;
  invR.y = 1.f / r.d.y;
  invR.z = 1.f / r.d.z;

  float3 T1;
  T1 = invR * (boxMin - r.o);
  float3 T2;
  T2 = invR * (boxMax - r.o);

  // Re-order intersections to find smallest and largest on each axis
  float3 tmin;
  tmin = min(T2, T1);
  float3 tmax;
  tmax = max(T2, T1);

  // Find the largest tmin and the smallest tmax
  float largest_tmin = max(max(tmin.x, tmin.y), max(tmin.x, tmin.z));
  float smallest_tmax = min(min(tmax.x, tmax.y), min(tmax.x, tmax.z));

  *tnear = largest_tmin;
  *tfar = smallest_tmax;

  return smallest_tmax > largest_tmin;
}

inline float3 matrix_multiply(float3 a, __constant float *matrix) {
  return (float3)(
      matrix[0] * a.x + matrix[1] * a.y + matrix[2] * a.z + matrix[3],
      matrix[4] * a.x + matrix[5] * a.y + matrix[6] * a.z + matrix[7],
      matrix[8] * a.x + matrix[9] * a.y + matrix[10] * a.z + matrix[11]);
}

// CONSTANTS // MUST be compiled in with defines
__constant int3 c_projSize = (int3){C_PROJ_SIZE};
__constant float3 c_boxMin = (float3){C_BOX_MIN};
__constant float3 c_boxMax = (float3){C_BOX_MAX};
__constant float3 c_spacing = (float3){C_SPACING};
__constant int3 c_volSize = (int3){C_VOL_SIZE};
__constant float c_tStep = C_TSTEP;

#if !C_RADIUS_ZERO
__constant float c_radius = C_RADIUS;
__constant float c_translatedVolumeTransformMatrices[SLAB_SIZE * 12] = {
    C_TRANS_VOL_TRN_MAT};
#endif

__constant sampler_t projectionSampler =
    CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP | CLK_FILTER_LINEAR;

// Can process stacks of at most SLAB_SIZE projections
// VECTOR_LENGTH must be compiled in

__kernel void kernel_forwardProject_test(__global float *dev_proj_out) {

  unsigned int i = get_global_id(0);
  unsigned int j = get_global_id(1);
  unsigned int numThread = j * c_projSize.x + i;

  if (i >= c_projSize.x || j >= c_projSize.y) {
    return;
  }

  for (unsigned int proj = 0; proj < c_projSize.z; proj++) {
    int projOffset = numThread + proj * c_projSize.x * c_projSize.y;
    for (unsigned int c = 0; c < VECTOR_LENGTH; c++) {
      dev_proj_out[projOffset * VECTOR_LENGTH + c] = 1.f;
      //      dev_proj_in[projOffset * VECTOR_LENGTH + c];
    }
  }
}

// KERNEL kernel_forwardProject
__kernel void kernel_forwardProject(
    __global float *dev_proj_out, __global float *dev_proj_in,
    __global float *dev_vol, __read_only image3d_t tex_vol,
    __constant float
        *c_translatedProjectionIndexTransformMatrices, // SLAB_SIZE * 12
    __constant float *c_sourcePos) {                   // SLAB_SIZE * 3
  unsigned int i = get_global_id(0);
  unsigned int j = get_global_id(1);
  unsigned int numThread = j * c_projSize.x + i;

  if (i >= c_projSize.x || j >= c_projSize.y) {
    return;
  }

  for (unsigned int proj = 0; proj < c_projSize.z; proj++) {
    // Setting ray origin
    struct Ray ray;
    ray.o.x = c_sourcePos[3 * proj];
    ray.o.y = c_sourcePos[3 * proj + 1];
    ray.o.z = c_sourcePos[3 * proj + 2];

#if C_RADIUS_ZERO
    float3 pixelPos = matrix_multiply(
        (float3)(i, j, 0),
        &(c_translatedProjectionIndexTransformMatrices[12 * proj]));
#else
    float3 posProj = matrix_multiply(
        (float3)(i, j, 0),
        &(c_translatedProjectionIndexTransformMatrices[12 * proj]));
    double a = posProj.x / c_radius;
    posProj.x = sin(a) * c_radius;
    posProj.z += (1. - cos(a)) * c_radius;
    float3 pixelPos = matrix_multiply(
        posProj, &(c_translatedVolumeTransformMatrices[12 * proj]));
#endif

    ray.d = pixelPos - ray.o;
    // rsqrt -> inverse square root, native_ available for only float and floatN
    ray.d = ray.d * native_rsqrt(dot(ray.d, ray.d));

    int projOffset = numThread + proj * c_projSize.x * c_projSize.y;

    // Declare variables used in the loop
    float tnear, tfar;

    // Detect intersection with box
    if (!intersectBox(ray, &tnear, &tfar, c_boxMin, c_boxMax) || tfar < 0.f) {
      for (unsigned int c = 0; c < VECTOR_LENGTH; c++) {
        dev_proj_out[projOffset * VECTOR_LENGTH + c] =
            dev_proj_in[projOffset * VECTOR_LENGTH + c];
      }
    } else {
      if (tnear < 0.f) {
        tnear = 0.f; // clamp to near plane
      }

      // Step length in mm
      float3 dirInMM = c_spacing * ray.d;

      float vStep = c_tStep * native_rsqrt(dot(dirInMM, dirInMM));
      float3 step_d = vStep * ray.d;

      // First position in the box
      float halfVStep = 0.5f * vStep;
      tnear = tnear + halfVStep;
      float3 pos = ray.o + tnear * ray.d;

      float sum[VECTOR_LENGTH];
      float sample[VECTOR_LENGTH];
      for (unsigned int c = 0; c < VECTOR_LENGTH; c++) {
        sum[c] = 0.0f;
        sample[c] = 0.0f;
      }

      float t;
      for (t = tnear; t <= tfar; t += vStep) {
        // Read from 3D texture from volume(s)
        float4 tmp_sample = read_imagef(tex_vol, projectionSampler,
                                        (float4)(pos.x, pos.y, pos.z, 0));

        // Accumulate
        sample[0] = tmp_sample.x;
#if VECTOR_LENGTH == 3
        sample[1] = tmp_sample.y;
        sample[2] = tmp_sample.z;
#endif
        for (unsigned int c = 0; c < VECTOR_LENGTH; c++) {
          sum[c] += sample[c];
        }

        // Step forward
        pos += step_d;
      }

      // Update the output projection pixels
      for (unsigned int c = 0; c < VECTOR_LENGTH; c++) {
        dev_proj_out[projOffset * VECTOR_LENGTH + c] =
            dev_proj_in[projOffset * VECTOR_LENGTH + c] +
            (sum[c] + (tfar - t + halfVStep) / vStep * sample[c]) * c_tStep;
      }
    }
  }
}

