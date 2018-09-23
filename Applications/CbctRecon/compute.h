#ifndef COMPUTE_H
#define COMPUTE_H

#include "cbctrecon.h"

#include <QString>

#include "rtkThreeDCircularProjectionGeometry.h"

template <typename ImageType>
double
GetFOVRadius(const rtk::ThreeDCircularProjectionGeometry::Pointer &geometry,
             const typename ImageType::Pointer &ProjStack);


template <typename T, typename ImageType>
bool GetOutputResolutionFromFOV(
    typename T::SizeType &sizeOutput, typename T::SpacingType &spacing,
    const rtk::ThreeDCircularProjectionGeometry::Pointer &geometry,
    const typename ImageType::Pointer &ProjStack,
    const QString &outputFilePath);

#endif // COMPUTE_H
