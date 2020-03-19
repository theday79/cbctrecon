echo "Running Windows CI script"

mkdir build
cd build

cmake -G"Visual Studio 16 2019" -A"x64" .. ^
  -DCMAKE_BUILD_TYPE=RelWithDebInfo ^
  -DCMAKE_PREFIX_PATH="C:/Qt/5.14.1/msvc2017_64/" ^
  -DHUNTER_ENABLED=ON ^
  -DRTK_USE_OPENCL=ON ^
  -DUSE_CUDA=ON ^
  -DITK_DIR="C:/Program Files (x86)/ITK/lib/cmake/ITK-5.1" ^
  -DUSE_SYSTEM_DCMTK=OFF ^
  -DUSE_SYSTEM_Plastimatch=OFF ^
  -DUSE_SYSTEM_ZLIB=OFF ^
  -DUSE_SYSTEM_dlib=OFF

cmake --build . --config RelWithDebInfo -j 12
REM Where N is the number of CPU cores you want to assign to compiling

cmake --build . --config RelWithDebInfo --target INSTALL

ctest -VV -C RelWithDebInfo
