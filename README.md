# manthan-preprocess


As a preprocessing step [Manthan](https://github.com/meelgroup/manthan) finds constant functions, called as unate. This is a framework to find 
positive and negative unates. The framework are build on code of [Arjun](https://github.com/meelgroup/arjun), 
and uses [CryptoMiniSAT](https://github.com/msoos/cryptominisat) as underlying SAT solver. 


## How to Build
To build on Linux, you will need the following:
```
sudo apt-get install build-essential cmake
sudo apt-get install zlib1g-dev libboost-program-options-dev libboost-serialization-dev
```

Then, build CryptoMiniSat, Louvain-Community, and Arjun:
```
git clone https://github.com/priyanka-golia/cryptominisat
cd cryptominisat
git checkout for-manthan-preprocess
mkdir build && cd build
cmake ..
make
sudo make install
sudo ldconfig

cd ../..
git clone https://github.com/meelgroup/louvain-community
cd louvain-community
mkdir build && cd build
cmake ..
make -j4
sudo make install
sudo ldconfig

cd ../..
git clone https://github.com/meelgroup/arjun
cd arjun
mkdir build && cd build
cmake ..
make
sudo make install
sudo ldconfig
```
