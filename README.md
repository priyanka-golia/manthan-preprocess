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

Then, build CryptoMiniSat, Louvain-Community, and manthan-preprocess:
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
mkdir build && cd build
cmake ..
make
sudo make install
sudo ldconfig
```
Command-line usage
-----

Let's take the file:
```
p cnf 4 7
c ret 1 3 0
c ind 2 4 0
1 3 2 0
-1 -3 -2 0
-1 3 2 0
1 -3 -2 0
-2 4 0
2 -4 0
4 0
```

The file has 4 variables and 5 clauses, this is reflected in the header
`p cnf 4 5` which gives the number of variables as the first number and the number of clauses as the second.
Every clause is ended by '0'. `c ret` represents the universally quantified variables, and `c ind` presents existentially quantified variables. We need to find the positive and negative unates in existentially quantified variables.

```
c executed with command line: ./build/preprocess test.qdimacs
[mis] using seed: 0
[mis] Y set : 2, 4, 
[mis] Orig size Y vars   : 2
[mis] X var set: 1, 3, 
[mis] Orig size X vars  : 2
[mis] Setup time: 0
no positive unates 2
no. negative_unate 0
[mis] Total time: 0.003478
``

Here both `2` and `4` variable is positive unate.
