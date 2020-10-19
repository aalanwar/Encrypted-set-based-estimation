# Privacy Preserving Set-Based Estimation Using Partially Homomorphic Encryption

This repo cotains the code and data for our paper  <br />
Amr Alanwar, Victor Gassmann, Xingkang He, Hazem Said, Henrik Sandberg, Karl Henrik Johansson, and Matthias Althoff "Privacy Preserving Set-Based Estimation Using Partially Homomorphic Encryption"<br />
We propose privacy-preserving set-based estimation protocols using partially homomorphic encryption. Set-based estimation constructs a set that guarantees the inclusion of the system state. We represent sets by zonotopes and constrained zonotopes as they can compactly represent high-dimensional sets and are closed under linear maps and Minkowski addition. By selectively encrypting some parameters of the used set representations, we are able to intersect sets in the encrypted domain, which enables guaranteed state estimation while ensuring the privacy goals. In particular, we show that our protocols achieve computational privacy using formal cryptographic definitions of computational indistinguishability. We demonstrate the efficiency of our approach by localizing a mobile quadcopter using custom ultra-wideband wireless devices.

We consider two problem setups:

1- Distributed sensor
Aggregator collects encrypted strips from each sensor and intersects them with the previous estimated set to obtain a new corrected set.  

<br /> <br />
<p align="center">
<img
src="output/meas2.png"
raw=true
alt="Subject Pronouns"
width=500
/>
</p>

2- Distributed sensor groups
Each sensor group manager collects a set of strips from its sensors and intersects them with previous sets. Then share with the aggregator an encrypted set. The aggregator collects encrypted sets from each sensor group and intersects them with the previous estimated set to obtain a new corrected set.  


Please refer to the paper for more technical details. 


## Building


install Visual Studio Code 
add extention c/c++
clone this repo
open Visual Studio Code
file--> open folder (open the repo folder)
Make sure that your matlab path are correct in the makefile CXXFLAGS line
view --> terminal 
sudo apt install libeigen3-dev
sudo apt-get install libboost-dev
sudo apt-get install libntl-dev
sudo apt-get install libgmp-dev

#build the project
make all


make debug



copy the tageted file to amr_test.pp
make all
./HW
This will generate file under "MATLAB/CMatFiles/FILENAME"
run the correspoding file in matlab
This will generate file under "cache/FILENAME"
run 
A-plot_ZS_cppAndMat
B-plot_ZE_cppAndMat

