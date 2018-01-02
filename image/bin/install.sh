#!/usr/bin/env bash

apt update
apt-get install -y wget libssl-dev openssl build-essential 

# install python3
wget https://www.python.org/ftp/python/3.5.0/Python-3.5.0.tgz 
tar xzvf Python-3.5.0.tgz && cd Python-3.5.0 
./configure 
make
make install
ln -s /usr/local/bin/python3 /usr/local/bin/python 

# install pip
apt install -y python3-pip python3-dev pkg-config libfreetype6-dev
pip3 install --upgrade pip

# remove make, configure and wget
apt remove -y build-essential wget
