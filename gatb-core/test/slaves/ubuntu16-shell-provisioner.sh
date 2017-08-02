#!/bin/bash
# Shell script for provisioning an Ubuntu 16.04 LTS slave on Inria cloudstack to compile GATB-CORE
# (for use with Jenkins, ci.inria.fr)

set -xv
set -e

# Configure hostname
# ------------------

#HOST_NAME=gatb-core-ubuntu16-docker
HOST_NAME=$1

[ -z "$HOST_NAME"] && { echo "Please give a HOST_NAME argument to this script..."; exit 1; } 

hostnamectl set-hostname $HOST_NAME

# Install necessary packages
# --------------------------

apt-get -y update

apt-get install -y --no-install-recommends \
        vim git wget make zlib1g-dev hwloc \
        doxygen graphviz \
        valgrind libxml2-utils cmake

# Install gcc-4.7 instead of gcc-5 
# --------------------------------
# Note : gcc 5.4.0 is already installed

apt-get install -y gcc-4.7 g++-4.7 gcc-4.7-base

update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.7 60 --slave /usr/bin/g++ g++ /usr/bin/g++-4.7
update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5   40 --slave /usr/bin/g++ g++ /usr/bin/g++-5
update-alternatives  --set gcc /usr/bin/gcc-4.7


# Install cppunit-1.12
# --------------------
# Note: the libcppunit-dev ubuntu package corresponds to version 1.13, which is not compatible with gatb-core

cd ~/    # now in /builds
git clone git://anongit.freedesktop.org/git/libreoffice/cppunit/ cppunit_gcc47
cd cppunit_gcc47
git checkout cppunit-1.12.1

./autogen.sh

./configure LDFLAGS=-Wl,--no-as-needed

make 

make check

make install

# mount point for external hard drive
[ -d /scratchdir ] || { mkdir /scratchdir; chown -R ci:ci /scratchdir; }


# Install Docker
# --------------
# Note: see https://docs.docker.com/engine/installation/linux/docker-ce/ubuntu/ 

# ... Install Docker, part 1

apt-get remove docker docker-engine docker.io

apt-get install \
     apt-transport-https \
     ca-certificates \
     curl \
     software-properties-common

curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -

apt-key fingerprint 0EBFCD88

sudo add-apt-repository \
     "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
     $(lsb_release -cs) \
     stable"

apt-get -y update
apt-get -y install docker-ce

# just to check
docker run hello-world

# ... Install Docker, part 2

getent group docker || groupadd docker
usermod -aG docker ci

# Some cleaning
# -------------

apt-get clean
