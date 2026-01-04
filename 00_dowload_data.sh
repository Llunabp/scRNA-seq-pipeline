#!/bin/bash 

set -e

wget --show-progress -O data.zip "https://www.dropbox.com/scl/fi/o0m64ck3y6x11nd6xpc1v/data.zip?rlkey=sw8g5avhog9vspv8muz73hmsq&st=gyb5jfqu&dl=0"
unzip data.zip
rm data.zip
