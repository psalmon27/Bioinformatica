#!/bin/bash

var=$(ls *.lock | wc -l 2>> error.log)
echo "$var"

if [[ $var -eq 0 ]];
then
touch archivo.lock
chmod 777 archivo.lock

flag=0
flag2=0


dpkg -s emboss 1>> packages.log 2>> error.log
if [[ $? -ne 0 ]];
then
read -p "emboss was not found. Do you want to install emboss (1 for yes, 0 for no):" RTA2
if [[ RTA2 -eq 1 ]];
then
pip install emboss 2>> error.log 1>> success.log
else
echo "This code can't run without emboss, please make sure to install the library before running"
flag=1
fi
fi

dpkg -s muscle 2>> error.log 1>> packages.log
if [[ $? -ne 0 ]];
then
read -p "Muscle was not found. Do you want to install Muscle (1 for yes, 0 for no):" RTA
if [[ RTA -eq 1 ]];
then
pip install muscle 2>> error.log 1>> success.log
else
echo "This code can't run without Muscle, please make sure to install the library before running"
flag2=1
fi
fi

if [[ $flag -eq 0 && $flag2 -eq 0 && $flag3 -eq 0 ]];
then
python3 script.py 2>> error.log
fi

rm *.lock
echo "Done running, check for the creation of new files"


else
echo "Can't run the code because it's already running in another terminal"
fi
