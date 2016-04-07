#!/bin/bash
echo  "#!/bin/bash" >> script.sh
echo  "mkdir tem" >> script.sh
xc="\"tcp\""
i=0
j=1
jj=0
a=-1.5001
b=0.0002
c=0.025
d=0.02
while [ "${j}" != "41" ]
do
  bd=$(echo "$b+$d" | bc)
  while [ "${i}" != "100" ] 
  do
    ac=$(echo "$a+$c" | bc)
    echo  "~/work/pro/iso/hsms " >> script.sh
    echo  "sed -i 's/$a/$ac/g' input" >> script.sh
    echo  "grep \"s1s3\" ~/work/pro/iso/corr | awk '{print \$3}' >> tem/data$j" >> script.sh
    sed -i 's/tcp/\\t/g' script.sh
    a=$(echo "$a+$c" | bc)
    i=$(($i+1))
  done
  echo  "sed -i 's/.9999/-1.5001/g' input" >> script.sh
  echo  "sed -i 's/$b/$bd/g' input" >> script.sh
  jj=$(($j-1))
  a=-1.5001
  i=0
  echo "paste tem/data$jj tem/data$j > tem/aaaa" >> script.sh
  echo "rm tem/data$j" >> script.sh
  echo "mv tem/aaaa tem/data$j" >> script.sh
  j=$(($j+1))
  b=$(echo "$b+$d" | bc)
done
