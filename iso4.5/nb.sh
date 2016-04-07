#!/bin/bash
echo  "#!/bin/bash" >> script.sh
xc="\"tcp\""
i=0
a=-1.0001
c=0.05
while [ "${i}" != "30" ] 
do
          d=$(echo "$a+$c" | bc)
          echo  "~/work/pro/iso/hsms " >> script.sh
          echo  "sed -i 's/$a/$d/g' input" >> script.sh
#          echo -e "grep ""CORRFUNC_ AW"" ./out.txt | awk '{print \$1 $xc \$2}' >> data.data" >> script.sh
          echo  "grep \"s1s3\" ~/work/pro/iso/corr | awk '{print \$2 $xc \$3}' >> data1.data" >> script.sh
#          echo  "grep \"s1s3\" ~/work/pro/iso111/three/corr/corr | awk '{print \$2 $xc \$3}' >> data2.data" >> script.sh
#          echo  "grep \"s2s3\" ~/work/pro/iso111/three/corr/corr | awk '{print \$2 $xc \$3}' >> data3.data" >> script.sh
#          echo -e "sed -i 's/correlation1cdot3/$a/g' data.data" >> script.sh
          sed -i 's/tcp/\\t/g' script.sh
          a=$(echo "$a+$c" | bc)
          i=$(($i+1))
done

