#!/bin/csh -f 

set POSCAR = $1
set cycle = $2

cp ${POSCAR} tmp_${POSCAR}
sed -i 's/[ \t]*$//' tmp_${POSCAR}
set nl_all = `grep . tmp_${POSCAR} | wc -l`
set nl = `echo "${nl_all}-8" | bc`

touch config_all
#sed 's/restart=false/restart=true/g' xcontrol > xcontrol_restart
#
xtb tmp_${POSCAR} --periodic --md --input xcontrol
set n = 1
while ( $n <= ${cycle} )
  xtb tmp_${POSCAR} --periodic --grad
  ./xtb2force tmp_${POSCAR}
  cat config >> config_all
  xtb tmp_${POSCAR} --periodic --md --restart --input xcontrol_restart
  cat tmp_${POSCAR} | head -8 > tmp.poscar
  sed -i 's/[ \t]*$//' xtb.trj
  cat xtb.trj | tail -n ${nl} | cut -c 3- >> tmp.poscar
  cat tmp.poscar > tmp_${POSCAR}
  @ n = $n + 1
end
rm tmp_${POSCAR} tmp.poscar
