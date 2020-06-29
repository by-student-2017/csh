#!/bin/csh -f 

set POSCAR = $1
set cycle = $2

cp ${POSCAR} tmp_${POSCAR}
sed -i 's/[ \t]*$//' tmp_${POSCAR}
set nl_all = `grep . tmp_${POSCAR} | wc -l`
set nl = `echo "${nl_all}-8" | bc`

touch config_all
sed 's/restart=false/restart=true/g' xcontrol > xcontrol_restart
#
xtb tmp_${POSCAR} --gfn 0 --opt
set n = 1
while ( $n <= ${cycle} )
  xtb tmp_${POSCAR} --gfn 0 --grad
  ./xtb2force tmp_${POSCAR}
  cat config >> config_all
  if ( $n >= 2 ) then
    xtb tmp_${POSCAR} --gfn 0 --md --restart --input xcontrol_restart
  else
    xtb tmp_${POSCAR} --gfn 0 --md --restart --input xcontrol
  endif
  cat tmp_${POSCAR} | head -8 > tmp.poscar
  sed -i 's/[ \t]*$//' xtb.trj
  cat xtb.trj | tail -n ${nl} | cut -c 3- >> tmp.poscar
  cat tmp.poscar > tmp_${POSCAR}
  @ n = $n + 1
end
rm tmp_${POSCAR} tmp.poscar
