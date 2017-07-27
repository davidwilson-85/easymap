<?php


$command = './easymap.sh project-name snp sim se genome.fa n/p n/p n/p n/p chr1.gff n/p 150+e 0,14-1,31-2,33-3,15-4,5-5,2+1,10000000+r+10 2+200,40+0,0+1+100 n/p n/p n/p oc ref mut f2wt & echo $!';

$output = shell_exec($command);

// Example to run a script from an inner directory
//shell_exec('cd ../..; '. $command);
shell_exec($command);

echo 'output: '. $output;

echo 'Done.'



?>
