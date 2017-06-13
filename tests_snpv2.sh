echo CASO 1___________________________________________________________________________________________
#_________________________________Case 1: Mutant in ref background, backcross, mutant parental control (ref bc mut)_______________________________________________________________
./easymap -n caso1 -w snp -sim -r at -g chr1+4.gff -rs ref -cr bc -co par_mut -sm 200+e -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,2010324+r+100 -ss 25+100,0+500,100+1+50+se


echo CASO 2/5___________________________________________________________________________________________
#_________________________________Case 2 and case 5: Backcross, f2wt control (bc mut)_________________________________________________________________________________________
./easymap -n caso5 -w snp -sim -r at -g chr1+4.gff -rs noref -cr bc -co f2wt -sm 200+e -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,2010324+r+100 -ss 25+100,0+500,100+1+50+se
./easymap -n caso2 -w snp -sim -r at -g chr1+4.gff -rs ref -cr bc -co f2wt -sm 200+e -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,2010324+r+100 -ss 25+100,0+500,100+1+50+se

echo CASO 3___________________________________________________________________________________________
#_________________________________Case 3: Mutant in ref background, outcross, mutant parental control (ref oc mut)________________________________________________________________
./easymap -n caso3 -w snp -sim -r at -g chr1+4.gff -rs ref -cr oc -co par_mut -sm 200+e -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,2010324+r+100 -ss 25+100,0+500,100+1+50+se

echo CASO 4___________________________________________________________________________________________
#_________________________________Case 4: Mutant in ref background, outcross, wt parental control (ref oc nomut)________________________________________________________________
./easymap -n caso4 -w snp -sim -r at -g chr1+4.gff -rs ref -cr oc -co par_nomut -sm 200+e -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,2010324+r+100 -ss 25+100,0+500,100+1+50+se


echo CASO 6___________________________________________________________________________________________
#_________________________________Case 6: Mutant in noref background, outcross, mutant parental control (noref oc mut)________________________________________________________________
./easymap -n caso6 -w snp -sim -r at -g chr1+4.gff -rs noref -cr oc -co par_mut -sm 200+e -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,2010324+r+100 -ss 25+100,0+500,100+1+50+se