#data
#data
file=$1
train=$file/train
test=$file/test
eval=$file/evaluation

# train=~/Project/Benckmark/data/real_data/rarefied/train
# test=~/Project/Benckmark/data/real_data/rarefied/test
# real=~/Project/Benckmark/data/real_data/rarefied/evaluation

OUT=~/Project/Benckmark/output/
filename=`basename $file`
trainname=`basename $train`
testname=`basename $test`
evalname=`basename $eval`


#change each parameter 
for i in 0 50 62
do
    for j in  `seq 0 1 1`
    do
	for k in `seq 0 1 1`
	do
	    for l in `seq 0 1 1`
	    do
		for m in `seq 0 1 1`
		do
		    ##modify input format
		    trainfile=$OUT${trainname}_${i}_${j}_${k}_${l}_${m}
		    testfile=$OUT${testname}_${i}_${j}_${k}_${l}_${m}
		    evalfile=$OUT${evalname}_${i}_${j}_${k}_${l}_${m}

                    seq2inp_new -bl $i -se $j -cp $k -si $l -ch $m $train |grep -v "#" > $trainfile
                    seq2inp_new -bl $i -se $j -cp $k -si $l -ch $m $test |grep -v "#" > $testfile
                    seq2inp_new -bl $i -se $j -cp $k -si $l -ch $m $eval |grep -v "#" > $evalfile

                    ##training(2 neurons)
                    nnbackprop -nh 2 -tf $testfile -syn ${trainfile}.dat $trainfile

		    #evaluating
		    printf `basename $evalfile`"\t" >> ${filename}_eval
		    nnforward -s ${trainfile}.dat $evalfile | grep -v "#" | gawk '{print $1,$3}' |evaluate >> ${filename}_eval
		done
	    done
	done
    done
done

wait

for f in ./output/*.dat
do
head -n 1 $f
done > info

for f in ./output/*.dat
do
awk 'NR==2{print $1}' $f
done > info2

wait

paste info info2 > $filename
cat $filename |tr ' ' '\t' >tmp
mv -f tmp $filename

paste $filename ${filename}_eval > tmp 
mv -f tmp $filename

rm info info2 ${filename}_eval
