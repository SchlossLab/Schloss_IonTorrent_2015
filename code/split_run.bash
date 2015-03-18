VARIABLE=$1
FILES=$(make print-$1 | sed -E 's/.*=//')

rm -f temp

for F in $FILES
do
	echo "make $F" >> temp
done


split -l 1 temp

for X in x??
do
	cat code/head.batch $X code/tail.batch > $X.qsub
	qsub $X.qsub
	rm $X.qsub $X
done

rm temp
