rm -f out.txt
mass=1
v=5
i=0
a=0
while [ $a -le 20 ]
do
mass=1
while [ $mass -le 5 ]
do
v=5
while [ $v -le 30 ]
	do
	python3 thbody.py $mass $v $i $a
	v=$(($v+5))
	((i++))
	echo $v
	done
((mass++))
echo $mass
done
a=$(($a+5))
echo $a
done
cp /home/mkap/Files/Astro/hardbinaries/*png /home/mkapy/Files/Astro/hardbinaries/plots/
rm -f *png
