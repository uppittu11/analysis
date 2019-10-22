outfile=$1
total=$2
n_old=0
diffr1=0
diffr2=0
while true; do
    echo -en "\r"
    n=$( grep -c imaframe $1 )

    if [ $n -eq 0 ]
    then
        echo -en "CALM DOWN! its still loading!"
        sleep 5
        continue
    fi

    if [ $n -eq $total ]
    then
        break
    fi

    if [ $n -ne $n_old ]
    then
        (( diffr = ( $total - $n ) * 5 / ( $n - $n_old ) / 60  ))
        (( pct = 100 * $n / $total ))
    fi

    echo -en "we're on frame "
    echo -en $n
    echo -en " of "
    echo -en $total
    echo -en " ("
    echo -en $pct
    echo -en "%)"

    (( avg = ( $diffr + $diffr1 + $diffr2 ) / 3 ))

    echo -en ". It will take "
    echo -en $avg
    echo -en " min to finish   "

    diffr2=$diffr1
    diffr1=$diffr

    n_old=$n
    sleep 5
done

sleep 10
tail -2 progress.log
