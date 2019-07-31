for dir in ../../mixed* ../../pure*; do
    echo $dir
    mv $dir/system.info $dir/aa_system.info
    python get_system_info.py $dir/branch1/5_npt/cg/start.hoomdxml $dir/cg_system.info
done

