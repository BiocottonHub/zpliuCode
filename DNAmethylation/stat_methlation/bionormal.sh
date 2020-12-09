for i in `ls ${1}`
do
bsub -q high -n 1 -J ${i} -R span[hosts=1] " python /data/cotton/zhenpingliu/QingxinSong_GB_DNAmethlation/stat_methlation/extract_read_count.py  ${1}/${i} ${1}/${i}_out"
sleep 1
done
