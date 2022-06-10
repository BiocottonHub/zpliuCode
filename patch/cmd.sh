#! 1. login to get cookie file 
python login.py 

#! 2. get transctipt motif

python patch test.fa test.out 

#! LSF 
#! work dir: /public/home/jqyou/data/PromoteFasta
bsub -q smp -n 1 -R span[hosts=1] -J Ga -e patch.err -o patch.out "python patch/patch.py Ga_pro_2K.fa  Ga_patch_out.txt "
bsub -q smp -n 1 -R span[hosts=1] -J Gr -e patch.err -o patch.out "python patch/patch.py Gr_pro_2K.fa  Gr_patch_out.txt "