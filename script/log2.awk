# @Author: zpliu 
# @Date: 2020-05-05 21:44:08 
# @Last Modified by:   zpliu 
# @Last Modified time: 2020-05-05 21:44:08 
# Usage:
# awk  -v fild1=7 -v fild2=8 -f log2.awk geneExpress.txt >out 



{

if($fild1<=0.5&&$fild2<=0.5)
    {
        print $0"\t"0
    }

if($fild1>0.5 && $fild2<=0.5)
    {	
		print $0"\t"log($fild1/0.5)/log(2)
    }
if($fild1<=0.5 && $fild2>0.5)
    {
        print $0"\t"log(0.5/$fild2)/log(2)
    }

if($fild1>0.5 && $fild2>0.5)
    {
		print $0"\t"log($fild1/$fild2)/log(2)
    }
}
