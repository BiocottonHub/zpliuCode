# @Author: zpliu 
# @Date: 2020-05-05 21:44:08 
# @Last Modified by:   zpliu 
# @Last Modified time: 2020-05-05 21:44:08 
# Usage:
# awk  -v fild1=7 -v fild2=8 -f log2.awk geneExpress.txt >out 



{

if($fild1<=1&&$fild2<=1)
    {
        print $0"\t"0
    }

if($fild1>1 && $fild2<=1)
    {	
		print $0"\t"log($fild1/1)/log(2)
    }
if($fild1<=1 && $fild2>1)
    {
        print $0"\t"log(1/$fild2)/log(2)
    }

if($fild1>1 && $fild2>1)
    {
		print $0"\t"log($fild1/$fild2)/log(2)
    }
}
